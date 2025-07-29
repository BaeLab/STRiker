import statistics
from typing import List, Tuple, KeysView, Dict
from collections import defaultdict
import pysam
import re
from openpyxl import Workbook
from openpyxl.styles import Font
from config import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import gaussian_kde
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator


def reverse_complement(seq: str) -> str:
    """
    This function returns the reverse complement of the input sequence
    """
    complement_table = str.maketrans("ATCGNatgc", "TAGCNtacg")
    return seq.translate(complement_table)[::-1]



def is_sense(FLAG:str) -> bool:
    """
    FLAG를 받아서 sense인지 antisense인지 판단한다.
    """
    if FLAG & 16 == 0:
        return True
    else:
        return False


def find_inserted_seq2(cigar:str, read_seq:str, aln_start:int, patho_start:int, patho_end: int) -> list:
    """
    Patho_start와 Patho_end 사이에 있는 insertion들을 찾아서  리스트 형태로 반환한다.
    """
    ins_seq_list = []
    length =""
    base_count = 0 # 이게 중요
    for i in cigar:
        if aln_start >= patho_end:
            break
        if i.isdigit():
            length += i
            continue
        length = int(length)
        # if i == "M":
        if i in ("M", "=", "X"):
            aln_start += length
            base_count += length
        
        elif i == "I":
            if aln_start < patho_start:
                base_count += length
            elif aln_start >= patho_start:
                inserted_seq = read_seq[base_count:base_count+length]
                ins_seq_list.append((inserted_seq, aln_start))
                base_count += length
        
        elif i == "D":
            aln_start += length
        
        elif i == "S":
            base_count += length
            pass

        
        length = ""
    return ins_seq_list



def line_validity_check(line: str) -> bool:
    if line.startswith("@"):
        return False
    line = line.split("\t")
    if line[2] == "*": # Remove unmapped read
        return False
    return True


def is_over_patho_range(patho_range:str, repeat_count:int) -> bool:
    """
    repeat count가 patho_range에 해당하는지 알아본다.
    """
    patho_range = patho_range.strip()
    if "-" in patho_range:
        start, end = map(int, patho_range.split("-"))
        if start <= repeat_count <= end:
            return True
        else:
            return False
    
    else:
        if ">" in patho_range:
            if repeat_count > int(patho_range.split(">")[1].strip()):
                return True
            else:
                return False
        elif "<" in patho_range:
            if repeat_count < int(patho_range.split("<")[1].strip()):
                return True
            else:
                return False
        else:
            if repeat_count == int(patho_range):
                return True
            else:
                return False
            

def process_patho_range(patho_range:str) -> Tuple[int]:
    """
    patho_range를 받아서 start, end로 나눠서 반환한다.
    """
    patho_range = patho_range.strip()
    if "-" in patho_range:
        start, end = map(int, patho_range.split("-"))
        return start, end
    else:
        if ">" in patho_range:
            return int(patho_range.split(">")[1].strip()), None
        elif "<" in patho_range:
            return None, int(patho_range.split("<")[1].strip())
        else:
            return int(patho_range), None


def analyze_pattern_occurrencesV6_option1(long_str:str, motifs:list) -> List[Tuple[str,int]]:
    """
    옵션 1: 현재 위치에서 매칭 가능한 모든 motif 중 가장 긴 것을 무조건 선택
    """
    i = 0
    result = []
    
    # motifs를 길이가 긴 순으로 정렬(가장 긴 motif부터 시도)
    motifs_sorted = sorted(motifs, key=lambda x: (-len(x), x))
    
    # '연속 미매칭' 문자들을 임시로 저장할 버퍼
    unmatched_buffer = []
    
    while i < len(long_str):
        matched = False
        matched_motif = None
        
        # 현재 위치에서 매칭 가능한 모든 motif 찾기
        possible_motifs = []
        for m in motifs_sorted:
            if long_str.startswith(m, i):
                possible_motifs.append(m)
        
        # 가장 긴 motif 선택 (이미 길이순으로 정렬되어 있으므로 첫 번째가 가장 긴 것)
        if possible_motifs:
            matched_motif = possible_motifs[0]
            matched = True
        
        if matched:
            # 미매칭 버퍼 처리
            if unmatched_buffer:
                unmatched_str = "".join(unmatched_buffer)
                if result and result[-1][0] == unmatched_str:
                    result[-1] = (unmatched_str, result[-1][1] + 1)
                else:
                    result.append((unmatched_str, 1))
                unmatched_buffer = []
            
            # 매칭된 motif 추가
            if result and result[-1][0] == matched_motif:
                result[-1] = (matched_motif, result[-1][1] + 1)
            else:
                result.append((matched_motif, 1))
            
            i += len(matched_motif)
        else:
            # 매칭 실패시 unmatched_buffer에 추가
            unmatched_buffer.append(long_str[i])
            i += 1
    
    # 마지막 unmatched_buffer 처리
    if unmatched_buffer:
        unmatched_str = "".join(unmatched_buffer)
        if result and result[-1][0] == unmatched_str:
            result[-1] = (unmatched_str, result[-1][1] + 1)
        else:
            result.append((unmatched_str, 1))
    
    return result


def analyze_pattern_occurrencesV6_option2(long_str:str, motifs:list) -> List[Tuple[str,int]]:
    """
    옵션 2: 앞뒤 컨텍스트를 보고 반복 패턴이 가장 잘 맞는 motif 선택
    """
    i = 0
    result = []
    
    # motifs를 길이가 긴 순으로 정렬
    motifs_sorted = sorted(motifs, key=lambda x: (-len(x), x))
    
    # '연속 미매칭' 문자들을 임시로 저장할 버퍼
    unmatched_buffer = []
    
    def get_context_score(motif, position, window_size=20):
        """주변 컨텍스트에서 해당 motif의 반복 점수 계산"""
        if position < 0 or position >= len(long_str):
            return 0
        
        # 앞뒤 window_size만큼의 컨텍스트 추출
        start = max(0, position - window_size)
        end = min(len(long_str), position + window_size + len(motif))
        context = long_str[start:end]
        
        # 컨텍스트에서 해당 motif가 몇 번 연속으로 등장하는지 계산
        motif_len = len(motif)
        max_consecutive = 0
        current_consecutive = 0
        
        for j in range(0, len(context) - motif_len + 1, motif_len):
            if context[j:j+motif_len] == motif:
                current_consecutive += 1
                max_consecutive = max(max_consecutive, current_consecutive)
            else:
                current_consecutive = 0
        
        return max_consecutive
    
    while i < len(long_str):
        matched = False
        matched_motif = None
        best_score = 0
        
        # 현재 위치에서 매칭 가능한 모든 motif 찾기
        possible_motifs = []
        for m in motifs_sorted:
            if long_str.startswith(m, i):
                possible_motifs.append(m)
        
        if possible_motifs:
            # 각 motif에 대해 컨텍스트 점수 계산
            for motif in possible_motifs:
                score = get_context_score(motif, i)
                # 점수가 같으면 더 긴 motif 우선, 점수가 다르면 점수가 높은 것 우선
                if score > best_score or (score == best_score and (matched_motif is None or len(motif) > len(matched_motif))):
                    best_score = score
                    matched_motif = motif
            
            matched = True
        
        if matched:
            # 미매칭 버퍼 처리
            if unmatched_buffer:
                unmatched_str = "".join(unmatched_buffer)
                if result and result[-1][0] == unmatched_str:
                    result[-1] = (unmatched_str, result[-1][1] + 1)
                else:
                    result.append((unmatched_str, 1))
                unmatched_buffer = []
            
            # 매칭된 motif 추가
            if result and result[-1][0] == matched_motif:
                result[-1] = (matched_motif, result[-1][1] + 1)
            else:
                result.append((matched_motif, 1))
            
            i += len(matched_motif)
        else:
            # 매칭 실패시 unmatched_buffer에 추가
            unmatched_buffer.append(long_str[i])
            i += 1
    
    # 마지막 unmatched_buffer 처리
    if unmatched_buffer:
        unmatched_str = "".join(unmatched_buffer)
        if result and result[-1][0] == unmatched_str:
            result[-1] = (unmatched_str, result[-1][1] + 1)
        else:
            result.append((unmatched_str, 1))
    
    return result


def analyze_pattern_occurrencesV6(long_str:str, motifs:list) -> List[Tuple[str,int]]:
    """
    기본 함수 - 옵션 1을 기본으로 사용
    """
    return analyze_pattern_occurrencesV6_option1(long_str, motifs)


def samfile_process(line:str) -> Tuple[str]:
    """
    sam file의 한 줄을 받아서 필요한 정보들을 추출한다.
    """
    line = line.split("\t")
    flag = int(line[1])
    chrom = line[2]
    cigar = line[5]
    aln_start = int(line[3])
    return (chrom, cigar, aln_start, flag)



def find_efficient_repeated_substrings(
    string, 
    min_length=2, 
    max_length=10, 
    repeat_threshold=10  # 새로 추가된 반복 횟수 threshold
):
    """
    문자열에서 반복되는 부분 문자열을 효율적으로 찾는 함수
    
    Args:
    string (str): 입력 문자열
    min_length (int): 최소 부분 문자열 길이
    max_length (int): 최대 부분 문자열 길이
    repeat_threshold (int): 최소 반복 횟수 threshold
    
    Returns:
    dict: 고유한 반복 부분 문자열과 반복 횟수
    """
    # 부분 문자열 위치를 저장할 해시 맵
    substring_positions = {}
    
    # 고유 문자 구성을 추적할 집합
    unique_char_compositions = set()
    
    # 결과 저장 딕셔너리
    repeated_substrings = {}
    
    # 부분 문자열의 길이별 순회
    for length in range(min_length, max_length + 1):
        # 부분 문자열 위치 초기화
        substring_positions.clear()
        
        # 문자열의 모든 부분 문자열 탐색
        for i in range(len(string) - length + 1):
            substring = string[i:i+length]
            
            # 부분 문자열의 고유한 문자 구성 생성
            char_composition = tuple(sorted((char, substring.count(char)) for char in set(substring)))
            
            # 이 부분 문자열의 모든 위치 저장
            if substring not in substring_positions:
                substring_positions[substring] = []
            substring_positions[substring].append(i)
        
        # threshold 이상 반복되는 부분 문자열 찾기
        for substring, positions in substring_positions.items():
            # threshold 이상 반복되는지 확인
            if len(positions) >= repeat_threshold:
                # 새로운 문자 구성인지 확인
                char_composition = tuple(sorted((char, substring.count(char)) for char in set(substring)))
                
                # 이미 본 문자 구성이 아니면 추가
                if char_composition not in unique_char_compositions:
                    # 다른 부분 문자열에 완전히 포함되지 않도록 확인
                    is_unique = True
                    for existing in repeated_substrings.keys():
                        if substring != existing and substring in existing:
                            is_unique = False
                            break
                    
                    # 고유한 부분 문자열이면 추가
                    if is_unique:
                        repeated_substrings[substring] = len(positions)
                        unique_char_compositions.add(char_composition)
    
    return repeated_substrings





def rotate_string(s):
    """
    문자열의 모든 가능한 회전 버전을 반환
    예: "ABCD" → ["ABCD", "BCDA", "CDAB", "DABC"]
    """
    return [s[i:] + s[:i] for i in range(len(s))]



def find_most_frequent_rotation(s, string, length, start, consecutive_threshold):
    """
    모든 회전 버전 중에서 가장 자주 등장한 버전을 반환
    """
    rotations = rotate_string(s)
    rotation_counts = {rotation: 0 for rotation in rotations}

    for i in range(consecutive_threshold):
        substring = string[start + i*length : start + (i+1)*length]
        if substring in rotations:
            rotation_counts[substring] += 1

    # 등장 횟수가 가장 많은 회전 버전 반환
    return max(rotation_counts, key=rotation_counts.get)


def is_repeated_by(key, candidate):
    """
    key가 candidate의 반복으로 이루어져 있는지 확인
    """
    if len(key) % len(candidate) != 0:
        return False  # 길이가 나누어 떨어지지 않으면 반복 불가능
    repeated = candidate * (len(key) // len(candidate))
    return repeated == key

def filter_keys_by_repetition(input_dict):
    """
    딕셔너리의 key 중 반복으로 구성된 key를 제거하고 작은 key만 남김
    """
    keys = sorted(input_dict.keys(), key=len)  # 길이순으로 정렬 (짧은 것부터 검사)
    valid_keys = set(keys)  # 제거된 키를 추적하기 위한 집합

    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if keys[j] in valid_keys and is_repeated_by(keys[j], keys[i]):
                valid_keys.discard(keys[j])  # 반복으로 구성된 키 제거

    # 필터링된 딕셔너리 생성
    return {key: input_dict[key] for key in valid_keys}





def write_coverage_percent(file_name: str, depth_dict: str, threshold: int = 30):
    """
    depth_dict를 받아서 coverage percentage를 파일로 저장한다.
    """

    def is_over_threshold(depth, threshold):
        if depth >= threshold:
            return True
        else:
            return False
    
    # with open(file_name.replace(".sam", "_coverage.txt"), "w") as cov_f:
    with open(file_name, "w") as cov_f:
        passed_gene = 0
        cov_f.write(f"Gene\tDepth\tOverThreshold({threshold})\n")
        for gene, depth in depth_dict.items():
            if is_over_threshold(depth, threshold):
                cov_f.write(f"{gene}\t{depth}\tPASS\n")
                passed_gene += 1
            else:
                cov_f.write(f"{gene}\t{depth}\tFAIL\n")
        cov_f.write(f"{passed_gene} genes of {len(depth_dict)} passed threshold({threshold})\n")
    print(f"Coverage percentage saved to {file_name}")
    return None



def canonical_rotation(s):
    """
    회전 불변성을 유지하되, 알파벳 순서가 아닌 다른 기준을 사용
    최소 circular string 알고리즘을 사용하여 일관된 대표 형태를 선택
    """
    if not s:
        return s
        
    # 가장 작은 회전 형태를 찾되, 알파벳 순서가 아닌 순환 문자열의 최소 표현을 사용
    rotations = rotate_string(s)
    
    # 원본 문자열 그대로 반환 (회전 불변성은 나중에 그룹화에서 처리)
    # 이렇게 하면 긴 motif도 제대로 처리되고, 원본 형태가 유지됨
    return s

def find_consecutive_repeated_substrings_with_rotation_optimized(string, min_length=2, max_length=10, consecutive_threshold=10) -> dict:
    """
    회전 불변성을 고려하여 연속으로 반복되는 부분 문자열을 찾는 함수의 최적화 버전.
    
    기존에는 각 후보 부분 문자열마다 전체 문자열에 대해 정규표현식 검색을 수행했으나,
    여기서는 슬라이딩 윈도우 방식으로 한 번의 전체 스캔으로 후보에 해당하는 연속 구간을 탐색합니다.
    
    Args:
        string (str): 검색할 전체 문자열.
        min_length (int): 후보 부분 문자열의 최소 길이.
        max_length (int): 후보 부분 문자열의 최대 길이.
        consecutive_threshold (int): 연속 반복 최소 횟수.
        
    Returns:
        dict: {canonical_rotation: 최대 반복 횟수} 형태의 결과.
    """
    consecutive_substrings = {}
    checked = set()  # 이미 canonical rotation을 처리한 후보 저장
    n = len(string)
    
    for L in range(min_length, max_length + 1):
        # 시작 인덱스: 최소 consecutive_threshold개의 블록을 포함할 수 있는 범위
        for start in range(n - L * consecutive_threshold + 1):
            candidate = string[start:start+L]
            can = canonical_rotation(candidate)
            if can in checked:
                continue
            checked.add(can)
            double_can = can + can  # 회전 여부 판별을 위한 미리 계산된 문자열
            max_repeats = 0
            i = 0
            while i <= n - L:
                current_block = string[i:i+L]
                # current_block이 can의 회전인지 판별 (회전이면 double_can 안에 포함됨)
                if current_block in double_can:
                    count = 0
                    # 연속 반복 구간 계산
                    while i <= n - L and string[i:i+L] in double_can:
                        count += 1
                        i += L
                    max_repeats = max(max_repeats, count)
                else:
                    i += 1
            if max_repeats >= consecutive_threshold:
                consecutive_substrings[can] = max_repeats
                
    return consecutive_substrings




def make_motif_dict(bam_file, STR_regions_dict, depth_dict, reference_motif_dict) -> Dict[str, Dict[str, List[int]]]:
    """
    This function takes Sam file handle and returns motif_dict
    motif_dict를 만들 때는 range 전체를 span하는 리드만 고려
    수정된 로직: CONSECUTIVE_THRESHOLD를 만족하는 read가 MINIMUM_MOTIF_COVERAGE 이상이어야 motif로 인정
    """
    # 임시 저장소: 각 motif별로 조건을 만족하는 read 카운트 추적
    temp_motif_candidates = defaultdict(lambda: defaultdict(list))
    motif_dict = defaultdict(lambda: defaultdict(list))
    bam_hdl = pysam.AlignmentFile(bam_file, "rb")
    for gene in STR_regions_dict.keys():
        chrom = STR_regions_dict[gene]["chrom"]
        start = STR_regions_dict[gene]["start"]
        end = STR_regions_dict[gene]["end"]

        # 조금 넓은 부분 span해서 본다
        # start = start - LEFT_TRIM
        # end = end + RIGHT_TRIM
        reference_motifs = reference_motif_dict.get(gene, {})

        # 원하는 부위의 sequence만을 가져온 뒤 insertion된 시퀀스에 대해 motif finding
        for read in bam_hdl.fetch(chrom, start, end):
            if not is_primary_and_mapped(read):
                continue
            if not(read.reference_start <= start and read.reference_end >= end):
                continue
            else:
                depth_dict[gene] += 1
            cigartuples = read.cigartuples # List of tuples [(operation, length), ...]
            # operation codes: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, ...
            query_pos = 0
            ref_pos = read.reference_start


            for op, length in cigartuples:
                if op == 0: # match (M)
                    query_pos += length
                    ref_pos += length
                elif op == 1: # insertion (I)
                    # ref_pos는 증가하지 않고, query_pos만 증가
                    if start <= ref_pos <= end and read.query_sequence is not None:
                        inserted_seq = read.query_sequence[query_pos:query_pos + length]
                        # motif finding
                        de_novo_motif = find_consecutive_base_motifs(inserted_seq, min_length=MINIMUM_MOTIF_LENGTH, max_length=MAXIMUM_MOTIF_LENGTH, consecutive_threshold=CONSECUTIVE_THRESHOLD)
                        if de_novo_motif:
                            de_novo_motif = filter_keys_by_repetition(de_novo_motif)
                            # 임시 저장소에 후보 motif들 저장
                            # 길이 우선 정렬 (긴 motif 우선 처리)
                            # for motif, count in sorted(de_novo_motif.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
                            for motif, count in sorted(de_novo_motif.items(), key=lambda x: (-len(x[0]))):
                                temp_motif_candidates[gene][motif].append(count)
                        else:
                            continue
                    query_pos += length
                elif op == 2: # deletion (D)
                    ref_pos += length
                elif op == 3: # skip (N)
                    ref_pos += length
                elif op == 4: # soft clipping (S)
                    query_pos += length
                elif op == 5: # hard clipping (H)
                    pass
                else:
                    pass
            
            # 현재 gene의 범위에 대해 reference motif를 찾음
            for reference_motif in reference_motifs:
                if read.query_sequence is not None:  # None 체크 추가
                    consecutive_count = count_consecutive_occurrences_optimized(read.query_sequence, reference_motif)
                    if consecutive_count > 0:  # 0보다 큰 경우만 추가
                        motif_dict[gene][reference_motif].append(consecutive_count)
    
    # BAM 파일 처리 완료 후 coverage 필터링 수행
    # print(f"[COVERAGE FILTERING] MINIMUM_MOTIF_COVERAGE = {MINIMUM_MOTIF_COVERAGE}")
    # 재현성을 위해 정렬된 순서로 처리
    # for gene in sorted(temp_motif_candidates.keys()):
    for gene in temp_motif_candidates.keys():
        motif_candidates = temp_motif_candidates[gene]
        # for motif in sorted(motif_candidates.keys()):
        for motif in motif_candidates.keys():
            count_list = motif_candidates[motif]
            # CONSECUTIVE_THRESHOLD를 만족하는 read가 MINIMUM_MOTIF_COVERAGE 이상인지 확인
            if len(count_list) >= MINIMUM_MOTIF_COVERAGE:
                motif_dict[gene][motif] = count_list
                # print(f"[COVERAGE FILTERING] {gene} - {motif}: {len(count_list)} reads ≥ {MINIMUM_MOTIF_COVERAGE} → ACCEPTED")
            # else:
            #     print(f"[COVERAGE FILTERING] {gene} - {motif}: {len(count_list)} reads < {MINIMUM_MOTIF_COVERAGE} → REJECTED")
    
    bam_hdl.close()
    return motif_dict



def make_pattern_dict(bam_file, STR_regions_dict, motif_dict):
    """
    Deletion gap을 *로 채우고, soft clipping을 제거한 뒤 pattern_dict생성
    """
    bam_hdl = pysam.AlignmentFile(bam_file, "rb")
    pattern_dict = defaultdict(list)
    consecutive_repeat_results = defaultdict(lambda: defaultdict(list))
    total_repeat_results = defaultdict(lambda: defaultdict(list))

    for gene in STR_regions_dict.keys():
        chrom = STR_regions_dict[gene]["chrom"]
        start = STR_regions_dict[gene]["start"]
        end = STR_regions_dict[gene]["end"]

        # 조금 넓은 부분 span해서 본다
        # start = start - LEFT_TRIM
        # end = end + RIGHT_TRIM
        for read in bam_hdl.fetch(chrom, start, end):
            if not is_primary_and_mapped(read):
                continue
            if not(read.reference_start <= start and read.reference_end >= end):
                continue
            # if read.reference_start <= start and read.reference_start + read.query_length >= end:
            read_seq = extract_sequence_with_deletionV2(read, start - LEFT_TRIM, end + RIGHT_TRIM)
            motifs = list(motif_dict.get(gene, {}).keys())
            # 옵션 1 방식 사용 (기본)
            pattern_list = analyze_pattern_occurrencesV6_option1(read_seq, motifs)

            pattern_dict[gene].append(pattern_list)
            # 범위내에 있는 read의 motif total갯수를 셈.
            # 재현성을 위해 정렬된 순서로 처리
            for motif, count in sorted(analyze_pattern_total_occurences(read_seq, motifs).items(), key=lambda x: x[0]):
                total_repeat_results[gene][motif].append(count)

            # consecutive repeat에 대해 각 리드별 motif의 연속 등장 최댓값 저장
            # 재현성을 위해 정렬된 순서로 처리
            for motif, max_consec_repeat_num in sorted(get_maximum_consecutive_repeats(pattern_list, motifs).items(), key=lambda x: x[0]):
                consecutive_repeat_results[gene][motif].append(max_consec_repeat_num)
        
    bam_hdl.close()
    # sorting
    get_length = lambda x: sum(len(seq) * count for seq, count in x)
    pattern_dict = dict(sorted(pattern_dict.items(), key=lambda x: x[0], reverse=False))
    for gene in pattern_dict.keys():
        pattern_dict[gene].sort(key=get_length, reverse=True)
    for gene in total_repeat_results.keys():
        total_repeat_results[gene] = dict(sorted(total_repeat_results[gene].items(), key=lambda x: x[0], reverse=False))
        consecutive_repeat_results[gene] = dict(sorted(consecutive_repeat_results[gene].items(), key=lambda x: x[0], reverse=False))
    total_repeat_results = dict(sorted(total_repeat_results.items(), key=lambda x: x[0], reverse=False))
    consecutive_repeat_results = dict(sorted(consecutive_repeat_results.items(), key=lambda x: x[0], reverse=False))

    return (pattern_dict, consecutive_repeat_results, total_repeat_results)


def process_linesV2_option2(bam_file, STR_regions_dict, motif_dict):
    """
    옵션 2: 컨텍스트 기반 패턴 매칭을 사용한 버전
    """
    bam_hdl = pysam.AlignmentFile(bam_file, "rb")
    pattern_dict = defaultdict(list)
    consecutive_repeat_results = defaultdict(lambda: defaultdict(list))
    total_repeat_results = defaultdict(lambda: defaultdict(list))

    for gene in STR_regions_dict.keys():
        chrom = STR_regions_dict[gene]["chrom"]
        start = STR_regions_dict[gene]["start"]
        end = STR_regions_dict[gene]["end"]

        # 조금 넓은 부분 span해서 본다
        start = start - LEFT_TRIM
        end = end + RIGHT_TRIM
        
        for read in bam_hdl.fetch(chrom, start, end):
            if not is_primary_and_mapped(read):
                continue
            read_seq = extract_sequence_with_deletionV2(read, start, end)
            motifs = list(motif_dict.get(gene, {}).keys())
            
            # 옵션 2 방식 사용 (컨텍스트 기반)
            pattern_list = analyze_pattern_occurrencesV6_option2(read_seq, motifs)

            pattern_dict[gene].append(pattern_list)
            # 범위내에 있는 read의 motif total갯수를 셈.
            # 재현성을 위해 정렬된 순서로 처리
            for motif, count in sorted(analyze_pattern_total_occurences(read_seq, motifs).items(), key=lambda x: x[0]):
                total_repeat_results[gene][motif].append(count)

            # consecutive repeat에 대해 각 리드별 motif의 연속 등장 최댓값 저장
            # 재현성을 위해 정렬된 순서로 처리
            for motif, max_consec_repeat_num in sorted(get_maximum_consecutive_repeats(pattern_list, motifs).items(), key=lambda x: x[0]):
                consecutive_repeat_results[gene][motif].append(max_consec_repeat_num)
        
    bam_hdl.close()
    # sorting
    get_length = lambda x: sum(len(seq) * count for seq, count in x)
    pattern_dict = dict(sorted(pattern_dict.items(), key=lambda x: x[0], reverse=False))
    for gene in pattern_dict.keys():
        pattern_dict[gene].sort(key=get_length, reverse=True)
    
    return (pattern_dict, consecutive_repeat_results, total_repeat_results)



def filter_motif_dict(motif_dict, minimum_motif_coverage, reference_motif_dict):
    """
    motif_dict를 받아서 각 motif의 최댓값으로 변환한다.
    새로운 로직: make_motif_dict에서 이미 coverage 필터링이 완료되었으므로 여기서는 최댓값 변환만 수행
    reference motif_dict에 있는 motif은 무조건 유지
    """
    # 재현성을 위해 정렬된 순서로 처리
    for gene in sorted(motif_dict.keys()):
        motif = motif_dict[gene]
        for k in sorted(motif.keys()):
            v = motif[k]
            # 리스트를 최댓값으로 변환
            if isinstance(v, list) and v:
                motif_dict[gene][k] = max(v)
            elif not v:  # 빈 리스트인 경우
                motif_dict[gene][k] = 0
        
        # reference motif는 무조건 유지하되, 없으면 0으로 설정
        for ref_motif in reference_motif_dict.get(gene, {}):
            if ref_motif not in motif_dict[gene]:
                motif_dict[gene][ref_motif] = 0
                print(f"[FILTER] {gene} - Reference motif {ref_motif} added with count 0")
    
    motif_dict = {k: v for k, v in motif_dict.items() if v}  # 빈 딕셔너리 제거
    return motif_dict



def apply_known_motif(motif_dict, STR_regions_dict):
    """
    motif_dict[gene]이 {motif: value} 형태일 때,
    motif가 known_motif의 회전형이면 해당 known_motif로 key를 바꿔서 값을 합산합니다.
    """
    new_motif_dict = {}

    # 재현성을 위해 정렬된 순서로 처리
    for gene in sorted(motif_dict.keys()):
        motif_data = motif_dict[gene]
        known_motif_list = STR_regions_dict[gene]["known_motif"].split("/")
        new_motifs = {}

        # 길이 우선 정렬 (긴 motif 우선 처리)
        for motif, value in sorted(motif_data.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
            replaced = False
            for known_motif in known_motif_list:
                if motif in rotate_string_set(known_motif):
                    new_motifs[known_motif] = new_motifs.get(known_motif, 0) + value
                    replaced = True
                    break
            if not replaced:
                new_motifs[motif] = new_motifs.get(motif, 0) + value
        
        new_motif_dict[gene] = new_motifs

    return new_motif_dict

def apply_known_motif_v2(motif_dict, STR_regions_dict):
    """
    통합 회전 그룹핑: known motif 매핑과 회전 관계 통합을 함께 수행
    1. known motif와 회전 관계에 있는 motif들을 known motif로 통합
    2. 나머지 motif들도 회전 관계에 있으면 하나로 그룹핑
    """
    # 재현성을 위해 정렬된 순서로 처리
    for gene in sorted(motif_dict.keys()):
        motif_data = motif_dict[gene]
        if not motif_data:
            continue
            
        # print(f"[GROUPING {gene}] 처리 전 motifs: {sorted(motif_data.keys(), key=len)}")
            
        known_motif_list = []
        if (gene in STR_regions_dict and 
            STR_regions_dict[gene]["known_motif"] != "None"):
            known_motif_list = STR_regions_dict[gene]["known_motif"].split("/")

        new_motif_data = {}
        processed_motifs = set()

        # 1단계: known motif 우선 처리
        for known_motif in known_motif_list:
            total_count = 0
            mapped_motifs = []
            known_rotations = rotate_string_set(known_motif)
            
            # 길이 우선 정렬 (긴 motif 우선 처리)
            for motif, count in sorted(motif_data.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
                if motif in known_rotations and motif not in processed_motifs:
                    total_count += count
                    mapped_motifs.append(motif)
                    processed_motifs.add(motif)
            
            if total_count > 0:
                new_motif_data[known_motif] = total_count
                # print(f"[GROUPING {gene}] Known motif '{known_motif}' ← {mapped_motifs} (total: {total_count})")

        # 2단계: 남은 motif들의 회전 그룹핑
        remaining_motifs = {k: v for k, v in motif_data.items() if k not in processed_motifs}
        if remaining_motifs:
            # print(f"[GROUPING {gene}] 회전 그룹핑 대상: {sorted(remaining_motifs.keys(), key=len)}")
            grouped_remaining = group_rotational_motifs(remaining_motifs, known_motifs=known_motif_list)
            # print(f"[GROUPING {gene}] 회전 그룹핑 결과: {sorted(grouped_remaining.keys(), key=len)}")
            new_motif_data.update(grouped_remaining)
        
        # print(f"[GROUPING {gene}] 최종 결과: {sorted(new_motif_data.keys(), key=len)}")
        motif_dict[gene] = new_motif_data

    return motif_dict




def rotate_string_set(s):
    """
    문자열의 모든 회전 버전을 set으로 반환 (중복 제거)
    """
    return set(s[i:] + s[:i] for i in range(len(s)))


def fill_deletion_gap(read_seq, cigar):
    """
    CIGAR string을 보고 deletion이 있으면 해당 gap을 *로 채워줌.
    추가로 soft clipping이 있으면 해당 부분을 제거함.
    soft clipping은 CIGAR의 첫 부분이나 마지막에만 있다고 가정함.
    """
    # CIGAR 문자열을 (숫자, 연산) 튜플의 리스트로 파싱
    tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    # 첫 토큰이 soft clipping ('S')이면, 해당 길이만큼 시퀀스 앞부분 제거
    if tokens and tokens[0][1] == 'S':
        soft_length = int(tokens[0][0])
        read_seq = read_seq[soft_length:]
        tokens = tokens[1:]  # 첫 토큰 제거

    # 마지막 토큰이 soft clipping ('S')이면, 해당 길이만큼 시퀀스 뒷부분 제거
    if tokens and tokens[-1][1] == 'S':
        soft_length = int(tokens[-1][0])
        read_seq = read_seq[:-soft_length]
        tokens = tokens[:-1]  # 마지막 토큰 제거

    base_count = 0  # read_seq 내에서 현재 위치 (pointer)
    # 나머지 토큰들을 순회하며 deletion gap을 채움
    for length_str, op in tokens:
        length = int(length_str)
        if op == "M" or op == "I":
            base_count += length
        elif op == "D":
            # 현재 위치에서 deletion 길이만큼 '*' 삽입
            read_seq = read_seq[:base_count] + "*" * length + read_seq[base_count:]
            base_count += length
        # 다른 연산자(op)가 있다면 필요한 경우 처리할 수 있음.
    return read_seq

def get_sequence_from_fasta(fasta_path: str, chrom: str, start: int, end: int) -> str:
    """
    FASTA 파일과 chromosome, 시작 및 끝 위치가 주어졌을 때 해당 범위의 시퀀스를 반환합니다.
    
    주의: pysam은 0-base coordinate를 사용합니다. 
           (예를 들어, 1-base coordinate를 사용하는 경우 start-1로 조정)
    """
    # pysam.FastaFile을 사용하면 인덱스(.fai)가 있으면 해당 부분만 읽습니다.
    with pysam.FastaFile(fasta_path) as fasta:
        sequence = fasta.fetch(chrom, start, end)
    return sequence


            

def minimal_repeating_unit(s):
    """
    문자열 s가 반복 패턴으로 구성되어 있다면, 가장 단순한(최소) 반복 단위를 반환합니다.
    예: "AAGAAG" → "AAG"
    """
    L = len(s)
    for d in range(1, L + 1):
        if L % d == 0:
            unit = s[:d]
            if unit * (L // d) == s:
                return unit
    return s



def find_consecutive_base_motifs(string, min_length=3, max_length=30, consecutive_threshold=3) -> dict:
    """
    문자열 내에서 연속 반복 영역을 탐색할 때,
    후보 블록이 완전한 반복(즉, 한 번에 단순 모티프 하나로 구성됨)인 경우에만
    base motif (회전 불변성을 적용하여 canonical한 형태)와 그 반복 횟수를 반환합니다.
    
    [조건]
      - 후보 블록은 string[i:i+m]를 사용하며, m은 min_length부터 max_length까지 시도합니다.
      - 후보 블록의 minimal_repeating_unit를 구한 후, 만약 후보 자체가 단순 모티프(즉, m == len(minimal_repeating_unit(candidate)))가 아니라면
        해당 후보는 건너뜁니다.
      - 또한, 같은 반복 구간의 중복 계산을 피하기 위해, 해당 모티프의 반복 구간이 시작하는 위치(반복 run의 시작)에서만 계산합니다.
    
    예를 들어, 
        seq = "ATTTT" * 19
    인 경우, 모티프 "ATTTT"가 19번 연속 반복되는 영역만을 검출하여
        {"ATTTT": 19}
    를 반환하게 됩니다.
    """
    result = {}
    n = len(string)
    
    # 모든 가능한 시작 위치와 모티프 길이에 대해 검사
    for i in range(n):
        for m in range(min_length, max_length+1):
            if i + m > n:
                break
            candidate = string[i:i+m]
            # base: candidate가 완전한 반복 구조라면 그 단순 모티프가 됨
            base = minimal_repeating_unit(candidate)
            # 만약 candidate가 이미 여러 복제된 형태라면(예: m != len(base)) 건너뜁니다.
            # 즉, 한 블록은 반드시 단순 모티프 하나여야 함.
            if m != len(base):
                continue
            
            # 반복 영역의 시작점에서만 계산하기 위해,
            # 바로 앞 m 길이 블록과 동일하다면 이미 해당 영역을 처리한 것으로 간주.
            if i - m >= 0 and string[i-m:i] == candidate:
                continue

            # candidate와 동일한 블록이 몇 번 연속 등장하는지 검사합니다.
            count = 0
            pos = i
            while pos + m <= n and string[pos:pos+m] == candidate:
                count += 1
                pos += m

            if count >= consecutive_threshold:
                # 원본 motif를 그대로 저장 (회전 불변성은 후처리에서 처리)
                result[candidate] = max(result.get(candidate, 0), count)
    return result


def count_consecutive_occurrences_optimized(long_string: str, pattern: str) -> int:
    """
    긴 문자열(long_string)에서 특정 문자열(pattern)이 연속해서 몇 번 등장하는지 찾는 최적화된 코드.

    Args:
        long_string (str): 검색할 긴 문자열.
        pattern (str): 찾을 패턴.

    Returns:
        int: pattern이 연속으로 반복된 최대 횟수.
    """
    # None 체크와 빈 문자열 체크 추가
    if long_string is None or pattern is None or not long_string or not pattern:
        return 0
        
    max_count = 0
    i = 0
    n = len(long_string)
    m = len(pattern)

    while i <= n - m:
        count = 0
        while i <= n - m and long_string[i:i+m] == pattern:
            count += 1
            i += m  # 패턴 길이만큼 건너뜀 (연속된 패턴 찾기)

        max_count = max(max_count, count)
        
        # 연속된 패턴이 없었으면 한 칸만 이동, 아니면 다음 검사 진행
        if count == 0:
            i += 1

    return max_count

def save_motif_as_xlsx(
                        output_file,
                        motif_dict,
                        reference_motif_dict_consc,
                        reference_motif_dict_total,
                        consecutive_repeat_result,
                        total_repeat_results,
                        depth_dict):
    """
    motif정보를 xlsx로 저장한다.
    reference에서 연속으로 등장한 motif의 횟수.
    reference에서 총 등장한 motif의 횟수.
    Sample에서 연속으로 등장한 motif의 통계값.
    Sample에서 총 등장한 motif의 통계값.
    """
    wb = Workbook()
    ws = wb.active
    header = ["Gene",
            "Motif type", 
            "Motif", 
            "Total count in reference",
            "Total count in Sample (mean)",
            "Consecutive count in Reference",
            "Consecutive count in Sample (mean)",
            "Read Depth",
            "Total count in Sample (max)",
            "Total count in Sample (min)",
            "Total count in Sample (median)",
            "Total count in Sample (variance)",
            "Total count in Sample (standard deviation)",
            "Consecutive count in Sample (max)",
            "Consecutive count in Sample (min)",
            "Consecutive count in Sample (median)",
            "Consecutive count in Sample (variance)",
            "Consecutive count in Sample (standard deviation)"]
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)

    # 재현성을 위해 정렬된 순서로 처리
    for gene in sorted(motif_dict.keys()):
        motifs = motif_dict[gene]
        for motif in sorted(motifs.keys()):
            count = motifs[motif]
            if motif in reference_motif_dict_consc.get(gene, {}):
                motif_type = "Reference"
                ref_total = reference_motif_dict_total.get(gene, {}).get(motif, 0)
                
                # 안전한 데이터 접근
                total_data = total_repeat_results.get(gene, {}).get(motif, [0])
                consc_data = consecutive_repeat_result.get(gene, {}).get(motif, [0])
                
                # 리스트가 비어있으면 [0]으로 설정
                if not total_data:
                    total_data = [0]
                if not consc_data:
                    consc_data = [0]
                
                read_total_mean = statistics.mean(total_data)
                ref_consc = reference_motif_dict_consc.get(gene, {}).get(motif, 0)
                read_consc_mean = statistics.mean(consc_data)
                depth = depth_dict.get(gene, 0)
                read_total_max = max(total_data)
                read_total_min = min(total_data)
                read_total_median = statistics.median(total_data)
                read_total_variance = statistics.variance(total_data) if len(total_data) > 1 else 0
                read_total_stdev = statistics.stdev(total_data) if len(total_data) > 1 else 0
                read_consc_max = max(consc_data)
                read_consc_min = min(consc_data)
                read_consc_median = statistics.median(consc_data)
                read_consc_variance = statistics.variance(consc_data) if len(consc_data) > 1 else 0
                read_consc_stdev = statistics.stdev(consc_data) if len(consc_data) > 1 else 0
                ws.append([gene,
                           motif_type,
                           motif,
                           ref_total,
                            read_total_mean,
                           ref_consc,
                            read_consc_mean,
                            depth,
                           read_total_max,
                            read_total_min,
                            read_total_median,
                            read_total_variance,
                            read_total_stdev,
                           read_consc_max,
                            read_consc_min,
                            read_consc_median,
                            read_consc_variance,
                            read_consc_stdev,
                           ])
            else:
                motif_type = "De novo"
                ref_total = "Under threshold" 
                
                # 안전한 데이터 접근 (De novo motif용)
                total_data = total_repeat_results.get(gene, {}).get(motif, [0])
                consc_data = consecutive_repeat_result.get(gene, {}).get(motif, [0])
                
                # 리스트가 비어있으면 [0]으로 설정
                if not total_data:
                    total_data = [0]
                if not consc_data:
                    consc_data = [0]
                
                read_total_mean = statistics.mean(total_data)
                ref_consc = "Underthreshold" 
                read_consc_mean = statistics.mean(consc_data)
                depth = depth_dict.get(gene, 0)
                read_total_max = max(total_data)
                read_total_min = min(total_data)
                read_total_median = statistics.median(total_data)
                read_total_variance = statistics.variance(total_data) if len(total_data) > 1 else 0
                read_total_stdev = statistics.stdev(total_data) if len(total_data) > 1 else 0
                read_consc_max = max(consc_data)
                read_consc_min = min(consc_data)
                read_consc_median = statistics.median(consc_data)
                read_consc_variance = statistics.variance(consc_data) if len(consc_data) > 1 else 0
                read_consc_stdev = statistics.stdev(consc_data) if len(consc_data) > 1 else 0
                ws.append([gene,
                           motif_type,
                           motif,
                           ref_total,
                            read_total_mean,
                           ref_consc,
                            read_consc_mean,
                            depth,
                           read_total_max,
                            read_total_min,
                            read_total_median,
                            read_total_variance,
                            read_total_stdev,
                           read_consc_max,
                            read_consc_min,
                            read_consc_median,
                            read_consc_variance,
                            read_consc_stdev,
                           ])
    wb.save(output_file)
    print("Saved as", output_file)
    
    # TXT 파일로도 저장 (디버깅용)
    # txt_file = output_file.replace('.xlsx', '_debug.txt')
    # with open(txt_file, 'w') as f:
    #     f.write('\t'.join(header) + '\n')
    #     # 재현성을 위해 정렬된 순서로 처리
    #     for gene in sorted(motif_dict.keys()):
    #         motifs = motif_dict[gene]
    #         for motif in sorted(motifs.keys()):
    #             count = motifs[motif]
    #             if motif in reference_motif_dict_consc.get(gene, {}):
    #                 motif_type = "Reference"
    #                 ref_total = reference_motif_dict_total.get(gene, {}).get(motif, 0)
                    
    #                 # 안전한 데이터 접근
    #                 total_data = total_repeat_results.get(gene, {}).get(motif, [0])
    #                 consc_data = consecutive_repeat_result.get(gene, {}).get(motif, [0])
                    
    #                 # 리스트가 비어있으면 [0]으로 설정
    #                 if not total_data:
    #                     total_data = [0]
    #                 if not consc_data:
    #                     consc_data = [0]
                    
    #                 read_total_mean = statistics.mean(total_data)
    #                 ref_consc = reference_motif_dict_consc.get(gene, {}).get(motif, 0)
    #                 read_consc_mean = statistics.mean(consc_data)
    #                 depth = depth_dict.get(gene, 0)
    #                 read_total_max = max(total_data)
    #                 read_total_min = min(total_data)
    #                 read_total_median = statistics.median(total_data)
    #                 read_total_variance = statistics.variance(total_data) if len(total_data) > 1 else 0
    #                 read_total_stdev = statistics.stdev(total_data) if len(total_data) > 1 else 0
    #                 read_consc_max = max(consc_data)
    #                 read_consc_min = min(consc_data)
    #                 read_consc_median = statistics.median(consc_data)
    #                 read_consc_variance = statistics.variance(consc_data) if len(consc_data) > 1 else 0
    #                 read_consc_stdev = statistics.stdev(consc_data) if len(consc_data) > 1 else 0
                    
    #                 row = [gene, motif_type, motif, ref_total, read_total_mean, ref_consc, read_consc_mean, depth,
    #                        read_total_max, read_total_min, read_total_median, read_total_variance, read_total_stdev,
    #                        read_consc_max, read_consc_min, read_consc_median, read_consc_variance, read_consc_stdev]
    #                 f.write('\t'.join(map(str, row)) + '\n')
    #             else:
    #                 motif_type = "De novo"
    #                 ref_total = "Under threshold" 
                    
    #                 # 안전한 데이터 접근 (De novo motif용)
    #                 total_data = total_repeat_results.get(gene, {}).get(motif, [0])
    #                 consc_data = consecutive_repeat_result.get(gene, {}).get(motif, [0])
                    
    #                 # 리스트가 비어있으면 [0]으로 설정
    #                 if not total_data:
    #                     total_data = [0]
    #                 if not consc_data:
    #                     consc_data = [0]
                    
    #                 read_total_mean = statistics.mean(total_data)
    #                 ref_consc = "Underthreshold" 
    #                 read_consc_mean = statistics.mean(consc_data)
    #                 depth = depth_dict.get(gene, 0)
    #                 read_total_max = max(total_data)
    #                 read_total_min = min(total_data)
    #                 read_total_median = statistics.median(total_data)
    #                 read_total_variance = statistics.variance(total_data) if len(total_data) > 1 else 0
    #                 read_total_stdev = statistics.stdev(total_data) if len(total_data) > 1 else 0
    #                 read_consc_max = max(consc_data)
    #                 read_consc_min = min(consc_data)
    #                 read_consc_median = statistics.median(consc_data)
    #                 read_consc_variance = statistics.variance(consc_data) if len(consc_data) > 1 else 0
    #                 read_consc_stdev = statistics.stdev(consc_data) if len(consc_data) > 1 else 0
                    
    #                 row = [gene, motif_type, motif, ref_total, read_total_mean, ref_consc, read_consc_mean, depth,
    #                        read_total_max, read_total_min, read_total_median, read_total_variance, read_total_stdev,
    #                        read_consc_max, read_consc_min, read_consc_median, read_consc_variance, read_consc_stdev]
    #                 f.write('\t'.join(map(str, row)) + '\n')
    
    # print("Debug TXT file saved as", txt_file)
    # return None


def save_debug_info_to_file(output_folder, bam_file, motif_dict, pattern_dict):
    """
    디버깅 정보를 파일로 저장하는 함수
    """
    import os
    import json
    from datetime import datetime
    
    # 디버그 폴더 생성
    debug_folder = os.path.join(output_folder, "debug_output")
    os.makedirs(debug_folder, exist_ok=True)
    
    sample_name = os.path.basename(bam_file).replace(".bam", "")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 1. Motif 정보를 파일로 저장
    motif_debug_file = os.path.join(debug_folder, f"{sample_name}_motif_debug.txt")
    with open(motif_debug_file, 'w', encoding='utf-8') as f:
        f.write(f"=== Motif Debug Information ===\n")
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for gene in sorted(motif_dict.keys()):
            f.write(f"[GENE: {gene}]\n")
            f.write(f"Total motifs found: {len(motif_dict[gene])}\n")
            
            # motif별 정보 출력
            for motif in sorted(motif_dict[gene].keys(), key=lambda x: (len(x), x)):
                count = motif_dict[gene][motif]
                f.write(f"  - {motif} (length={len(motif)}): max_repeat={count}\n")
            f.write("\n")
    
    
    # 3. Pattern dictionary를 JSON으로 저장 (읽기 쉽게)
    pattern_json_file = os.path.join(debug_folder, f"{sample_name}_pattern_dict.json")
    
    # JSON serializable 형태로 변환
    pattern_dict_serializable = {}
    for gene, reads in pattern_dict.items():
        pattern_dict_serializable[gene] = []
        for read in reads:
            pattern_dict_serializable[gene].append(read)  # read는 이미 [(pattern, count), ...] 형태
    
    with open(pattern_json_file, 'w', encoding='utf-8') as f:
        json.dump(pattern_dict_serializable, f, indent=2, ensure_ascii=False)
    
    # print(f"Debug files saved:")
    # print(f"  - Motif debug: {motif_debug_file}")
    # print(f"  - Pattern debug: {pattern_debug_file}")
    # print(f"  - Pattern JSON: {pattern_json_file}")


def analyze_pattern_total_occurences(read_seq: str, motifs: list) -> List[Tuple[str, int]]:
    """
    read_seq에서 motifs에 있는 motif들이 몇 번 등장하는지 반환
    """
    # 재현성을 위해 길이가 같을 때는 알파벳 순으로 정렬
    motifs = sorted(motifs, key=lambda x: (-len(x), x))
    result = {motif: 0 for motif in motifs}
    for motif in motifs:
        count = read_seq.count(motif)
        result[motif] = count
    return result

def extract_sequence_with_deletion(read, target_start, target_end):
    """
    주어진 read에서 target 범위에 해당하는 시퀀스를 추출하되,
    reference에서 deletion인 부분은 -로 채워줌.
    align시작 위치가 target_start보다 크면 앞부분에 *로 채워줌
    """
    result = ""
    query_pos = 0
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op  in (0, 7, 8):  # M , =, X (match or mismatch)
            for i in range(length):
                if target_start <= ref_pos < target_end:
                    result += read.query_sequence[query_pos]
                query_pos += 1
                ref_pos += 1
        elif op == 1:  # I (insertion, reference에는 없음)
            query_pos += length  # 그냥 건너뜀
        elif op == 2:  # D (deletion, read에는 없음)
            for i in range(length):
                if target_start <= ref_pos < target_end:
                    result += "-"
                ref_pos += 1
        elif op == 4:  # soft clipping
            query_pos += length
        elif op == 5:  # hard clipping (read에 없음)
            continue
        else:
            # 필요하면 다른 CIGAR 연산 처리
            pass

        # 종료 조건: reference 범위를 초과하면 중단
        if ref_pos >= target_end:
            break
    # 만약 read의 시작 위치가 target_start보다 크면 앞부분에 *로 채움
    if read.reference_start > target_start:
        result = "*" * (read.reference_start - target_start) + result

    return result

def get_maximum_consecutive_repeats(pattern_list: list, motifs: list) -> Dict[str, int]:
    """
    pattern_list에서 motif가 몇 번 연속으로 등장하는지 찾는 함수.
    """
    result_dict = {}
    # result_dict 초기화
    for motif in motifs:
        result_dict[motif] = 0

    for motif, count in pattern_list:
        if motif in motifs:
            if count > result_dict[motif]:
                result_dict[motif] = count
    return result_dict
                

def is_primary_and_mapped(read):
    return (
        not read.is_secondary and
        not read.is_supplementary and
        not read.is_unmapped
    )



def gaussian(x, a, mu, sigma):
    return a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def show_kde_v2(ax, bam_file, STR_regions_dict, gene, read_threshold):
    bam_hdl = pysam.AlignmentFile(bam_file, "rb")
    chrom = STR_regions_dict[gene]["chrom"]
    start = STR_regions_dict[gene]["start"]
    end = STR_regions_dict[gene]["end"]
    expansion_number = STR_regions_dict[gene]["expansion_number"]

    motif_size = len(STR_regions_dict[gene]["known_motif"].split("/")[0])
    pathogenic_size = motif_size * expansion_number
    start -= LEFT_TRIM
    end += RIGHT_TRIM
    reference_length = end - start
    read_length_list = []

    for read in bam_hdl.fetch(chrom, start, end):
        if not is_primary_and_mapped(read):
            continue
        if read.reference_start < start and read.reference_end > end:
            seq = extract_sequence_with_deletionV2(read, start, end)
            seq = seq.replace("-", "").replace("*", "")
            read_length_list.append(len(seq))
    bam_hdl.close()

    if len(read_length_list) < read_threshold or np.std(read_length_list) < 1e-3:
        ax.text(0.5, 0.5, "Insufficient or flat data", ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.set_title(f"{gene} - KDE skipped", fontsize=5)
        return


    # Reference length을 0으로 맞추기
    # read_length_list = [x - reference_length for x in read_length_list]
    # reference_length = 0

    kde = gaussian_kde(read_length_list, bw_method=0.5)
    x_vals = np.linspace(min(read_length_list), max(read_length_list), 1000)
    y_vals = kde(x_vals)

    ax.fill_between(x_vals, y_vals, color='skyblue', alpha=0.5)
    ax.plot(x_vals, y_vals, color='skyblue')

    peaks, _ = find_peaks(y_vals, distance=10)
    # colors = ["darkblue", "darkorange", "green", "purple", "brown"]
    colors = ["darkblue", "limegreen", "cadetblue", "darkolivegreen", "brown"]
    legend_elements = []

    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1, label=f"Reference Length={reference_length}"),
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label=f"Pathogenic Length={reference_length + pathogenic_size}")
    ])

    if len(peaks) > 0:
        top_indices = np.argsort(y_vals[peaks])[::-1][:3]
        for i, idx in enumerate(top_indices):
            p = peaks[idx]
            color = colors[i % len(colors)]
            peak_x = x_vals[p]
            prominence = y_vals[p]  # 해당 peak의 높이
            window = int((prominence / max(y_vals)) * len(x_vals) * 0.3)
            window = max(100, min(window, 500))  # 너무 작거나 크지 않게 제한k
            left = max(p - window, 0)
            right = min(p + window, len(x_vals))
            x_fit = x_vals[left:right]
            y_fit = y_vals[left:right]

            try:
                popt, _ = curve_fit(lambda x, a, mu, sigma: a * np.exp(-(x - mu)**2 / (2 * sigma**2)),
                                    x_fit, y_fit,
                                    p0=[y_vals[p], peak_x, 5])
                a_fit, mu_fit, sigma_fit = popt
                fwhm = 2.355 * sigma_fit
                repeat_num = (mu_fit - reference_length) / motif_size

                ax.plot(x_fit, gaussian(x_fit, *popt), '--', color=color, linewidth=0.8, alpha=0.5)
                ax.axvline(mu_fit, color=color, linestyle="--", linewidth=1)
                # ax.fill_betweenx([0, a_fit], mu_fit - fwhm/2, mu_fit + fwhm/2, color=color, alpha=0.2)

                ax.text(mu_fit, a_fit + 0.01 * max(y_vals),
                        f"{mu_fit:.1f} (±{fwhm/2:.1f})", color=color,
                        fontsize=4, ha='center', va='bottom')

                legend_elements.append(Line2D(
                    [0], [0],
                    color=color,
                    linestyle="--",
                    label=f"μ={mu_fit:.1f}, σ={sigma_fit:.1f}, repeats={repeat_num:.1f}"
                ))

            except RuntimeError:
                ax.text(peak_x, y_vals[p], "Fit failed", fontsize=4, color="red")

    else:
        ax.text(0.5, 0.5, "No peaks found", ha='center', va='center', transform=ax.transAxes, fontsize=6)

    ax.axvline(x=reference_length, color="gray", linestyle="--", linewidth=1)
    ax.axvline(x=reference_length + pathogenic_size, color="red", linestyle="--", linewidth=2)


    ax.set_xlabel("Read Length", fontsize=4)
    ax.set_ylabel("Frequency", fontsize=4)
    ax.set_title(f"Read Length Distribution for {gene}", fontsize=7)
    ax.tick_params(axis='both', which='major', labelsize=4)
    ax.legend(handles=legend_elements, fontsize=4, loc='upper right', title_fontsize=4)


def plot_repeat_number_distribution(ax, bam_file, STR_regions_dict, gene, read_threshold):
    """
    각 리드들에 대해서 motif의 등장횟수를 bar로 나타내는 함수.
    - x축: motif 반복 횟수 (repeat number)
    - y축: 해당 repeat number를 가진 read의 개수
    """
    bam_hdl = pysam.AlignmentFile(bam_file, "rb")
    chrom = STR_regions_dict[gene]["chrom"]
    start = STR_regions_dict[gene]["start"]
    end = STR_regions_dict[gene]["end"]
    expansion_number = STR_regions_dict[gene]["expansion_number"]

    # repeat number를 셀 때는 Spanning 필요없음.
    # start -= LEFT_TRIM
    # end += RIGHT_TRIM
    repeat_number_dict = defaultdict(int)
    motif = STR_regions_dict[gene]["known_motif"].split("/")[0]

    for read in bam_hdl.fetch(chrom, start, end):
        if not is_primary_and_mapped(read):
            continue
        if read.reference_start < start and read.reference_end > end:
            seq = extract_sequence_with_deletionV2(read, start, end)
            count = seq.count(motif)
            repeat_number_dict[count] += 1
    bam_hdl.close()
    if sum(repeat_number_dict.values()) < read_threshold:
        ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', transform=ax.transAxes, fontsize=6)
        ax.set_title(f"{gene} - Insufficient data", fontsize=5)
        return

    # 정렬된 repeat 숫자 기준으로 x, y 준비
    x_vals = sorted(repeat_number_dict.keys())
    y_vals = [repeat_number_dict[x] for x in x_vals]
    legend_elements = [
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label=f"Pathogenic {motif} repeat Number={expansion_number}"),
    ]

    if x_vals:
        min_x = min(x_vals + [expansion_number])
        max_x = max(x_vals + [expansion_number])
        ax.bar(x_vals, y_vals, width=0.8, color='skyblue', edgecolor='black')
        ax.set_xlim(min_x - 1, max_x + 1)  # x축 범위 설정
    else: # 데이터가 없는 경우
        ax.bar([expansion_number], [0], width=0.8, color='lightgray', edgecolor='black')
        ax.set_xlim(expansion_number - 2, expansion_number + 2)
        ax.set_ylim(0, 1)
        ax.text(expansion_number, 0.5, "No valid reads", ha='center', fontsize=6, color='gray')




    # 기준선 (확장 기준)
    ax.axvline(x=expansion_number, color="red", linestyle="--", linewidth=2)

    # 시각적 요소
    ax.set_xlabel("Repeat Number", fontsize=8)
    ax.set_ylabel("Read Count", fontsize=8)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))  # x축 눈금 정수로 설정
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # y축 눈금 정수로 설정
    ax.set_title(f"Repeat Number Distribution for {gene}", fontsize=7)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.legend(handles=legend_elements, fontsize=6, loc='upper right', title_fontsize=4)
    

    





def extract_sequence_with_deletionV2(read, target_start, target_end):
    """
    주어진 read에서 target 범위([target_start, target_end))에 해당하는 시퀀스를 추출합니다.
    
    - M, =, X: target 영역 내의 query base를 그대로 사용.
    - I: insertion은 reference 좌표를 이동하지 않으므로, 
         현재 ref_pos가 target 내(또는 target_end에 딱 붙어 있는 경우)면 삽입된 query bases를 그대로 추가.
    - D, N: deletion은 target 내에서는 '-'로 채움.
    - S: soft clipping은 target 내에서는 '-'로 채움.
    - 만약 read의 alignment가 target_start보다 늦게 시작하면 그 전은 '*'로 패딩.
    - target 범위는 reference 좌표 기준이므로, deletion은 길이에 영향을 주지 않으나, insertion은 read의 길이를 늘림.
    """
    result = ""
    query_pos = 0
    ref_pos = read.reference_start

    # 왼쪽 패딩: read가 target_start보다 늦게 시작하면, 그 차이만큼 '*'를 채움.
    if ref_pos > target_start:
        result += "*" * (ref_pos - target_start)

    # CIGAR 튜플 순회
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X: match/mismatch → ref와 query 모두 소비
            # base 단위로 처리
            for i in range(length):
                # 해당 base가 target 영역 내라면 추가
                if target_start <= ref_pos < target_end:
                    result += read.query_sequence[query_pos]
                query_pos += 1
                ref_pos += 1
            # reference 좌표를 소모하는 op이므로, target을 완전히 넘겼으면 break
            if ref_pos > target_end:
                break

        elif op == 1:  # I: insertion → query만 소비, ref_pos 변화 없음
            # 삽입 op는 현재 ref_pos에 “붙어 있음”
            # target 범위 내(또는 ref_pos가 target_end에 딱 맞는 경우)라면 삽입된 시퀀스를 추가
            if target_start <= ref_pos <= target_end:
                result += read.query_sequence[query_pos: query_pos + length]
            query_pos += length

        elif op in (2, 3):  # D, N: deletion/skipped → ref만 소비, query 없음
            for i in range(length):
                if target_start <= ref_pos < target_end:
                    result += "-"  # deletion 부분은 '-'로 채움
                ref_pos += 1
            if ref_pos > target_end:
                break

        elif op == 4:  # S: soft clipping → query만 소비, ref_pos 변화 없음
            # soft clipping이 target 내라면 '-'로 채움
            if target_start <= ref_pos <= target_end:
                result += "-" * length
            query_pos += length

        elif op == 5:  # H: hard clipping → query에도 없으므로 무시
            continue

        else:
            # 기타 op는 무시
            pass

    # 만약 read의 alignment가 target_end 전에 끝났다면,
    # target 끝까지 '*'로 채워줌.
    if ref_pos < target_end:
        result += "*" * (target_end - ref_pos)

    return result


def simplify_pattern_dict(pattern_dict, motif_dict):
    """
    motif_dict에 없는 motif를 가진 tuple들을 하나로 합쳐서 pattern_dict를 단순화함.
    plotly에서의 trace 수를 줄이기 위함.
    """
    simplified_pattern_dict = defaultdict(list)
    for gene in pattern_dict.keys():
        pattern_lists = pattern_dict[gene]
        for i in range(len(pattern_lists)):
            pattern_list = pattern_lists[i]
            simplified_pattern_list = []
            simplified_pattern = []
            for pattern, count in pattern_list:
                if pattern not in motif_dict.get(gene, {}):
                    simplified_pattern.append(pattern * count)
                else:
                    if simplified_pattern:
                        simplified_pattern_list.append(("".join(simplified_pattern), 1))
                        simplified_pattern = []
                    # motif_dict에 있는 motif는 그대로 사용
                    simplified_pattern_list.append((pattern, count))
            if simplified_pattern:
                simplified_pattern_list.append(("".join(simplified_pattern), 1))
                simplified_pattern = []
            simplified_pattern_dict[gene].append(simplified_pattern_list)
    return simplified_pattern_dict



def find_motif_sequence(text, motifs):
    """문자열에서 알려진 motif들의 순서와 위치를 찾는 함수"""
    # motif를 길이 순으로 정렬 (긴 것부터)
    sorted_motifs = sorted(motifs, key=len, reverse=True)
    
    result = []
    position = 0
    
    while position < len(text):
        found = False
        
        # 가장 긴 motif부터 확인
        for motif in sorted_motifs:
            if text[position:position + len(motif)] == motif:
                result.append((motif, position))
                position += len(motif)
                found = True
                break
        
        if not found:
            # 매칭되지 않는 문자 처리
            result.append((text[position], position))
            position += 1
    
    return result



def analyze_consecutive_motifs(sequence):
    """연속된 motif들을 그룹화해서 (motif, 연속횟수) 형태로 반환"""
    if not sequence:
        return []
    
    consecutive_groups = []
    current_motif = sequence[0][0]  # 첫 번째 motif
    current_count = 1
    
    for i in range(1, len(sequence)):
        motif = sequence[i][0]
        if motif == current_motif:
            current_count += 1
        else:
            consecutive_groups.append((current_motif, current_count))
            current_motif = motif
            current_count = 1
    
    # 마지막 그룹 추가
    consecutive_groups.append((current_motif, current_count))
    
    return consecutive_groups



def group_rotational_motifs(motif_dict, known_motifs=None):
    """
    회전 관계에 있는 motif들을 그룹화하여 하나로 통합
    선택 기준: 1) known_motif로 시작하는 것 우선, 2) 카운트, 3) 길이, 4) 알파벳 순서
    
    Args:
        motif_dict: {motif: count} 형태의 딕셔너리
        known_motifs: 알려진 motif 리스트 (예: ['GAA', 'GAAGAAGCAGAA'])
    """
    if not motif_dict:
        return motif_dict
    
    if known_motifs is None:
        known_motifs = []
        
    grouped_result = {}
    processed_motifs = set()

    for motif, count in motif_dict.items():
        if motif in processed_motifs:
            continue
            
        # 현재 motif의 모든 회전 형태를 찾음
        rotations = rotate_string_set(motif)
        
        # 회전 그룹 내에서 가장 적합한 motif 찾기
        best_motif = motif
        best_count = count
        total_count = count
        group_members = [(motif, count)]
        
        # 회전 관계에 있는 다른 motif들 찾기
        for other_motif, other_count in motif_dict.items():
            if other_motif != motif and other_motif in rotations and other_motif not in processed_motifs:
                group_members.append((other_motif, other_count))
                total_count += other_count
                processed_motifs.add(other_motif)
        
        # 모든 회전 그룹 멤버들을 processed로 마킹
        for member_motif, _ in group_members:
            processed_motifs.add(member_motif)
            
        # 대표 motif 선택 로직
        best_motif = select_best_motif(group_members, known_motifs)
        
        # 대표 motif와 총 카운트로 저장
        grouped_result[best_motif] = total_count
        
        # 디버깅용 출력
        if len(group_members) > 1:
            member_info = [(m, len(m), c) for m, c in group_members]
            # print(f"[GROUPING] 회전 그룹: {member_info} → 선택: {best_motif}")
    
    return grouped_result


def select_best_motif(group_members, known_motifs):
    """
    회전 그룹에서 최적의 대표 motif를 선택
    
    선택 기준:
    1. known_motif로 시작하는 것 우선
    2. known_motif가 여러 개면 카운트가 더 높은 것
    3. known_motif가 없으면 카운트가 높은 것
    4. 카운트가 같으면 길이가 긴 것
    5. 길이도 같으면 알파벳 순서
    """
    if not group_members:
        return None
    
    # known_motif로 시작하는 것들 찾기
    known_starting_candidates = []
    for motif, count in group_members:
        for known_motif in known_motifs:
            if motif.startswith(known_motif):
                known_starting_candidates.append((motif, count, known_motif))
                break
    
    # 1. known_motif로 시작하는 것이 있으면 그 중에서 선택
    if known_starting_candidates:
        # 카운트 내림차순, 길이 내림차순, 알파벳 오름차순으로 정렬
        best_candidate = sorted(known_starting_candidates, 
                              key=lambda x: (-x[1], -len(x[0]), x[0]))[0]
        return best_candidate[0]
    
    # 2. known_motif로 시작하는 것이 없으면 일반 기준으로 선택
    # 카운트 내림차순, 길이 내림차순, 알파벳 오름차순으로 정렬
    best_motif = sorted(group_members, 
                       key=lambda x: (-x[1], -len(x[0]), x[0]))[0]
    return best_motif[0]