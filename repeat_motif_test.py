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

def analyze_motif_usage(sequence):
    """전체 Motif 사용 통계 분석"""
    from collections import Counter
    
    # 각 motif의 총 사용 횟수
    motif_counts = Counter([motif for motif, _ in sequence])
    
    return motif_counts

def visualize_motif_sequence(text, sequence):
    """Motif 분해 결과를 시각적으로 표시"""
    print("=== 원본 문자열 ===")
    print(text)
    print(f"길이: {len(text)}")
    
    print("\n=== 분해 결과 ===")
    visual = ""
    position_info = []
    
    for i, (motif, start_pos) in enumerate(sequence):
        if len(motif) == 1:  # 단일 문자는 괄호로 표시
            visual += f"({motif})"
            position_info.append(f"{i+1}. '({motif})' - 위치 {start_pos}")
        else:  # 알려진 motif는 대괄호로 표시
            visual += f"[{motif}]"
            position_info.append(f"{i+1}. '[{motif}]' - 위치 {start_pos}")
    
    print(visual)
    print("\n=== 상세 분해 정보 ===")
    for info in position_info:
        print(info)

def find_unknown_segments(sequence, known_motifs):
    """알려진 motif가 아닌 부분들을 찾기"""
    unknown_segments = []
    
    for motif, pos in sequence:
        if motif not in known_motifs:
            unknown_segments.append((motif, pos))
    
    return unknown_segments

# 분석 실행
text = "GAAGAAGCAGAAGAAGAAGCAGAAGAAGAAGAAGCAGAACCGAAGGAAGAAGCAGAA"
known_motifs = ["GAA", "GAAGAAGCAGAA"]

print("=== MOTIF 분석 결과 ===")
print("Sequence")
sequence = find_motif_sequence(text, known_motifs)
print(sequence)
# visualize_motif_sequence(text, sequence)

# 연속 그룹 분석
consecutive_groups = analyze_consecutive_motifs(sequence)
print(consecutive_groups)
# print(f"\n=== 연속 그룹 분석 ===")
# for i, (motif, count) in enumerate(consecutive_groups, 1):
#     if len(motif) == 1:
#         print(f"{i}. ({motif}, {count})")
#     else:
#         print(f"{i}. ({motif}, {count})")

# # 전체 사용 빈도
# total_counts = analyze_motif_usage(sequence)
# print(f"\n=== 전체 사용 빈도 ===")
# for motif, count in total_counts.items():
#     if len(motif) == 1:
#         print(f"'{motif}' (단일 문자): {count}번")
#     else:
#         print(f"'{motif}': {count}번")

# print(f"\n=== 알려진 MOTIF 이외의 구성 요소 ===")
# unknown_segments = find_unknown_segments(sequence, known_motifs)
# if unknown_segments:
#     unique_unknowns = list(set([seg[0] for seg in unknown_segments]))
#     for segment in unique_unknowns:
#         positions = [pos for seg, pos in unknown_segments if seg == segment]
#         print(f"'{segment}' - 위치: {positions}")
# else:
#     print("모든 부분이 알려진 motif로 구성됨")

# print(f"\n=== 전체 구성 비율 ===")
# total_length = len(text)
# for motif, count in total_counts.items():
#     coverage = (len(motif) * count / total_length) * 100
#     print(f"'{motif}': {coverage:.1f}% ({len(motif) * count}/{total_length} 문자)")