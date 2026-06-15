# RESOLVED (2026-06-15) reference total count 중복 계산 수정:
#   기존 seq.count(motif)는 서열 공유 motif(예: TTTC ⊂ TTTCC)를 중복으로 셌음.
#   process_gene_reference_motif에서 greedy partition 기반으로 변경 (sample 경로와 일관).
from collections import defaultdict
import analyze_func_pysam
from config import *
import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import multiprocessing as mp
from functools import partial


def load_csv_data(csv_file):
    STR_regions_dict = defaultdict(dict)
    depth_dict = defaultdict(int)
    with open(csv_file, "r") as f:
        next(f)  # remove header
        for line in f:
            gene, chrom, start, end, known_motif, expansion_number = line.strip().split(",")
            STR_regions_dict[gene] = {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "known_motif": known_motif,
                "expansion_number": int(expansion_number)
            }
            depth_dict[gene] = 0
    return STR_regions_dict, depth_dict




def process_gene_reference_motif(gene_data, fasta_file):
    """
    단일 gene에 대해 reference motif 찾기를 수행하는 함수
    """
    gene, gene_info = gene_data
    patho_start = gene_info["start"]
    patho_end = gene_info["end"]
    chrom = gene_info["chrom"]
    known_motif_list = gene_info["known_motif"].split("/")

    seq = analyze_func_pysam.get_sequence_from_fasta(fasta_file, chrom, patho_start-REFERENCE_LEFT_TRIM, patho_end + REFERENCE_RIGHT_TRIM)
    consecutive_substrings = analyze_func_pysam.find_consecutive_base_motifs(seq,
                                                                                min_length=REFERENCE_MINIMUM_MOTIF_LENGTH,
                                                                                max_length=REFERENCE_MAXIMUM_MOTIF_LENGTH, 
                                                                                consecutive_threshold=REFERENCE_CONSECUTIVE_THRESHOLD)
    
    # Reference motif 처리: 원본 motif를 그대로 유지하면서 known motif와 매핑
    new_consecutive_substrings = {}
    
    # 길이 우선 정렬 (긴 motif 우선 처리)
    for motif, count in sorted(consecutive_substrings.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
        final_motif = motif  # 기본적으로 원본 motif 유지
        
        # known motif와 회전해서 일치하는 경우에만 매핑
        for known_motif in known_motif_list:
            if motif in analyze_func_pysam.rotate_string_set(known_motif):
                final_motif = known_motif
                break
        
        # 매핑된 motif로 카운트 합산
        if final_motif in new_consecutive_substrings:
            new_consecutive_substrings[final_motif] += count
        else:
            new_consecutive_substrings[final_motif] = count
    
    consecutive_substrings = new_consecutive_substrings

    # Reference motif dictionary에 저장
    reference_motif_dict_consc = {}
    reference_motif_dict_total = {}

    # 길이 우선 정렬 (긴 motif 우선 처리)
    # NOTE: 여기의 total(seq.count)은 서열 공유 motif를 중복 계산한다. 최종 보고 전
    #       motif_dict 확정 후 recompute_reference_total_no_overlap로 덮어쓴다(main_*).
    for motif, count in sorted(consecutive_substrings.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
        reference_motif_dict_consc[motif] = count
        reference_motif_dict_total[motif] = seq.count(motif)

    # reference 윈도우 서열을 함께 반환 → motif_dict 확정 후 중복 없는 total 재계산에 사용
    return gene, reference_motif_dict_consc, reference_motif_dict_total, seq


def main_parallel(bam_file, csv_file, fasta_file, num_processes=None, output_dir=None):
    """
    multiprocessing을 사용하는 메인 함수
    """
    if num_processes is None:
        num_processes = min(mp.cpu_count(), 8)  # 최대 8개 프로세스 사용

    # 결과물은 output_dir 아래에 생성 (raw data 폴더 오염 방지).
    # output_dir 미지정 시 기존 동작(BAM 폴더) 유지.
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(bam_file))
    motif_results_folder = os.path.join(output_dir, "motif_results")
    os.makedirs(motif_results_folder, exist_ok=True)
    
    # CSV 파일에서 데이터를 읽어옴
    STR_regions_dict, depth_dict = load_csv_data(csv_file)

    motif_ref_denovo_file_name = bam_file.replace(".bam", "_motif_ref_denovo.xlsx")
    motif_ref_denovo_file_name = os.path.join(motif_results_folder, os.path.basename(motif_ref_denovo_file_name))
    coverage_file = bam_file.replace(".bam", "_coverage.txt")
    coverage_file = os.path.join(motif_results_folder, os.path.basename(coverage_file))

    reference_motif_dict_consc = defaultdict(dict)
    reference_motif_dict_total = defaultdict(dict)
    reference_seq_dict = {}

    # GRCh38 reference genome에서 tandem repeat motif를 찾음. (병렬처리)
    print(f"Processing reference motifs with {num_processes} processes...")
    with mp.Pool(processes=num_processes) as pool:
        gene_data_list = list(STR_regions_dict.items())
        process_gene_func = partial(process_gene_reference_motif, fasta_file=fasta_file)
        results = pool.map(process_gene_func, gene_data_list)

    # 결과를 딕셔너리로 변환
    for gene, ref_consc, ref_total, seq in results:
        reference_motif_dict_consc[gene] = ref_consc
        reference_motif_dict_total[gene] = ref_total
        reference_seq_dict[gene] = seq

    # motif dict만들기, 여기서의 motif_dict에 reference motif도 포함되어 있음
    print("Creating motif dictionary...")
    motif_dict = analyze_func_pysam.make_motif_dict_parallel(bam_file, STR_regions_dict, depth_dict, reference_motif_dict_consc, num_processes)

    # motif_dict는 {gene: {motif: [count, ...], ...}, ...} 형태
    motif_dict = analyze_func_pysam.filter_motif_dict(motif_dict, reference_motif_dict_consc)
    # motif중에 known_motif와 회전해서 일치하는 motif가 있다면 해당 motif를 사용한다.
    motif_dict = analyze_func_pysam.apply_known_motif_v2(motif_dict, STR_regions_dict)

    # reference total을 중복 없이 재계산 (서열 공유 motif의 이중 카운트 제거). 보고되는
    # reference motif set 기준 greedy partition → sample 경로와 동일한 척도.
    reference_motif_dict_total = analyze_func_pysam.recompute_reference_total_no_overlap(
        reference_seq_dict, motif_dict, reference_motif_dict_consc)

    print("Creating pattern dictionary...")
    pattern_dict, consecutive_repeat_results, total_repeat_results = analyze_func_pysam.make_pattern_dict_parallel(bam_file, STR_regions_dict, motif_dict, num_processes)
    pattern_dict = analyze_func_pysam.simplify_pattern_dict(pattern_dict, motif_dict)

    # reference motif과 de novo motif을 비교하기 위한 파일 생성 
    print("Saving results...")
    analyze_func_pysam.save_motif_as_xlsx(motif_ref_denovo_file_name, motif_dict, reference_motif_dict_consc, reference_motif_dict_total, consecutive_repeat_results, total_repeat_results, depth_dict)

    # # Percent coverage 계산 및 저장    
    analyze_func_pysam.write_coverage_percent(coverage_file, depth_dict, threshold=COVERAGE_THRESHOLD)

    return (pattern_dict, motif_dict, consecutive_repeat_results, total_repeat_results)


def main_sequential(bam_file, csv_file, fasta_file, output_dir=None):
    # 결과물은 output_dir 아래에 생성 (raw data 폴더 오염 방지).
    # output_dir 미지정 시 기존 동작(BAM 폴더) 유지.
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(bam_file))
    motif_results_folder = os.path.join(output_dir, "motif_results")
    os.makedirs(motif_results_folder, exist_ok=True)
    motif_dict = defaultdict(lambda: defaultdict(list))
    # CSV 파일에서 데이터를 읽어옴
    STR_regions_dict, depth_dict = load_csv_data(csv_file)

    motif_ref_denovo_file_name = bam_file.replace(".bam", "_motif_ref_denovo.xlsx")
    motif_ref_denovo_file_name = os.path.join(motif_results_folder, os.path.basename(motif_ref_denovo_file_name))
    coverage_file = bam_file.replace(".bam", "_coverage.txt")
    coverage_file = os.path.join(motif_results_folder, os.path.basename(coverage_file))

    reference_motif_dict_consc = defaultdict(dict)
    reference_motif_dict_total = defaultdict(dict)
    reference_seq_dict = {}

    # GRCh38 reference genome에서 tandem repeat motif를 찾음.
    for gene in STR_regions_dict.keys():
        patho_start = STR_regions_dict[gene]["start"]
        patho_end = STR_regions_dict[gene]["end"]
        chrom = STR_regions_dict[gene]["chrom"]

        seq = analyze_func_pysam.get_sequence_from_fasta(fasta_file, chrom, patho_start-REFERENCE_LEFT_TRIM, patho_end + REFERENCE_RIGHT_TRIM)
        reference_seq_dict[gene] = seq
        # consecutive_substrings = analyze_func.find_consecutive_repeated_substrings_with_rotation_optimized_v2(seq, min_length=3, max_length=30, consecutive_threshold=3)
        consecutive_substrings = analyze_func_pysam.find_consecutive_base_motifs(seq,
                                                                                    min_length=REFERENCE_MINIMUM_MOTIF_LENGTH,
                                                                                    max_length=REFERENCE_MAXIMUM_MOTIF_LENGTH, 
                                                                                    consecutive_threshold=REFERENCE_CONSECUTIVE_THRESHOLD)
        

        # Reference motif 처리: 원본 motif를 그대로 유지하면서 known motif와 매핑
        new_consecutive_substrings = {}
        known_motif_list = STR_regions_dict[gene]["known_motif"].split("/")
        
        # 길이 우선 정렬 (긴 motif 우선 처리)
        for motif, count in sorted(consecutive_substrings.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
            final_motif = motif  # 기본적으로 원본 motif 유지
            
            # known motif와 회전해서 일치하는 경우에만 매핑
            for known_motif in known_motif_list:
                if motif in analyze_func_pysam.rotate_string_set(known_motif):
                    final_motif = known_motif
                    break
            
            # 매핑된 motif로 카운트 합산
            if final_motif in new_consecutive_substrings:
                new_consecutive_substrings[final_motif] += count
            else:
                new_consecutive_substrings[final_motif] = count
        consecutive_substrings = new_consecutive_substrings

        # Reference motif dictionary에 저장
        # 길이 우선 정렬 (긴 motif 우선 처리)
        for motif, count in sorted(consecutive_substrings.items(), key=lambda x: (-len(x[0]), -x[1], x[0])):
            reference_motif_dict_consc[gene][motif] = count
            reference_motif_dict_total[gene][motif] = seq.count(motif)
        
    # motif dict만들기, 여기서의 motif_dict에 reference motif도 포함되어 있음
    motif_dict = analyze_func_pysam.make_motif_dict(bam_file, STR_regions_dict, depth_dict, reference_motif_dict_consc)
    # motif_dict는 {gene: {motif: [count, ...], ...}, ...} 형태
    motif_dict = analyze_func_pysam.filter_motif_dict(motif_dict, reference_motif_dict_consc)
    # motif중에 known_motif와 회전해서 일치하는 motif가 있다면 해당 motif를 사용한다.
    motif_dict = analyze_func_pysam.apply_known_motif_v2(motif_dict, STR_regions_dict)

    # reference total을 중복 없이 재계산 (서열 공유 motif의 이중 카운트 제거).
    reference_motif_dict_total = analyze_func_pysam.recompute_reference_total_no_overlap(
        reference_seq_dict, motif_dict, reference_motif_dict_consc)

    # custom motif 처리

    pattern_dict, consecutive_repeat_results, total_repeat_results = analyze_func_pysam.make_pattern_dict(bam_file, STR_regions_dict, motif_dict)
    pattern_dict = analyze_func_pysam.simplify_pattern_dict(pattern_dict, motif_dict)


    # reference motif과 de novo motif을 비교하기 위한 파일 생성 
    analyze_func_pysam.save_motif_as_xlsx(motif_ref_denovo_file_name, motif_dict, reference_motif_dict_consc, reference_motif_dict_total, consecutive_repeat_results, total_repeat_results, depth_dict)

    # # Percent coverage 계산 및 저장    
    analyze_func_pysam.write_coverage_percent(coverage_file, depth_dict, threshold=COVERAGE_THRESHOLD)

    return (pattern_dict, motif_dict, consecutive_repeat_results, total_repeat_results)

def nice_tick_step(n, target_tick_count=7):
    """
    예쁜 숫자의 step 값을 반환함.
    전체 길이 n과 원하는 tick 갯수를 기준으로 적당한 step을 결정.
    """
    raw_step = n / target_tick_count
    # 1, 2, 5, 10, 20, 50, 100... 같은 기본 단위
    base_steps = [1, 2, 5]
    
    # 10의 거듭제곱 단위를 기반으로 예쁜 step 계산
    magnitude = 10 ** int(np.floor(np.log10(raw_step)))
    for base in base_steps:
        step = base * magnitude
        if step >= raw_step:
            return step
    return 10 * magnitude  # fallback

def plot_pattern(ax, pattern_dict, motif_dict, gene, read_threshold):
    from matplotlib.patches import Patch

    # pattern_dict에 gene이 없거나 read 수가 부족한 경우
    if gene not in pattern_dict or len(pattern_dict[gene]) < read_threshold:
        ax.text(0.5, 0.5, f"Not enough reads for {gene}", fontsize=10, ha='center', va='center')
        ax.axis('off')
        return

    # motif가 없는 경우 처리
    if gene not in motif_dict or len(motif_dict[gene]) == 0:
        # motif 없이 회색으로만 표시
        all_patterns = []
        pattern_colors = {}
    else:
        # 재현성을 위해 길이와 알파벳 순으로 정렬
        all_patterns = list(sorted(motif_dict[gene].keys(), key=lambda x: (len(x), x)))

        # colormap = "Paired"
        colormap = "Accent"
        if len(all_patterns) > 12:
            colormap = "tab20"
        elif len(all_patterns) > 18:
            raise ValueError("Too many patterns to display. Please reduce the number of patterns.")

        if colormap == "Accent":
            pattern_colors = dict(zip(all_patterns, sns.color_palette(colormap, len(all_patterns))))
        else:
            tab20 = sns.color_palette("tab20", 20)
            filtered_tab20 = [color for i, color in enumerate(tab20) if i not in [14,15]]
            pattern_colors = dict(zip(all_patterns, filtered_tab20))

    gray = sns.color_palette("gray", 1)[0]

    rows = []
    for read in pattern_dict[gene]:
        row_colors = []
        for pattern, count in read:
            row_colors.extend([pattern_colors.get(pattern, gray)] * len(pattern)  * count)
        rows.append(row_colors)

    max_len = max(len(row) for row in rows)
    for i in range(len(rows)):
        if len(rows[i]) < max_len:
            rows[i].extend([(1,1,1)] * (max_len - len(rows[i])))

    rgb_array = np.array(rows)
    num_cols = rgb_array.shape[1]
    num_reads = rgb_array.shape[0]

    xtick_step = nice_tick_step(num_cols)
    ytick_step = nice_tick_step(num_reads)

    ax.imshow(rgb_array, aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(0, num_cols, xtick_step))
    ax.set_xticklabels(np.arange(0, num_cols, xtick_step), fontsize=8)
    ax.set_yticks(np.arange(0, num_reads, ytick_step))
    ax.set_yticklabels(np.arange(0, num_reads, ytick_step), fontsize=8)
    ax.set_xlabel("Sequence Length", fontsize=7)
    ax.set_ylabel("Read Count", fontsize=7)
    ax.set_title(f"{gene} - Pattern Visualization", fontsize=7)

    # Legend 추가 (motif가 있을 때만)
    if len(pattern_colors) > 0:
        legend_elements = [
            Patch(facecolor=color, label=pattern) for pattern, color in pattern_colors.items()
        ]
        ax.legend(handles=legend_elements, title="Pattern",bbox_to_anchor=(1.0,1.0), loc='upper left', ncol=1, frameon=True,
                  fontsize=4, title_fontsize=5)

    


def save_gene_plots_with_heatmap_v2(pattern_dict, motif_dict, bam_file, STR_regions_dict, input_file, read_threshold=5, output_dir=None):
    """
    repeat number histogram을 추가하여 4행 3열로 구성된 PDF 파일을 생성합니다.

    PDF는 output_dir/gene_panel_output/ 아래에 저장된다. output_dir 미지정 시
    기존 동작(input_file 폴더)을 유지한다. 파일명은 input_file의 basename을 사용.
    """
    # motif가 없는 gene도 포함하기 위해 STR_regions_dict 사용
    gene_list = list(STR_regions_dict.keys())

    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(input_file))
    output_folder = os.path.join(output_dir, "gene_panel_output")
    os.makedirs(output_folder, exist_ok=True)

    rows, cols = 4, 2
    plots_per_page = rows  # gene per page
    figsize = (11.69, 11.69)  # A4 landscape
    filename = os.path.join(output_folder, f"{os.path.basename(input_file)}_gene_panel_output.pdf")

    with PdfPages(filename) as pdf:
        for page_start in range(0, len(gene_list), plots_per_page):
            fig, axes = plt.subplots(rows, cols, figsize=figsize, constrained_layout=True)
            axes = axes.reshape(rows, cols)

            for row in range(rows):
                gene_index = page_start + row
                if gene_index >= len(gene_list):
                    for col in range(cols):
                        axes[row, col].axis('off')
                    continue

                gene = gene_list[gene_index]

                # 왼쪽: heatmap
                plot_pattern(axes[row, 0], pattern_dict, motif_dict, gene, read_threshold=read_threshold)

                # 가운데: KDE plot
                analyze_func_pysam.show_kde_v2(axes[row, 1], bam_file, STR_regions_dict, gene, read_threshold=read_threshold)

                # 오른쪽: Repeat Number Histogram
                # analyze_func_pysam.plot_repeat_number_distribution(axes[row, 2], bam_file, STR_regions_dict, gene, read_threshold=read_threshold)

            pdf.savefig(fig)
            plt.close(fig)



VERSION = "v1.2.0"


def parse_args(argv=None):
    """Parse STRiker command-line arguments.

    Returns:
        argparse.Namespace with: csv_file, fasta_file, bam_file (paths),
        output (output directory), process (int or None for sequential mode).
    """
    class _Formatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
        # show argument defaults, but keep the epilog examples verbatim
        pass

    parser = argparse.ArgumentParser(
        prog="STRiker.py",
        description=(
            "STRiker - Short Tandem Repeat analyzer for nanopore BAMs.\n"
            "Finds reference + de novo motifs and calls per-allele repeat length\n"
            "from read-length KDE; writes a motif xlsx, a coverage report, and a\n"
            "per-gene panel PDF."
        ),
        formatter_class=_Formatter,
        epilog="Examples:\n"
               "  # sequential, results into ./striker_out\n"
               "  python STRiker.py regions.csv ref.fa input.bam -o striker_out\n\n"
               "  # 8 processes\n"
               "  python STRiker.py regions.csv ref.fa input.bam -o striker_out -p 8",
    )
    parser.add_argument("csv_file",
                        help="CSV of STR regions (chr,start,end,motif).")
    parser.add_argument("fasta_file",
                        help="Reference genome FASTA (indexed; uppercase recommended).")
    parser.add_argument("bam_file",
                        help="Input BAM aligned to the reference.")
    parser.add_argument("-o", "--output", default=".",
                        help="Output directory for all results "
                             "(motif_results/ and gene_panel_output/ are created "
                             "inside). Keeps results out of the raw-data folder.")
    parser.add_argument("-p", "--process", type=int, default=None, metavar="N",
                        help="Number of processes for multiprocessing. "
                             "If omitted, runs sequentially.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"STRiker {VERSION}")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    csv_file, fasta_file, bam_file = args.csv_file, args.fasta_file, args.bam_file
    output_dir = args.output
    num_processes = args.process

    # File existence checks
    for label, path in (("BAM", bam_file), ("FASTA", fasta_file), ("CSV", csv_file)):
        if not os.path.exists(path):
            print(f"Error: {label} file {path} does not exist.")
            sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Run analysis
    if num_processes is not None:
        print(f"Running analysis with multiprocessing ({num_processes} processes)...")
        pattern_dict, motif_dict, consecutive_repeat_results, total_repeat_results = main_parallel(
            bam_file, csv_file, fasta_file, num_processes, output_dir=output_dir)
    else:
        print("Running analysis sequentially...")
        pattern_dict, motif_dict, consecutive_repeat_results, total_repeat_results = main_sequential(
            bam_file, csv_file, fasta_file, output_dir=output_dir)

    STR_regions_dict, _ = load_csv_data(csv_file)

    print("Generating plots...")
    save_gene_plots_with_heatmap_v2(
        pattern_dict,
        motif_dict,
        bam_file=bam_file,
        STR_regions_dict=STR_regions_dict,
        input_file=bam_file,
        output_dir=output_dir
    )
    print(f"Analysis completed! Results in: {os.path.abspath(output_dir)}")