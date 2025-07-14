from collections import defaultdict
import analyze_func_pysam
from config import *
import sys
import os 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


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




def main_sequential(bam_file, csv_file, fasta_file):
    # motif_results라는 폴더 없으면 생성
    motif_results_folder = os.path.join(os.path.dirname(bam_file), "motif_results")
    os.makedirs(motif_results_folder, exist_ok=True)
    motif_dict = defaultdict(lambda: defaultdict(list))
    # CSV 파일에서 데이터를 읽어옴
    STR_regions_dict, depth_dict = load_csv_data(csv_file)

    motif_ref_denovo_file_name = bam_file.replace(".bam", "_motif_ref_denovo.xlsx")
    motif_ref_denovo_file_name = os.path.join(motif_results_folder, os.path.basename(motif_ref_denovo_file_name))
    coverage_file = bam_file.replace(".bam", "_coverage.txt")
    coverage_file = os.path.join(motif_results_folder, os.path.basename(coverage_file))

    reference_motif_dict_consc = defaultdict(lambda: defaultdict(int))
    reference_motif_dict_total = defaultdict(lambda: defaultdict(int))
    
    # GRCh38 reference genome에서 tandem repeat motif를 찾음.
    for gene in STR_regions_dict.keys():
        patho_start = STR_regions_dict[gene]["start"]
        patho_end = STR_regions_dict[gene]["end"]
        chrom = STR_regions_dict[gene]["chrom"]

        seq = analyze_func_pysam.get_sequence_from_fasta(fasta_file, chrom, patho_start-REFERENCE_LEFT_TRIM, patho_end + REFERENCE_RIGHT_TRIM)
        # consecutive_substrings = analyze_func.find_consecutive_repeated_substrings_with_rotation_optimized_v2(seq, min_length=3, max_length=30, consecutive_threshold=3)
        consecutive_substrings = analyze_func_pysam.find_consecutive_base_motifs(seq,
                                                                                    min_length=REFERENCE_MINIMUM_MOTIF_LENGTH,
                                                                                    max_length=REFERENCE_MAXIMUM_MOTIF_LENGTH, 
                                                                                    consecutive_threshold=REFERENCE_CONSECUTIVE_THRESHOLD)
        
        # 회전 그룹핑은 나중에 apply_known_motif_v2에서 통합 처리하므로 여기서는 제거

        # Reference motif 처리: 원본 motif를 그대로 유지하면서 known motif와 매핑
        new_consecutive_substrings = {}
        known_motif_list = STR_regions_dict[gene]["known_motif"].split("/")
        
        for motif, count in consecutive_substrings.items():
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
        for motif, count in consecutive_substrings.items():
            reference_motif_dict_consc[gene][motif] = count
            reference_motif_dict_total[gene][motif] = seq.count(motif)
        

    # motif dict만들기, 여기서의 motif_dict에 reference motif도 포함되어 있음
    motif_dict = analyze_func_pysam.make_motif_dict(bam_file, STR_regions_dict, depth_dict, reference_motif_dict_consc)

    # coverage에 못미치는 motif들 제거, reference motif는 무조건 살리기, 그리고 각 motif별 최댓값을 저장함.
    # motif_dict는 {gene: {motif: [count, ...], ...}, ...} 형태
    motif_dict = analyze_func_pysam.filter_motif_dict(motif_dict, MINIMUM_MOTIF_COVERAGE, reference_motif_dict_consc)
    # motif중에 known_motif와 회전해서 일치하는 motif가 있다면 해당 motif를 사용한다.
    motif_dict = analyze_func_pysam.apply_known_motif_v2(motif_dict, STR_regions_dict)
    
    pattern_dict, consecutive_repeat_results, total_repeat_results = analyze_func_pysam.process_linesV2(bam_file, STR_regions_dict, motif_dict)
    pattern_dict = analyze_func_pysam.simplify_pattern_dict(pattern_dict, motif_dict)



    # reference motif과 de novo motif을 비교하기 위한 파일 생성 
    analyze_func_pysam.save_motif_as_xlsx(motif_ref_denovo_file_name, motif_dict, reference_motif_dict_consc, reference_motif_dict_total,consecutive_repeat_results,total_repeat_results, depth_dict)

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

    if len(pattern_dict[gene]) < read_threshold:
        ax.text(0.5, 0.5, f"Not enough reads for {gene}", fontsize=10, ha='center', va='center')
        ax.axis('off')
        return

    all_patterns = list(sorted(motif_dict[gene].keys(), key=lambda x: (len(x), x)))
    
    # 디버깅: motif_dict에 있는 motif들과 실제 pattern에서 사용되는 motif들 비교
    pattern_motifs = set()
    for read in pattern_dict[gene]:
        for pattern, count in read:
            pattern_motifs.add(pattern)
    
    print(f"[DEBUG {gene}] motif_dict에 있는 motifs: {sorted(motif_dict[gene].keys(), key=len)}")
    print(f"[DEBUG {gene}] pattern에서 사용된 motifs: {sorted(pattern_motifs, key=len)}")
    
    missing_in_pattern = set(motif_dict[gene].keys()) - pattern_motifs
    only_in_pattern = pattern_motifs - set(motif_dict[gene].keys())
    
    if missing_in_pattern:
        print(f"[DEBUG {gene}] ⚠️  motif_dict에만 있는 motifs (plot에서 누락): {sorted(missing_in_pattern, key=len)}")
    if only_in_pattern:
        print(f"[DEBUG {gene}] ⚠️  pattern에만 있는 motifs: {sorted(only_in_pattern, key=len)}")
    if not missing_in_pattern and not only_in_pattern:
        print(f"[DEBUG {gene}] ✅ motif_dict와 pattern이 완전히 일치!")
    colormap = "Paired"
    if len(all_patterns) > 12:
        colormap = "tab20"
    elif len(all_patterns) > 18:
        raise ValueError("Too many patterns to display. Please reduce the number of patterns.")

    if colormap == "Paired":
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
    ax.set_xlabel("Sequence Length", fontsize=4)
    ax.set_ylabel("Read Index", fontsize=4)
    ax.set_title(f"{gene} - Pattern Visualization", fontsize=7)

    # Legend 추가
    legend_elements = [
        Patch(facecolor=color, label=pattern) for pattern, color in pattern_colors.items()
    ]
    ax.legend(handles=legend_elements, title="Pattern",bbox_to_anchor=(1.0,1.0), loc='upper left', ncol=1, frameon=True,
              fontsize=4, title_fontsize=5)

    # 범례는 따로 처리 (원하면 첫 페이지만 넣을 수도 있음)
def save_gene_plots_with_heatmap(pattern_dict, motif_dict, bam_file, STR_regions_dict, input_file, read_threshold=5):
    gene_list = list(motif_dict.keys())
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    output_folder = os.path.join(os.path.dirname(input_file), "gene_panel_output")
    os.makedirs(output_folder, exist_ok=True)
    plots_per_page = 4
    rows, cols = 4, 2
    figsize = (8.27, 11.69)
    filename = os.path.join(output_folder, f"{os.path.basename(input_file)}_gene_panel_output.pdf")
    with PdfPages(filename) as pdf:
        for page_start in range(0, len(gene_list), plots_per_page):
            fig, axes = plt.subplots(rows, cols, figsize=figsize, constrained_layout=True)
            axes = axes.flatten()

            for row in range(plots_per_page):
                gene_index = page_start + row
                if gene_index >= len(gene_list):
                    axes[row*2].axis('off')
                    axes[row*2+1].axis('off')
                    continue

                gene = gene_list[gene_index]

                # 왼쪽: heatmap
                plot_pattern(axes[row*2], pattern_dict, motif_dict, gene, read_threshold=read_threshold)
                # 오른쪽: read length KDE
                analyze_func_pysam.show_kde_v2(axes[row*2+1], bam_file, STR_regions_dict, gene, read_threshold=read_threshold)

            pdf.savefig(fig)
            plt.close(fig)



def save_gene_plots_with_heatmap_v2(pattern_dict, motif_dict, bam_file, STR_regions_dict, input_file, read_threshold=5):
    """
    repeat number histogram을 추가하여 4행 3열로 구성된 PDF 파일을 생성합니다.
    """
    from matplotlib.lines import Line2D
    gene_list = list(motif_dict.keys())

    output_folder = os.path.join(os.path.dirname(input_file), "gene_panel_output")
    os.makedirs(output_folder, exist_ok=True)

    rows, cols = 4, 3
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
                analyze_func_pysam.plot_repeat_number_distribution(axes[row, 2], bam_file, STR_regions_dict, gene, read_threshold=read_threshold)

            pdf.savefig(fig)
            plt.close(fig)



if __name__ == "__main__":

    # check arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <csv_file> <fasta_file> <bam_file> ")
        sys.exit(1)
    csv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    bam_file = sys.argv[3]
    if not os.path.exists(bam_file):
        print(f"Error: BAM file {bam_file} does not exist.")
        sys.exit(1)
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file {fasta_file} does not exist.")
        sys.exit(1)
    if not os.path.exists(csv_file):
        print(f"Error: CSV file {csv_file} does not exist.")
        sys.exit(1)



    pattern_dict, motif_dict, consecutive_repeat_results, total_repeat_results= main_sequential(bam_file, csv_file, fasta_file=fasta_file)
    STR_regions_dict, depth_dict = load_csv_data(csv_file)
    save_gene_plots_with_heatmap_v2(
        pattern_dict,
        motif_dict,
        bam_file=bam_file,
        STR_regions_dict=STR_regions_dict,
        input_file=bam_file
    )
    # print("Sample name:", os.path.basename(bam_file))
    # print("Gene panel output saved to:", os.path.join(os.path.dirname(bam_file), "gene_panel_output"))