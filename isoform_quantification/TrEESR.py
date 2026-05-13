from collections import defaultdict
import numpy as np
from TransELS import infer_read_len
from construct_feature_matrix import calculate_all_condition_number,calculate_all_condition_number_long_read
from parse_annotation_main import parse_reference_annotation
from generate_output import generate_TrEESR_output
from construct_feature_matrix import generate_all_feature_matrix_short_read
from construct_long_reads_feature_matrix import generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment,filter_regions_read_length
from parse_alignment_main import parse_alignment,is_transcriptome_aligned_sam,build_transcript_genome_map,parse_transcriptome_sr_alignment
from identifiability import compute_and_output_identifiability
import config
import datetime
from pathlib import Path
import shutil
def TrEESR(ref_file_path,output_path,short_read_alignment_file_path,long_read_alignment_file_path,sr_region_selection,filtering,threads,lr_region_selection='real_data',READ_LEN=150,READ_JUNC_MIN_MAP_LEN=0):
    # long_read_alignment_file_path 可以是单个路径字符串或路径列表
    if isinstance(long_read_alignment_file_path, list):
        lr_sam_list = [p for p in long_read_alignment_file_path if p is not None]
    elif long_read_alignment_file_path is not None:
        lr_sam_list = [long_read_alignment_file_path]
    else:
        lr_sam_list = []
    # 只清理 cal_K_value 自身会重建的子目录，保留 quantify 的输出文件
    Path(output_path).mkdir(parents=True, exist_ok=True)
    temp_dir = Path(f'{output_path}/temp')
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    Path(f'{output_path}/temp/machine_learning/').mkdir(parents=True, exist_ok=True)
    # 推断 SR 读长：只读取第一条读段长度，不做完整解析
    if short_read_alignment_file_path is not None:
        _sr_path_for_len = short_read_alignment_file_path[0] if isinstance(short_read_alignment_file_path, list) else short_read_alignment_file_path
        READ_LEN = infer_read_len(_sr_path_for_len)
        print(f'[INFO] SR read length inferred from SAM: {READ_LEN}', flush=True)
    start_time = datetime.datetime.now()
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict,same_structure_isoform_dict,removed_gene_isoform_dict = \
        parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,sr_region_selection)
    end_time_1 = datetime.datetime.now()
    print('[INFO] Done in %.3f s'%((end_time_1-start_time).total_seconds()),flush=True)
    print('[INFO] Calculating the condition number...',flush=True)
    gene_regions_points_list,gene_range,gene_interval_tree_dict = process_annotation_for_alignment(gene_exons_dict,gene_points_dict)
    # SR A 矩阵：对每个 SR SAM 各推断读长、各建一个矩阵
    if isinstance(short_read_alignment_file_path, list):
        sr_sam_list = [p for p in short_read_alignment_file_path if p is not None]
    elif short_read_alignment_file_path is not None:
        sr_sam_list = [short_read_alignment_file_path]
    else:
        sr_sam_list = []
    sr_theoretical_matrix_dict = None  # 理论 A 矩阵（与 README 一致），仅在 real_data 模式下额外计算
    if sr_sam_list:
        sr_matrix_list = []
        for sr_sam in sr_sam_list:
            sr_read_len = infer_read_len(sr_sam)
            print(f'[INFO] SR read length inferred from {sr_sam}: {sr_read_len}', flush=True)
            if sr_region_selection == 'real_data':
                sr_read_count, sr_read_len, num_SRs = parse_alignment(sr_sam, READ_JUNC_MIN_MAP_LEN, gene_points_dict,
                    gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, False, False, threads)
                config.READ_LEN = sr_read_len
                sr_matrix = generate_all_feature_matrix_short_read(gene_isoforms_dict, SR_gene_regions_dict, sr_read_count, sr_read_len, SR_genes_regions_len_dict, num_SRs, False)
                # 额外计算理论 A 矩阵（与 read_length 模式相同的区域过滤）用于可识别性分析
                # 注意：eff_length 计算需要 inner region 长度，因此传入完整的 SR_genes_regions_len_dict
                if sr_theoretical_matrix_dict is None:
                    theo_SR_regions_dict, _ = filter_regions_read_length(
                        SR_gene_regions_dict, gene_points_dict, SR_genes_regions_len_dict,
                        READ_JUNC_MIN_MAP_LEN, sr_read_len, sr_read_len)
                    sr_theoretical_matrix_dict = calculate_all_condition_number(gene_isoforms_dict, theo_SR_regions_dict, SR_genes_regions_len_dict, sr_read_len, allow_multi_exons=False)
            else:
                sr_matrix = calculate_all_condition_number(gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, sr_read_len, allow_multi_exons=False)
            sr_matrix_list.append(sr_matrix)
        short_read_gene_matrix_dict = sr_matrix_list if len(sr_matrix_list) > 1 else sr_matrix_list[0]
    else:
        SR_read_len = READ_LEN
        short_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,SR_read_len,allow_multi_exons=False)
    # LR A 矩阵：每个 LR SAM 各自构建独立矩阵
    if lr_sam_list:
        lr_matrix_list = []
        for k, lr_sam in enumerate(lr_sam_list):
            rc, rl, total_len, num_LRs, _, _ = parse_alignment(lr_sam, READ_JUNC_MIN_MAP_LEN, gene_points_dict, gene_range, gene_interval_tree_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, True, filtering, threads)
            filtered_LR_regions_dict = None
            # 始终从过滤前的 rl 计算全局中位读长，保证与区域过滤所用阈值一致
            all_lr_read_lengths = []
            for chr_name in rl:
                for gene_name in rl[chr_name]:
                    for region in rl[chr_name][gene_name]:
                        all_lr_read_lengths.extend(rl[chr_name][gene_name][region])
            global_median_lr_read_len = int(np.median(all_lr_read_lengths)) if all_lr_read_lengths else 1000
            if lr_region_selection == 'read_length':
                print(f'[INFO] LR SAM {k} global median read length for region filtering: {global_median_lr_read_len} bp', flush=True)
                filtered_LR_regions_dict, _ = filter_regions_read_length(
                    LR_gene_regions_dict, gene_points_dict, LR_genes_regions_len_dict,
                    READ_JUNC_MIN_MAP_LEN, READ_JUNC_MIN_MAP_LEN, global_median_lr_read_len)
            lr_regions_to_use = filtered_LR_regions_dict if filtered_LR_regions_dict is not None else LR_gene_regions_dict
            if lr_region_selection == 'read_length' and filtered_LR_regions_dict is not None:
                for chr_name in rc:
                    for gene_name in rc[chr_name]:
                        gene_filtered = filtered_LR_regions_dict.get(chr_name, {}).get(gene_name, {})
                        for region in list(rc[chr_name][gene_name].keys()):
                            if region not in gene_filtered:
                                del rc[chr_name][gene_name][region]
                        if chr_name in rl and gene_name in rl[chr_name]:
                            for region in list(rl[chr_name][gene_name].keys()):
                                if region not in gene_filtered:
                                    del rl[chr_name][gene_name][region]
            lr_matrix = generate_all_feature_matrix_long_read(gene_isoforms_dict, lr_regions_to_use, rc, rl, LR_genes_regions_len_dict, raw_isoform_exons_dict, keep_all_regions=(lr_region_selection == 'read_length'), global_median_lr_read_len=global_median_lr_read_len)
            lr_matrix_list.append(lr_matrix)
        long_read_gene_matrix_dict = lr_matrix_list if len(lr_matrix_list) > 1 else lr_matrix_list[0]
    else:
        long_read_gene_matrix_dict = calculate_all_condition_number_long_read(gene_isoforms_dict, LR_gene_regions_dict, allow_multi_exons=True)
    end_time_2 = datetime.datetime.now()
    print('[INFO] Done in %.3f s'%((end_time_2-end_time_1).total_seconds()),flush=True)
    raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict = defaultdict(dict),defaultdict(dict),defaultdict(dict)
    raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            raw_gene_num_exon_dict[chr_name][gene_name] = len(raw_gene_exons_dict[chr_name][gene_name])
            gene_num_exon_dict[chr_name][gene_name] = len(gene_exons_dict[chr_name][gene_name])
            gene_num_isoform_dict[chr_name][gene_name] = len(gene_isoforms_dict[chr_name][gene_name])
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name]:
                raw_isoform_num_exon_dict[isoform_name] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                isoform_length_dict[isoform_name] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                num_isoforms_dict[isoform_name] =  len(raw_isoform_exons_dict[chr_name][gene_name])
    info_dict_list = [raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict,raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict]
    # SR 和 LR 均使用注释结构 + 读长过滤的理论矩阵，不依赖实际读段信息
    sr_output_matrix_dict = sr_theoretical_matrix_dict if sr_theoretical_matrix_dict is not None else short_read_gene_matrix_dict
    gene_feature_dict = generate_TrEESR_output(output_path,sr_output_matrix_dict,long_read_gene_matrix_dict,info_dict_list,same_structure_isoform_dict,removed_gene_isoform_dict,gene_points_dict,has_lr_data=bool(lr_sam_list),has_sr_data=bool(sr_sam_list))
    # For identifiability: if SR SAM is transcriptome-aligned, recompute SR b vector via coordinate conversion
    sr_matrix_for_identifiability = short_read_gene_matrix_dict if sr_sam_list else None
    if sr_sam_list and sr_region_selection == 'real_data' and is_transcriptome_aligned_sam(sr_sam_list[0]):
        print('[INFO] Transcriptome-aligned SR SAM detected. Converting read positions for identifiability analysis...', flush=True)
        tx_map = build_transcript_genome_map(ref_file_path)
        tx_sr_matrices = []
        for tx_sam in sr_sam_list:
            tx_sr_read_count, tx_sr_read_len, tx_num_SRs = parse_transcriptome_sr_alignment(
                tx_sam, tx_map, READ_JUNC_MIN_MAP_LEN,
                gene_points_dict, gene_range, gene_interval_tree_dict,
                SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict)
            tx_sr_matrices.append(generate_all_feature_matrix_short_read(
                gene_isoforms_dict, SR_gene_regions_dict, tx_sr_read_count,
                tx_sr_read_len, SR_genes_regions_len_dict, tx_num_SRs, False))
        sr_matrix_for_identifiability = tx_sr_matrices if len(tx_sr_matrices) > 1 else tx_sr_matrices[0]
    lr_matrix_for_identifiability = long_read_gene_matrix_dict if lr_sam_list else None
    compute_and_output_identifiability(output_path, sr_matrix_for_identifiability, lr_matrix_for_identifiability, sr_theoretical_matrix_dict)
    ref_annotation_dict_list = [gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict]
    return gene_feature_dict
def get_kvalues_dict(ref_file_path,threads,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=10):
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=True)
    return long_read_gene_matrix_dict


