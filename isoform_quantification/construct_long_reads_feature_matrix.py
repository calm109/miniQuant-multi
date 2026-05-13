from pathlib import Path
import numpy as np
from collections import defaultdict
from construct_feature_matrix import (
    calculate_condition_number, is_multi_isoform_region, get_condition_number,
    calculate_eff_length, cal_weight_region
)
import config

def check_full_rank(isoform_region_matrix):
    if isoform_region_matrix.size == 0:
        return False
    return np.linalg.matrix_rank(isoform_region_matrix) == isoform_region_matrix.shape[1]

def get_full_length_region_set(gene_regions_dict, rname, gname):
    regions_set = set()
    isoform_region_dict = defaultdict(lambda: set())
    for region in gene_regions_dict[rname][gname]:
        for isoform in gene_regions_dict[rname][gname][region]:
            isoform_region_dict[isoform].add(region)
    for isoform in isoform_region_dict:
        max_region_exon_num = 0
        for region in isoform_region_dict[isoform]:
            region_exon_num = region.count(':')
            if max_region_exon_num < region_exon_num:
                max_region_exon_num = region_exon_num
        for region in isoform_region_dict[isoform]:
            region_exon_num = region.count(':')
            if region_exon_num == max_region_exon_num:
                regions_set.add(region)
    return regions_set

def construct_region_abundance_matrix_long_read(region_read_count_dict, region_names_indics):
    region_read_count_matrix = np.zeros((len(region_names_indics)))
    for region_name in region_read_count_dict:
        region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name]
    return region_read_count_matrix

def generate_all_feature_matrix_long_read(gene_isoforms_dict, gene_regions_dict, gene_regions_read_count, gene_regions_read_length, gene_region_len_dict, raw_isoform_exons_dict, keep_all_regions=False, global_median_lr_read_len=1000):
    gene_matrix_dict = dict()
    for chr_name in gene_regions_read_count:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_regions_read_count[chr_name]:
            num_region_with_read_counts = 0
            isoform_names = gene_isoforms_dict[chr_name][gene_name]
            full_length_region_set = get_full_length_region_set(gene_regions_dict, chr_name, gene_name)

            region_isoform_dict = {}
            for region in gene_regions_read_count[chr_name][gene_name].copy():
                if keep_all_regions:
                    if gene_regions_read_count[chr_name][gene_name][region] > 0:
                        num_region_with_read_counts += 1
                elif config.add_full_length_region == 'all':
                    if gene_regions_read_count[chr_name][gene_name][region] == 0:
                        if region not in full_length_region_set:
                            del gene_regions_read_length[chr_name][gene_name][region]
                            del gene_regions_read_count[chr_name][gene_name][region]
                            continue
                    else:
                        num_region_with_read_counts += 1
                elif config.add_full_length_region in ['nonfullrank', 'none']:
                    if gene_regions_read_count[chr_name][gene_name][region] == 0:
                        if config.add_full_length_region == 'none':
                            del gene_regions_read_length[chr_name][gene_name][region]
                            del gene_regions_read_count[chr_name][gene_name][region]
                        continue
                    else:
                        num_region_with_read_counts += 1
                region_isoform_dict[region] = gene_regions_dict[chr_name][gene_name][region]

            matrix_dict = calculate_condition_number(region_isoform_dict, isoform_names, config.normalize_lr_A)

            if not keep_all_regions and config.add_full_length_region == 'nonfullrank':
                region_isoform_dict = {}
                for region in gene_regions_read_count[chr_name][gene_name].copy():
                    if gene_regions_read_count[chr_name][gene_name][region] == 0:
                        if not check_full_rank(matrix_dict['isoform_region_matrix']):
                            if region not in full_length_region_set:
                                del gene_regions_read_length[chr_name][gene_name][region]
                                del gene_regions_read_count[chr_name][gene_name][region]
                                continue
                        else:
                            del gene_regions_read_length[chr_name][gene_name][region]
                            del gene_regions_read_count[chr_name][gene_name][region]
                            continue
                    region_isoform_dict[region] = gene_regions_dict[chr_name][gene_name][region]
                matrix_dict = calculate_condition_number(region_isoform_dict, isoform_names, config.normalize_lr_A)

            region_read_count_dict = gene_regions_read_count[chr_name][gene_name]
            region_len_dict = gene_region_len_dict[chr_name][gene_name]
            region_read_length = gene_regions_read_length[chr_name][gene_name]

            if keep_all_regions:
                # 二值 0/1 矩阵，不做有效长度加权，直接列归一化
                if config.normalize_lr_A:
                    sum_A = matrix_dict['isoform_region_matrix'].sum(axis=0)
                    sum_A[sum_A == 0] = 1
                    matrix_dict['isoform_region_matrix'] = matrix_dict['isoform_region_matrix'] / sum_A
            else:
                # Apply SR-style effective length weighting with global median LR read length
                matrix_dict['region_eff_length_dict'] = calculate_eff_length(region_len_dict, global_median_lr_read_len)
                matrix_dict['isoform_region_matrix'] = cal_weight_region(matrix_dict)
                if config.normalize_lr_A:
                    sum_A = matrix_dict['isoform_region_matrix'].sum(axis=0)
                    sum_A[sum_A == 0] = 1
                    matrix_dict['isoform_region_matrix'] = matrix_dict['isoform_region_matrix'] / sum_A

            matrix_dict['region_abund_matrix'] = construct_region_abundance_matrix_long_read(
                region_read_count_dict, matrix_dict['region_names_indics'])

            num_LRs_mapped_gene = 0
            for region in region_read_count_dict:
                num_LRs_mapped_gene += region_read_count_dict[region]
            matrix_dict['num_LRs_mapped_gene'] = num_LRs_mapped_gene
            matrix_dict['num_exons'] = {}
            for isoform in isoform_names:
                matrix_dict['num_exons'][isoform] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform]['start_pos'])
            matrix_dict['multi_isoforms_count'], matrix_dict['unique_isoforms_count'] = 0, 0
            for region in matrix_dict['region_names_indics']:
                index = matrix_dict['region_names_indics'][region]
                count = matrix_dict['region_abund_matrix'][index]
                if is_multi_isoform_region(matrix_dict, region):
                    matrix_dict['multi_isoforms_count'] += count
                else:
                    matrix_dict['unique_isoforms_count'] += count

            try:
                matrix_dict['condition_number'] = get_condition_number(matrix_dict['isoform_region_matrix'])
            except:
                matrix_dict['condition_number'] = (float('nan'), float('nan'), float('nan'), float('nan'))

            gene_matrix_dict[chr_name][gene_name] = matrix_dict

    return gene_matrix_dict
