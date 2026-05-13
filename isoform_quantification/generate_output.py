from pathlib import Path
import numpy as np
import dill as pickle
import io
from util import output_matrix_info
from construct_feature_matrix import get_condition_number
import config
import scipy.stats

def get_stats(arr):
    if len(arr) == 0:
        return np.array([float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')])
    s = scipy.stats.describe(arr)
    return np.array([s.minmax[0],s.minmax[1],s.mean,s.variance,s.skewness,s.kurtosis])
def generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list,same_structure_isoform_dict,removed_gene_isoform_dict,gene_points_dict,has_lr_data=True,has_sr_data=True,sr_actual_matrix_for_info=None):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    [raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict,raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict] = info_dict_list
    # out_dict = short_read_gene_matrix_dict.copy()
    # bio = io.BytesIO()
    # for chr in out_dict:
    #     for gene in out_dict[chr]:
    #         bio.write(str.encode('{}\n'.format(gene)))
    #         np.savetxt(bio, out_dict[chr][gene]['isoform_region_matrix'],fmt='%.d',delimiter=',')
    
    # mystr = bio.getvalue().decode('latin1')
    # with open(output_path+'/sr_A.out','w') as f:
    #     f.write(mystr)
    # bio = io.BytesIO()
    # for chr in long_read_gene_matrix_dict:
    #     for gene in long_read_gene_matrix_dict[chr]:
    #         bio.write(str.encode('{}\n'.format(gene)))
    #         np.savetxt(bio, long_read_gene_matrix_dict[chr][gene]['isoform_region_matrix'],fmt='%.d',delimiter=',')
    
    # mystr = bio.getvalue().decode('latin1')
    # with open(output_path+'/lr_A.out','w') as f:
    #     f.write(mystr)
    # 支持多 SR/LR 矩阵（列表）或单个矩阵；结构性操作（K值/奇异值）统一用第一个
    sr_matrix_list = short_read_gene_matrix_dict if isinstance(short_read_gene_matrix_dict, list) else [short_read_gene_matrix_dict]
    sr_matrix_first = sr_matrix_list[0]
    lr_matrix_list = long_read_gene_matrix_dict if isinstance(long_read_gene_matrix_dict, list) else [long_read_gene_matrix_dict]
    lr_matrix_first = lr_matrix_list[0]
    list_of_all_genes_chrs = []
    for chr_name in lr_matrix_first:
        if chr_name in sr_matrix_first:
            for gene_name in lr_matrix_first[chr_name]:
                if gene_name in sr_matrix_first[chr_name]:
                    list_of_all_genes_chrs.append((gene_name,chr_name))
    if config.output_matrix_info:
        # matrix_info 显示实际 per-sample 矩阵（含 b 向量）；k 值计算用理论矩阵
        _sr_info_raw = sr_actual_matrix_for_info if sr_actual_matrix_for_info is not None else short_read_gene_matrix_dict
        _sr_info_list = _sr_info_raw if isinstance(_sr_info_raw, list) else [_sr_info_raw]
        sr_for_matrix = (_sr_info_list if len(_sr_info_list) > 1 else _sr_info_list[0]) if has_sr_data else None
        lr_for_matrix = (lr_matrix_list if len(lr_matrix_list) > 1 else lr_matrix_first) if has_lr_data else None
        output_matrix_info(sr_for_matrix, lr_for_matrix, list_of_all_genes_chrs, gene_points_dict, output_path)
    n_sr = len(sr_matrix_list)
    n_lr = len(lr_matrix_list)
    sr_prefixes = ['SR'] if n_sr == 1 else [f'SR_{i}' for i in range(n_sr)]
    lr_prefixes = ['LR'] if n_lr == 1 else [f'LR_{i}' for i in range(n_lr)]

    # 预计算 Hybrid 条件数（仅用实际有数据的平台的矩阵堆叠）
    n_actual_lr = n_lr if has_lr_data else 0
    n_actual_sr = n_sr if has_sr_data else 0
    need_hybrid = (n_actual_lr + n_actual_sr) > 1
    hybrid_cond_dict = {}  # (gene_name, chr_name) -> (kvalue, reg_cond, gen_cond, svp)
    if need_hybrid:
        hybrid_dicts = (lr_matrix_list if has_lr_data else []) + (sr_matrix_list if has_sr_data else [])
        for (gname, rname) in list_of_all_genes_chrs:
            A_list = [d[rname][gname]['isoform_region_matrix']
                      for d in hybrid_dicts
                      if rname in d and gname in d[rname]]
            if len(A_list) > 1:
                try:
                    hybrid_cond_dict[(gname, rname)] = get_condition_number(np.vstack(A_list))
                except Exception:
                    hybrid_cond_dict[(gname, rname)] = (float('nan'),) * 4
            else:
                hybrid_cond_dict[(gname, rname)] = (float('nan'),) * 4

    # 奇异值文件：每个矩阵单独输出（多矩阵时加编号后缀）
    for i, sr_dict in enumerate(sr_matrix_list):
        fname = 'SR_singular_values.out' if n_sr == 1 else f'SR_{i}_singular_values.out'
        with open(output_path + '/' + fname, 'w') as f:
            f.write('Gene\tSingular_values\n')
            for (gname, rname) in list_of_all_genes_chrs:
                if rname in sr_dict and gname in sr_dict[rname]:
                    svalues = ','.join([str(v) for v in sr_dict[rname][gname]['singular_values']])
                else:
                    svalues = 'NA'
                f.write('{}\t{}\n'.format(gname, svalues))
            for chr_name in removed_gene_isoform_dict:
                for gene_name in removed_gene_isoform_dict[chr_name]:
                    f.write('{}\tNA\n'.format(gene_name))
    if has_lr_data:
        for i, lr_dict in enumerate(lr_matrix_list):
            fname = 'LR_singular_values.out' if n_lr == 1 else f'LR_{i}_singular_values.out'
            with open(output_path + '/' + fname, 'w') as f:
                f.write('Gene\tSingular_values\n')
                for (gname, rname) in list_of_all_genes_chrs:
                    if rname in lr_dict and gname in lr_dict[rname]:
                        svalues = ','.join([str(v) for v in lr_dict[rname][gname]['singular_values']])
                    else:
                        svalues = 'NA'
                    f.write('{}\t{}\n'.format(gname, svalues))
                for chr_name in removed_gene_isoform_dict:
                    for gene_name in removed_gene_isoform_dict[chr_name]:
                        f.write('{}\tNA\n'.format(gene_name))

    gene_feature_dict = {}

    # kvalues_gene.out：每个 SR/LR 矩阵各输出一组 K 值列
    sr_gene_header = '\t'.join(['{0}_k_value\t{0}_regular_condition_number\t{0}_generalized_condition_number\t{0}_A_dim'.format(p) for p in sr_prefixes])
    lr_gene_header = '\t'.join(['{0}_k_value\t{0}_regular_condition_number\t{0}_generalized_condition_number\t{0}_A_dim'.format(p) for p in lr_prefixes])
    hybrid_gene_header = '\tHybrid_k_value\tHybrid_regular_condition_number\tHybrid_generalized_condition_number' if need_hybrid else ''
    with open(output_path+"/kvalues_gene.out",'w') as f:
        f.write('Gene\tChr\tNum_isoforms\tNum_exons\tNum_split_exons\t{}\t{}{}\n'.format(sr_gene_header, lr_gene_header, hybrid_gene_header))
        for (gene_name,chr_name) in list_of_all_genes_chrs:
            num_isoforms,num_exons,num_split_exons = gene_num_isoform_dict[chr_name][gene_name],raw_gene_num_exon_dict[chr_name][gene_name],gene_num_exon_dict[chr_name][gene_name]
            row_vals = [gene_name, chr_name, num_isoforms, num_exons, num_split_exons]
            for sr_dict in sr_matrix_list:
                SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,_ = sr_dict[chr_name][gene_name]['condition_number']
                SR_A_dim = sr_dict[chr_name][gene_name]['isoform_region_matrix'].shape
                row_vals += [SR_kvalue, SR_regular_condition_number, SR_generalized_condition_number, SR_A_dim]
            for lr_dict in lr_matrix_list:
                LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,_ = lr_dict[chr_name][gene_name]['condition_number']
                LR_A_dim = lr_dict[chr_name][gene_name]['isoform_region_matrix'].shape
                row_vals += [LR_kvalue, LR_regular_condition_number, LR_generalized_condition_number, LR_A_dim]
            if need_hybrid:
                H_kvalue, H_reg, H_gen, _ = hybrid_cond_dict.get((gene_name, chr_name), (float('nan'),)*4)
                row_vals += [H_kvalue, H_reg, H_gen]
            f.write('\t'.join([str(v) for v in row_vals]) + '\n')
            gene_feature_dict[gene_name] = [sr_matrix_first[chr_name][gene_name]['condition_number'][2], lr_matrix_first[chr_name][gene_name]['condition_number'][2]]
        for chr_name in removed_gene_isoform_dict:
            for gene_name in removed_gene_isoform_dict[chr_name]:
                info_dict = removed_gene_isoform_dict[chr_name][gene_name]['info']
                num_isoforms,num_exons,num_split_exons = info_dict['num_isoforms'],info_dict['num_exons'],info_dict['num_split_exons']
                na_cols = ['NA'] * (4 * (n_sr + n_lr) + (3 if need_hybrid else 0))
                f.write('\t'.join([str(v) for v in [gene_name, chr_name, num_isoforms, num_exons, num_split_exons] + na_cols]) + '\n')

    # kvalues_isoform.out：isoform 自身属性，基因级别指标请 join kvalues_gene.out
    with open(output_path+"/kvalues_isoform.out",'w') as f:
        f.write('Isoform\tGene\tChr\tNum_exons\tIsoform_length\tNum_isoforms\n')
        for (gene_name,chr_name) in list_of_all_genes_chrs:
            for isoform_name in sr_matrix_first[chr_name][gene_name]['isoform_names_indics']:
                num_exons,isoform_length,num_isoforms = raw_isoform_num_exon_dict[isoform_name],isoform_length_dict[isoform_name],num_isoforms_dict[isoform_name]
                f.write('\t'.join([str(v) for v in [isoform_name, gene_name, chr_name, num_exons, isoform_length, num_isoforms]]) + '\n')
        for chr_name in removed_gene_isoform_dict:
            for gene_name in removed_gene_isoform_dict[chr_name]:
                isoform_info_dict = removed_gene_isoform_dict[chr_name][gene_name]['isoform_info']
                for isoform_name in isoform_info_dict:
                    num_exons,isoform_length,num_isoforms = isoform_info_dict[isoform_name]['num_exons'],isoform_info_dict[isoform_name]['isoform_length'],removed_gene_isoform_dict[chr_name][gene_name]['info']['num_isoforms']
                    f.write('\t'.join([str(v) for v in [isoform_name, gene_name, chr_name, num_exons, isoform_length, num_isoforms]]) + '\n')
    return gene_feature_dict
def generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict,same_structure_isoform_dict,removed_gene_isoform_dict,gene_points_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    if config.output_matrix_info:
        output_matrix_info(short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_points_dict,output_path)
    # with open(output_path+'lr.pkl','wb') as f:
    #     pickle.dump(long_read_gene_matrix_dict,f)
    with open(output_path+"/expression_gene.out",'w') as f_gene:
        with open(output_path+"/expression_isoform.out",'w') as f_isoform:
            f_gene.write('Gene\tChr\tTPM\n')
            f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tAlpha\n')
            # f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\n')
            for gene_name,chr_name in list_of_all_genes_chrs:
                tpm_sum = 0
                # sr_expected_counts_sum = 0
                # lr_expected_counts_sum = 0
                for isoform_name in sr_matrix_first[chr_name][gene_name]['isoform_names_indics']:
                    start_pos = min(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                    end_pos = max(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'])
                    isoform_len = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                    isoform_index = sr_matrix_first[chr_name][gene_name]['isoform_names_indics'][isoform_name]
                    tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'][isoform_index]
                    # sr_expected_counts = gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_expected_counts'][isoform_index]
                    # lr_expected_counts = gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_expected_counts'][isoform_index]
                    alpha = gene_isoform_tpm_expression_dict[chr_name][gene_name]['alpha']
                    tpm_sum += tpm
                    # sr_expected_counts_sum += sr_expected_counts
                    # lr_expected_counts_sum += lr_expected_counts
                    # if chr_name in same_structure_isoform_dict:
                    #     if gene_name in same_structure_isoform_dict[chr_name]:
                    #         if isoform_name in same_structure_isoform_dict[chr_name][gene_name]:
                    #             num_same_structure_isoforms = len(same_structure_isoform_dict[chr_name][gene_name][isoform_name])
                    #             tpm = tpm/(num_same_structure_isoforms+1)
                                # sr_expected_counts = sr_expected_counts/(num_same_structure_isoforms+1)
                                # lr_expected_counts = lr_expected_counts/(num_same_structure_isoforms+1)
                                # for same_structure_isoform in same_structure_isoform_dict[chr_name][gene_name][isoform_name]:
                                #     f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(same_structure_isoform,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,sr_expected_counts,lr_expected_counts,alpha))
                    f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,alpha))
                    # SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number = sr_matrix_first[chr_name][gene_name]['condition_number']
                    # LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number = lr_matrix_first[chr_name][gene_name]['condition_number']

                    # f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
                f_gene.write('{}\t{}\t{}\n'.format(gene_name,chr_name,tpm_sum))
                
            # censored gene output  
            for chr_name in removed_gene_isoform_dict:
                for gene_name in removed_gene_isoform_dict[chr_name]:
                    f_gene.write('{}\t{}\t0\t0\t0\n'.format(gene_name,chr_name))
                    for isoform_name in removed_gene_isoform_dict[chr_name][gene_name]['isoforms_length_dict']:
                        start_pos = min(removed_gene_isoform_dict[chr_name][gene_name]['raw_isoform_exons_dict'][isoform_name]['start_pos'])
                        end_pos = max(removed_gene_isoform_dict[chr_name][gene_name]['raw_isoform_exons_dict'][isoform_name]['end_pos'])
                        isoform_len = removed_gene_isoform_dict[chr_name][gene_name]['isoforms_length_dict'][isoform_name]
                        tpm = 0
                        sr_expected_counts = 0
                        lr_expected_counts = 0
                        f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,sr_expected_counts,lr_expected_counts))
            
    # with open(output_path+"/expression_isoform_result1.out",'w') as f:
    #     f.write('Isoform\tGene\tChr\tTPM\tTPM_SR\tTPM_LR\n')
    #     for chr_name in short_read_gene_matrix_dict:
    #             for gene_name in short_read_gene_matrix_dict[chr_name]:
    #                 for isoform_name in sr_matrix_first[chr_name][gene_name]['isoform_names_indics']:
    #                     isoform_index = sr_matrix_first[chr_name][gene_name]['isoform_names_indics'][isoform_name]
    #                     tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_by_gene']
    #                     perfect_tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['perfect_tpm_by_gene']
    #                     tpm_sr = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_sr']
    #                     tpm_lr = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_lr']