import argparse
from TrEESR import TrEESR
from TransELS import TransELS
from EM import EM,EM_SR,EM_hybrid,EM_hybrid_multi
import config
import os
# import os
# os.system("taskset -p 0xfffff %d" % os.getpid())
# affinity_mask = os.sched_getaffinity(0)
# os.sched_setaffinity(0, affinity_mask)
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_TrEESR = subparsers.add_parser('cal_identifiability', aliases=['cal_K_value','TrEESR'],help='Calculate isoform identifiability (K values)')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for calculation of K value')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-srsam','--short_read_sam_path', type=str, nargs='+', help="The path(s) of short read sam file(s); read length inferred per file",required=False)
    optional_TrEESR.add_argument('-lrsam','--long_read_sam_path', type=str, nargs='+', help="The path(s) of long read sam file(s)",required=False)
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_TrEESR.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:read_length][read_length,num_exons,real_data]")
    optional_TrEESR.add_argument('--lr_region_selection',type=str, default='read_length',help="LR region selection methods [default:read_length][read_length,real_data]")
    optional_TrEESR.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")
    optional_TrEESR.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:True][True,False]")
    optional_TrEESR.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_TrEESR.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_TrEESR.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_TrEESR.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_TrEESR.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_TrEESR.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_TrEESR.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--add_full_length_region',type=str, default='nonfullrank',help="Whether add full length region[default:nonfullrank] [all,nonfullrank,none]")
    optional_TrEESR.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    weight_path = os.path.dirname(os.path.realpath(__file__))+'/weights/nanosim_weight_dict.pkl'
    # assert os.path.exists(weight_path)
    optional_TrEESR.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")
    optional_TrEESR.add_argument('--kde_lr_model_path',type=str, default=None,help="Path to a pre-trained KDE model (joblib format) for LR region weight calculation")
    optional_TrEESR.add_argument('--kde_lr',action='store_true',default=True,help="Train KDE model from input LR data to weight the LR A matrix")
    optional_TrEESR.add_argument('--keep_kde_lr',action='store_true',default=False,help="Keep the trained KDE model file after the run (default: delete after use)")
    optional_TrEESR.add_argument('--identi_data', type=str, default=None,
        help="Path to feature_data.tsv for computing identifiability WITHOUT SAM files. "
             "Columns (tab-separated, header optional): "
             "(1) transcript_name, (2) gene_name, (3) gene_count, "
             "(4) transcript_abundance, (5) read_length. "
             "When this option is set, -srsam / -lrsam are not required.")

    parser_EM = subparsers.add_parser('quantify', aliases=['EM'],help='Isoform quantification by EM algorithm')
    requiredNamed_EM = parser_EM.add_argument_group('required named arguments for isoform quantification')
    requiredNamed_EM.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_EM.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_EM = parser_EM.add_argument_group('optional arguments')
    optional_EM.add_argument('-lrsam','--long_read_sam_path', type=str, nargs='+', help="The path(s) of long read sam file(s)",required=False,default=None)
    optional_EM.add_argument('-srsam','--short_read_sam_path', type=str, nargs='+', help="The path(s) of short read sam file(s)",default=None)
    optional_EM.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    optional_EM.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    optional_EM.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    optional_EM.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    optional_EM.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, kallisto,Salmon, RSEM] [default:kallisto]",default='kallisto')
    # optional_EM.add_argument('--kallisto_index', type=str, help="kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/kallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    optional_EM.add_argument('--alpha',type=str,default='0.5', help="Alpha[default:0.5]: SR and LR balance parameter (0~1, or 'adaptive' to auto-predict per gene community)")
    optional_EM.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    optional_EM.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    optional_EM.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    optional_EM.add_argument('--training',type=str,default='False', help="Generate training dict")
    optional_EM.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    optional_EM.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    optional_EM.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_EM.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_EM.add_argument('--use_weight_matrix',type=str, default='False',help="Whether use weight matrix[default:True][True,False]")
    optional_EM.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    # optional_EM.add_argument('--same_struc_isoform_handling',type=str, default='keep',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_EM.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_EM.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_EM.add_argument('--sr_design_matrix',type=str, default='weight',help="How to calculate design matrix [default:weight][weight,binary]")
    optional_EM.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_EM.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    optional_EM.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:read_length][read_length,num_exons,real_data]")
    optional_EM.add_argument('--keep_sr_exon_region',type=str, default='nonfullrank',help="Keep exon region for SR if using real data to filter region nonfullrank: only keep zero count exon region in non fulll rank gene [default:nonfullrank][nonfullrank,all,none]")
    optional_EM.add_argument('--region_weight_path',type=str, default=None,help="Mili LR region weight path")
    optional_EM.add_argument('--kde_lr_model_path',type=str, default=None,help="Path to a pre-trained KDE model (joblib format) for LR region weight calculation")
    optional_EM.add_argument('--kde_lr',action='store_true',default=False,help="Train KDE model from input LR data to weight the LR A matrix")
    optional_EM.add_argument('--keep_kde_lr',action='store_true',default=False,help="Keep the trained KDE model file after the run (default: delete after use)")
    optional_EM.add_argument('--EM_choice',type=str, default='LR',help="EM_choice[SR,LR,hybrid]")
    optional_EM.add_argument('--iter_theta',type=str, default='False',help="Whether use updated theta to re-calculate conditional prob [True,False]")
    optional_EM.add_argument('--kde_path',type=str, default='/fs/project/PCON0009/Au-scratch2/haoran/_projects/long_reads_rna_seq_simulator/models/kde_H1-hESC_dRNA',help="KDE model path")
    optional_EM.add_argument('--eff_len_option',type=str, default='kallisto',help="Calculation of effective length option [kallisto,RSEM]")
    optional_EM.add_argument('--EM_SR_num_iters',type=int, default=200,help="Number of EM SR iterations")
    optional_EM.add_argument('--EM_output_frequency',type=int, default=200,help="Frequency(in itertations) of outputting EM results")
    optional_EM.add_argument('--pretrained_model_path',type=str, default='cDNA-ONT',help="The pretrained model path to identify the alpha")
    optional_EM.add_argument('--alpha_df_path',type=str, default=None,help="Alpha df path")
    optional_EM.add_argument('--inital_theta','--initial_theta',type=str, default='uniform',help="initial_theta [LR,SR,LR_unique,SR_unique,uniform,hybrid,hybrid_unique,random]")
    optional_EM.add_argument('--inital_theta_eps','--initial_theta_eps',type=float, default=0.0,help="initial_theta eps [float]")
    optional_EM.add_argument('--eps_strategy',type=str, default='add_eps_small',help="how to add initial_theta eps [add_eps_all,add_eps_small]. (add_eps_small: add isoform with theta < eps with eps. add_eps: add eps to all isoforms)")
    optional_EM.add_argument('--isoform_start_end_site_tolerance',type=int, default=20,help="Isoform Start and end site tolerance for mapping long reads")
    optional_EM.add_argument('--junction_site_tolerance',type=int, default=5,help="Junction site tolerance for mapping long reads")
    optional_EM.add_argument('--read_len_dist_sm_dict_path',type=str, default=None,help="The path of read length distribution for long reads")
    optional_EM.add_argument('--LR_cond_prob_calc',type=str, default='form_2',help="How to calculate LR length distribution [form_1,form_2]")
    optional_EM.add_argument('--singular_values_tol',type=float,default=0,help="Singular value tolerence")
    optional_EM.add_argument('--lr_weights', type=float, nargs='+', default=None,
        help="Per-sample weights for LR files (same order as -lrsam, will be normalized)")
    optional_EM.add_argument('--sr_weights', type=float, nargs='+', default=None,
        help="Per-sample weights for SR files (same order as -srsam, will be normalized)")
    optional_EM.add_argument('--use_quality_weights', action='store_true', default=False,
        help="Use unique mapping rate to automatically weight samples (default: equal weights)")
    optional_EM.add_argument('--normalize_q', action='store_true', default=False,
        help="Normalize each sample's q by its sum before weighting in multi-sample E-step (default: False)")

    args = parser.parse_args()
    if args.filtering == 'True':
        args.filtering = True
    else:
        args.filtering = False

    config.output_path = args.output_path
    config.threads = args.threads
    config.same_struc_isoform_handling = 'keep'
    config.READ_JUNC_MIN_MAP_LEN = args.READ_JUNC_MIN_MAP_LEN
    config.multi_exon_region_weight = args.multi_exon_region_weight
    config.sr_region_selection = args.sr_region_selection
    config.region_weight_path = args.region_weight_path
    config.kde_lr_model_path = args.kde_lr_model_path
    config.use_kde_lr = args.kde_lr
    config.keep_kde_lr = args.keep_kde_lr
    config.sr_design_matrix = args.sr_design_matrix
    if args.output_matrix_info == 'True':
        config.output_matrix_info = True
    else:
        config.output_matrix_info = False
    config.keep_sr_exon_region = args.keep_sr_exon_region
    if args.normalize_sr_A == 'True':
        config.normalize_sr_A = True
    else:
        config.normalize_sr_A = False
    if args.normalize_lr_A == 'True':
        config.normalize_lr_A = True
    else:
        config.normalize_lr_A = False
    if args.use_weight_matrix == 'True':
        config.use_weight_matrix = True
    else:
        config.use_weight_matrix = False
    config.add_full_length_region = args.add_full_length_region
    config.singular_values_tol = args.singular_values_tol

    if args.subparser_name in ['cal_identifiability','cal_K_value','TrEESR']:
        print('[INFO] Computing identifiability metrics')
        identi_data_path = getattr(args, 'identi_data', None)
        if identi_data_path is not None:
            # --identi_data mode: no SAM files needed
            from TrEESR import TrEESR_identi_data
            print(f'[INFO] --identi_data mode: {identi_data_path}', flush=True)
            TrEESR_identi_data(
                ref_file_path=args.gtf_annotation_path,
                output_path=args.output_path,
                feature_data_path=identi_data_path,
                sr_region_selection=args.sr_region_selection,
                threads=args.threads,
                READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,
            )
        else:
            TrEESR(args.gtf_annotation_path,args.output_path,args.short_read_sam_path,args.long_read_sam_path,args.sr_region_selection,args.filtering,args.threads,lr_region_selection=args.lr_region_selection,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    elif args.subparser_name in ['quantify','EM']:
        config.kde_path = args.kde_path
        if args.training == 'True':
            args.training = True
        else:
            args.training = False
        print('[INFO] Isoform quantification by miniQuant',flush=True)
        # 判断是否为多平台模式：有任何输入数据均走 multi 路径
        lr_sam_list = args.long_read_sam_path  # None 或 list
        sr_sam_list = args.short_read_sam_path  # None 或 list
        n_lr = len(lr_sam_list) if lr_sam_list else 0
        n_sr = len(sr_sam_list) if sr_sam_list else 0
        is_multi = (n_lr + n_sr) >= 1
        if is_multi:
            # 多平台模式
            config.EM_SR_num_iters = args.EM_SR_num_iters
            config.inital_theta_eps = args.inital_theta_eps
            config.EM_output_frequency = args.EM_output_frequency
            config.isoform_start_end_site_tolerance = args.isoform_start_end_site_tolerance
            config.junction_site_tolerance = args.junction_site_tolerance
            config.eps_strategy = args.eps_strategy
            config.read_len_dist_sm_dict_path = args.read_len_dist_sm_dict_path
            config.LR_cond_prob_calc = args.LR_cond_prob_calc
            if args.pretrained_model_path in ['cDNA-ONT','dRNA-ONT','cDNA-PacBio']:
                args.pretrained_model_path = os.path.dirname(os.path.realpath(__file__))+'/pretrained_models/' + args.pretrained_model_path +'/'
            config.pretrained_model_path = args.pretrained_model_path
            if sr_sam_list:
                # hybrid 或 LR+SR 多平台
                # alpha 初始值：未指定则暂设0.5，自适应权重模式下会被 EM_algo_hybrid_multi 覆盖
                if args.alpha == 'adaptive':
                    config.alpha = 0.5  # 占位，自适应权重时由程序根据质量分推导并覆盖
                else:
                    config.alpha = float(args.alpha)
                config.alpha_df_path = None
                config.inital_theta = args.inital_theta
                em_choice = args.EM_choice if args.EM_choice != 'LR' else 'LIQA_modified'
            else:
                # 纯LR多平台
                config.alpha = 1
                config.alpha_df_path = None
                config.inital_theta = 'LR'
                em_choice = args.EM_choice if args.EM_choice != 'LR' else 'LIQA_modified'
            # 用户指定权重则归一化后传入；否则传None，由EM_algo_hybrid_multi根据唯一比对率自适应计算
            # 直接传原始值，由 EM_algo_hybrid_multi 统一归一化（保留LR/SR组间比例）
            _lr_weights = list(args.lr_weights) if args.lr_weights is not None else None
            _sr_weights = list(args.sr_weights) if args.sr_weights is not None else None
            config.use_quality_weights = args.use_quality_weights
            config.normalize_q = args.normalize_q
 
            EM_hybrid_multi(args.gtf_annotation_path,sr_sam_list or [],lr_sam_list or [],args.output_path,multi_mapping_filtering=args.multi_mapping_filtering,threads=args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=em_choice,lr_weights=_lr_weights,sr_weights=_sr_weights,use_quality_weights=args.use_quality_weights,alpha=args.alpha)
            return
        # 单文件模式：提取第一个元素以保持向后兼容
        if lr_sam_list is not None:
            args.long_read_sam_path = lr_sam_list[0]
        if sr_sam_list is not None:
            args.short_read_sam_path = sr_sam_list[0]
        if (args.short_read_sam_path is None) or (args.alpha == 1.0):
            args.alpha = 1.0
            args.SR_quantification_option = 'Mini'
        if (args.alpha == 'adaptive'):
            alpha = 'adaptive'
        else:
            try:
                alpha = float(args.alpha)
            except:
                raise Exception('Alpha given is not numeric')
        if (args.beta == 'adaptive'):
            beta = 'adaptive'
        else:
            try:
                beta = float(args.beta)
            except:
                raise Exception('Beta given is not numeric')
        # if args.SR_quantification_option not in ['Mili','kallisto','Salmon','RSEM']:
        #     raise Exception('SR_quantification_option is not valid.Options: [Mili, kallisto,Salmon, RSEM]')
        if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
            args.multi_mapping_filtering = 'no_filtering'
        SR_fastq_list = []
        if args.short_read_fastq is not None:
            SR_fastq_list = [args.short_read_fastq]
        elif args.short_read_mate1_fastq is not None:
            SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
        if args.DL_model is None:
            args.DL_model = args.SR_quantification_option + '.pt'
        config.EM_SR_num_iters = args.EM_SR_num_iters
        config.inital_theta_eps = args.inital_theta_eps
        config.EM_output_frequency = args.EM_output_frequency
        config.isoform_start_end_site_tolerance = args.isoform_start_end_site_tolerance
        config.junction_site_tolerance = args.junction_site_tolerance
        config.eps_strategy = args.eps_strategy
        config.read_len_dist_sm_dict_path = args.read_len_dist_sm_dict_path
        config.LR_cond_prob_calc = args.LR_cond_prob_calc
        if args.pretrained_model_path in ['cDNA-ONT','dRNA-ONT','cDNA-PacBio']:
            args.pretrained_model_path = os.path.dirname(os.path.realpath(__file__))+'/pretrained_models/' + args.pretrained_model_path +'/'
        config.pretrained_model_path = args.pretrained_model_path
        if args.EM_choice == 'SR':
            config.eff_len_option = args.eff_len_option
            args.long_read_sam_path = None
            args.alpha = 0
            args.inital_theta = 'SR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        elif args.EM_choice == 'hybrid':
            # args.alpha = 0.5
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            if args.alpha_df_path is None:
                config.alpha_df_path = args.output_path +'/hybrid_alpha.tsv'
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        else:
            if args.EM_choice == 'LR':
                args.EM_choice = 'LIQA_modified'
            args.short_read_sam_path = None
            args.alpha = 1
            args.inital_theta = 'LR'
            config.alpha = args.alpha
            config.alpha_df_path = args.alpha_df_path
            config.inital_theta = args.inital_theta
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)    
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()