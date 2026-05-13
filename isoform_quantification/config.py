same_struc_isoform_handling = 'keep'
add_full_length_region = 'all'
READ_JUNC_MIN_MAP_LEN = 1
multi_exon_region_weight = 'regular'
output_matrix_info = False
normalize_sr_A = False
sr_region_selection = 'real_data'
keep_sr_exon_region = 'nonfullrank'
normalize_lr_A = True
kallisto_index = ''
eff_len_option = 'kallisto'
EM_SR_num_iters = 1000
alpha_df_path = None
alpha = 0.5
EM_output_frequency = 50
isoform_start_end_site_tolerance = 20
junction_site_tolerance = 5
eps_strategy = 'add_eps_small'
pseudo_count_SR = None
pseudo_count_LR = 1
std_f_len = None
mean_f_len = None
read_len_dist_sm_dict_path = None
LR_cond_prob_calc = 'form_1'
sr_design_matrix = 'regular'
output_path = None
pretrained_model_path = None
threads = None
singular_values_tol = 0
lr_global_weights = None   # None or list of global weights (all LR+SR sum=1)
sr_global_weights = None   # None or list of global weights (all LR+SR sum=1)
lr_within_weights = None   # None or list of within-group weights (sum=1 within LR group)
sr_within_weights = None   # None or list of within-group weights (sum=1 within SR group)
use_quality_weights = False  # Whether to use unique mapping rate as sample weights (default: equal weights)
normalize_q = False          # Whether to normalize q by its sum before weighting in multi-sample E-step (default: False)