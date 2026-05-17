# hybrid
import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from pathlib import Path
import glob
import os
import gc
from EM_hybrid.EM_LR import prepare_LR
from EM_hybrid.EM_SR import prepare_hits
from machine_learning.predict_alpha import predict_alpha
import scipy
import config
import datetime
import scipy.sparse
import scipy.stats
import time
import warnings
from EM_hybrid.util import sp_unique
warnings.filterwarnings('ignore')
def E_step_SR(ANT,theta_eff_len_product_arr):
    q = ANT.multiply(theta_eff_len_product_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def E_step_LR(cond_prob,theta_arr):
    q = cond_prob.multiply(theta_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def _weighted_E_step_LR(cond_prob_list, theta_arr, lr_weights):
    q = np.zeros(theta_arr.shape[0])
    for w, cp in zip(lr_weights, cond_prob_list):
        if cp.shape[0] == 0:
            continue
        q_i = np.asarray(E_step_LR(cp, theta_arr)).flatten()
        n_i = q_i.sum()
        if n_i > 0:
            q += w * q_i / n_i
    return q
def _weighted_E_step_SR(ANT_list, theta_eff_len_product_arr, sr_weights):
    q = np.zeros(theta_eff_len_product_arr.shape[0])
    for w, ant in zip(sr_weights, ANT_list):
        if ant.shape[0] == 0:
            continue
        q_j = np.asarray(E_step_SR(ant, theta_eff_len_product_arr)).flatten()
        n_j = q_j.sum()
        if n_j > 0:
            q += w * q_j / n_j
    return q

def _compute_sr_quality(hits_dir, threads):
    """从ANT矩阵计算SR样本质量分（唯一比对率）。唯一比对率越高越可靠，权重越大。"""
    num_unique = 0
    num_total = 0
    for worker_id in range(threads):
        fpath = f'{hits_dir}/{worker_id}_ANT.npz'
        if Path(fpath).exists():
            ANT = scipy.sparse.load_npz(fpath).tocsr()
            row_nnz = np.diff(ANT.indptr)
            num_unique += int((row_nnz == 1).sum())
            num_total += ANT.shape[0]
    return num_unique / max(num_total, 1)

def _compute_lr_quality(cond_prob_dir, threads):
    """从cond_prob矩阵计算LR样本质量分（唯一比对率）。"""
    num_unique = 0
    num_total = 0
    for worker_id in range(threads):
        fpath = f'{cond_prob_dir}/{worker_id}_cond_prob.npz'
        if Path(fpath).exists():
            cp = scipy.sparse.load_npz(fpath).tocsr()
            row_nnz = np.diff(cp.indptr)
            num_unique += int((row_nnz == 1).sum())
            num_total += cp.shape[0]
    return num_unique / max(num_total, 1)
def M_step(isoform_q_arr_SR_all,isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha):
    ss = eff_len_arr / ((theta_arr * eff_len_arr).sum())
    # if alpha_df is None:
    new_theta_arr = ((1-alpha) * isoform_q_arr_SR_all + alpha * isoform_q_arr_LR_all) / ((1-alpha) * isoform_q_arr_SR_all.sum() * ss + alpha * isoform_q_arr_LR_all.sum())
    new_theta_arr = np.nan_to_num(new_theta_arr,0)
    # else:
    #     new_theta_arr = ((1-alpha_df) * isoform_q_arr_SR_all + alpha_df * isoform_q_arr_LR_all) / ((1-alpha_df) * isoform_q_df_SR.sum() * ss + alpha_df * isoform_q_arr_LR_all.sum())
    if new_theta_arr.sum() != 0:
        new_theta_arr = new_theta_arr/new_theta_arr.sum()
    new_theta_arr[new_theta_arr<1e-100] = 0
    return new_theta_arr
def M_step_weighted(q_SR, q_LR, theta_arr, eff_len_arr):
    """M-step using global weights directly (no alpha parameter needed).
    q_SR and q_LR are already globally weighted, so their sums encode the
    effective SR/LR balance (equivalent to 1-alpha and alpha respectively).
    """
    ss = eff_len_arr / ((theta_arr * eff_len_arr).sum())
    sr_sum = float(np.asarray(q_SR).sum())
    lr_sum = float(np.asarray(q_LR).sum())
    new_theta_arr = (np.asarray(q_SR).flatten() + np.asarray(q_LR).flatten()) / (sr_sum * ss + lr_sum)
    new_theta_arr = np.nan_to_num(new_theta_arr, 0)
    if new_theta_arr.sum() != 0:
        new_theta_arr = new_theta_arr / new_theta_arr.sum()
    new_theta_arr[new_theta_arr < 1e-100] = 0
    return new_theta_arr
def build_dummy_gene_reads(isoform_gene_dict,isoform_index_dict,num_isoforms):
    gene_index_dict = {}
    for isoform,gene in isoform_gene_dict.items():
        if gene not in gene_index_dict:
            gene_index_dict[gene] = set()
        gene_index_dict[gene].add(isoform_index_dict[isoform])
    num_genes = 0
    row = []
    col = []
    data = []
    gene_list = []
    for gene,isoform_indices in gene_index_dict.items():
        for isoform_index in isoform_indices:
            row.append(num_genes)
            col.append(isoform_index)
            data.append(1)
        num_genes += 1
        gene_list.append(gene)
    dummy_gene = scipy.sparse.coo_matrix((data, (row, col)), shape=(num_genes, num_isoforms))
    return dummy_gene,gene_list
def get_connected_components(cond_prob_matrix):
    biadjacency = scipy.sparse.csr_matrix(cond_prob_matrix)
    adjacency = scipy.sparse.bmat([[None, biadjacency], [biadjacency.T, None]], format='csr')
    adjacency.sort_indices()
    n_components,labels = scipy.sparse.csgraph.connected_components(adjacency,return_labels=True)
    return labels
def serialize_community_worker(worker_id,worker_label,isoform_labels,output_path,sr_dirs=None,lr_dirs=None):
    community_batch_index = worker_id
    worker_isoform = []
    worker_community = []
    for community_index in worker_label:
        community_isoforms_index = (isoform_labels == community_index)
        community_isoform = np.where(community_isoforms_index)[0]
        worker_isoform.append(community_isoform)
        worker_community.append(community_index)
    with open(f'{output_path}/temp/community/{community_batch_index}_isoform.pkl','wb') as f:
        pickle.dump(worker_isoform,f, protocol=4)
    with open(f'{output_path}/temp/community/{community_batch_index}_community.pkl','wb') as f:
        pickle.dump(worker_community,f, protocol=4)
    for reads_batch_id in range(config.threads):
        if sr_dirs is not None:
            for s_idx, sr_dir in enumerate(sr_dirs):
                fpath = f'{sr_dir}/{reads_batch_id}_ANT.npz'
                worker_ANT_s = []
                if Path(fpath).exists():
                    ANT_s = scipy.sparse.csc_matrix(scipy.sparse.load_npz(fpath))
                    for community_index in worker_label:
                        community_isoforms_index = (isoform_labels == community_index)
                        comm_ANT = scipy.sparse.csr_matrix(ANT_s[:, community_isoforms_index])
                        comm_ANT = comm_ANT[comm_ANT.getnnz(1) > 0]
                        worker_ANT_s.append(comm_ANT)
                else:
                    for community_index in worker_label:
                        community_isoforms_index = (isoform_labels == community_index)
                        worker_ANT_s.append(scipy.sparse.csr_matrix((0, int(community_isoforms_index.sum()))))
                with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_ANT_sr{s_idx}.pkl', 'wb') as f:
                    pickle.dump(worker_ANT_s, f, protocol=4)
            for l_idx, lr_dir in enumerate(lr_dirs):
                fpath = f'{lr_dir}/{reads_batch_id}_cond_prob.npz'
                worker_cond_prob_l = []
                if Path(fpath).exists():
                    cond_prob_l = scipy.sparse.csc_matrix(scipy.sparse.load_npz(fpath))
                    for community_index in worker_label:
                        community_isoforms_index = (isoform_labels == community_index)
                        comm_cp = scipy.sparse.csr_matrix(cond_prob_l[:, community_isoforms_index])
                        comm_cp = comm_cp[comm_cp.getnnz(1) > 0]
                        worker_cond_prob_l.append(comm_cp)
                else:
                    for community_index in worker_label:
                        community_isoforms_index = (isoform_labels == community_index)
                        worker_cond_prob_l.append(scipy.sparse.csr_matrix((0, int(community_isoforms_index.sum()))))
                with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_cond_prob_lr{l_idx}.pkl', 'wb') as f:
                    pickle.dump(worker_cond_prob_l, f, protocol=4)
        else:
            worker_ANT = []
            worker_cond_prob = []
            ANT = scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{reads_batch_id}_ANT.npz')
            ANT = scipy.sparse.csc_matrix(ANT)
            cond_prob = scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{reads_batch_id}_cond_prob.npz')
            cond_prob = scipy.sparse.csc_matrix(cond_prob)
            for community_index in worker_label:
                community_isoforms_index = (isoform_labels == community_index)
                community_ANT = scipy.sparse.csr_matrix(ANT[:,community_isoforms_index])
                community_ANT = community_ANT[community_ANT.getnnz(1)>0]
                community_cond_prob = scipy.sparse.csr_matrix(cond_prob[:,community_isoforms_index])
                community_cond_prob = community_cond_prob[community_cond_prob.getnnz(1)>0]
                worker_ANT.append(community_ANT)
                worker_cond_prob.append(community_cond_prob)
            with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_ANT.pkl','wb') as f:
                pickle.dump(worker_ANT,f, protocol=4)
            with open(f'{output_path}/temp/community/{community_batch_index}_{reads_batch_id}_cond_prob.pkl','wb') as f:
                pickle.dump(worker_cond_prob,f, protocol=4)
    # for worker_label,community_batch_index in zip(all_worker_labels,range(len(all_worker_labels))):
    #     worker_ANT = []
    #     worker_cond_prob = []
    #     if worker_id == 0:
    #         worker_isoform = []
    #         worker_community = []
    #     for community_index in worker_label:
    #         community_SR_index = (SR_labels == community_index).nonzero()[0]
    #         community_LR_index = (LR_labels == community_index).nonzero()[0]
    #         community_isoforms_index = (isoform_labels == community_index)
    #         community_isoform = np.where(community_isoforms_index)[0]
    #         worker_community_SR_index = community_SR_index[(community_SR_index < num_processed_SR + ANT.shape[0]) & (community_SR_index >= num_processed_SR)].copy()
    #         worker_community_LR_index = community_LR_index[(community_LR_index < num_processed_LR + cond_prob.shape[0]) & (community_LR_index >= num_processed_LR)].copy()
    #         worker_community_SR_index -= num_processed_SR
    #         worker_community_LR_index -= num_processed_LR
    #         community_ANT = ANT[worker_community_SR_index,:][:,community_isoforms_index]
    #         community_cond_prob = cond_prob[worker_community_LR_index,:][:,community_isoforms_index]
    #         worker_ANT.append(community_ANT)
    #         worker_cond_prob.append(community_cond_prob)
    #         if worker_id == 0:
    #             worker_isoform.append(community_isoform)
    #             worker_community.append(community_index)
    #     with open(f'{output_path}/temp/community/{community_batch_index}_{worker_id}_ANT.pkl','wb') as f:
    #         pickle.dump(worker_ANT,f, protocol=4)
    #     with open(f'{output_path}/temp/community/{community_batch_index}_{worker_id}_cond_prob.pkl','wb') as f:
    #         pickle.dump(worker_cond_prob,f, protocol=4)
    #     if worker_id == 0:
    #         with open(f'{output_path}/temp/community/{community_batch_index}_isoform.pkl','wb') as f:
    #             pickle.dump(worker_isoform,f, protocol=4)
    #         with open(f'{output_path}/temp/community/{community_batch_index}_community.pkl','wb') as f:
    #             pickle.dump(worker_community,f, protocol=4)
def serialize_community(isoform_labels,output_path,threads,sr_dirs=None,lr_dirs=None):
    Path(f'{output_path}/temp/community/').mkdir(exist_ok=True,parents=True)
    pool = mp.Pool(threads)
    unique_labels = np.unique(isoform_labels)
    np.random.shuffle(unique_labels)
    all_worker_labels = np.array_split(unique_labels, threads)
    futures = []
    for worker_id in range(threads):
        futures.append(pool.apply_async(serialize_community_worker,
            (worker_id,all_worker_labels[worker_id],isoform_labels,output_path,sr_dirs,lr_dirs),
            error_callback=callback_error))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
def construct_community_helper(isoform_gene_dict,isoform_index_dict,output_path,threads):
    st = time.time()
    ANT = []
    cond_prob = []
    # num_processed_SRs = []
    # num_processed_LRs = []
    for worker_id in range(threads):
        worker_ANT = scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{worker_id}_ANT.npz')
        worker_cond_prob = scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{worker_id}_cond_prob.npz')
        ANT.append(worker_ANT)
        cond_prob.append(worker_cond_prob)
        # num_processed_SRs.append(worker_ANT.shape[0])
        # num_processed_LRs.append(worker_cond_prob.shape[0])
    del worker_ANT
    del worker_cond_prob
    ANT = scipy.sparse.vstack(ANT)
    cond_prob = scipy.sparse.vstack(cond_prob)
    ANT = scipy.sparse.csr_matrix(ANT)
    cond_prob = scipy.sparse.csr_matrix(cond_prob)
    time1 = time.time()
    # if threads > 1:
    #     num_processed_SRs = [0]+[sum(np.array(num_processed_SRs)[:i]) for i in range(threads-1)]
    #     num_processed_LRs = [0]+[sum(np.array(num_processed_LRs)[:i]) for i in range(threads-1)]
    # else:
    #     num_processed_SRs = [0]
    #     num_processed_LRs = [0]
    num_SRs,num_LRs,num_isoforms = ANT.shape[0], cond_prob.shape[0],ANT.shape[1]
#     print(f'Number of SRs:{num_SRs}')
#     print(f'Number of LRs:{num_LRs}')
    dummy_gene,gene_list = build_dummy_gene_reads(isoform_gene_dict,isoform_index_dict,num_isoforms)
    time2 = time.time()
    reads_matrix_uniq = sp_unique(scipy.sparse.vstack([ANT.sign(),cond_prob.sign()]), axis=0)
    time3 = time.time()
    del ANT
    del cond_prob
    cond_prob_matrix = scipy.sparse.vstack([reads_matrix_uniq,dummy_gene])
    labels = get_connected_components(cond_prob_matrix)
    time4 = time.time()
    isoform_labels = labels[cond_prob_matrix.shape[0]:]
    gene_community_id_dict = {}
    for label,gene in zip(labels[reads_matrix_uniq.shape[0]:cond_prob_matrix.shape[0]],gene_list):
        gene_community_id_dict[gene] = label
    with open(f'{output_path}/temp/machine_learning/gene_community_id_dict.pkl','wb') as f:
        pickle.dump(gene_community_id_dict,f)
    # print(time1-st)
    # print(time2-time1)
    # print(time3-time2)
    # print(time4-time3)
    return isoform_labels,num_SRs,num_LRs
def construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads,sr_dirs=None,lr_dirs=None):
    isoform_labels,num_SRs,num_LRs = construct_community_helper(isoform_gene_dict,isoform_index_dict,output_path,threads)
    st = time.time()
    serialize_community(isoform_labels,output_path,threads,sr_dirs=sr_dirs,lr_dirs=lr_dirs)
    time5 = time.time()
    # print(time5-st)
    return num_SRs,num_LRs
def EM_worker(worker_id,output_df,output_path,eff_len_arr,num_SRs,num_LRs):
    num_iters = config.EM_SR_num_iters
    min_diff = 1e-6
    with open(f'{output_path}/temp/community/{worker_id}_isoform.pkl','rb') as f:
        worker_isoform = pickle.load(f)
    with open(f'{output_path}/temp/community/{worker_id}_community.pkl','rb') as f:
        worker_community = pickle.load(f)

    lr_weights = config.lr_global_weights
    sr_weights = config.sr_global_weights
    is_weighted = (lr_weights is not None) or (sr_weights is not None)

    if is_weighted:
        num_lr = len(lr_weights) if lr_weights is not None else 0
        num_sr = len(sr_weights) if sr_weights is not None else 0
        worker_cond_prob_per_sample = [[] for _ in range(num_lr)]
        worker_ANT_per_sample = [[] for _ in range(num_sr)]
        for i in range(config.threads):
            for l in range(num_lr):
                with open(f'{output_path}/temp/community/{worker_id}_{i}_cond_prob_lr{l}.pkl','rb') as f:
                    temp = pickle.load(f)
                if i == 0:
                    worker_cond_prob_per_sample[l] = temp
                else:
                    for j in range(len(temp)):
                        worker_cond_prob_per_sample[l][j] = scipy.sparse.vstack(
                            [worker_cond_prob_per_sample[l][j], temp[j]])
            for s in range(num_sr):
                with open(f'{output_path}/temp/community/{worker_id}_{i}_ANT_sr{s}.pkl','rb') as f:
                    temp = pickle.load(f)
                if i == 0:
                    worker_ANT_per_sample[s] = temp
                else:
                    for j in range(len(temp)):
                        worker_ANT_per_sample[s][j] = scipy.sparse.vstack(
                            [worker_ANT_per_sample[s][j], temp[j]])
    else:
        worker_ANT, worker_cond_prob = [], []
        for i in range(config.threads):
            with open(f'{output_path}/temp/community/{worker_id}_{i}_ANT.pkl','rb') as f:
                temp_ANT = pickle.load(f)
            with open(f'{output_path}/temp/community/{worker_id}_{i}_cond_prob.pkl','rb') as f:
                temp_cond_prob = pickle.load(f)
            if i == 0:
                worker_ANT, worker_cond_prob = temp_ANT, temp_cond_prob
            else:
                for j in range(len(temp_ANT)):
                    worker_ANT[j] = scipy.sparse.vstack([worker_ANT[j], temp_ANT[j]])
                    worker_cond_prob[j] = scipy.sparse.vstack([worker_cond_prob[j], temp_cond_prob[j]])

    all_LR_expression_df = []
    all_SR_expression_df = []
    all_community_iteration_df = []

    for community_idx, (community_isoform, community_id) in enumerate(zip(worker_isoform, worker_community)):
        if config.alpha_df_path is None:
            alpha = float(config.alpha)
        else:
            alpha_df = pd.read_csv(config.alpha_df_path,sep='\t')
            if len(alpha_df.columns) == 2:
                alpha_df['community_id'] = alpha_df['community_id'].astype(int)
                alpha_df = alpha_df.set_index('community_id')
                alpha = alpha_df.loc[community_id].values[0]
            elif len(alpha_df.columns) == 3:
                alpha_df['id'] = alpha_df['id'].astype(int)
                alpha_df = alpha_df.set_index('id')
                alpha = alpha_df.loc[community_id]['alpha']

        if is_weighted:
            community_cond_prob_list = [worker_cond_prob_per_sample[l][community_idx] for l in range(num_lr)]
            community_ANT_list = [worker_ANT_per_sample[s][community_idx] for s in range(num_sr)]
            community_num_LRs = sum(cp.shape[0] for cp in community_cond_prob_list)
            community_num_SRs = sum(ant.shape[0] for ant in community_ANT_list)
        else:
            community_cond_prob = worker_cond_prob[community_idx]
            community_ANT = worker_ANT[community_idx]
            community_num_LRs = community_cond_prob.shape[0]
            community_num_SRs = community_ANT.shape[0]

        theta_arr = np.ones(shape=(community_isoform.shape[0]))
        theta_arr = theta_arr/theta_arr.sum()
        community_eff_len_arr = np.array(eff_len_arr)[community_isoform]
        community_iteration_df = []
        for i in range(num_iters):
            theta_eff_len_product_arr = theta_arr * community_eff_len_arr /((theta_arr * community_eff_len_arr).sum())
            if is_weighted:
                isoform_q_arr_LR = _weighted_E_step_LR(community_cond_prob_list, theta_arr, lr_weights if lr_weights is not None else [])
                isoform_q_arr_SR = _weighted_E_step_SR(community_ANT_list, theta_eff_len_product_arr, sr_weights if sr_weights is not None else [])
                new_theta_arr = M_step(isoform_q_arr_SR, isoform_q_arr_LR, theta_arr, community_eff_len_arr, alpha)
            else:
                isoform_q_arr_LR = E_step_LR(community_cond_prob,theta_arr)
                isoform_q_arr_SR = E_step_SR(community_ANT,theta_eff_len_product_arr)
                new_theta_arr = M_step(isoform_q_arr_SR,isoform_q_arr_LR,theta_arr,community_eff_len_arr,alpha)
            new_theta_arr = np.array(new_theta_arr).flatten()
            diff = np.abs(theta_arr[new_theta_arr>1e-7] - new_theta_arr[new_theta_arr>1e-7])/new_theta_arr[new_theta_arr>1e-7]
            iteration_df = pd.DataFrame({'theta':theta_arr},index=community_isoform)
            iteration_df['iteration'] = i
            community_iteration_df.append(iteration_df)
            theta_arr = new_theta_arr
            if diff[diff > min_diff].shape[0] == 0:
                break
        LR_expression = community_num_LRs/num_LRs * theta_arr * 1e6
        LR_expression_df = pd.DataFrame({'TPM':LR_expression,'theta':theta_arr},index=community_isoform)
        LR_expression_df['community'] = community_id
        LR_expression_df['community_num_LRs'] = community_num_LRs
        all_LR_expression_df.append(LR_expression_df)
        theta_eff_len_product_arr = theta_arr /((theta_arr * community_eff_len_arr).sum())
        theta_eff_len_product_arr = np.nan_to_num(theta_eff_len_product_arr,0)
        community_expression = (community_num_SRs * theta_eff_len_product_arr).sum()
        transcript_expression = community_num_SRs * theta_eff_len_product_arr
        SR_expression_df = pd.DataFrame({'transcript_expression':transcript_expression,'theta':theta_arr},index=community_isoform)
        SR_expression_df['community_expression'] = community_expression
        SR_expression_df['community'] = community_id
        SR_expression_df['community_num_SRs'] = community_num_SRs
        all_SR_expression_df.append(SR_expression_df)
        community_iteration_df = pd.concat(community_iteration_df)
        community_iteration_df = output_df.join(community_iteration_df,on='Index',how='inner').sort_values(['iteration','Isoform'])
        community_iteration_df['community'] = community_id
        all_community_iteration_df.append(community_iteration_df)
    if len(all_LR_expression_df) == 0:
        LR_TPM_df = None
    else:
        all_LR_expression_df = pd.concat(all_LR_expression_df)
        LR_TPM_df = output_df.join(all_LR_expression_df,on='Index',how='inner')
    if len(all_SR_expression_df) == 0:
        SR_TPM_df = None
    else:
        all_SR_expression_df = pd.concat(all_SR_expression_df)
        SR_TPM_df = output_df.join(all_SR_expression_df,on='Index',how='inner')
    if len(all_community_iteration_df) == 0:
        all_community_iteration_df = None
    else:
        all_community_iteration_df = pd.concat(all_community_iteration_df)
    return LR_TPM_df,SR_TPM_df,all_community_iteration_df
def callback_error(result):
    print('ERR:', result,flush=True)
def EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads,num_SRs,num_LRs):
    print('[INFO] Start quantification...')  
    st = time.time()
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        futures.append(pool.apply_async(EM_worker,(worker_id,output_df,output_path,eff_len_arr,num_SRs,num_LRs,),error_callback=callback_error))
    all_LR_TPM_df = []
    all_SR_TPM_df = []
    all_iteration_df = []
    for future in futures:
        LR_TPM_df,SR_TPM_df,iteration_df = future.get()
        all_LR_TPM_df.append(LR_TPM_df)
        all_SR_TPM_df.append(SR_TPM_df)
        all_iteration_df.append(iteration_df)
    # all_LR_TPM_df['TPM'] = all_LR_TPM_df['TPM']/all_LR_TPM_df['TPM'].sum() * 1e6
    pool.close()
    pool.join()
    all_LR_TPM_df = pd.concat(all_LR_TPM_df)
    all_LR_TPM_df['num_expected_LRs'] = all_LR_TPM_df['community_num_LRs'] * all_LR_TPM_df['theta']
    all_LR_TPM_df[['Isoform','Gene','TPM','theta','community','community_num_LRs']].to_csv(f'{output_path}/LR_EM_expression.out',sep='\t',index=False)
    all_SR_TPM_df = pd.concat(all_SR_TPM_df)
    all_SR_TPM_df['TPM'] = all_SR_TPM_df['community_expression']/(all_SR_TPM_df['transcript_expression'].sum()) * all_SR_TPM_df['theta'] * 1e6
    all_SR_TPM_df['num_expected_SRs'] = all_SR_TPM_df['community_num_SRs'] * all_SR_TPM_df['theta']
    all_SR_TPM_df[['Isoform','Gene','TPM','Effective length','theta','community','community_num_SRs']].to_csv(f'{output_path}/SR_EM_expression.out',sep='\t',index=False)
    all_SR_TPM_df = all_SR_TPM_df.set_index('Isoform').join(all_LR_TPM_df.set_index('Isoform')[['num_expected_LRs']]).reset_index()
    all_SR_TPM_df[['Isoform','Gene','Effective length','TPM','num_expected_SRs','num_expected_LRs']].sort_values(by=['Gene','Isoform']).to_csv(f'{output_path}/Isoform_abundance.out',sep='\t',index=False)
    all_iteration_df = pd.concat(all_iteration_df)
    all_iteration_df.to_csv(f'{output_path}/EM_iterations.tsv',sep='\t',index=False)
    duration = (time.time() - st)
    print('[INFO] Done in {} seconds at {}!'.format(duration,str(datetime.datetime.now())),flush=True)
def merge_ant_files(hits_dirs,output_path,threads):
    """将多个SR数据集的ANT文件合并到规范目录 temp/hits_dict/。"""
    canonical_dir = f'{output_path}/temp/hits_dict'
    Path(canonical_dir).mkdir(exist_ok=True,parents=True)
    for worker_id in range(threads):
        matrices = []
        for d in hits_dirs:
            fpath = f'{d}/{worker_id}_ANT.npz'
            if Path(fpath).exists():
                matrices.append(scipy.sparse.load_npz(fpath))
        if matrices:
            merged = scipy.sparse.vstack(matrices)
            scipy.sparse.save_npz(f'{canonical_dir}/{worker_id}_ANT.npz',merged)

def EM_algo_hybrid_multi(isoform_len_dict,isoform_gene_dict,gene_isoforms_dict,sr_sam_list,lr_align_dirs,output_path,threads,EM_choice,lr_weights=None,sr_weights=None):
    """多平台hybrid EM：多个LR+SR数据集，支持per-sample加权。"""
    from EM_hybrid.EM_LR import prepare_LR,merge_cond_prob_files
    # 准备isoform索引（与EM_algo_hybrid相同）
    isoform_len_df = pd.Series(isoform_len_dict)
    isoform_list = sorted(isoform_len_dict.keys())
    isoform_index_dict = {}
    isoform_len_arr = []
    for i,isoform in enumerate(isoform_list):
        isoform_index_dict[isoform] = i
        isoform_len_arr.append(isoform_len_dict[isoform])
    isoform_index_series = pd.Series(isoform_index_dict)
    isoform_index_series.name = 'Index'
    isoform_index_series.index.name = 'Isoform'
    output_df = isoform_index_series.to_frame().reset_index().sort_values(by="Index").set_index('Index')
    isoform_gene_df = pd.Series(isoform_gene_dict)
    isoform_gene_df.name = 'Gene'
    isoform_gene_df.index.name = 'Isoform'
    output_df = output_df.join(isoform_gene_df,on='Isoform')
    isoform_len_arr = np.array(isoform_len_arr)
    eff_len_arr = isoform_len_arr.copy()
    gene_isoform_index = {}
    for rname in gene_isoforms_dict:
        for gname in gene_isoforms_dict[rname]:
            gene_isoform_index[gname] = []
            for isoform in gene_isoforms_dict[rname][gname]:
                gene_isoform_index[gname].append(isoform_index_dict[isoform])
    # 逐个处理SR数据集，同时收集质量分用于自适应权重
    hits_dirs = []
    sr_quality_scores = []
    for i,sr_sam in enumerate(sr_sam_list):
        hits_dir = f'{output_path}/temp/hits_dict_sr{i}'
        print(f'[INFO] Start preparing short reads data {i}... at {datetime.datetime.now()}',flush=True)
        theta_SR_i,eff_len_i,_ = prepare_hits(sr_sam,output_path,isoform_index_dict,gene_isoform_index,threads,hits_dir=hits_dir)
        if i == 0:
            eff_len_arr = eff_len_i  # 使用第一个SR的有效长度
        hits_dirs.append(hits_dir)
        # Rename SR feature dicts with sample index so extract_features accumulates all samples
        for worker_id in range(threads):
            src = f'{output_path}/temp/machine_learning/SR_feature_dict_{worker_id}'
            if os.path.exists(src):
                os.rename(src, f'{output_path}/temp/machine_learning/SR_feature_dict_sr{i}_{worker_id}')
        if sr_weights is None:
            q = _compute_sr_quality(hits_dir, threads)
            sr_quality_scores.append(q)
            print(f'[INFO] SR sample {i} unique mapping rate (quality score): {q:.4f}',flush=True)
    output_df['Effective length'] = eff_len_arr
    print(f'[INFO] Prepare short reads data done at {datetime.datetime.now()}',flush=True)
    # 逐个处理LR数据集，同时收集质量分
    cond_prob_dirs = []
    lr_quality_scores = []
    for i,lr_align_dir in enumerate(lr_align_dirs):
        cond_prob_dir = f'{output_path}/temp/cond_prob_lr{i}'
        print(f'[INFO] Start preparing long reads data {i}...',flush=True)
        prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path,lr_align_dir=lr_align_dir,cond_prob_dir=cond_prob_dir)
        cond_prob_dirs.append(cond_prob_dir)
        # Rename LR feature dicts with sample index so extract_features accumulates all samples
        for fpath in glob.glob(f'{output_path}/temp/machine_learning/LR_feature_dict_*'):
            fname = os.path.basename(fpath)
            if not fname.startswith('LR_feature_dict_lr'):
                new_fname = fname.replace('LR_feature_dict_', f'LR_feature_dict_lr{i}_')
                os.rename(fpath, f'{output_path}/temp/machine_learning/{new_fname}')
        if lr_weights is None:
            q = _compute_lr_quality(cond_prob_dir, threads)
            lr_quality_scores.append(q)
            print(f'[INFO] LR sample {i} unique mapping rate (quality score): {q:.4f}',flush=True)
    print(f'[INFO] Prepare long reads data done at {datetime.datetime.now()}',flush=True)
    # 合并用于community结构检测
    merge_cond_prob_files(cond_prob_dirs,output_path,threads)
    merge_ant_files(hits_dirs,output_path,threads)
    # 设置per-sample权重：所有样本权重（LR+SR）统一归一化，总和=1
    n_lr = len(lr_align_dirs)
    n_sr = len(sr_sam_list) if sr_sam_list else 0
    if sr_sam_list:
        if lr_weights is None and sr_weights is None:
            # 自适应模式：按唯一比对率计算各组内权重（各组内和=1）
            lr_q = lr_quality_scores
            sr_q = sr_quality_scores
        else:
            lr_q = lr_weights if lr_weights is not None else lr_quality_scores
            sr_q = sr_weights if sr_weights is not None else sr_quality_scores
        all_q = list(lr_q) + list(sr_q)
        total = sum(all_q)
        g = [q / total for q in all_q] if total > 0 else [1.0 / (n_lr + n_sr)] * (n_lr + n_sr)
        config.lr_global_weights = g[:n_lr]
        config.sr_global_weights = g[n_lr:]
        print(f'[INFO] Global weights (all sum=1): LR={[round(w,4) for w in config.lr_global_weights]}, SR={[round(w,4) for w in config.sr_global_weights]}')
        # 使用预训练模型为每个 gene community 预测最优 alpha（与原始单文件方式一致）
        config.alpha = 'adaptive'
        config.alpha_df_path = f'{output_path}/adaptive_alpha_df.tsv'
    else:
        # 纯LR模式：只有LR，alpha=1，不需要预测
        lr_q = lr_weights if lr_weights is not None else lr_quality_scores
        total_w = sum(lr_q)
        config.lr_global_weights = [w / total_w for w in lr_q] if total_w > 0 else [1.0 / n_lr] * n_lr
        config.sr_global_weights = None
        config.alpha = 1.0
        config.alpha_df_path = None
        print(f'[INFO] LR-only mode. LR global weights (sum=1): {[round(w,4) for w in config.lr_global_weights]}')
    # 构建community并运行EM
    print('[INFO] Start constructing the community...',flush=True)
    st = time.time()
    num_SRs, num_LRs = construct_community(isoform_gene_dict, isoform_index_dict, output_path, threads,
                                            sr_dirs=hits_dirs if sr_sam_list else None,
                                            lr_dirs=cond_prob_dirs)
    duration = (time.time() - st)
    print('[INFO] Done in {} seconds at {}!'.format(duration,str(datetime.datetime.now())),flush=True)
    if sr_sam_list and config.alpha == 'adaptive':
        model_files = glob.glob(config.pretrained_model_path + '/*.pkl') if config.pretrained_model_path else []
        if not model_files:
            print(f'[WARNING] No pretrained model files found in: {config.pretrained_model_path}')
            print('[WARNING] Falling back to fixed alpha=0.5')
            config.alpha = 0.5
            config.alpha_df_path = None
        elif not Path(config.alpha_df_path).exists():
            predict_alpha(output_path, num_SRs, num_LRs)
        else:
            print('[INFO] Using existing alpha from ' + str(config.alpha_df_path))
    print(f'[INFO] LR global weights = {config.lr_global_weights}')
    print(f'[INFO] SR global weights = {config.sr_global_weights}')
    print(f'[INFO] Alpha mode: {config.alpha}')
    EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads,num_SRs,num_LRs)

def EM_algo_hybrid(isoform_len_dict,isoform_gene_dict,gene_isoforms_dict,SR_sam,output_path,threads,EM_choice):
   # prepare arr
    isoform_len_df = pd.Series(isoform_len_dict)
    isoform_list = sorted(isoform_len_dict.keys())
    isoform_index_dict = {}
    isoform_len_arr = []
    for i,isoform in enumerate(isoform_list):
        isoform_index_dict[isoform] = i
        isoform_len_arr.append(isoform_len_dict[isoform])
    isoform_index_series = pd.Series(isoform_index_dict)
    isoform_index_series.name = 'Index'
    isoform_index_series.index.name = 'Isoform'
    output_df = isoform_index_series.to_frame().reset_index().sort_values(by="Index").set_index('Index')
    isoform_gene_df = pd.Series(isoform_gene_dict)
    isoform_gene_df.name = 'Gene'
    isoform_gene_df.index.name = 'Isoform'
    output_df = output_df.join(isoform_gene_df,on='Isoform')
    isoform_len_arr = np.array(isoform_len_arr)
    eff_len_arr = isoform_len_arr.copy()
    # prepare SR
    gene_isoform_index = {}
    for rname in gene_isoforms_dict:
        for gname in gene_isoforms_dict[rname]:
            gene_isoform_index[gname] = []
            for isoform in gene_isoforms_dict[rname][gname]:
                gene_isoform_index[gname].append(isoform_index_dict[isoform])
    print('[INFO] Start preparing short reads data... at {}'.format(str(datetime.datetime.now())),flush=True)
    theta_SR_arr,eff_len_arr,SR_num_batches_dict = prepare_hits(SR_sam,output_path,isoform_index_dict,gene_isoform_index,threads)
    output_df['Effective length'] = eff_len_arr
    print('[INFO] Prepare short reads data done at {}'.format(str(datetime.datetime.now())),flush=True)
    print('[INFO] Start preparing long reads data...',flush=True)
    theta_LR_arr,_,LR_num_batches_dict = prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path)
    num_SRs = theta_SR_arr.sum()
    num_LRs = theta_LR_arr.sum()
    print('[INFO] Prepare long reads data done at {}'.format(str(datetime.datetime.now())),flush=True)
    # print(f'Number of SRs/eff_len:{num_SRs}')
    # print(f'Number of LRs:{num_LRs}')
    # print(f'Pseudo_count_SR:'+str(config.pseudo_count_SR))
    # print(f'Pseudo_count_LR:'+str(config.pseudo_count_LR),flush=True)
    # write_result_to_tsv(f'{output_path}/SR_count.tsv',output_df,theta_SR_arr.flatten())
    # write_result_to_tsv(f'{output_path}/LR_count.tsv',output_df,theta_LR_arr.flatten())
    # np.savez_compressed(f'{output_path}/initial_theta',theta=theta_arr)
    # Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    # print('Using {} as initial theta'.format(config.inital_theta))
    
    # theta_arr = theta_arr/theta_arr.sum()
    print('[INFO] Start constructing the community...',flush=True)
    st = time.time()
    num_SRs,num_LRs = construct_community(isoform_gene_dict,isoform_index_dict,output_path,threads)
    duration = (time.time() - st)
    print('[INFO] Done in {} seconds at {}!'.format(duration,str(datetime.datetime.now())),flush=True)
    print('[INFO] Extract features and predict best alpha...',flush=True)
    if config.alpha == 'adaptive':
        if not Path(config.alpha_df_path).exists():
            predict_alpha(output_path,num_SRs,num_LRs)
        else:
            print('[INFO] Using alpha from '+str(config.alpha_df_path))
    else:
        config.alpha_df_path = None
        print('[INFO] Using fixed alpha = '+str(config.alpha))
    EM_manager(isoform_gene_dict,isoform_index_dict,eff_len_arr,output_df,output_path,threads,num_SRs,num_LRs)
    


