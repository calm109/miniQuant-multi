from parse_alignment import map_read,parse_read_line,map_read_to_region
from patch_mp import patch_mp_connection_bpo_17560
from parse_annotation_main import check_valid_region
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import numpy as np
import time
import random
import multiprocessing as mp
import os
import bisect
import re as _re
from util import check_region_type,sync_reference_name
import config

# from memory_profiler import profile
# def parse_alignment_iteration(alignment_file_path,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
#                     start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#                     READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,map_f,line_nums):
def debuginfoStr(info):
    print(info,flush=True)
    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]
    mem = int(memusage.strip())/1024
    print('Mem consumption: '+str(mem),flush=True)
def parse_alignment_iteration(alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,long_read,aln_line_marker):
    os.nice(10)
    start_file_pos,num_lines = aln_line_marker
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_count = {}
        local_gene_regions_read_length = {}
        local_gene_regions_read_pos = {}
        aln_file.seek(start_file_pos)
        line_num_ct = 0
        max_buffer_size = 1e1
        buffer_size = 0
        for line in aln_file: # 这里改成逐行读取sam文件，之前是一次性读取整个sam文件到内存中，现在改成逐行读取，减少内存占用
            line_num_ct += 1
            if line_num_ct > num_lines:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'): # 如果该行的RNAME字段为*，说明该read没有比对上任何位置，直接跳过
                    continue
                aln_line = parse_read_line(line) # 解析sam文件中的一行，提取出该read的比对信息，包括比对位置、比对质量、比对的基因和区间等
                mapping = map_f(points_dict,interval_tree_dict, filtered_gene_regions_dict, # 根据该read的比对信息，判断该read是否比对到已知的基因-区间上，如果比对上了，记录该read比对到的基因-区间的信息，包括染色体、基因名称、区间名称等
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
                if (mapping['read_mapped']):
                    random.seed(mapping['read_name'])
                    for mapping_area in [random.choice(mapping['mapping_area'])]: # 如果该read比对上了多个基因-区间，则随机选择一个基因-区间作为该read的比对区域，之前是把所有比对区域都考虑进来，现在改成随机选择一个，减少后续统计的复杂度
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name'] # 这里提取出该read比对到的基因-区间的信息，包括染色体、基因名称、区间名称等
                        if rname not in local_gene_regions_read_count:
                            local_gene_regions_read_count[rname],local_gene_regions_read_length[rname] = {},{}
                            local_gene_regions_read_pos[rname] = {}
                        if gname not in local_gene_regions_read_count[rname]:
                            local_gene_regions_read_count[rname][gname],local_gene_regions_read_length[rname][gname] = {},{}
                            local_gene_regions_read_pos[rname][gname] = {}
                        if region_name not in local_gene_regions_read_count[rname][gname]:
                            local_gene_regions_read_count[rname][gname][region_name],local_gene_regions_read_length[rname][gname][region_name] = 0,[]
                            local_gene_regions_read_pos[rname][gname][region_name] = []
                        local_gene_regions_read_count[rname][gname][region_name] += 1 
                        # if long_read:
                        local_gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
                        local_gene_regions_read_pos[rname][gname][region_name].append(mapping)
                    buffer_size += 1
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
            if buffer_size > max_buffer_size:
                temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos))
                local_gene_regions_read_count,local_gene_regions_read_length = {},{}
                local_gene_regions_read_pos = {}
                buffer_size = 0
        if buffer_size > 0:
            temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos))
    return # 这里改成多进程版本，之前是单线程版本直接返回统计结果，现在改成多进程版本，使用 temp_queue 将统计结果传回主进程，由 mapping_listener 函数进行汇总统计
def mapping_listener(temp_queue,gene_regions_read_count,gene_regions_read_length,gene_regions_read_pos):
    num_mapped_lines = 0
    num_lines = 0
    while True:
        msg = temp_queue.get()
        if msg == 'kill':
            break
        else:
            local_gene_regions_read_count,local_gene_regions_read_length,local_gene_regions_read_pos = msg
            for rname in local_gene_regions_read_count:
                for gname in local_gene_regions_read_count[rname]:
                    for region in local_gene_regions_read_count[rname][gname]:
                        num_mapped_lines += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_count[rname][gname][region] += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_length[rname][gname][region] += local_gene_regions_read_length[rname][gname][region]
                        gene_regions_read_pos[rname][gname][region] += local_gene_regions_read_pos[rname][gname][region]

            # for mapping in local_all_mappings:
            #     num_lines += 1
            #     if len(mapping['gene_candidates'])>0:
            #         num_mapped_to_gene += 1
            #     if (mapping['read_mapped']):
            #         num_mapped_lines += 1
            #         for mapping_area in [random.choice(mapping['mapping_area'])]:
            #             rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
            #             if region_name in gene_regions_read_count[rname][gname]:
            #                 gene_regions_read_count[rname][gname][region_name] += 1 
            #                 gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
            #         read_lens.append(mapping['read_length'])
            #         read_names.update(local_read_names)
    return gene_regions_read_count,gene_regions_read_length,num_mapped_lines,gene_regions_read_pos

# @profile
def get_aln_line_marker(alignment_file_path,threads):
    with open(alignment_file_path, 'r') as aln_file:
        line_offset = []
        offset = 0
        for line in aln_file:
            if line[0] != '@':
                line_offset.append(offset)
            offset += len(line)
    num_aln_lines = len(line_offset)
    if num_aln_lines == 0:
        return []
    actual_threads = min(threads, num_aln_lines)
    chunksize, extra = divmod(num_aln_lines, actual_threads)
    if extra:
        chunksize += 1
    aln_line_marker = []
    for i in range(actual_threads):
        aln_line_marker.append((line_offset[i*chunksize],chunksize))
    return aln_line_marker

def parse_alignment(alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,long_read,filtering,threads):
    patch_mp_connection_bpo_17560()
    start_t = time.time()
    manager = mp.Manager()
    gene_regions_read_count = {}
    gene_regions_read_length ={}
    gene_regions_read_pos = {}
    global filtered_gene_regions_dict
    # 初始化基因-区间的读数统计字典，结构为 gene_regions_read_count[chr][gene][region] = count，同时初始化一个过滤后的基因-区间字典 filtered_gene_regions_dict 用于后续快速判断某个区间是否需要考虑
    filtered_gene_regions_dict = defaultdict(lambda: defaultdict(dict)) 
    for rname in gene_regions_dict:
        gene_regions_read_pos[rname] = {}
        gene_regions_read_count[rname],gene_regions_read_length[rname] = {},{}
        for gname in gene_regions_dict[rname]:
            gene_regions_read_count[rname][gname],gene_regions_read_length[rname][gname] = {},{}
            gene_regions_read_pos[rname][gname] = {}
            # if (not long_read):
            #     per_gene_regions_dict = filter_regions(gene_regions_dict[rname][gname],long_read=False)
            # else:
            #     per_gene_regions_dict =  filter_regions(gene_regions_dict[rname][gname],long_read=True)
            per_gene_regions_dict =  gene_regions_dict[rname][gname]
            # 这里提前把所有已知区域都初始化为0，这样即使没有read比对到某区域，该区域仍然存在于字典中（对SR模式的矩阵构建有意义）。
            for region in per_gene_regions_dict:
                gene_regions_read_count[rname][gname][region] = 0
                gene_regions_read_length[rname][gname][region] = []
                gene_regions_read_pos[rname][gname][region] = []
                filtered_gene_regions_dict[rname][gname][region] = True
    if alignment_file_path == None:
        return gene_regions_read_count,150,0
    # Create sorted end and start positions
    global start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST
    start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST = dict(),dict(),dict(),dict(),list(gene_range.keys())
    CHR_LIST = list(gene_range.keys())
    # 对于每条染色体，分别对基因的起始位置和结束位置进行排序，并记录排序后的位置信息和对应的基因名称，存储在 start_pos_list、start_gname_list、end_pos_list 和 end_gname_list 中，以便后续快速定位比对区域所属的基因
    for rname in CHR_LIST:     
        # Sort based on start position
        temp_list = sorted(gene_range[rname], key=itemgetter(1))
        start_pos_list[rname] = [temp_list[j][1] for j in range(len(temp_list))]
        start_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
        # Sort based on end position
        temp_list = sorted(gene_range[rname], key=itemgetter(2))
        end_pos_list[rname] = [temp_list[j][2] for j in range(len(temp_list))]
        end_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
    global points_dict,interval_tree_dict
    points_dict,interval_tree_dict = gene_points_dict,gene_interval_tree_dict
    map_f = map_read
    debuginfoStr('Before MP mem usage')
    pool = mp.Pool(threads+1)
    # partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path)
    temp_queue = manager.Queue()    
    watcher = pool.apply_async(mapping_listener, args=(temp_queue,gene_regions_read_count,gene_regions_read_length,gene_regions_read_pos))
    # parse_alignment_iteration函数读取并处理sam文件
    partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,long_read) # 这里改成多进程版本，之前是单线程版本直接调用 parse_alignment_iteration 函数，现在改成多进程版本，使用 pool.apply_async 来异步执行 parse_alignment_iteration 函数，并将结果通过 temp_queue 传回主进程，由 mapping_listener 函数进行汇总统计
    futures = [] 
    aln_line_marker = get_aln_line_marker(alignment_file_path,threads)
    for marker in aln_line_marker:
        futures.append(pool.apply_async(partial_read_alignment,(marker,)))
    # with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #     futures = [executor.submit(partial_read_alignment,marker) for marker in aln_line_marker]
    #     concurrent.futures.wait(futures)
    for future in futures:
        future.get()
    temp_queue.put('kill')
    gene_regions_read_count,gene_regions_read_length,num_mapped_lines,gene_regions_read_pos = watcher.get()
    pool.close()
    pool.join()
    read_lens = []
    for rname in gene_regions_read_length:
        for gname in gene_regions_read_length[rname]:
            for region in gene_regions_read_length[rname][gname]:
                read_lens += gene_regions_read_length[rname][gname][region]
        
    if (long_read):
        num_long_reads = 0
        gene_full_length_region_dict = defaultdict(lambda:{})
        isoform_max_num_exons_dict = {}
        for rname in gene_regions_dict:
            for gname in gene_regions_dict[rname]:
                regions_set = set()
                isoform_region_dict = defaultdict(lambda:set())
                for region in gene_regions_dict[rname][gname]:
                    for isoform in gene_regions_dict[rname][gname][region]:
                        isoform_region_dict[isoform].add(region)
                for isoform in isoform_region_dict:
                    max_region_exon_num = 0
                    longest_region = ''
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if max_region_exon_num < region_exon_num:
                            max_region_exon_num = region_exon_num
                            longest_region = region
                    isoform_max_num_exons_dict[isoform] = max_region_exon_num
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if region_exon_num == max_region_exon_num:
                            regions_set.add(region)
                # 对于每个基因，记录其所有全长区域（即外显子数最多的那些区域），后续在过滤比对到该基因的 read 时，如果该 read 比对到的区域不是全长区域且该 read 的长度占对应 isoform 长度的比例不高，则认为该 read 可能是部分比对到该区域的，过滤掉该 read
                gene_full_length_region_dict[rname][gname] = regions_set
        filtered_gene_regions_read_length = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
        for rname in gene_regions_read_length.copy():
            for gname in gene_regions_read_length[rname].copy():    
                # region_lens = []
                for region in gene_regions_read_length[rname][gname].copy():
                    region_exon_num = region.count(':')
                    read_length_list = []
                    for read_length in gene_regions_read_length[rname][gname][region]:
                        is_valid_read = True
                        if (filtering):
                            for isoform in gene_regions_dict[rname][gname][region]:
                                if (isoform_max_num_exons_dict[isoform] - region_exon_num  > 7) and (read_length / gene_isoforms_length_dict[rname][gname][isoform] <= 0.2):
                                        # 如果该区域对应的 isoform 的最大外显子数比该区域的外显子数多超过7，并且该 read 的长度占该 isoform 长度的比例不超过20%，则认为该 read 可能是部分比对到该区域的，过滤掉该 read
                                        is_valid_read = False 
                                        filtered_gene_regions_read_length[rname][gname][region].append(read_length)
                                        break
                        if (is_valid_read):
                            read_length_list.append(read_length)
                    gene_regions_read_length[rname][gname][region] = read_length_list
                    gene_regions_read_count[rname][gname][region] = len(read_length_list)
                    num_long_reads += len(read_length_list)
                    # 保留全部零读段区域（包括非full-length区域），使可识别性矩阵A与SR对称（纯注释驱动）
                    # EM定量路径中零读段区域会被 add_full_length_region 逻辑的 continue 跳过，不影响结果
        return gene_regions_read_count,gene_regions_read_length,sum(read_lens),num_long_reads,filtered_gene_regions_read_length,gene_regions_read_pos
    # 对于 SR 数据，统计所有比对到已知区域的 read 的长度，并计算平均 read 长度，后续在构建基因-区间矩阵时可以根据该平均 read 长度来判断某个区间是否有足够的覆盖度（即是否有足够多的 read 比对到该区间且这些 read 的长度足够长），从而决定是否保留该区间
    else:
        all_read_len = []
        for rname in gene_regions_read_count.copy():
            for gname in gene_regions_read_count[rname].copy():    
                # region_lens = []
                for region in gene_regions_read_count[rname][gname].copy():
                    if gene_regions_read_count[rname][gname][region] == 0:
                        if config.sr_region_selection == 'real_data':
                            if config.keep_sr_exon_region == 'nonfullrank':
                                pass
                            elif config.keep_sr_exon_region == 'all':
                                if check_region_type(region) not in ['one_exon','two_exons','exons']:
                                    del gene_regions_read_count[rname][gname][region]
                            else:
                                del gene_regions_read_count[rname][gname][region]
                    else:
                        all_read_len += gene_regions_read_length[rname][gname][region]
        SR_read_len = sum(all_read_len) / len(all_read_len) if all_read_len else getattr(config, 'READ_LEN', 150)
        return gene_regions_read_count,SR_read_len,num_mapped_lines


# ---------------------------------------------------------------------------
# Transcriptome-aligned SR SAM support
# ---------------------------------------------------------------------------

def is_transcriptome_aligned_sam(sam_path):
    """Return True if @SQ SN headers look like transcript IDs rather than chromosomes."""
    chr_pat = _re.compile(r'^(chr[^_\s]+|\d{1,2}|X|Y|MT|M|W|Z)$', _re.IGNORECASE)
    with open(sam_path) as f:
        for line in f:
            if not line.startswith('@SQ'):
                if not line.startswith('@'):
                    break
                continue
            for field in line.split('\t'):
                if field.startswith('SN:'):
                    sn = field[3:].strip()
                    return not chr_pat.match(sn)
    return False


def build_transcript_genome_map(ref_file_path):
    """Parse GTF to build {transcript_id: {'chr', 'strand', 'exons': [(tx_s, tx_e, g_s, g_e), ...]}}
    where exons are sorted in transcript order (tx_s ascending) and g_s/g_e are 1-based genomic coords.
    """
    tx_exons = defaultdict(list)
    with open(ref_file_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            chr_name = fields[0]
            strand = fields[6]
            try:
                g_start = int(fields[3])
                g_end = int(fields[4])
            except ValueError:
                continue
            m = _re.search('transcript_id "([^"]*)"', fields[8])
            if not m:
                continue
            tx_id = m.group(1)
            tx_exons[tx_id].append((g_start, g_end, chr_name, strand))

    tx_map = {}
    for tx_id, exons in tx_exons.items():
        chr_name = exons[0][2]
        strand = exons[0][3]
        # Sort exons in transcript order: ascending genomic for +, descending for -
        if strand == '+':
            sorted_exons = sorted(exons, key=lambda x: x[0])
        else:
            sorted_exons = sorted(exons, key=lambda x: -x[1])
        exon_list = []
        tx_offset = 0
        for (g_s, g_e, _, _) in sorted_exons:
            exon_len = g_e - g_s + 1
            exon_list.append((tx_offset, tx_offset + exon_len - 1, g_s, g_e))
            tx_offset += exon_len
        entry = {'chr': chr_name, 'strand': strand, 'exons': exon_list}
        tx_map[tx_id] = entry
        # Also index by unversioned ID (e.g. ENST00000456328.2 -> ENST00000456328)
        # so SAMs without version numbers can still be matched
        tx_id_base = tx_id.rsplit('.', 1)[0]
        if tx_id_base != tx_id:
            tx_map.setdefault(tx_id_base, entry)
    return tx_map


def _tx_to_genome_segments(tx_id, tx_pos_1based, ref_len, tx_map):
    """Convert a transcript-coordinate interval to genomic segments.
    tx_pos_1based: 1-based start on transcript; ref_len: number of reference bases consumed.
    Returns (chr_name, genome_start_1based, read_len_list) or None.
    read_len_list is [M, N, M, ...] alternating for map_read_to_region.
    """
    if tx_id not in tx_map:
        # Fallback: SAM has versioned ID but tx_map has unversioned (or vice versa)
        tx_id_base = tx_id.rsplit('.', 1)[0]
        if tx_id_base != tx_id and tx_id_base in tx_map:
            tx_id = tx_id_base
        else:
            return None
    info = tx_map[tx_id]
    strand = info['strand']
    exons = info['exons']  # (tx_s, tx_e, g_s, g_e) in transcript order

    tx_s0 = tx_pos_1based - 1          # 0-based inclusive start on transcript
    tx_e0 = tx_s0 + ref_len - 1        # 0-based inclusive end on transcript

    genomic_segs = []
    for (tx_exon_s, tx_exon_e, g_s, g_e) in exons:
        if tx_exon_s > tx_e0:
            break
        if tx_exon_e < tx_s0:
            continue
        ovlp_tx_s = max(tx_s0, tx_exon_s)
        ovlp_tx_e = min(tx_e0, tx_exon_e)
        off_s = ovlp_tx_s - tx_exon_s
        off_e = ovlp_tx_e - tx_exon_s
        if strand == '+':
            geno_s = g_s + off_s
            geno_e = g_s + off_e
        else:
            # transcript goes high→low in genome: tx_exon_s maps to g_e, tx_exon_e maps to g_s
            geno_s = g_e - off_e
            geno_e = g_e - off_s
        genomic_segs.append((geno_s, geno_e))

    if not genomic_segs:
        return None

    genomic_segs.sort(key=lambda x: x[0])

    genome_start = genomic_segs[0][0]
    read_len_list = []
    for i, (seg_s, seg_e) in enumerate(genomic_segs):
        read_len_list.append(seg_e - seg_s + 1)
        if i < len(genomic_segs) - 1:
            gap = genomic_segs[i + 1][0] - seg_e - 1
            read_len_list.append(gap)

    chr_name = info['chr']
    converted = sync_reference_name(chr_name)
    if converted.isnumeric():
        chr_name = converted
    return (chr_name, genome_start, read_len_list)


def parse_transcriptome_sr_alignment(alignment_file_path, tx_map, READ_JUNC_MIN_MAP_LEN,
                                      gene_points_dict, gene_range, gene_interval_tree_dict,
                                      gene_regions_dict, genes_regions_len_dict,
                                      gene_isoforms_length_dict):
    """Parse transcriptome-aligned SR SAM by converting reads to genomic coordinates.
    Returns (gene_regions_read_count, SR_read_len, num_mapped) — same format as parse_alignment for SR.
    Only primary alignments (no FLAG 256/2048) are counted.
    """
    # Initialize read count dicts with all regions = 0
    gene_regions_read_count = {}
    gene_regions_read_length = {}
    for rname in gene_regions_dict:
        gene_regions_read_count[rname] = {}
        gene_regions_read_length[rname] = {}
        for gname in gene_regions_dict[rname]:
            gene_regions_read_count[rname][gname] = {}
            gene_regions_read_length[rname][gname] = {}
            for region in gene_regions_dict[rname][gname]:
                gene_regions_read_count[rname][gname][region] = 0
                gene_regions_read_length[rname][gname][region] = []

    # Build sorted gene range lists for region lookup
    CHR_LIST = list(gene_range.keys())
    start_pos_list, end_pos_list = {}, {}
    start_gname_list, end_gname_list = {}, {}
    for rname in CHR_LIST:
        temp = sorted(gene_range[rname], key=lambda x: x[1])
        start_pos_list[rname] = [t[1] for t in temp]
        start_gname_list[rname] = [t[0] for t in temp]
        temp = sorted(gene_range[rname], key=lambda x: x[2])
        end_pos_list[rname] = [t[2] for t in temp]
        end_gname_list[rname] = [t[0] for t in temp]

    all_read_len = []
    num_mapped = 0
    cigar_op_pat = _re.compile(r'(\d+)([MIDNSHP=X])')

    with open(alignment_file_path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.split('\t')
            if len(fields) < 6 or fields[2] == '*':
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            if flag & 256 or flag & 2048:   # skip secondary / supplementary
                continue

            tx_id = fields[2]
            # Handle pipe-delimited reference names (e.g. ENST00000508313.3|ENSG...|...)
            if '|' in tx_id:
                tx_id = tx_id.split('|')[0]
            try:
                tx_pos = int(fields[3])     # 1-based position on transcript
            except ValueError:
                continue
            cigar_str = fields[5]
            if cigar_str == '*':
                continue

            # Count reference-consuming bases (M/=/ X/D) for transcript extent
            ref_len = sum(int(n) for n, op in cigar_op_pat.findall(cigar_str) if op in ('M', '=', 'X', 'D'))
            if ref_len == 0:
                continue
            read_m_len = sum(int(n) for n, op in cigar_op_pat.findall(cigar_str) if op in ('M', '=', 'X'))

            result = _tx_to_genome_segments(tx_id, tx_pos, ref_len, tx_map)
            if result is None:
                continue
            chr_name, genome_start, geno_read_len_list = result
            if len(geno_read_len_list) % 2 == 0:   # must start/end with M segment
                continue
            if chr_name not in CHR_LIST:
                continue

            # Compute genomic read end position
            read_end = genome_start - 1
            for seg in geno_read_len_list:
                read_end += seg

            tolerance = 20
            if read_end - genome_start > 2 * tolerance:
                s_idx = bisect.bisect_right(start_pos_list[chr_name], genome_start + tolerance)
                e_idx = bisect.bisect_left(end_pos_list[chr_name], read_end - tolerance)
            else:
                s_idx = bisect.bisect_right(start_pos_list[chr_name], genome_start + 2)
                e_idx = bisect.bisect_left(end_pos_list[chr_name], read_end - 2)

            gene_cands = (set(end_gname_list[chr_name][e_idx:]) &
                          set(start_gname_list[chr_name][:s_idx]))

            best_region, best_overlap, best_map_len, best_gene = '', 0, 1, None
            for gname in sorted(gene_cands):
                region, overlap, map_len = map_read_to_region(
                    genome_start, geno_read_len_list,
                    gene_points_dict[chr_name][gname],
                    gene_interval_tree_dict[chr_name][gname],
                    gene_regions_dict[chr_name][gname],
                    '', READ_JUNC_MIN_MAP_LEN)
                if region == '':
                    continue
                if overlap > best_overlap or (abs(overlap - best_overlap) <= 2 and map_len < best_map_len):
                    best_region, best_overlap, best_map_len, best_gene = region, overlap, map_len, gname

            if best_gene and best_region:
                gene_regions_read_count[chr_name][best_gene][best_region] += 1
                gene_regions_read_length[chr_name][best_gene][best_region].append(read_m_len)
                all_read_len.append(read_m_len)
                num_mapped += 1

    SR_read_len = sum(all_read_len) / len(all_read_len) if all_read_len else getattr(config, 'READ_LEN', 150)
    print(f'[INFO] Transcriptome SR: mapped {num_mapped} reads, avg read len {SR_read_len:.1f}', flush=True)
    return gene_regions_read_count, SR_read_len, num_mapped