import time
import os
import multiprocessing as mp
import platform
import operator
import numpy as np
########################################
# from __future__ import print_function
# from sys import argv
import pandas as pd
from collections import defaultdict
########################################

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

# type
# TYPE = 'mouse'
TYPE = 'human'

if SYSTEM_NM == 'Linux':
    # REAL
    pass
else:
    # DEV
    WORK_DIR = "D:/000_WORK/SeoSangYeon/20201126_off_target_finder__SaCas9_off-target_activity/WORK_DIR/"

IN = 'input/'
OU = 'output/'

CAS_OFF_RESULT = 'FirstResult/'

FLTD_CDS_INFO = 'all_ccds_filtered_201130_CCDS_' + TYPE + '_current.txt'
OFF_TRGT_FORM = 'cleavage_pos_in_shortest_cds_w_whole_' + TYPE + '_gene_SaCas9.txt_*_result.txt'

LEN_SEQ = 21

BIN_LABEL = [  # [1.1, 1.0 ... range of Soff-target]
            [1.1, 1.0]
            , [1.0, 0.2]
            , [0.2, 0.05]
            , [0.05, 0.0]
            ]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################


def predict_activity_single(guide_seq, target_seq):
    model_dir = 'pairwise-library-screen-master/models/'

    num_mismatches = range(20, 25)

    # Check input
    guide_seq = guide_seq.replace('U', 'T')

    if len(target_seq) > len(guide_seq):
        target_seq = target_seq[0:len(guide_seq)]  # This deals with target sequences that include the PAM

    pair_key = (guide_seq, target_seq)

    if len(guide_seq) not in num_mismatches:
        raise Exception('Model not found for this guide length.')

    if '-' in guide_seq or '-' in target_seq:
        raise Exception('This model does not support gaps/bulges.')

    # Define mismatches

    base_pairings = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
    mismatch_ids = ['dT:rC', 'dT:rG', 'dT:rU', 'dG:rA', 'dG:rG', 'dG:rU', 'dC:rA', 'dC:rC', 'dC:rU', 'dA:rA', 'dA:rC',
                    'dA:rG']

    mm_dict = {_: __ for _, __ in zip(base_pairings, mismatch_ids)}

    # Read models

    par_dict = defaultdict(dict)
    for guide_len in num_mismatches:

        model_file = model_dir + 'rstan_mult_%s_coeffs.txt' % guide_len
        model_df = pd.read_table(model_file)

        par_dict[guide_len] = {}

        for i in range(1, guide_len + 1):
            par_dict[guide_len][i] = defaultdict(dict)

        for row in model_df.iterrows():
            row_values = row[1]
            par_dict[guide_len][guide_len - row_values['Position'] + 1][row_values['Mismatch']] = row_values['Median']

    # Add up penalties

    pred_activity = 1
    num_mm = 0

    for i, [dna, rna] in enumerate(zip(target_seq, guide_seq)):
        mm_pos = i + 1

        if mm_pos == 1:  # Mismatches at at the 5' position are not truly determined by the model
            continue
        if dna != rna:
            if dna + rna in mm_dict:
                mm_coeff = par_dict[guide_len][mm_pos][mm_dict[dna + rna]]
                num_mm += 1
                pred_activity *= mm_coeff

    return guide_seq, target_seq, pred_activity, num_mm


def process(path_arr):
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    cds_info = util.read_tsv_ignore_N_line(WORK_DIR + IN + FLTD_CDS_INFO)

    """
    3-3-1) Tier: off-target position이 CDS (Tier I)인지 non-CDS (Tier II)인지 분류
    seq_idx list by key as chromosome
    """
    cds_dict = logic_prep.get_cds_idx_list_by_chr(cds_info)
    cds_info.clear()

    result_dict = {}
    for path in path_arr:
        with open(path) as f:
            while True:
                tmp_line = f.readline().replace('\n', '')
                if tmp_line == '':
                    """
                    GGTGTAGAGCTGGAACATGAANNNNNN	X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF	67064737	GGTaTAGtGaTGGAAgATGAAATGAAT	-	4
                    GAGGTTGAGCAGGGCGCGCTGNNNNNN	X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF	16787030	GcaGgTGcGCAGGGCGCGCTGCGGGGT	+	4
                                there is empty line like this
                    TCTCCTCTGCCTCTTCTTTCGNNNNNN	X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF	39576853	TCTCCTCTGCCTCTcCTTTgaGGGGGT	+	3
                    TCTCCTCTGCCTCTTCTTTCGNNNNNN	X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF	49292670	TgTCtcCTGCCTCcTCTTTCGAGGGAT	-	4
                    """
                    # check if there is another line
                    tmp_line = f.readline().replace('\n', '')
                    if tmp_line == '':
                        break

                tmp_arr = tmp_line.split('\t')
                g_seq = tmp_arr[0][:LEN_SEQ].upper()
                tmp_chr = 'chr' + tmp_arr[1].split(' ')[0]
                tmp_pos = int(tmp_arr[2])
                t_seq = tmp_arr[3][:LEN_SEQ].upper()
                tmp_strnd = tmp_arr[4].replace(' ', '')
                if tmp_strnd == '+':
                    tmp_pos += 19
                else:
                    tmp_pos += 10
                num_mismatch = int(tmp_arr[-1])
                guide_seq, _, scr_off_trgt, _ = predict_activity_single(g_seq, t_seq)

                # check bin I (in cds) or bin II(not in any cds)
                in_cds_flag = logic.check_seq_in_cds(cds_dict[tmp_chr], tmp_pos)
                # check Soff-trgt Tier [1, 0.2, 0.05, 0]
                idx_bin_label = logic.check_which_bin_label(BIN_LABEL, scr_off_trgt)

                idx_bin_label *= 2
                if not in_cds_flag:
                    idx_bin_label += 1

                # TODO debugging
                if guide_seq == 'GGTACACCTTCGCTGGTTCAC':
                    print(in_cds_flag, ': in_cds_flag')
                    print(idx_bin_label, ': idx_bin_label')
                    print(tmp_chr, ': tmp_chr')
                    print(tmp_pos, ': tmp_pos')
                    print()
                    with open(WORK_DIR + OU + guide_seq + '.txt' + str(tmp_pos), 'a') as f:
                        f.write(str(cds_dict[tmp_chr]) + '\n')

                if guide_seq in result_dict:
                    result_dict[guide_seq][idx_bin_label] += 1
                    result_dict[guide_seq][-1] += scr_off_trgt
                else:
                    init_arr = []
                    for i in range(len(BIN_LABEL) * 2 + 1):
                        init_arr.append(0)
                    result_dict.update({guide_seq: init_arr})
                    result_dict[guide_seq][idx_bin_label] += 1
                    result_dict[guide_seq][-1] += scr_off_trgt
        cds_dict.clear()
        print('DONE :', str(path))
    return result_dict


def multi_processing():
    util = Util.Utils()
    logic = Logic.Logics()

    sources = util.get_files_from_dir(WORK_DIR + CAS_OFF_RESULT + OFF_TRGT_FORM)

    cnt_cpu = len(sources)
    if MULTI_CNT < cnt_cpu:
        cnt_cpu = MULTI_CNT
    splited_sources = np.array_split(sources, cnt_cpu)

    pool = mp.Pool(processes=cnt_cpu)
    pool_list = pool.map(process, splited_sources)

    # merge dicts
    tot_dict = {}
    for tmp_dict in pool_list:
        for guide_key, val_arr in tmp_dict.items():
            if guide_key in tot_dict:
                for i in range(len(val_arr)):
                    tot_dict[guide_key][i] += val_arr[i]
            else:
                tot_dict.update({guide_key: val_arr})
    pool_list.clear()

    # dict_to_list
    dict_to_list = []
    for guide_key, val_arr in tot_dict.items():
        tmp_arr = [guide_key]
        tmp_arr += val_arr
        dict_to_list.append(tmp_arr)

    dict_to_list.sort(key=operator.itemgetter(1, 2, 3, 4, 5, 6, 7, 8))

    result_list = []
    for i in range(len(dict_to_list)):
        tmp_arr = []
        for val in dict_to_list[i][:-1]:
            tmp_arr.append(val)
        tmp_arr.append(i + 1)
        tmp_arr.append(logic.cal_sguide_score(dict_to_list[i][-1]))
        result_list.append(tmp_arr)

    header = ['guide_seq', 'Tier_I_bin_I', 'Tier_II_bin_I', 'Tier_I_bin_II', 'Tier_II_bin_II', 'Tier_I_bin_III',
              'Tier_II_bin_III', 'Tier_I_bin_IV', 'Tier_II_bin_IV', 'off_trgt_rank', 'Sguide']
    util.make_tsv(WORK_DIR + OU + 'result_' + FLTD_CDS_INFO, header, result_list)
    util.make_excel(WORK_DIR + OU + 'result_' + FLTD_CDS_INFO.replace('.txt', ''), header, result_list)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    multi_processing()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))