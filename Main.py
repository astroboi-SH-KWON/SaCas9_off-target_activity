import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform
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

if SYSTEM_NM == 'Linux':
    # REAL
    pass
else:
    # DEV
    WORK_DIR = "D:/000_WORK/SeoSangYeon/20201126_off_target_finder__SaCas9_off-target_activity/WORK_DIR/"

IN = 'input/'
OU = 'output/'

FLTD_CDS_INFO = 'all_ccds_filtered_201130_CCDS_mouse_current.txt'  # mouse
# FLTD_CDS_INFO = 'all_ccds_filtered_201130_CCDS_human_current.txt'  # human

FILE_NM_FORM = 'cleavage_pos_in_shortest_cds_w_whole_mouse_gene_SaCas9.txt_*_result.txt'  # mouse
# FILE_NM_FORM = 'cleavage_pos_in_shortest_cds_w_whole_human_gene_SaCas9.txt_*_result.txt'  # human

LEN_SEQ = 21

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


def analize_individual(off_trgt, result_dict):
    print('st : analize_individual\n', off_trgt)
    logic = Logic.Logics()

    with open(off_trgt) as f:
        while True:
            tmp_line = f.readline()
            if tmp_line == '': break
            tmp_arr = tmp_line.split('\t')
            g_seq = tmp_arr[0][:LEN_SEQ].upper()
            t_seq = tmp_arr[3][:LEN_SEQ].upper()
            num_mismatch = int(tmp_arr[-1])
            guide_seq, _, pred_activity, _ = predict_activity_single(g_seq, t_seq)

            if guide_seq in result_dict:
                result_dict[guide_seq][0][num_mismatch] += 1
                result_dict[guide_seq][1] += pred_activity
            else:
                # e.g. guide RNA sequence, # of 0-bp mismatched targets, ...,  # of 6-bp mismatched targets
                result_dict.update({guide_seq: [[0, 0, 0, 0, 0, 0, 0], pred_activity]})
                result_dict[guide_seq][0][num_mismatch] += 1

    result_fn, _ = os.path.splitext(off_trgt)
    with open(result_fn.replace(IN[:-1], OU[:-1]) + '_mis_match', 'w') as mis_mat_ou_f:
        with open(result_fn.replace(IN[:-1], OU[:-1]) + '_guide_score', 'w') as scores_ou_f:
            for guide_seq, val_arr in result_dict.items():

                # make file for e.g. guide seq, # of 0-bp mismatched targets, ...,  # of 6-bp mismatched targets
                tmp_row = ''
                tmp_row += guide_seq + '\t'
                for num_mismat in val_arr[0]:
                    tmp_row += str(num_mismat) + '\t'
                mis_mat_ou_f.write(tmp_row[:-1] + '\n')

                # guide score
                num_0_mis_mat = val_arr[0][0]
                if num_0_mis_mat == 1:
                    sum_pred_act = val_arr[1]
                    sguide_score = logic.cal_sguide_score(sum_pred_act)
                    scores_ou_f.write(guide_seq + '\t' + str(sguide_score) + '\n')

    print('en : analize_individual\n', off_trgt)


def multi_processing_individual():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + IN + FILE_NM_FORM)

    for off_trgt in sources:
        proc = mp.Process(target=analize_individual, args=(off_trgt, {}))
        proc.start()


def multi_processing_all_in_one():
    util = Util.Utils()

    sources = util.get_files_from_dir(WORK_DIR + IN + FILE_NM_FORM)

    manager = mp.Manager()
    result_dict = manager.dict()
    jobs = []
    for off_trgt in sources:
        proc = mp.Process(target=analize_individual, args=(off_trgt, result_dict))
        jobs.append(proc)
        proc.start()

    for ret_proc in jobs:
        ret_proc.join()


def test():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()
    cds_info = util.read_tsv_ignore_N_line(WORK_DIR + IN + FLTD_CDS_INFO)

    """
    3-3-1) Tier: off-target position이 CDS (Tier I)인지 non-CDS (Tier II)인지 분류
    seq_idx list by key as chromosome
    """
    trgt_pos = 97104886
    trgt_pos += 15
    cds_idx_list = logic_prep.get_cds_idx_list_by_chr(cds_info)['chr1']
    print(not logic.check_seq_in_cds(cds_idx_list, trgt_pos))


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # multi_processing_individual()
    test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))