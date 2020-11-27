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
# import Logic
# import LogicPrep
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

    print(guide_seq, target_seq, pred_activity, num_mm, sep='\t')
    # return guide_seq, target_seq, pred_activity, num_mm


def test():
    util = Util.Utils()
    df = util.read_excel_to_df(WORK_DIR + IN + '201127_For KSH.xlsx', header=None)

    len_df = len(df[df.columns[0]])

    for i in range(len_df):
        g_seq = df.loc[i][0][:LEN_SEQ].upper()
        t_seq = df.loc[i][3][:LEN_SEQ].upper()
        predict_activity_single(g_seq, t_seq)



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))