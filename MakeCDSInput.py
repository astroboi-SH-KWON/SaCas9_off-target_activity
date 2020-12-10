import time
import os
import platform

import Util
import LogicPrep
#################### st env ####################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "/media/backup/ref/Ensemble_GRCm38_p6/"  # 20201130 mouse
    # REF_DIR = "/media/backup/ref/Ensemble_GRCh38_p13/"  # 20201130 human
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    WORK_DIR = 'D:/000_WORK/SeoSangYeon/20201126_off_target_finder__SaCas9_off-target_activity/WORK_DIR/'
PROJECT_NAME = WORK_DIR.split("/")[-2]

IN = 'input/'
OU = 'output/'

# NON_FLT_CDS_INFO = '201130_CCDS_mouse_current.txt'  # mouse
NON_FLT_CDS_INFO = '201130_CCDS_human_current.txt'  # human

# FLTD_CDS_INFO = 'filtered_201130_CCDS_mouse_current.txt'  # mouse
FLTD_CDS_INFO = 'filtered_201130_CCDS_human_current.txt'  # human

#################### en env ####################


def make_filtered_mouse_ccds_current_file():
    print('make_filtered_mouse_ccds_current_file')
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    ccds_list = []
    if SYSTEM_NM == 'Linux':
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))
    else:
        # ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0)[:3000])
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))

    # st plan A : filter out non Public, non Identical
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Public', 5)
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Identical', -1)
    # get the highest num of ccds_id in each gene
    ccds_list = logic_prep.get_highest_ccds_id_among_same_gen_id(ccds_list)
    # en plan A

    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list']
    util.make_tsv(WORK_DIR + IN + 'highest_ccds_' + FLTD_CDS_INFO, header, ccds_hg38_form_list)


def make_filtered_ccds_current_file_by_shortest_cdn():
    print('make_filtered_ccds_current_file_by_shortest_cdn')
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    ccds_list = []
    if SYSTEM_NM == 'Linux':
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))
    else:
        # ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0)[:3000])
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))

    # st plan A : filter out non Public, non Identical
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Public', 5)
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Identical', -1)

    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)

    filted_ccds_list = logic_prep.get_shortest_cdn_among_same_gen_id(ccds_hg38_form_list)  # 20201201
    # en plan A

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list']
    util.make_tsv(WORK_DIR + IN + 'shortest_cdn_' + FLTD_CDS_INFO, header, filted_ccds_list)


def make_filtered_ccds_current_file_by_all_ccds_id():
    print('make_filtered_ccds_current_file_by_all_ccds_id')
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    ccds_list = []
    if SYSTEM_NM == 'Linux':
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))
    else:
        # ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0)[:3000])
        ccds_list.extend(util.read_tsv_ignore_N_line(WORK_DIR + IN + NON_FLT_CDS_INFO, n_line=0))

    # st plan A : filter out non Public, non Identical
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Public', 5)
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Identical', -1)

    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list']
    util.make_tsv(WORK_DIR + IN + 'all_ccds_' + FLTD_CDS_INFO, header, ccds_hg38_form_list)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # make_filtered_mouse_ccds_current_file()
    # make_filtered_ccds_current_file_by_shortest_cdn()  # 20201201
    make_filtered_ccds_current_file_by_all_ccds_id()  # 20201209
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
