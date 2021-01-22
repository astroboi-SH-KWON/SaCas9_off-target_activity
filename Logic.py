
from astroboi_bio_tools.ToolLogic import ToolLogics
class Logics(ToolLogics):
    def cal_sguide_score(self, sum_scores):
        return 100.0 / (1.0 + sum_scores)

    # def check_seq_in_cds(self, cds_idx_list, trgt_pos):
    #     for cds_idx_arr in cds_idx_list:
    #         for idx in cds_idx_arr:
    #             if trgt_pos == idx:
    #                 return True
    #     return False

    def check_seq_in_cds(self, cds_idx_list, trgt_pos):
        for cds_idx_arr in cds_idx_list:
            if trgt_pos in cds_idx_arr:
                return True
        return False

    def check_which_bin_label(self, range_arr, score):
        for idx in range(len(range_arr)):
            if range_arr[idx][0] > score >= range_arr[idx][1]:
                return idx
