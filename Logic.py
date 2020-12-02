
from astroboi_bio_tools.ToolLogic import ToolLogics
class Logics(ToolLogics):
    def cal_sguide_score(self, data_list):
        result_list = []
        for val_arr in data_list:
            guide_seq = val_arr[0]
            sum_scores = val_arr[1]
            # 3-4) from README.md
            s_guide = 100.0 / (1.0 + sum_scores)
            result_list.append([guide_seq, s_guide])
        return result_list
