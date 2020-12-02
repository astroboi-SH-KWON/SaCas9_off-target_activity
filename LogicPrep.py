
from astroboi_bio_tools.ToolLogicPrep import ToolLogicPreps
class LogicPreps(ToolLogicPreps):
    def make_dict_to_list(self, result_dict):
        result_list = []
        for guide_seq, sum_scores in result_dict.items():
            result_list.append([guide_seq, sum_scores])
        return result_list
