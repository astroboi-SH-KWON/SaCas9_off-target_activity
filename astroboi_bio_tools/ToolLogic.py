from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62


class ToolLogics:
    def __init__(self):
        pass

    def complement_char(self, ch):
        complement_char_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        try:
            return complement_char_dict[ch]
        except:
            print("complement_char : [" + ch + "]")
            raise Exception

    def make_complement_string(self, trgt_seq):
        comp_seq = ""
        for ch in trgt_seq:
            try:
                comp_seq += self.complement_char(ch)
            except:
                raise Exception
        return comp_seq

    """
    checkSeqByChar : match sequences by char
    :param
        seq_char :
        rule_char : 
    :return
        boolean
    """
    def checkSeqByChar(self, seq_char, rule_char):
        flag = False
        if rule_char == 'N':
            return True
        elif rule_char in 'ACGTU':
            if seq_char == rule_char:
                return True
        elif rule_char == 'R':
            if seq_char in 'AG':
                return True
        elif rule_char == 'Y':
            if seq_char in 'CT':
                return True
        """
        add more rules of "ACGTU"
        """

        return flag

    """
    match : match sequence with "same length" strings
    :param
        i : index of seq
        seq_str : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """
    def match(self, i, seq_str, rule_str):
        if len(seq_str) == i:
            return True
        if self.checkSeqByChar(seq_str[i], rule_str[i]):
            return self.match(i + 1, seq_str, rule_str)
        else:
            return False

    """
    get_matched_seq_idx_list : matched sequence index from long strings
    :param
        i : index of rule_str in long_seq
        long_seq : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        index list
    """
    def get_matched_seq_idx_list(self, long_seq, rule_str):
        idx_list = []
        len_rule = len(rule_str)
        for i in range(len(long_seq) - len_rule):
            if self.match(0, long_seq[i: i + len_rule], rule_str):
                idx_list.append(i)
        return idx_list

    """
    by using the BLOSUM62 matrix, together with a gap open penalty of 10 and a gap extension penalty of 0.5 (using globalds)
    """
    def get_pairwise2_globalds_result(self, asequence, bsequence, gap_open_penalty=10, extension_penalty=0.5, matrx=blosum62):
        alignments = pairwise2.align.globalds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                              matrx, -gap_open_penalty, -extension_penalty)
        alignments_result = pairwise2.format_alignment(*alignments[0])
        align_arr = alignments_result.split("\n")
        return align_arr[0], align_arr[1], align_arr[2], alignments_result

    """
    LSPADKTNVKAA
      |..|..|
    --PEEKSAV---
    Score=16
    <BLANKLINE>  
    """
    def get_pairwise2_localds_result(self, asequence, bsequence, gap_open_penalty=10, extension_penalty=1, matrx=blosum62):
        alignments = pairwise2.align.localds(asequence.upper().replace(" ", ""), bsequence.upper().replace(" ", ""),
                                              matrx, -gap_open_penalty, -extension_penalty)
        alignments_result = pairwise2.format_alignment(*alignments[0])
        align_arr = alignments_result.split("\n")
        return ''.join([i for i in align_arr[0] if not i.isdigit()]), ''.join(
            [i for i in align_arr[1] if not i.isdigit()]), ''.join(
            [i for i in align_arr[2] if not i.isdigit()]), alignments_result