import re

from astroboi_bio_tools.ToolLogicPrep import ToolLogicPreps
class LogicPreps(ToolLogicPreps):
    def make_dict_to_list(self, result_dict):
        result_list = []
        for guide_seq, sum_scores in result_dict.items():
            result_list.append([guide_seq, sum_scores])
        return result_list

    def make_dict_to_list2(self, result_dict):
        result_list = []
        for guide_seq, val_arr in result_dict.items():
            tmp_arr = [guide_seq]
            tmp_arr.extend(val_arr)
            result_list.append(tmp_arr)
        return result_list

    def get_data_with_trgt_strng(self, ccds_list, trgt_str, idx):
        return [cds_arr for cds_arr in ccds_list if cds_arr[idx] == trgt_str]

    def get_highest_ccds_id_among_same_gen_id(self, ccds_list):
        gen_id_dict = {}
        result_list = []
        for ccds_arr in ccds_list:
            gene = ccds_arr[2]
            ccds_id = ''.join(x for x in ccds_arr[4] if x.isdigit())
            ccds_arr.append(ccds_id)
            if gene in gen_id_dict:
                gen_id_dict[gene].append(ccds_arr)
            else:
                gen_id_dict.update({gene: [ccds_arr]})

        for ccds_list in gen_id_dict.values():
            """
            sangyeon  오후 6:16 20201123
            그리고 CCDS id의 순선는 큰 숫자가 더 최근에 나온 data네요!
            """
            sorted_ccds_list = self.sort_list_by_ele(ccds_list, -1, True)
            result_list.append(sorted_ccds_list[0][:-1])
        return result_list

    """
    mouse_ccds form : 
        #chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type
    hg38_refFlat form : 
        GeneSym   NMID    Chrom   Strand  Transcript_Start   End ORFStart    End #Exon   ExonS_list  ExonE_list

    GeneSyn = gene
    NMID = ccds_id (이 부분은 human을 filter하실 때, gene 하나당 NMID가 여러개였던걸 줄이신거면 동일한 것으로 생각됩니다)
    Chrom = #chromosome
    Strand = cds_strand
    Transcript_Start =
    End =
    ORF Start = cds_from
    End = cds_to
    #Exon, ExonS_list, ExonE_list = cds_locations
    """
    def transform_mouse_ccds_form_to_hg38_refFlat(self, ccds_list):
        result_list = []
        for ccds_arr in ccds_list:
            gene_sym = ccds_arr[2]
            nm_id = ccds_arr[4]
            chrom = 'chr' + str(ccds_arr[0])
            strand = ccds_arr[6]

            transcript_st = ccds_arr[7]
            transcript_en = ccds_arr[8]
            orf_st = ccds_arr[7]
            orf_en = ccds_arr[8]

            cds_locations_arr = re.findall('\d+', ccds_arr[9])
            num_exon = len(cds_locations_arr) // 2
            exon_s_list = ''
            for cds_loc in [cds_locations_arr[i] for i in range(len(cds_locations_arr)) if i % 2 == 0]:
                exon_s_list = exon_s_list + str(cds_loc) + ','
            exon_e_list = ''
            for cds_loc in [cds_locations_arr[i] for i in range(len(cds_locations_arr)) if i % 2 != 0]:
                exon_e_list = exon_e_list + str(cds_loc) + ','

            result_list.append(
                [gene_sym, nm_id, chrom, strand, transcript_st, transcript_en, orf_st, orf_en, num_exon, exon_s_list,
                 exon_e_list])

        return result_list

    """
    CCDS_id가 여러개인 gene에서 하나씩만 가져오는 기준이 
    현재는 가장 최근의 CCDS_id만을 가져오는 방식으로 filter를 했었는데, ==> self.get_highest_ccds_id_among_same_gen_id
    혹시 이 부분을 가장 짧은 cds를 가지는 CCDS_id로 가져오게 변경해서 다시 뽑아주실 수 있나요? 
    """
    def get_shortest_cdn_among_same_gen_id(self, ccds_hg38_form_list):
        gen_id_dict = {}
        result_list = []
        for ccds_arr in ccds_hg38_form_list:
            gene = ccds_arr[0]
            start_idx_arr, end_idx_arr = self.get_orf_strt_end_idx_arr(ccds_arr)
            idx_list = self.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)
            ccds_arr.append(len(idx_list))
            if gene in gen_id_dict:
                gen_id_dict[gene].append(ccds_arr)
            else:
                gen_id_dict.update({gene: [ccds_arr]})

        for ccds_list in gen_id_dict.values():
            sorted_ccds_list = self.sort_list_by_ele(ccds_list, -1, False)
            result_list.append(sorted_ccds_list[0][:-1])
        return result_list

    def get_orf_strt_end_idx_arr(self, cds_arr):
        orf_strt_idx = int(cds_arr[6])
        orf_end_idx = int(cds_arr[7])

        start_idx_arr = [orf_strt_idx]
        start_idx_arr.extend([int(cds_idx) for cds_idx in re.findall('\d+', cds_arr[9])[1:]])
        end_idx_arr = [int(cds_idx) for cds_idx in re.findall('\d+', cds_arr[10])[:-1]]
        end_idx_arr.append(orf_end_idx)

        # check UTR in first and last codon
        sorted_start_idx_arr = sorted(start_idx_arr)
        sorted_end_idx_arr = sorted(end_idx_arr)
        strt_idx = sorted_start_idx_arr.index(orf_strt_idx)
        end_idx = sorted_end_idx_arr.index(orf_end_idx)
        return sorted_start_idx_arr[strt_idx: end_idx + 1], sorted_end_idx_arr[strt_idx: end_idx + 1]

    def get_idx_num_frm_strt_to_end_list(self, start_idx_arr, end_idx_arr):
        result_list = []
        for idx in range(len(start_idx_arr)):
            nxt_idx = int(start_idx_arr[idx])
            end_idx = int(end_idx_arr[idx])
            result_list.extend([tmp_idx for tmp_idx in range(nxt_idx, end_idx)])
        return result_list

    def get_cds_idx_list_by_chr(self, cds_list):
        result_dict = {}
        for val_arr in cds_list:
            chr_key = val_arr[2]
            st_cdn_arr, en_cdn_arr = self.get_orf_strt_end_idx_arr(val_arr)
            idx_list = self.get_idx_num_frm_strt_to_end_list(st_cdn_arr, en_cdn_arr)
            if chr_key in result_dict:
                result_dict[chr_key].append(idx_list)
            else:
                result_dict.update({chr_key: [idx_list]})
        return result_dict
