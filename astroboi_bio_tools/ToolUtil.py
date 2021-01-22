import glob
from Bio import SeqIO
import openpyxl
import os
import pandas as pd


class ToolUtils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    """
    get file lists in target dir by target ext
    :param
        path : target dir + "*." + target ext
    :return
        ['target dir/file_name.target ext', 'target dir/file_name.target ext' ...]
    """
    def get_files_from_dir(self, path):
        return glob.glob(path)

    def split_big_file_to_files(self, big_f, num_split, max_row):
        # filter out unapproved chromosome
        file_nm_arr = ['chrX', 'chrY']
        for f_num in range(1, 23):
            file_nm_arr.append("chr" + str(f_num))

        with open(big_f) as input_f:
            for num in range(num_split):
                with open(big_f + str(num), 'w') as out_f:
                    cnt = 0
                    for tmp_line in input_f:

                        # filter out unapproved chromosome
                        if tmp_line.split('\t')[0] not in file_nm_arr:
                            continue

                        cnt += 1
                        out_f.write(tmp_line)
                        if cnt == max_row:
                            break

    def read_tsv_ignore_N_line(self, path, n_line=1, deli_str="\t"):
        result_list = []
        with open(path, "r") as f:
            for ignr_line in range(n_line):
                header = f.readline()
                print(header)
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == '':
                    break

                result_list.append(tmp_line.split(deli_str))
        return result_list

    def make_excel_row(self, sheet, row, data_arr, col=1):
        for idx in range(len(data_arr)):
            sheet.cell(row=row, column=(col + idx), value=data_arr[idx])

    def make_excel(self, path, header, data_list, strt_idx=0):
        print("start make_excel :", path)
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        self.make_excel_row(sheet, row, header[strt_idx:])

        for data_arr in data_list:
            row += 1
            self.make_excel_row(sheet, row, data_arr[strt_idx:])

        workbook.save(filename=path + self.ext_xlsx)
        print("end make_excel :", path, "\n")

    def make_tsv(self, path, header, data_list, strt_idx=0, deli='\t'):
        print("start make_tsv :", path)
        with open(path, 'w') as f:
            tmp_head = ''
            for head in header[strt_idx:]:
                tmp_head += (head + deli)
            f.write(tmp_head[:-1] + "\n")

            for data_arr in data_list:
                tmp_row = ''
                for row_val in data_arr[strt_idx:]:
                    tmp_row += (str(row_val) + deli)
                f.write(tmp_row[:-1] + "\n")
        print("end make_tsv :", path, "\n")

    """
    :param
        path : file path with ext
        f_format : file format (ex : fasta, genbank...)
    """
    def read_file_by_biopython(self, path, f_format):
        seq_record = SeqIO.read(path, f_format)
        return str(seq_record.seq).upper(), str(seq_record.seq.complement()).upper()

    """
    :param
        init_split_file = {'big_file_path': WORK_DIR + FASTQ + str(fn_nm) + FASTQ_EXT
                                , 'num_row': 4000000
                                , 'splited_files_dir': WORK_DIR + FASTQ + str(fn_nm) + "/"
                                , 'output_file_nm': only_file_name_without_ext
                                , 'output_file_ext': '.fastq'
                               }
    """
    def split_big_file_by_row(self, init):
        big_file_path = init['big_file_path']
        num_row = init['num_row']
        splited_files_dir = init['splited_files_dir']
        output_file_nm = init['output_file_nm']
        output_file_ext = init['output_file_ext']

        os.makedirs(splited_files_dir, exist_ok=True)

        with open(big_file_path) as fin:
            fout = open('{}/{}_{}{}'.format(splited_files_dir, output_file_nm, '0', output_file_ext), "w")
            for i, line in enumerate(fin):
                fout.write(line)
                if (i + 1) % num_row == 0:
                    fout.close()
                    fout = open(
                        '{}/{}_{}{}'.format(splited_files_dir, output_file_nm, str(i // num_row + 1), output_file_ext),
                        "w")
            fout.close()

    # conda install -c anaconda xlrd
    def get_sheet_names(self, path):
        df = pd.read_excel(path, None)
        return [k for k in df.keys()]

    """
    :return df
    
    len_df = len(df[df.columns[0]])

        for i in range(len_df):
            rg_seq = df.loc[i][2]  # row: i, column: 2
    """
    def read_excel_to_df(self, path, sheet_name='Sheet1', header=0):
        return pd.read_excel(path, sheet_name=sheet_name, header=header)
