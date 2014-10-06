__author__ = 'dtgillis'
import os


class BlastN:

    def __init__(self, query_file, subject_file, num_jobs):

        self.query_file = query_file
        self.subject_file =subject_file
        self.num_jobs = num_jobs

    def break_down_query(self):

        # read in the subject file
        query_file = [line for line in open(self.query_file, 'r')]
        total_lines = len(query_file)
        lines_per_chunk = total_lines / self.num_jobs
        chunk_list = []
        start_line = 0

        for i in range(self.num_jobs):
            if i != self.num_jobs - 1:
                tmp_line_num = start_line + lines_per_chunk
                tmp_line = query_file[tmp_line_num]
                while ">" not in tmp_line:
                    tmp_line_num -= 1
                    tmp_line = query_file[tmp_line_num]
            else:
                tmp_line_num = -1

            chunk_list.append((start_line, tmp_line_num-1))
            out_file = open(self.query_file + '-' + str(i) + '.fasta', 'w')
            out_file.writelines(query_file[start_line:tmp_line_num])
            start_line = tmp_line_num


    def create_subject_db(self):
        cmd = 'makeblastdb -in ' + self.subject_file + ' -hash_index -dbtype nucl'
        os.system(cmd)

















