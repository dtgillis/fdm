__author__ = 'dtgillis'


class BlastN:

    def __init__(self, query_file, subject_file, num_jobs):

        self.query_file = query_file
        self.subject_file =subject_file
        self.num_jobs = num_jobs

    def break_down_subject(self):

        # read in the subject file
        subject_file = [line for line in open(self.subject_file, 'r')]
        total_lines = len(subject_file)
        lines_per_chunk = total_lines / self.num_jobs
        chunk_list = []
        start_line = 0

        for i in range(self.num_jobs):
            tmp_line_num = start_line + lines_per_chunk
            tmp_line = subject_file[tmp_line_num]
            while ">" not in tmp_line:
                tmp_line_num -= 1
                tmp_line = subject_file[tmp_line_num]

            chunk_list.append((start_line, tmp_line_num-1))
            out_file = open(self.subject_file + '-' + str(i) + '.fasta', 'w')
            out_file.writelines(subject_file[start_line:tmp_line_num])
            start_line = tmp_line_num













