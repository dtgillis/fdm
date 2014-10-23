__author__ = 'dtgillis'
import os
from textwrap import dedent
import stat


class TBlastXBig:

    def __init__(self, query_file, subject_file, num_jobs, work_dir):

        self.query_file = query_file
        self.subject_file =subject_file
        self.num_jobs = num_jobs
        self.work_dir = work_dir

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
                tmp_line_num = len(query_file)

            chunk_list.append((start_line, tmp_line_num-1))
            out_file = open(self.query_file + '-' + str(i) + '.fasta', 'w')
            out_file.writelines(query_file[start_line:tmp_line_num])
            start_line = tmp_line_num

    def create_subject_db(self):
        cmd = 'makeblastdb -in ' + self.subject_file + ' -hash_index -dbtype nucl'
        os.system(cmd)

    def write_blast_submit_script(self):

        script_file = self.work_dir + os.sep + 'blast_query.sh'
        script_writer = open(script_file, 'w')
        script_writer.write(dedent('''
        #!/bin/bash
        query=$1
        db=$2
        tblastx -num_threads 2 -query $query -db $db -outfmt 6 -out $query.blast.out -evalue .01
        '''))
        script_writer.close()
        os.chmod(script_file,)

        submit_script = self.work_dir + os.sep + 'submit_script.sh'
        script_writer = open(submit_script, 'w')
        script_writer.write(dedent('''
        #!/bin/bash
        QUERY_PREFIX=''' + self.query_file.split(os.sep)[-1] + '''
        SUBJECT_DB=''' + self.subject_file.split(os.sep)[-1] + '''
        for query in $QUERY_PREFIX-*.fasta
        do
            bsub -q week -n 2 -R "span[hosts=1]" -o $query.blast.o -e $query.blase.e ./blast_query.sh $query $SUBJECT_DB
        done;
        '''))
        script_writer.close()
        os.chmod(submit_script, stat.S_IXUSR)


class TBlastX:

    def __init__(self, query, subject, tmp_dir):

        self.query = query
        self.subject = subject
        self.out_file = tmp_dir + os.sep + query.split(os.sep)[-1] + '.tblastx.out'

    def run_tblastx(self):

        cmd = 'tblastx -query {0:s} -subject {1:s} -outfmt 6 -out {2:s} -num_alignments 1 -evalue 1'.format(
            self.query, self.subject, self.out_file
        )
        os.system(cmd)


















