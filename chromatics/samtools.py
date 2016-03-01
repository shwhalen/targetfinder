import chromatics
import io
import pandas as pd
import subprocess

def samtools(arguments):
    cmdline = 'samtools {}'.format(arguments)
    p = subprocess.Popen(cmdline, shell = True, stdout = subprocess.PIPE)
    stdout, _ = p.communicate()
    return chromatics.read_bed(io.StringIO(stdout.decode('utf-8'))).set_index(0)

if __name__ == '__main__':
    reads_df = samtools('view test_dataset1_2.hicup.bam')
    print(reads_df[reads_df[4] < 30])
    print(reads_df.loc['SRR071233.1357221'])

    reads_df = samtools('view -q 30 test_dataset1_2.hicup.bam')
    assert len(reads_df[reads_df[4] < 30]) == 0
    print(reads_df.loc['SRR071233.1357221'])
