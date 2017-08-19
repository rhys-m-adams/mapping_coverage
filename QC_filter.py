#!/opt/hpc/bin/python2.7
#/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
import sys
from itertools import izip

class file4(file):
    def next(self):
        out = []
        for l in range(4):
            temp = self.readline()
            out.append(temp.rstrip())
        if '' in out:
            raise StopIteration
        return out


def ave_phred(x):
    intx = map(ord, x)
    return sum(intx)/float(len(x)) - 33


def read_Illumina(file1, file2, out1, out2, QC_cutoff):
    f1 = file4(file1)
    f2 = file4(file2)
    fout1 = open(out1,'w')
    fout2 = open(out2,'w')
    for line_1, line_2 in izip(f1, f2):
        if (ave_phred(line_1[3]) >= QC_cutoff) and (ave_phred(line_2[3]) >= QC_cutoff):
            fout1.writelines("%s\n" % l for l in line_1)
            fout2.writelines("%s\n" % l for l in line_2)
    
    fout1.close()
    fout2.close()
    f1.close()
    f2.close()

file_name_1 = sys.argv[1]
file_name_2 = sys.argv[2]
file_name_3 = sys.argv[3]
file_name_4 = sys.argv[4]
QC_cutoff = int(sys.argv[5])
read_Illumina(file_name_1, file_name_2, file_name_3, file_name_4, QC_cutoff)