import pandas
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import pearsonr
from Bio import SeqIO
window_size = 100
out_name = ''
coverage_table = pandas.read_csv(sys.argv[1],sep='\t', names=['genomeID', 'start', 'end', 'count'])
fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
sequences = {fasta.id:str(fasta.seq) for fasta in fasta_sequences}
GC = np.array([(bp == 'C') or (bp == 'G') for bp in sequences.values()[0]])

pos = []
counts = []
for start_pos, end_pos, count in zip(coverage_table['start'],coverage_table['end'],coverage_table['count']):
    pos.extend(range(start_pos, end_pos))
    counts.extend([int(count)] * (end_pos - start_pos))

fig = plt.figure(figsize=(7.3, 2))
ax = fig.add_axes([0.12,0.25,0.82,0.6])
ax.plot(pos, counts, drawstyle='steps-pre')
ax.set_xlabel('position (bp)')
ax.set_ylabel('coverage (count)')
plt.savefig('coverage%s.pdf'%out_name)
plt.close()

ave_GC = []
ave_count = []
for ii in range(0, len(pos), window_size):
    ave_GC.append(np.mean(GC[ii:(ii+window_size)]))
    ave_count.append(np.mean(counts[ii:(ii+window_size)]))

fig = plt.figure(figsize=(3.5, 3.5))
ax = fig.add_axes([0.25,0.25,0.7,0.6])
ax.scatter(ave_GC, ave_count)
#ax.axvline(np.mean(GC),c=[0.3,0.3,0.3])
ax.set_xlabel('GC (%i bp window)'%window_size)
ax.set_ylabel('coverage (%i bp window)'%window_size)
ax.set_title('correlation:%.2f'%pearsonr(ave_GC,ave_count)[0])
xlims = ax.get_xlim()
ax.set_xlim(left=np.max([xlims[0], 0]), right=np.min([xlims[1],1]))
ylims = ax.get_ylim()
ax.set_ylim(bottom=np.max([ylims[0], 0])) 
plt.savefig('GC%s.pdf'%out_name)
plt.close()