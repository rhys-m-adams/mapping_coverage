import pandas
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import pearsonr
from Bio import SeqIO

def window_averaged(in_list, window_size):
    #Average and return every window in in_list, with width defined by window_size
    out = []
    for ii in range(0, len(in_list), window_size):
        out.append(np.mean(in_list[ii:(ii+window_size)]))
    return out

def plot_genome(pos, attribute, attribute_name, file_prefix):
    #plot genome wide plots.
    #pos is the position in the genome, should be sorted
    # attribute should be something like GC or coverage, for instance
    # attribute_name is y-axis label
    # file_prefix, name of output file
    fig = plt.figure(figsize=(7.3, 2))
    ax = fig.add_axes([0.12,0.25,0.82,0.6])
    ax.plot(pos, attribute, drawstyle='steps-pre')
    ax.set_xlabel('position (bp)')
    ax.set_ylabel(attribute_name)
    plt.savefig('%s.pdf'%file_prefix)
    plt.close()

def plot_pos_correlation(x,y,xname,yname,file_prefix):
    #scatter plot with correlation and p-value title. The x- and y- axes are not allowed to be less than 0
    #x - list of values to be plotted along the x-axis
    #y - list of values to be plotted along the y-axis
    #xname - name of the xlabel
    #yname - name of the ylabel
    #file_prefix - name of the output plot
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.add_axes([0.25,0.15,0.7,0.7])
    ax.scatter(x, y)
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    corr, pval = pearsonr(x, y)
    power = np.floor(np.log10(pval))
    ax.set_title('correlation:%.2f'%corr+'\n'+r'p-val:$%.1f \times 10^{%i}$'%(pval/10**power,power))
    xlims = ax.get_xlim()
    ax.set_xlim(left=np.max([xlims[0], 0]))
    ylims = ax.get_ylim()
    ax.set_ylim(bottom=np.max([ylims[0], 0])) 
    plt.savefig('%s.pdf'%file_prefix)
    plt.close()

if __name__ == '__main__':
    window_size = 100 #default window-size for averaging GC content and counts is 100
    out_name = '' #default out_name is empty
    bedtools_name = sys.argv[1]
    fasta_reference_name = sys.argv[2]
    coverage_table = pandas.read_csv(bedtools_name,sep='\t', names=['genomeID', 'position', 'count']) #get positions and coverage counts from bedtools genomic count file
    fasta_sequences = SeqIO.parse(open(fasta_reference_name),'fasta') #get the reference fasta_sequences
    sequences = {fasta.id:str(fasta.seq) for fasta in fasta_sequences} # genomic/reference fasta sequences are stored in a dictionary in case I want to analyze multiple genomes/chromosomes/reference sequences at a time
    if len(sys.argv)>3: #if 3 or more input argumetns, set the window size as the 3rd argument
        window_size = int(sys.argv[3])
    
    if len(sys.argv)>4: #if 4 or more input arguments, set the output name as the 4th argument
        out_name = sys.argv[4]

    for genome in sequences.keys(): # for each genome in the reference file, output the coverage and GC correlation
        GC = np.array([(bp == 'C') or (bp == 'G') for bp in sequences[genome]]) #make a boolean list of whether the base pairs in the first sequence are G or C, used to calculate GC content
        current_table_genome = coverage_table['genomeID'] == genome
        pos = np.array(coverage_table['position'].loc[current_table_genome]).flatten()
        counts = np.array(coverage_table['count'].loc[current_table_genome]).flatten()
        ave_GC = window_averaged(GC, window_size) # calculate average GC content
        ave_count = window_averaged(counts, window_size) # calculate average read coverage, paired with ave_GC
        
        if len(sequences.keys())==1: 
        #if there are multiply fasta entries, then output plots with fasta entry names, otherwise no fasta entry names. This is to satisfy the requirement
        #that no additionaly options will output the names coverage.pdf, and correlation.pdf
            out_prefix = out_name
        else:
            out_prefix = out_name + genome
        
        plot_genome(pos, counts, 'coverage (count)', out_prefix+'coverage')#output the genome coverage file


        plot_pos_correlation(ave_GC, ave_count, 'GC (%i bp window)'%window_size, 'coverage (%i bp window)'%window_size, out_prefix+'correlation') #scatter plot with correlation statistics of GC content versus RNA coverage
