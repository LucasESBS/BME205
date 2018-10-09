##################################################################
## MNM.py file: find underrepresented motifs in fasta sequence  ##
## by: Lucas Seninge; credits for help to: Rick Tao,            ##
## James McBroome                                               ##
## thanks to: BME cohort 2018                                   ##
##################################################################

#Import a few modules
from itertools import product
from scipy.stats import binom
from math import sqrt
import sys
import argparse

#Parser -- created my own
class Parser:
    """Class to parse commandline arguments to MNM.py"""
    
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='Find under represented motifs in a fasta file and rank them by z-score',
                                              formatter_class=argparse.RawDescriptionHelpFormatter,
                                              epilog='''Example usage: python MNM.py < path/to/genome.fa > path/to/report.out''')
        self.parser.add_argument('-m','--minMotif',nargs='?',
                            const=3,
                            default=3,
                            help='minimum motif size to evaluate. Default: 3', required=False, type=int)
        self.parser.add_argument('-M','--maxMotif',nargs='?',
                            const=8,
                            default=8,
                            help='maximum motif size to evaluate. Default: 8', required=False, type=int)
        self.parser.add_argument('-c','--cutoff',nargs='?',
                            const=-3,
                            default=-3,
                            help='cut-off for z-score. Default: -3', required=False, type=int)

        self.args = self.parser.parse_args()

#Class FastaReader (credit: Pr. David Bernick)
class FastAreader :
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
 
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header,sequence


#Reverse strand function - to be used anywhere
def reverse_strand(sequence):
    """This function takes a DNA sequence as input and convert it to its reverse strand
    param: sequence: A DNA sequence as a string
    return: Reverse strand as a string"""
    
    if set(sequence).issubset('ATCG'):
        rev_dict={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        return "".join(rev_dict.get(base, base) for base in reversed(sequence))
    else:
        return('reverse_strand only support canonical bases!')

#Class MissingMotif
class MissingMotif:
    """A class to find  under represented motif of DNA in bacterial genome
    Default values for minMotif, maxMotif and cuttoff are set but can be
    changed by specifying commandline arguments"""
    
    def __init__(self, minMotif, maxMotif, cutoff):
        self.seq=FastAreader().readFasta()
        self.minMotif=minMotif
        self.maxMotif=maxMotif
        self.cutoff=cutoff
        self.count_dict={}
        self.genome_size=0
    
    def initialize_count_dict(self):
        """This function prepopulates the count dictionary and link kmer count and 
        its reverse using mutable object (list)"""
        
        for k in range(1,self.maxMotif+1):
            kmer_size=k
            #From stackoverflow
            kmer_list = [''.join(i) for i in product(['A','T','C','G'], repeat = k)]
            for kmer in kmer_list:
                self.count_dict[kmer] = self.count_dict[reverse_strand(kmer)] = [0]
    
        
    def kmer_dict_count(self):
        """This function counts kmer over the input sequence, skipping Ns"""
        
        #Iterate on fasta sequences if several
        for head, sequence in self.seq:
            #Update genome_size for later statistics
            self.genome_size+=len(sequence)
    
            #Iterate over all k-mers > max_size because we need the smaller ones 
            for kmer_size in range(1,self.maxMotif+1):
                #Count k-mers in the sequence
                for kmer in range(len(sequence)-(self.maxMotif-1)):
                    #Check for N in the kmer
                    if 'N' not in sequence[kmer:kmer+kmer_size]:
                        #To deal with palindromes, we should add 2 when kmer=reverse;
                        #because by adding 1 we undercount palindromes
                        #To match the output of Pr. Bernick, we leave 1 here
                        self.count_dict[sequence[kmer:kmer+kmer_size]][0]+=1

        
    def kmer_expected_count(self):
        """This function computes expected kmer counts from existing kmers using the course
        formula to build our Null model. Expected counts are put in the value (a list) of self.count_dico"""
        #Compute expected kmer count from actual counts
        for kmer, count in self.count_dict.items():
            #We can only compute expected counts starting from min_size
            if len(kmer)>=self.minMotif and len(count)<2:
                #Avoid division by zero!
                if self.count_dict[kmer[1:-1]][0]==0:
                    self.count_dict[kmer].append(0)
                else:
                    self.count_dict[kmer].append((self.count_dict[kmer[:-1]][0]*self.count_dict[kmer[1:]][0])/self.count_dict[kmer[1:-1]][0])
                    
        
    def clean_reverse(self):
        """This function aims at removing the reverse strand kmers after computing the expected counts
        We counted the reverse kmers at the same time we counted the forward, so they can be removed now
        To avoid printing duplicates.
        Also removes kmers<min_size"""
        
        clean_dict={}
        for kmer,value in self.count_dict.items():
            if (kmer not in clean_dict and reverse_strand(kmer) not in clean_dict) and len(kmer)>=self.minMotif:
                clean_dict[kmer]=value
        
        #Update object
        self.clean_dict=clean_dict
     
    
    
    def compute_zscore_pval(self):
        """This function compute a zscore and a pvalue given the null distribution
        built from the expected number of counts. The 2 values are added to the list in the clean_dict
        in the following order: zscore, pval"""
        
        #Iterate over the kmer dictionary
        for kmer, value in self.clean_dict.items():
            #Only compute for length>min_size
            if len(kmer)>=self.minMotif:
                #Compute zscore from counts, expected counts and standard deviation
                score = (self.clean_dict[kmer][0]-self.clean_dict[kmer][1])/(sqrt(self.clean_dict[kmer][1]*(1-(self.clean_dict[kmer][1]/self.genome_size))))
                self.clean_dict[kmer].append(score)
                
                #Lower tailed test? 
                pval=binom.cdf(self.clean_dict[kmer][0], self.genome_size, (self.clean_dict[kmer][1]/self.genome_size))
                self.clean_dict[kmer].append(pval)
            
        #Remove kmers with zscore above threshold
        self.clean_dict = {kmer:value for kmer,value in self.clean_dict.items() if value[2] < self.cutoff}
        
    def output_result(self):
        """This function output a file with tab separated fields as a summary
        Kmers are ordered by length and zscore, and alphabetized. Reverse kmers is also printed."""
        
        #Insert sequence + reverse sequence in clean_dict value list at index 0
        #This way we can use the string formatting from Pr. Bernick
        for kmer in self.clean_dict.keys():
            self.clean_dict[kmer][0:0]=[kmer,reverse_strand(kmer)]
            #Sort first 2 elements (strand and reverse) alphabetically
            self.clean_dict[kmer]=sorted(self.clean_dict[kmer][:2])+self.clean_dict[kmer][2:]
        
        #Print headers
        print('sequence: reverse\tcount\tExpect\tZscore\tp-value')
        #Print by sorted by length & z-score; print seq, revseq, count, expect, z-score and pval
        #From stackoverflow
        for kmer in sorted(self.clean_dict.items(), key=lambda x: (-len(x[0]), x[1][4])):
            print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}\t{5:0.2e}'.format(kmer[1][0], kmer[1][1], kmer[1][2], kmer[1][3], kmer[1][4], kmer[1][5]))      
        

def main():
    """Main function to be run for MNM.py using the MissingMotif class"""
    #Parse arguments
    myParser=Parser()
    minMotif=myParser.args.minMotif
    maxMotif=myParser.args.maxMotif
    cutoff=myParser.args.cutoff
    #Raise error for some values
    if minMotif<3 or maxMotif>8 or cutoff>0:
        raise ValueError('minMotif should be greater or equal to 3, and maxMotif lower or equal to 8. Z-score cutoff should be lower than 0')
    else:
        #Count and produce output
        run=MissingMotif(minMotif=minMotif,maxMotif=maxMotif, cutoff=cutoff)
        run.initialize_count_dict()
        run.kmer_dict_count()
        run.kmer_expected_count()
        run.clean_reverse()
        run.compute_zscore_pval()
        run.output_result() 

if __name__=='__main__':
    main()   




