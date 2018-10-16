#!/usr/bin/env/ python3
#Lucas Seninge (lseninge)
#Group members: bme-2018, stack overflow

#Last updated: 15 October 2018
##################################################################
# File: randomizedMotifSearch.py
# executable: randomizedMotifSearch.py
# Purpose: find a consensus promoter motif given a fasta file 
# containing CRISPR array upstream sequences 
# Author: Lucas Seninge
# Credits for help to: Alex B.        
# thanks to: BME cohort 2018                               
##################################################################

""" This script can be used to find a consensus sequence of given length in input sequences
    Example usage:Example usage: python randomizedMotifSearch.py < path/to/genome.fa -i=10000 -p=1 -k=13 > path/to/report.out
"""


#Import a few useful modules
import sys
import argparse
from random import randint
from random import shuffle
import numpy as np
from scipy.stats import entropy

#Parser -- derived from assignment
class Parser:    
    def __init__(self):
        """Class to parse commandline arguments to randomizedMotifSearch.py
        attributes: arguments received by the commandline including
        number of iterations, number of pseudocounts, length of the motif. 
        Optional arguments are available to use Gibbs sampling, scramble the 
        input sequences and print the contributing sequences to the motifs.
        Defaults values are used if values are not specified by user for
        optional features.
        """
        
        self.parser = argparse.ArgumentParser(description='Find consensus promoter sequence from input upstream sequences',
                                              formatter_class=argparse.RawDescriptionHelpFormatter,
                                              epilog='''Example usage: python randomizedMotifSearch.py < path/to/genome.fa -i=10000 -p=1 -k=13 > path/to/report.out''')
        #Iterations -- mandatory
        self.parser.add_argument('-i','--iterations',nargs='?',
                            help='Number of iterations to run for convergence', required=True, type=int)
        #Pseudocounts -- mandatory
        self.parser.add_argument('-p','--pseudocounts',nargs='?',
                            help='Number of pseudocounts to add', required=True, type=int)
        #Length -- mandatory
        self.parser.add_argument('-k','--kmer',nargs='?',
                            help='Length of the motif to be found', required=True, type=int)
        #Scramble -- extra credit, optional (default to False)
        self.parser.add_argument('-r','--random',
                            action='store_true',
                            help='Shuffle the input sequences. Default to False', required=False)
        #Gibbs sampling -- extra credit, optional (default to False)
        self.parser.add_argument('-g','--gibbs',
                            action='store_true',
                            help='Use gibbs sampling to find consensus. Default to False', required=False)
        #Print contributing sequences -- extra credit, optional (default to False)
        self.parser.add_argument('-m','--mprint',
                            action='store_true',
                            help='Print contributing sequence to the motifs. Default to False', required=False)
        
        self.args = self.parser.parse_args()

#FastaReader from Pr. Bernick
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


#Class to find the consensus sequence using randomized motif search        
class ConsensusCRISPRPromoter:
    def __init__(self, sequence, iterations, pseudocounts, motif_length, scramble, gibbs, mprint):
        """
        Find a consensus sequence from input sequences using randomized search/gibbs sampling strategy.

        Args:
            sequence: sequences of dna to be analyzed
            iterations: number of iterations to perform
            pseudocounts: number of pseudocounts
            motif_length: length of the consensus sequence to build
            scramble: if True, scramble (shuffle) input sequences
            gibbs: if True, uses gibbs sampling to find best motifs
            mprint: if True, print each sequence of the best motifs and the name of the contributing sequence
        """
        
        self.sequence=sequence
        self.iterations=iterations
        self.pseudocounts=pseudocounts
        self.motif_length=motif_length
        self.nb_sequence=len(self.sequence)
        #This will be used to map back to our profile matrix when needed
        self.base_map=['A','T','C','G']
        
        #Initialize best_score to high value
        self.best_score=float('inf')
        
        #Extra credits
        self.scramble=scramble
        self.gibbs=gibbs
        self.mprint=mprint
        
    
    def randomized_motif(self):
        """
        This method randomly select a collection of kmers at each round and tries to converge towards
        an optimal collection of kmers for motif finding, by minimizing a scoring function.

        Returns:
            best_motifs: the convergent solution of randomized motif search as a list of kmers
        """
        
        #Initialize best_motifs(kmer) and motifs(kmer) 
        best_motifs=[]
        motifs=[]
        #Iterate over sequences in file
        for head,sequence in self.sequence.items():
            #Initialize random kmer of length motif_length in sequence
            rand_start=randint(0,(len(sequence)-self.motif_length))
            random_kmer=sequence[rand_start:(rand_start+self.motif_length)]
            motifs.append(random_kmer)
        
        best_motifs=motifs
        #Loop to find minimal score ie best motifs
        while True:
            #Compute profile and score associated with it
            profile=self.build_profile(motifs)
            #Compute new motifs based on most probable kmers
            motifs=self.probable_motif(profile)
            #Compare scores with best_motifs
            if self.score_profile(self.build_profile(motifs))<self.score_profile(self.build_profile(best_motifs)):
                best_motifs=motifs
            else:
                return best_motifs
            
    
    #Extra Credit - Gibbs Sampling
    def gibbs_sampler(self, N=1000):
        """
        This method uses the gibbs sampling strategy to find the best collection of kmers for motif finding.

        Args:
            N: number of iterations per trajectory

        Returns:
            best_motifs: the convergent solution of gibbs sampling as a list of kmers
        """
        
        #Initialize best_motifs(kmer) and motifs(kmer) -- same as randomizedMotifSearch
        best_motifs=[]
        motifs=[]
        #Iterate over sequences in file
        for head,sequence in self.sequence.items():
            #Initialize random kmer of length motif_length in sequence
            rand_start=randint(0,(len(sequence)-self.motif_length))
            random_kmer=sequence[rand_start:(rand_start+self.motif_length)]
            motifs.append(random_kmer)
        
        #N iterations
        for j in range(1,N):
            #Select random DNA sequence
            i=randint(0,self.nb_sequence-1)
            #Remove i-th DNA sequence kmer from motifs  to build profile
            #(kmer are ordered by sequence since we use append)
            profile=self.build_profile((motifs[:i]+motifs[i+1:]))
            #Randomly generate a kmer from i-th sequence based on profile - ith seq
            #Calling gibbs_random_kmer method
            motifs[i]=self.gibbs_random_kmer(seq=list(self.sequence.values())[i],profile=profile)
            #Compare score with best_motifs
            if self.score_profile(self.build_profile(motifs))<self.score_profile(self.build_profile(best_motifs)):
                best_motifs=motifs
        
        return best_motifs
            
            
    def gibbs_random_kmer(self, seq, profile):
        """
        This method generates a random kmer given a profile and a sequence of DNA.
        This is to be used in the gibbs_sampler function to generate a new, random kmer from the ith sequence.

        Args:
            seq: a string for DNA sequence selected
            profile : a np.array generated without the kmer of the ith sequence

        Returns:
            random_kmer: new kmer generated based on seq and profile probabilities
        """
        
        #Initialize dictionary of kmer:probability to pick one based on the 'random dice' later
        kmer_proba={}
        #Sliding over kmers on input seq
        for k in range(len(seq)-(self.motif_length-1)):
            #Probabilities based on profile
            #Initialize probability
            p=1
            kmer=seq[k:k+self.motif_length]

            #Compute probability for kmer given profile
            for pos, base in enumerate(kmer):
                p=p*profile[self.base_map.index(base)][pos]
            #Update dictionary with kmer and its probability given our profile
            kmer_proba[kmer]=p
        #Correct probabilities with the sum C of kmers probabilities
        #So we can use it as a random kmer generator
        C=sum(kmer_proba.values())
        kmer_proba={key:value/C for key,value in kmer_proba.items()}
        #Pick a kmer given new probabilities -- From Stackoverflow
        random_kmer=np.random.choice(a=list(kmer_proba.keys()), p=list(kmer_proba.values()))
        
        return random_kmer
        
        
        
    def build_profile(self, motifs):
        """
        Compute a profile based on a motifs list, 
        adding pseudocounts to soften the probability distribution.

        Args:
            motifs: a list of randomly selected kmers

        Returns:
            profile: a matrix with positional probabilities for each base
        """
        
        #Initialize profile with zeros. Order of rows A,T,C,G
        profile=np.zeros(shape=(4,self.motif_length))
        #Initialize dictionary to count bases// From stack overflow
        count_base={'A':[0]*self.motif_length, 'T':[0]*self.motif_length, 'C':[0]*self.motif_length,'G':[0]*self.motif_length}
        #Iterate over random kmers
        for seq in motifs:
            #Find base at given position
            for pos, base in enumerate(seq):
                count_base[base][pos]+=1
        
        #Fill profile matrix
        #This list helps us to map values in the matrix using dictionary keys
        for base, seq in count_base.items():
            for pos in range(len(seq)):
                profile[self.base_map.index(base)][pos]=count_base[base][pos]

                    
        #Add pseudocounts -- don't forget to add pseudocounts as 'sequences'
        profile+=self.pseudocounts
        profile=profile/(self.nb_sequence+self.pseudocounts*4)
        
        return profile    
    
    def probable_motif(self, profile):
        """
        This method computes the most probable motif of DNA sequences from a profile matrix.

        Args:
            profile: matrix containing positional probabilities

        Returns:
            new_motif: list of strings of most probable kmer for each DNA sequence
        """
        
        #Iterate on the DNA sequences
        new_motif=[]
        for seq in self.sequence.values():
            #Initialize most probable kmer and probability
            best_kmer=''
            best_p=0
            #Iterate over the sequence to get each kmers
            for k in range(len(seq)-(self.motif_length-1)):
                p=1
                kmer=seq[k:k+self.motif_length]
                #Compute probability for given kmer
                for pos, base in enumerate(kmer):
                    p=p*profile[self.base_map.index(base)][pos]
                #Compare with best so far
                if p>best_p:
                    best_kmer=kmer
                    best_p=p
            #Add the most probable kmer to the list
            new_motif.append(best_kmer)
        
        return new_motif
                    
                    
                
    def score_profile(self, profile):
        """
        This method scores a given profile using entropy

        Args:
            profile: a np.array representing a profile
            
        Returns:
            sum_entropy: The score of the given profile
        """
        
        sum_entropy=0
        for pos in range(self.motif_length):
            sum_entropy+=entropy(profile[:,pos], base=2)

        return sum_entropy
    
    
    def iterate_motif_search(self):
        """
        This method iterates i times the selected method to find a consensus sequence
        in an input fasta file. The method check for the score at each iteration and only
        keeps the best one.
        """
        
        #Before iterating, if scramble==True, scramble sequences
        if self.scramble==True:
            self.scramble_sequences()
        
        #Initialize stopper: if the best score doesn't change for 1000 iter: Break
        #This to avoid running useless iterations,especially for the gibbs sampling method
        stopper=0
        
        #Iterate i times
        for i in range(self.iterations):
            
            
            #Check if gibbs sampling should be used for getting best_motifs
            if self.gibbs==False:
                #Get best_motifs from the randomizedMotifSearch method
                best_motifs=self.randomized_motif()
            else:
                #Get best_motifs from the gibbs_sampler method
                best_motifs=self.gibbs_sampler()
                
            #Build profile for final best_motifs
            best_profile=self.build_profile(best_motifs)
            #Score it
            score=self.score_profile(best_profile)
            #Update class attributes if better (lower score) during the iterations
            if score<self.best_score:
                self.best_profile=best_profile
                self.best_score=score
                self.best_motifs=best_motifs
                #Reset stopper
                stopper=0
            else:
                stopper+=1
                if stopper>1000:
                    break
                

    
    def build_consensus(self):
        """
        This method build a consensus sequence given the best motifs/profile 
        found during the analysis steps.
        The consensus sequence is printed to a file using STDOUT.
        """
        
        #Use this list to build the sequence from the profile np.array
        consensus=[]
        #Iterate over the profile to build the consensus 
        for pos in range(self.motif_length):
            #New version -- in the unlikely case where several bases have the same max. probability
            #We use argwhere instead of argmax (which would only return the first occurence indice)
            #and randomize the indice choice to avoid overselecting a base (eg. As)
            #From Stack Overflow
            ind_max=np.argwhere(self.best_profile[:,pos]==np.amax(self.best_profile[:,pos])).flatten().tolist()
            if len(ind_max)>1:
                ind=np.random.choice(ind_max)
            else:
                ind=ind_max[0]
            consensus.append(self.base_map[ind])
        #Print consensus sequence and score
        print('Consensus sequence\tScore')
        print(''.join(base for base in consensus),'\t',self.best_score)
        
        #mprint option -- Extra credit
        #Print sequence in best_motifs and name of contributing sequence
        if self.mprint==True:
            print('Sequence in motifs\tContributing sequence')
            for head in list(self.sequence.keys()):
                print(self.best_motifs[list(self.sequence.keys()).index(head)],'\t',head.split(' ')[0])
            
    
    
    def scramble_sequences(self):
        """
        This method scrambles each DNA sequences in the set of sequences
        to establish a baseline score to be evaluate against non-scrambled run.
        """
        
        #Iterate over sequences in self.sequence dictionary
        for head, seq in self.sequence.items():
            #From Stack Overflow
            l_seq=list(seq)
            #Shuffle
            shuffle(l_seq)
            #Join and update
            self.sequence[head]=''.join(l_seq)
            
        

def main():
    """
    Main function to be run for randomizedMotifSearch.py
    using the ConsensusCRISPRPromoter class.
    """
    
    #Parse arguments
    myParser=Parser()
    iterations=myParser.args.iterations
    pseudocounts=myParser.args.pseudocounts
    kmer=myParser.args.kmer
    
    scramble=myParser.args.random
    gibbs=myParser.args.gibbs
    mprint=myParser.args.mprint
    #Open file(STDIN)
    open_file=FastAreader('')
    seq=open_file.readFasta()
    seq_dict={}
    #Sequence are not big so we store them into a dictionary instead of reloading the file at each iteration
    #This is to deal with the fact that we can't loop over a generator several time
    for head, seq in seq:
        seq_dict[head]=seq

    #Run method -- arguments using argparse
    run=ConsensusCRISPRPromoter(sequence=seq_dict, 
                                iterations=iterations, 
                                pseudocounts=pseudocounts, 
                                motif_length=kmer, 
                                scramble=scramble, 
                                gibbs=gibbs, 
                                mprint=mprint)
    #Iterate method
    run.iterate_motif_search()
    #Build and print consensus motif
    run.build_consensus()

        
if __name__=='__main__':
    main()
        
        
                
        
        
