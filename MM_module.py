#!/usr/local/bin/python3
"""
This is the documentation for MM_module. The functions below are used in the
execution of the different paths of MMID. The main dependencies of these module
are the sys, os and math packages but all of them should be present in a correct
python instalation.
This modeule only cotain documented functions used in the main script.
The authors of this module are Andreu Bofill and David Mas.
"""

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #       David Mas      #
                        ########################


########################
#        Modules       #
########################

import sys
import os
import math


def training_testing_sets_separation(filename, ratio, training_file, testing_file):
    """It separete the data set in training and testing"""
    bed50_file = open(filename,"r")
    training = open(training_file, 'w')
    testing = open(testing_file, 'w')
    count = 0
    for line in bed50_file:
        count += 1
        if count % ratio == 0:
            testing.write(line)
        else:
            training.write(line)
    bed50_file.close()
    training.close()
    testing.close()


def background_separation(filename):
    """It separate from the sequence the background from the signal"""
    bed50_file = open(filename,"r")
    for line in bed50_file:
        line = line.strip().upper().split("\t")
        yield line[10][:50]
        yield line[10][-50:]
    bed50_file.close()


def signal_separation(filename):
    """It separate from the sequence the signal from the background"""
    bed50_file = open(filename,"r")
    for line in bed50_file:
        line = line.strip().upper().split("\t")
        yield line[10][50:-50]
    bed50_file.close()

def none_separation(filename):
    """It yields sequences from a Interval file format"""
    bed50_file = open(filename,"r")
    for line in bed50_file:
        line = line.strip().upper().split("\t")
        yield line[10]
    bed50_file.close()


def generate_kmers(k,y=''):
    """
    This function is a generator of kmers for DNA bases.
    It uses a recursive aproach.
    The only required argument is the order (k).
    it also accepts prefix defined as y='prefix'
    the default value for y is ''.
    """
    if k==0:
        yield y
    else:
        for m in ['A','C','T','G']:
            kmer=y+m
            for item in generate_kmers(k-1,kmer):
                yield item


def dictionary_kmers(k):
    """
    It generates a dictionary for all possible kmers using the
    function generate_kmers. The value for all k-mers is 1 because is used
    to add pseudocounts
    """
    kmer_dict = {}
    sys.stderr.write("Gnerating all possible kmers\n")
    for i in generate_kmers(k):
        kmer_dict[i]=1
    sys.stderr.write("returning that initial dictionary\n")
    return kmer_dict

def get_kmers(sequence,k):
    """From a sequence extract all the possible k-mers in that sequence"""
    kmers_list=[]
    length=len(sequence)
    i=0
    while i<length:
        if len(sequence[i:i+k])==k:
            yield sequence[i:i+k]
        i=i+1

def build_hash(seq_list,k):
    """
    From a list of sequences it builds a dictionary of the k-mers that have
    appeared in the whole list. It generates  a empty dictionary and then add
    the k-mers when they start appearing in the sequences.
    This function uses get_kmers function to obtain the kmers represented
    in a sequence.
    """
    kmer_dict = {}
    total_length=0
    for i in seq_list:
        for j in get_kmers(i,k):
            total_length=total_length + 1
            if j in kmer_dict:
                kmer_dict[j]= 1 + kmer_dict[j]
            else:
                kmer_dict[j]= 1
    for x in kmer_dict:
        kmer_dict[x]=kmer_dict[x]/total_length
    return kmer_dict

def build_hash_pseudocount(seq_list,k):
    """
    From a list of sequences it builds a dictionary of the k-mers that have
    appeared in the whole list. It generates  a pseudocounts
    dictionary and then add the k-mers when they start appearing in the sequences.
    This function uses 2 more functions, dictionary_kmers to start the dictionary
    and get_kmers to obtain the kmers represented in a sequence.
    """
    length=len(seq_list)
    sys.stderr.write("Building hash from {}\n".format(length))
    kmer_dict=dictionary_kmers(k)

    total_length=0
    for i in seq_list:
        for j in get_kmers(i,k):
            kmer_dict[j]= 1 + kmer_dict[j]
            total_length=total_length + 1
    for x in kmer_dict:
        kmer_dict[x]=kmer_dict[x]/total_length
    return kmer_dict

def print_hash(dict_signal,dict_background,outfilename):
    """
    This function prints the MM in a output file. The format is a simple
    tabular format. The k-mer in the left and 7 digits of the probability after
    a tab. In some k-orders this would be too few decimal points.
    """
    count=0
    out_file_name=open(outfilename, "w")
    out_file_name.write("This file contains the information of the MM generated by in MMDI:\n")
    out_file_name.write("Background table\n")
    for i in dict_background:
        out_file_name.write("{}\t{:.7f}\n".format(i,dict_background[i]))
    out_file_name.write("Signal table\n")
    for j in dict_signal:
        out_file_name.write("{}\t{:.7f}\n".format(j,dict_signal[j]))

def windows_score(data, k, l, back_dict, sign_dict):
    """
    This function gets as input a list of sequences, the markov model represented
    in bg and fg dictionaries, the k-mer order and the window size. It yields the
    score of the top windows for the entire sequence for each sequence in the list.
    """
    for sentence in data:
        L = len(sentence)
        n = 0
        if (L-l < k):
            if (l >= L):
                next
            else:
                next
        top_windows = None
        bot_windows = None
        while (n < L-(l-1)):
            m=0
            windows = sentence[n:n+l]
            total_score = 0
            while m < (len(windows)-(k-1)):
                total_score += math.log((sign_dict[(windows[m:m+k])])/back_dict[(windows[m:m+k])])
                m +=1
            if top_windows is None:
                top_windows = total_score
            elif top_windows < total_score:
                top_windows = total_score

            n +=1
        if top_windows is None:
            pass
        else:
            yield top_windows
def windows_score_query(sequence, k, l, back_dict, sign_dict,t):
    """
    This function gets as input a list of sequences, the markov model represented
    in bg and fg dictionaries, the k-mer order, the window size and the threashold
    defined by the user. It yields the score and the sequence that have generated
    the hit. It is used in the query evaluation.
    """
    L = len(sequence)
    n = 0 #window iteration
    if (L-l < k):
        if (l >= L):
            next
        else:
            next
    while (n < L-(l-1)):
        total_score=0
        m=0 #kmer iteration
        windows = sequence[n:n+l]
        while m < (len(windows)-(k-1)):
            total_score += math.log((sign_dict[(windows[m:m+k])])/back_dict[(windows[m:m+k])])
            m +=1
        if total_score>t:
            yield windows,total_score
        n +=1

def identifier(t,sequence,l,k,back_dict,sign_dict):
    """
    This function takes as an argument a specified threshold, a sequence, the
    model dictionaries, the k-mer order and the window size desired by the user.
    Independently, it first scores the window by its differents kmers and sets
    a veredict. After that, another part of the function turn that values in FP
    and TP, if the veredict is positive, based on the positions of the window
    and in FN or TN if the veredict is negative.
    """
    L = len(sequence)
    x=50
    y=L-50
    n = 0
    TP=0
    FP=0
    TN=0
    FN=0
    if (L-l < k):
        if (l >= L):
            next
        else:
            next
    while (n < L-(l-1)):
        m=0
        windows = sequence[n:n+l]
        total_score = 0
        while m < (len(windows)-(k-1)):
            total_score += math.log((sign_dict[(windows[m:m+k])])/back_dict[(windows[m:m+k])])
            m +=1
        if total_score >= t:
            for i in range(n,n+l):
                if i >= x and i <= y:
                    TP=TP+1
                else:
                    FP=FP+1
        else:
            for i in range(n,n+l):
                if i < x or i > y:
                    TN=TN+1
                else:
                    FN=FN+1
        n +=1
    yield (TP,FP,TN,FN)

def generate_weblogo(dict_MM_fg,dict_MM_bg):
    """
    From a Markov Model this function prints in the Standard Output the most
    represented elements in the Signal dictionary. Basicaly, it sorts the
    dictionary by the values and prints the sequences. The web logo can be
    generated in weblogo.berkeley.edu"""
    length=len(dict_MM_fg)
    n_seq=0.01*length
    count=0
    for w in sorted(dict_MM_fg, key=dict_MM_fg.get, reverse=True):
        count+=1
        sys.stdout.write("{}\n".format(w))
        if n_seq<count:
            break

def read_MM(file_name):
    """
    This function takes as an argument an input file with the information from
    a Markov Model generated with MMID in a previous run. It is recomended for
    large order models.
    """
    MM_file=open(file_name, "r")
    next(MM_file)
    next(MM_file)
    back_dict={}
    sign_dict={}
    for line in MM_file:
        line=line.strip()
        if line=='Signal table':
            break
        else:
            fields=line.split()
            back_dict[fields[0]]=float(fields[1])
    for line in MM_file:
        line=line.strip()
        fields=line.split()
        sign_dict[fields[0]]=float(fields[1])
    return back_dict,sign_dict

def FASTA_iterator(fasta_filename):
    """Reads all the fasta sequences from a file and yields them."""
    fh=open(fasta_filename,"r")
    current_sequence=[]
    for line in fh:
        line=line.strip()
        if line.startswith(">"):
            if current_sequence:
                sequence="".join(current_sequence)
                yield sequence
            identifier=line
            current_sequence=[]
        else:
            current_sequence.append(line)
    sequence="".join(current_sequence)
    yield sequence
    fh.close()
