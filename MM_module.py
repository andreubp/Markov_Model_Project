#!/usr/local/bin/python3
import sys
import os
import math

def add_50(filename):
    bed_file = open(filename,'r')
    for line in bed_file:
        line = line.split("\t")
        line[1] = int(line[1]) - 50
        line[2] = int(line[2]) + 50
        print (line)

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
    """It separate from the sequence the background from the signal"""
    bed50_file = open(filename,"r")
    for line in bed50_file:
        line = line.strip().upper().split("\t")
        yield line[10][50:-50]
    bed50_file.close()

def none_separation(filename):
    """It separate from the sequence the background from the signal"""
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
    example:
        get_kmers(3)
    """
    if k==0:
        yield y
    else:
        for m in ['A','C','T','G']:
            kmer=y+m
            for item in generate_kmers(k-1,kmer):
                yield item


def dictionary_kmers(k):
    """It generates a dictionary for all possible kmers"""
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


### For simple iterators, yield from iterable is essentially
### just a shortened form of for item in iterable: yield item:

def build_hash(seq_list,k):
    """From a list of sequences it builds a dictionary without pseudocounts"""
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
    """From a list of sequences it builds a dictionary with pseudocounts"""
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
    """This function prints the MM in a output file"""
    count=0
    out_file_name=open(outfilename, "w")
    out_file_name.write("This file contains the information of the MM generated by in MMDI:\n")
    out_file_name.write("\tBackground table\n")
    out_file_name.write("bg\tA\tC\tG\tT\n")
    for i in dict_background:
        order=len(i)
        break
    for j in generate_kmers(order-1):
        out_file_name.write("{}\t".format(j))
        for u in  ['A','C','G','T']:
            if j+u in dict_background:
                out_file_name.write("{:.3f}\t".format(dict_background[j+u]))
            else:
                out_file_name.write("NA\t")
        out_file_name.write("\n")
    out_file_name.write("\tSignal table\n")
    out_file_name.write("fg\tA\tC\tG\tT\n")
    for j in generate_kmers(order-1):
        out_file_name.write("{}\t".format(j))
        for u in  ['A','C','G','T']:
            if j+u in dict_signal:
                out_file_name.write("{:.3f}\t".format(dict_signal[j+u]))
            else:
                out_file_name.write("NA\t")
        out_file_name.write("\n")

def windows_score(data, k, l, back_dict, sign_dict):
    for sentence in data:
        L = len(sentence)
        n = 0
        if (L-l < k):
            if (l >= L):
                #raise ValueError("The length of the windows: %s chosed is longer thant the sentence of length: %s" %(l,L))
                next
            else:
                #raise ValueError("The length of the windows: %s in the sentence of length: %s  is too long for the MMorder: %s " %(l,L,k))
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

if __name__ == '__main__':
    print("hello")
