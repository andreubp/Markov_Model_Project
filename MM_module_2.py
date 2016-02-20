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


def background_signal_separation(filename):
    """It separate from the sequence the background from the signal"""
    bed50_file = open(filename,"r")
    background = ()
    signal =()
    for line in bed50_file:
        line = line.strip().upper().split("\t")
        background += (line[10][:50],)
        background += (line[10][-50:],)
        signal += (line[10][50:-50],)
    return (background, signal)

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
    for i in generate_kmers(k):
        kmer_dict[i]=1
    return (kmer_dict)

def get_kmers(sequence,k):
    """From a sequence extract all the possible k-mers in that sequence"""
    kmers_list=[]
    length=len(sequence)
    if k<1 or length < k:
        sys.exit()
    i=0
    while i<length:
#        print(i)
        if len(sequence[i:i+k])==k:
            kmers_list.append(sequence[i:i+k])
        i=i+1
    return kmers_list

### For simple iterators, yield from iterable is essentially
### just a shortened form of for item in iterable: yield item:

def build_hash(seq_list,k):
    """From a list of sequences it builds a dictionary without pseudocounts"""
    total_length=0
    kmer_dict = {}
    for i in seq_list:
        kmers=get_kmers(i,k)
        for j in kmers:
            if j in kmer_dict: #dictionary_kmers(2) ?!
                kmer_dict[j]= i.count(j) + kmer_dict[j]
            else:
                kmer_dict[j]=i.count(j)
    for z in kmer_dict:
        total_length=total_length+kmer_dict[z]
    for x in kmer_dict:
        kmer_dict[x]=kmer_dict[x]/total_length
    return kmer_dict

def build_hash_pseudocount(seq_list,k):
    """From a list of sequences it builds a dictionary with pseudocounts"""
    kmer_dict=dictionary_kmers(k)
    total_length=4**(k)
    for i in seq_list:
        kmers=get_kmers(i,k)
        for j in kmers:
            kmer_dict[j]= i.count(j) + kmer_dict[j]
    for z in kmer_dict:
        total_length=total_length+kmer_dict[z]-1
    for x in kmer_dict:
        kmer_dict[x]=kmer_dict[x]/total_length
    return kmer_dict


def print_hash(dict_signal,dict_background,outfilename):
    """This function prints the MM in a output file"""
    count=0
    out_file_name=open(outfilename, "w")
    out_file_name.write("This file contains the information of the MM used in MMDI:\n")
    out_file_name.write("\tBackground respect to signal\n")
    out_file_name.write("Background a:\t\tSignal a:\n")
    for i in dict_background:
        out_file_name.write("{0} = {1:.3f}\t\t".format(i,dict_background[i]))
        count=count+dict_background[i]
        if i in dict_signal:
            out_file_name.write("{0} = {1:.3f}\n".format(i,dict_signal[i]))
        else:
            out_file_name.write(i+" = NA\n")
    out_file_name.write("\tThe sum of the at from Background is: {0}\n".format(count))
    count=0
    out_file_name.write("\tSignal respect to background\n")
    out_file_name.write("Signal a:\t\tBackground a:\n")
    for i in dict_signal:
        out_file_name.write("{0} = {1:.3f}\t\t".format(i,dict_signal[i]))
        count=count+dict_background[i]
        if i in dict_background:
            out_file_name.write("{0} = {1:.3f}\n".format(i,dict_background[i]))
        else:
            out_file_name.write(i+" = NA\n")
    out_file_name.write("\tThe sum of the at from Signal is: {0}\n".format(count))

def windows(data, k, l, back_dict, sign_dict):
    for sentence in data:
        L = len(sentence)
        n = 0

        if (L-l < k):
            if (l >= L):
                raise ValueError("The length of the windows: %s chosed is longer thant the sentence of length: %s" %(l,L))
            else:
                raise ValueError("The length of the windows: %s in the sentence of length: %s  is too long for the MMorder: %s " %(l,L,k))
        print ("\n")
        while (n < L-(l-1)):
            m=0
            windows = sentence[n:n+l]
            print (sentence)
            print (windows)
            total_score = 0
            while m < (len(windows)-(k-1)):
                print (n,m)
                print(windows[m:m+k])
                total_score += math.log((sign_dict[(windows[m:m+k])])/back_dict[(windows[m:m+k])])
                m +=1
            yield (total_score)
            n +=1
    #print (total_score)



if __name__ == '__main__':
    (background, signal) = background_signal_separation("example.bed")
    print("BACKGROUND KMER_DICT \n",build_hash(background,3))
    print ("\n\nSIGNAL KMER_DICT \n", build_hash(signal,3))
