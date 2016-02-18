#!/usr/local/bin/python3
import sys
import os

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
        kmer_dict[i]=0
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


if __name__ == '__main__':
    (background, signal) = background_signal_separation("example.bed")
    print("BACKGROUND KMER_DICT \n",build_hash(background,3))
    print ("\n\nSIGNAL KMER_DICT \n", build_hash(signal,3))
