#!/usr/local/bin/python3
import sys




def add_50(filename):
    bed_file = open(filename,'r')
    for line in bed_file:
        line = line.split("\t")
        line[1] = int(line[1]) - 50
        line[2] = int(line[2]) + 50
        print (line)

def training_testing_sets_separation(filename, ratio, training_file, testing_file):
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
    bed50_file = open(filename,"r")
    background = ()
    signal =()
    for line in bed50_file:
        line = line.strip().split("\t")
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
            #yield from get_kmers(k-1,kmer)
            for item in get_kmers(k-1,kmer):
                yield item

def get_kmers(sequence,k):
    """jhsdj"""
    kmers_list=[]
    length=len(sequence)
    print(length)
    if k<1 or length < k:
        sys.exit()
    i=0
    while i<length:
        print(i)
        if len(sequence[i:i+k])==k:
            kmers_list.append(sequence[i:i+k])
        i=i+1
    return kmers_list



if __name__ == '__main__':
	# print(get_kmers(3))
	# for i in get_kmers(3):
	# 	print(i)
    # print(get_kmers('AGTCGATGC',3))
    #print()



### For simple iterators, yield from iterable is essentially
### just a shortened form of for item in iterable: yield item:
