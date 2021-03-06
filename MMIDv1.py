#!/usr/local/bin/python3

########################
#      M.M.D.I. v1     #
########################

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
import MM_module as m
import argparse


####### INDEX ########
# 0) SYS ARGV
# 1) FROM a CLIP's BED (interval format) obatain a training and testing set.
# 2) FROM the training set separate the sequences in signal (+) and background (-)
# 3) Generate the kmer dictionary from the observed kmers in the sequences.
# 4) Store the output in a file.


######## 0. ARG ########

parser = argparse.ArgumentParser(description='MMDI')
parser.add_argument('-i','--input',dest="infile", action="store", default= None,
                    help="Input Interval formated file")
parser.add_argument('-o','--MM_output', dest='outfile', action='store',
                    default=None,
                    help='The output file where the MM is going to be written.')

args = parser.parse_args()


######## 1 TT SETS ######
training_filename='training_'+args.infile
testing_filename='testing_'+args.infile
m.training_testing_sets_separation(args.infile,3,training_filename,testing_filename)

######## 2 +/- ###########
(background, signal)=m.background_signal_separation(training_filename)

######## 3 kmer dict #####

back_dict=m.build_hash(background,3)
sign_dit=m.build_hash(signal,3)

######### 4 print ########

m.print_hash(sign_dit,back_dict,args.outfile)

sys.exit()
