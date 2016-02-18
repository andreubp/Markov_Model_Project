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
count=0
out_file_name=open(args.outfile, "w")
out_file_name.write("This file contains the information of the MM used in MMDI:\n")
out_file_name.write("\tBackground respect to signal\n")
out_file_name.write("Background a:\t\tSignal a:\n")
for i in back_dict:
    out_file_name.write("{0} = {1:.2f}\t\t".format(i,back_dict[i]))
    count=count+back_dict[i]
    if i in sign_dit:
        out_file_name.write("{0} = {1:.2f}\n".format(i,sign_dit[i]))
    else:
        out_file_name.write(i+" = NA\n")
out_file_name.write("\tThe sum of the at from Background is: {0}\n".format(count))
count=0
out_file_name.write("\tSignal respect to background\n")
out_file_name.write("Signal a:\t\tBackground a:\n")
for i in sign_dit:
    out_file_name.write("{0} = {1:.2f}\t\t".format(i,sign_dit[i]))
    count=count+back_dict[i]
    if i in back_dict:
        out_file_name.write("{0} = {1:.2f}\n".format(i,back_dict[i]))
    else:
        out_file_name.write(i+" = NA\n")
out_file_name.write("\tThe sum of the at from Signal is: {0}\n".format(count))
sys.exit()
