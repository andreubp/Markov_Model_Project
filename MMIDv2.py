#!/usr/local/bin/python3

########################
#    M.M.D.I. v2       #
########################

###############################################################################
#                                                                             #
###############################################################################


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
import MM_module_2 as m
import argparse


####### INDEX ########
# 0) SYS ARGV
# 1) FROM a CLIP's BED (interval format) obatain a training and testing set.
# 2) FROM the training set separate the sequences in signal (+) and background (-)
# 3A) Generate the kmer dictionary from the observed kmers in the sequences.
# 3B) Generate a all possible kmer dictionary and generate pseudocounts on the
#       fly for the not finded kmers. (if pseudocounts)
# 4) Store the output in a file.



######## 0. ARG ########

parser = argparse.ArgumentParser(description='MMDI')
parser.add_argument('-i','--input',dest="infile", action="store", default= None,
                    help="Input Interval formated file")
parser.add_argument('-o','--MM_output', dest='outfile', action='store',
                    default=None,
                    help='The output file where the MM is going to be written.')
parser.add_argument('-p','--pseudocounts', dest='pseudocounts',
                    action='store_true',
                    default=False,
                    help="Add pseudocounts in the generated MM"
                    )
parser.add_argument('-k','--order',dest='order',action='store',type=int,
                    default=3,
                    help='It sets the order of the Markov Model used.')

##############
parser.add_argument('-l','--windows',dest='wsize',action='store',type=int,
                    default=None,
                    help='It sets the windows size to analize the testing set.')
##############

args = parser.parse_args()

######## 1 TT SETS ######
training_filename='training_'+args.infile
testing_filename='testing_'+args.infile
m.training_testing_sets_separation(args.infile,3,training_filename,testing_filename)

######## 2 +/- ###########
(background, signal)=m.background_signal_separation(training_filename)
(background_testing,signal_testing)=m.background_signal_separation(testing_filename)

######## 3 kmer dict #####
if (args.pseudocounts):
    back_dict=m.build_hash_pseudocount(background,args.order)
    sign_dict=m.build_hash_pseudocount(signal,args.order)
else:
    back_dict=m.build_hash(background,args.order)
    sign_dict=m.build_hash(signal,args.order)

######### 4 print ########

m.print_hash(sign_dict,back_dict,args.outfile)

######## 5 windows #######
print(m.windows(background_testing,args.order,args.wsize, back_dict,sign_dict))
print(m.windows(signal_testing,args.order,args.wsize, back_dict,sign_dict))

sys.exit()
