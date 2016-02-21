#!/usr/local/bin/python3

########################
#    M.M.D.I. v2       #
########################

###############################################################################
#              version2: This version is including the production             #
# of the scores for the testing background and foreground sets to             #
# generate the ROC and PPV curves.                                            #
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
import MM_module as m
import argparse


####### INDEX ########
# 0) SYS ARGV
# 1) FROM a CLIP's BED (interval format) obatain a training and testing set.
# 2) FROM the training set separate the sequences in signal (+) and background (-)
# 3A) Generate the kmer dictionary from the observed kmers in the sequences.
# 3B) Generate a all possible kmer dictionary and generate pseudocounts on the
#       fly for the not finded kmers. (if pseudocounts)
# 4) Store the model in a file.
# 5) Produce the scores without threshold for plotting



######## 0. ARG ########

parser = argparse.ArgumentParser(description='MMDI')
parser.add_argument('-i','--input',dest="infile", action="store", default= None,
                    help="Input Interval formated file")
parser.add_argument('-o','--MM_output', dest='outfile', action='store',
                    default=None,
                    help='The root output file where the MM and score values is going to be written.')
parser.add_argument('-p','--pseudocounts', dest='pseudocounts',
                    action='store_true',
                    default=False,
                    help="Add pseudocounts in the generated MM"
                    )
parser.add_argument('-k','--order',dest='order',action='store',type=int,
                    default=3,
                    help='It sets the order of the Markov Model used.')


parser.add_argument('-l','--windows',dest='wsize',action='store',type=int,
                    default=None,
                    help='It sets the windows size to analize the testing set.')


args = parser.parse_args()

######## 1 TT SETS ######
training_filename='training_'+args.infile
testing_filename='testing_'+args.infile
sys.stderr.write("Starting the Log:\nThe input data is from {}\n".format(args.infile))
sys.stderr.write("Arguments:\nMM Order:{}\nWindows Length:{}\n".format(args.order,args.wsize))
if args.pseudocounts:
    sys.stderr.write("This run is using pseudocounts.\n\n\n")
sys.stderr.write("Starting the training and testing separation...\n\n")
m.training_testing_sets_separation(args.infile,3,training_filename,testing_filename)

######## 2 +/- ###########
sys.stderr.write("Starting the background and foreground separation...\n\n")
(background, signal)=m.background_signal_separation(training_filename)
(background_testing,signal_testing)=m.background_signal_separation(testing_filename)

######## 3 kmer dict #####
sys.stderr.write("Generating the Markov Model...\n\n")
if (args.pseudocounts):
    back_dict=m.build_hash_pseudocount(background,args.order)
    sign_dict=m.build_hash_pseudocount(signal,args.order)
else:
    back_dict=m.build_hash(background,args.order)
    sign_dict=m.build_hash(signal,args.order)

######### 4 print ########
output_filename=args.outfile+'.MM'
sys.stderr.write("Printing the model in {}...\n\n".format(output_filename))
m.print_hash(sign_dict,back_dict,output_filename)

######## 5 windows #######
sys.stderr.write("Starting the computation of the Scores...\n\n")
output_filename_fg=args.outfile+'_fg.th'
output_filename_bg=args.outfile+'_bg.th'
out_fg=open(output_filename_fg, "w")
out_bg=open(output_filename_bg, "w")
sys.stderr.write("Printing the background scores at {}...\n\n".format(output_filename_bg))
for i in m.windows_score(background_testing,args.order,args.wsize, back_dict,sign_dict):
    out_bg.write("{}\n".format(i))
sys.stderr.write("Printing the foreground scores at {}...\n\n".format(output_filename_fg))
for i in m.windows_score(signal_testing,args.order,args.wsize, back_dict,sign_dict):
    out_fg.write("{}\n".format(i))
