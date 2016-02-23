#!/usr/local/bin/python3

########################
#    M.M.D.I. v2.1     #
########################

###############################################################################
#              version2.1: This version is including the production           #
# of the scores for the testing background and foreground sets to             #
# generate the ROC and PPV curves.  the 2.1 uses iterations to get            #
# multiple values.                                                            #
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
                    default=True,
                    help="Add pseudocounts in the generated MM"
                    )

args = parser.parse_args()

k_list=[2,3,4,5,6,7,8]
l_list=[15,20,25,30,35,40,45]
for k in k_list:
    for l in l_list:
        if k<l:
            sys.stderr.write("The log begins, k-order = {} and windows length = {}".format(k,l))
            training_filename='training_'+args.infile
            testing_filename='testing_'+args.infile
            sys.stderr.write("Starting the Log:\nThe input data is from {}\n".format(args.infile))
            sys.stderr.write("Arguments:\nMM Order:{}\nWindows Length:{}\n".format(k,l))
            if args.pseudocounts:
                sys.stderr.write("This run is using pseudocounts.\n\n\n")
            sys.stderr.write("Starting the training and testing separation...\n\n")
            m.training_testing_sets_separation(args.infile,3,training_filename,testing_filename)

            ######## 2 +/- ###########


            ######## 3 kmer dict #####


            if (args.pseudocounts):
                sys.stderr.write("Starting the training background and foreground separation...\n\n")
                back_dict=m.build_hash_pseudocount([x for x in m.background_separation(training_filename)],k)
                sys.stderr.write("The MM for the bg is generated...\n\n")
                sign_dict=m.build_hash_pseudocount([x for x in m.signal_separation(training_filename)],k)
                sys.stderr.write("The MM for the bg is generated...\n\n")
            else:
                back_dict=m.build_hash([x for x in m.background_separation(training_filename)],k)
                sign_dict=m.build_hash([x for x in m.signal_separation(training_filename)],k)

            ######### 4 print ########
            output_filename=args.outfile+"_k_"+str(k)+"_l_"+str(l)+'.MM'
            sys.stderr.write("Printing the model in {}...\n\n".format(output_filename))
            m.print_hash(sign_dict,back_dict,output_filename)

            ######## 5 windows #######
            sys.stderr.write("Starting the testing background and foreground separation...\n\n")

            sys.stderr.write("Starting the computation of the Scores...\n\n")
            output_filename_fg=args.outfile+"_k_"+str(k)+"_l_"+str(l)+'_fg.th'
            output_filename_bg=args.outfile+"_k_"+str(k)+"_l_"+str(l)+'_bg.th'
            out_fg=open(output_filename_fg, "w")
            out_bg=open(output_filename_bg, "w")
            sys.stderr.write("Printing the background scores at {}...\n\n".format(output_filename_bg))
            for i in m.windows_score([x for x in m.background_separation(testing_filename)],k,l, back_dict,sign_dict):
                out_bg.write("{}\n".format(i))
            sys.stderr.write("Printing the foreground scores at {}...\n\n".format(output_filename_fg))
            for i in m.windows_score([x for x in m.signal_separation(testing_filename)],k,l, back_dict,sign_dict):
                out_fg.write("{}\n".format(i))
        else:
            pass

######## 1 TT SETS ######
