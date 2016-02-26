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

"""
Here we define the Documentation of the main Script
"""

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
count=0
k_list=[2,4,6,8,10,12]
k_list_david=[1,3,5,7,9,11]
l_list=[20,25,30,35,40,45,50]
## run

### python3 MMIDv3.1.py -i ENCF... -o results/v3.1 -p
for k in k_list:
    for l in l_list:
        count=1+count
        if k<l:
            sys.stderr.write("The log begins, k-order = {} and windows length = {}\n".format(k,l))
            sys.stderr.write("We are in the {} of the simulations\n".format((count//54)*100))
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

            output_filename= args.outfile+"_k_"+str(k)+"_l_"+str(l)+'.data'

            output=open(output_filename, "w")

            list_t=[-40]+[x for x in range(-10,10,1)]+[40]
            for i in list_t:
                TPR=0
                FPR=0
                PPV=0
                TP=0
                FP=0
                TN=0
                FN=0
                for x in  m.none_separation(testing_filename):
                    for y in m.identifier(i,x,l,k,back_dict,sign_dict):
                        TP+=y[0]
                        FP+=y[1]
                        TN+=y[2]
                        FN+=y[3]
                TPR=TP/(TP+FN)
                FPR=FP/(FP+TN)
                PPV=TP/(TP+FP)
                output.write("{0}\t{1}\t{2}\t{3}\n".format(i,TPR,FPR,PPV))
            #for i in m.windows_score([x for x in m.none_separation(testing_filename)], args.order, args.wsize, back_dict, sign_dict):
            #    output_all_top.write("{}\n".format(i))
        else:
            pass

######## 1 TT SETS ######
