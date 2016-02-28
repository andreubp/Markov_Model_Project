#!/usr/local/bin/python3

                        ########################
                        #    M.M.D.I. v4       #
                        ########################

###############################################################################
#              version4: This is the last stable version of MMID.             #
# It includes a userfriendly paths to drive the function of the software.     #
# Moreover, it includes the possibility of generating a Markov Model from     #
# a CLIP interval file format and print, evaluate and use that Markov Model.  #
# It also accept a multi FASTA file as query to identify DNA sequences        #
# related with the identified Motif.                                          #
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
# PATH 1 - Building and evaluating the model from a CLIP file in interval format.
# 0) SYS ARGV
# 1) FROM a CLIP's BED (interval format) obatain a training and testing set.
# 2) FROM the training set separate the sequences in signal (+) and background (-)
# 3A) Generate the kmer dictionary from the observed kmers in the sequences.
# 3B) Generate a all possible kmer dictionary and generate pseudocounts on the
#       fly for the not finded kmers. (if pseudocounts)
# 4) Store the model in a file.
# 5A) Evaluate the model with separated datasets. if classic
# 5B) Evaluate the model considering the boundary regions. if boundaries.
# PATH 2 - Using a previously build Markov Model to identify domains in a query multifasta file
# 0) SYS ARGV
# 1) Read the FASTA file
# 2) Score and identify the motifs using a defined threshold.

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
                    help="Add pseudocounts in the generated MM")
parser.add_argument('-k','--order',dest='order',action='store',type=int,
                    default=3,
                    help='It sets the k-mer order for the Markov Model.')
parser.add_argument('-l','--windows',dest='wsize',action='store',type=int,
                    default=None,
                    help='It sets the windows size to analize the testing set.')
parser.add_argument('-t','--threshold',dest='threshold',action='store',type=float,
                    default=None,
                    help='It sets the threshold for the motif identification.')
parser.add_argument('-b','--boundaries',dest='boundaries',action='store_true',
                    default=False,
                    help='The evaluation of the model will consider the boundaries.')
parser.add_argument('-q','--query',dest='query',action='store',
                    default=None,
                    help='It sets the multi FASTA input file to test.')
parser.add_argument('--train',dest='train',action='store_true',
                    default=False,
                    help='It defines the path of the main script.')
parser.add_argument('--identify',dest='identify',action='store_true',
                    default=False,
                    help='It defines the path of the main script.')
parser.add_argument('-m','--model',dest='model',action='store',
                    default=None,
                    help='It defines the file that contains the Markov Model to be used.')
parser.add_argument('-w','--web_logo',dest='weblogo',action='store_true',
                    default=False,
                    help='It generate the top represented k-mers in the MM.')

args = parser.parse_args()

if args.train:
    ######## 1 TT SETS ######
    training_filename='training_'+args.infile
    testing_filename='testing_'+args.infile
    sys.stderr.write("\t\tMMID v.4\n\t\tAdvanced Genome Bioinformatics\n\t\tAndreu Bofill & David Mas\n\nThe input CLIP data is from {}\n".format(args.infile))
    sys.stderr.write("Arguments:\nMM Order:{}\nWindows Length:{}\n".format(args.order,args.wsize))
    if args.pseudocounts:
        sys.stderr.write("This run is using pseudocounts.\n\n")
    sys.stderr.write("First Step: Training and Testing file separation:\n\n")
    m.training_testing_sets_separation(args.infile,3,training_filename,testing_filename)
    sys.stderr.write("DONE!\n\n")
    ##### 2 and 3 : Generating the fg and bg dictionaries ####
    if (args.pseudocounts):
        sys.stderr.write("Second Step: Generating the Markov Models\n")
        back_dict=m.build_hash_pseudocount([x for x in m.background_separation(training_filename)],args.order)
        sys.stderr.write("The MM for the background DONE!\n")
        sign_dict=m.build_hash_pseudocount([x for x in m.signal_separation(training_filename)],args.order)
        sys.stderr.write("The MM for the foreground DONE!\n\n")
    else:
        sys.stderr.write("Second Step: Generating the Markov Models\n")
        back_dict=m.build_hash([x for x in m.background_separation(training_filename)],args.order)
        sys.stderr.write("The MM for the background DONE!\n")
        sign_dict=m.build_hash([x for x in m.signal_separation(training_filename)],args.order)
        sys.stderr.write("The MM for the foreground DONE!\n\n")
    #### 4. printing the MM in a file ####
    output_mm=args.outfile+'.MM'
    sys.stderr.write("Printing the model in {}\n".format(output_mm))
    m.print_hash(sign_dict,back_dict,output_mm)
    #### 5. Evaluation of the model ####
    sys.stderr.write("Third Step, Evaluation of the model\n")
    if (args.boundaries):
        sys.stderr.write("Boundaries evaluation have been selected.\nThe threshold range goes from -40 to 40\n")
        output_filename= args.outfile+'.data'
        output=open(output_filename, "w")
        list_t=[-40]+[x for x in range(-10,10,1)]+[40]
        for i in list_t:
            sys.stderr.write("{}\n".format(i))
            TPR=0
            FPR=0
            PPV=0
            TP=0
            FP=0
            TN=0
            FN=0
            for x in  m.none_separation(testing_filename):
                for y in m.identifier(i,x,args.wsize,args.order,back_dict,sign_dict):
                    TP+=y[0]
                    FP+=y[1]
                    TN+=y[2]
                    FN+=y[3]
            if (TP+FN) != 0:
                TPR=TP/(TP+FN)
            else:
                TPR=0
            if (FP+TN) != 0:
                FPR=FP/(FP+TN)
            else:
                FPR=0
            if (TP+FP) != 0:
                PPV=TP/(TP+FP)
            else:
                PPV=FN/(FN+TN)
            output.write("{0}\t{1}\t{2}\t{3}\n".format(i,TPR,FPR,PPV))
    else:
        output_filename_fg=args.outfile+'_fg.th'
        output_filename_bg=args.outfile+'_bg.th'
        out_fg=open(output_filename_fg, "w")
        out_bg=open(output_filename_bg, "w")
        sys.stderr.write("Printing the background scores at {}\n".format(output_filename_bg))
        for i in m.windows_score([x for x in m.background_separation(testing_filename)],args.order,args.wsize, back_dict,sign_dict):
            out_bg.write("{}\n".format(i))
        sys.stderr.write("Printing the foreground scores at {}\n".format(output_filename_fg))
        for i in m.windows_score([x for x in m.signal_separation(testing_filename)],args.order,args.wsize, back_dict,sign_dict):
            out_fg.write("{}\n".format(i))
    sys.stderr.write("To plot the evaluation result, use the R script attached.\n")
    sys.stderr.write("MMID have finished the model building step succesfully.\n")

if args.identify:
    if args.model == None:
        model_file=output_mm
    else:
        model_file=args.model
    sys.stderr.write("\t\tMMID v.4\n\t\tAdvanced Genome Bioinformatics\n\t\tAndreu Bofill & David Mas\n\nThe input Markov Model data is from {}\n".format(model_file))
    sys.stderr.write("Identifying motifs in the FASTA file {}\n".format(args.query))
    sys.stdout.write("Motif Identifier.\nDefined Threshold:{0}\tFasta file input:{1}\tMarkov Model input:{2}\n".format(args.threshold,args.query,model_file))
    bg_MM,sg_MM=m.read_MM(model_file)
    for i in m.FASTA_iterator(args.query):
        for j,x in m.windows_score_query(i, args.order, args.wsize, bg_MM, sg_MM,args.threshold):
            sys.stdout.write("Sequence:{0}\tScore:{1:.3f}\n".format(j,x))

if args.weblogo:
    if args.model == None:
        model_file=output_mm
    else:
        model_file=args.model
    sys.stdout.write("Input MSA to generate a web logo:\n")
    bg_MM,sg_MM=m.read_MM(model_file)
    m.generate_weblogo(sg_MM,bg_MM)
    sys.stdout.write("In order to get the actual logo go to weblogo.berkeley.edu and paste the generated k-mers.\n")

if not args.identify and not args.train and not args.weblogo:
    sys.stderr.write("This Script needs some paths to run. Use the arguemnt -h to get help.\n")
