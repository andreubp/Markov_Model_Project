
                    ################################
                    #      Plotting Rscript        #
                    ################################

  ###########################################################################
  # To plot the results from the evaluation process the data-sets must be   #
  # previously executed by MMID. This Script is an example of how to plot   #
  # the evaluation results from the MMID --train option. File names used    #
  # in this scripts are examples and must be changed to the actual files.   #
  # Use this guide as a recomendation to plot the resuts.                   #
  ###########################################################################

                        ########################
                        #        Authors:      #
                        #     Andreu Bofill    #
                        #       David Mas      #
                        ########################

################################
#      IMPORTING LIBRARIES     #
################################

library(PRROC)
library(ggplot2)
setwd("/Users/davidmasp/Documents/results/")


#########################################
#      Separated datasets evaluation    #
#########################################


## STEP 1 Load the .th files

filename_fg="filenname_fg.th"
filename_bg="filenname_bg.th"

fg_data=read.table(filename_fg)
bg_data=read.table(filename_bg)

fg=sample(fg_data$V1,500)
bg=sample(bg_data$V1,500)

# STEP 2 Generate the histograms from both score datasets.

dat=data.frame(scores=c(fg,bg),yy=(c(rep("Foreground",500),rep("Background",500))))
print(ggplot(dat, aes(x=scores)) + ggtitle(paste("k=",k,"L=",l,sep=" ")) + geom_histogram(data=subset(dat,yy == 'Foreground'),fill = "red", alpha = 0.5, binwidth = 2) + geom_histogram(data=subset(dat,yy == 'Background'),fill = "blue", alpha = 0.5, binwidth = 2))

# STEP 3 computing the roc and pr curve using PRROC package.

roc=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
pr,pr.curve(scores.class0 = fg,scores.class1 = bg, curve = TRUE)

# STEP 4 Ploting both curves:

ggplot(data.frame(roc$curve),aes(x=X1,y=X2,color=X3)) + geom_line(size=1.5) + labs(x="Sensitivity",y="FPR",title=format(roc$auc,digits=3), colour="Threshold") + scale_colour_gradient2(low="red", mid="orange",high="yellow")
ggplot(data.frame(pr$curve),aes(x=X1,y=X2,color=X3)) + geom_line(size=1.5) + labs(x="Recall",y="Precision",title=format(pr$auc.integral,digits=3), colour="Threshold") + scale_colour_gradient2(low="red", mid="orange",high="yellow") + ylim(0,1)


#########################################
#         Boundaries Evaluation         #
#########################################

## STEP 1 Load the .data files

data=read.table("file.data")

## STEP 2 Compute the AUC computation

AUC=sum(data$V2)/length(data$V2)
AUC

# STEP 3 Ploting both curves:

ggplot(data,aes(y=V2,x=V3,colour = V1)) + geom_line(size=1.5) + geom_point() + ggtitle("K=9 l=35, AUC=0,7") +ylim(0, 1)+xlim(0, 1)
ggplot(data,aes(y=V4,x=V2,colour= V1)) + geom_point(size=5) + geom_line(size=1.5) + ylim(0, 1) + ggtitle("K=9 l=25 - AUC=")

##########################################################
#         k-mer order and window size evaluation         #
##########################################################

# STEP 0 Generating the data sets:

# Before ploting this results the diferent data sets have to be generated.

# STEP 1A Loading the results (Boundaries Evaluation):

AUC_vect=c()
k_vect=c()
l_vect=c()
k_list=c(3,5,7,9) #this should change depending on the generated datasets
l_list=c(15,25,35,45) #this should change depending on the generated datasets
for (k in k_list){
  for (l in l_list){
    #### here we open the results file from the python script.
    filename=paste("filename_k_",k,"_l_",l,'.data',sep="")
    data=read.table(filename)
    AUC=sum(data$V2)/length(data$V2)
    AUC_vect=c(AUC_vect,AUC)
    k_vect=c(k_vect,k)
    l_vect=c(l_vect,l)
  }
}
data_k_l=data.frame(k_vect,l_vect,AUC_vect)

# STEP 2A Ploting acuracy in function of k-mer order and l window size (Boundaries Evaluation):

ggplot(data_k_l) + geom_tile(aes(k_vect, l_vect, fill = AUC_vect)) + scale_fill_distiller(palette ="Spectral")
qplot(l_vect,AUC_vect,data=data_k_l) + geom_point() + geom_smooth(method = "lm",formula=y ~ poly(x, 2)) + ggtitle("AUC vs l-window")
qplot(k_vect,AUC_vect,data=data_k_l) + geom_point() + geom_smooth(method = "lm",formula=y ~ poly(x, 2)) + ggtitle("AUC vs K")

# STEP 1B Loading the results (Diferent Datasets Evaluation):

k_list=c(2,4,6,8,10,12)
l_list=c(15,20,25,30,35,40,45,50)
for (k in k_list){
  for (l in l_list){
    #### here we open the results file from the python script.
    filename_fg=paste("v2.1_k_",k,"_l_",l,'_fg.th',sep="")
    filename_bg=paste("v2.1_k_",k,"_l_",l,'_bg.th',sep="")
    fg_data=read.table(filename_fg)
    bg_data=read.table(filename_bg)
    #### here we take a sample of our scores to get a equal amount of background and foreground reads.
    #### Moreover, the script gets faster without loosing shape.
    fg=sample(fg_data$V1,500)
    bg=sample(bg_data$V1,500)
    #### Here we plot the histograms for each condition of the model.
    dat=data.frame(scores=c(fg,bg),yy=(c(rep("Foreground",500),rep("Background",500))))
    print(ggplot(dat, aes(x=scores)) + ggtitle(paste("k=",k,"L=",l,sep=" ")) + geom_histogram(data=subset(dat,yy == 'Foreground'),fill = "red", alpha = 0.5, binwidth = 2) + geom_histogram(data=subset(dat,yy == 'Background'),fill = "blue", alpha = 0.5, binwidth = 2))
    #### Finally we just compute the ROC and pressision curve using the PRROC package, available at CRAN.
    rocname=paste("roc_k_",k,"_l_",l,sep="")
    prname=paste("pr_k_",k,"_l_",l,sep="")
    assign(rocname,roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE))
    assign(prname,pr.curve(scores.class0 = fg,scores.class1 = bg, curve = TRUE))
  }
}

# STEP 2B Ploting roc and pr curves in terms of k-mer order.

# the datasets have to be generated a priori

plot(roc_k_2_l_45, color=1, auc.main = FALSE)
plot(roc_k_4_l_45, add = TRUE, color = 2)
plot(roc_k_6_l_45, add = TRUE, color = 3)
plot(roc_k_8_l_45, add = TRUE, color = 4)
plot(roc_k_10_l_45, add = TRUE, color = 5)
plot(roc_k_12_l_45, add = TRUE, color = 6)
legend("bottomright",c("k=2","k=4","k=6","k=8","k=10","k=12"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))

plot(pr_k_2_l_45, color=1, auc.main = FALSE)
plot(pr_k_4_l_45, add = TRUE, color = 2)
plot(pr_k_6_l_45, add = TRUE, color = 3)
plot(pr_k_8_l_45, add = TRUE, color = 4)
plot(pr_k_10_l_45, add = TRUE, color = 5)
plot(pr_k_12_l_45, add = TRUE, color = 6)
legend("bottomright",c("k=2","k=4","k=6","k=8","k=10","k=12"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))

# STEP 3B Ploting roc and pr curves in terms of l window size.

# the datasets have to be generated a priori

plot(roc1, color=1, auc.main = FALSE)
plot(roc2, add = TRUE, color = 2)
plot(roc3, add = TRUE, color = 3)
plot(roc4, add = TRUE, color = 4)
plot(roc5, add = TRUE, color = 5)
plot(roc6, add = TRUE, color = 6)
legend("bottomright",c("l=15","l=20","l=25","l=35","l=45","l=50"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))
plot(pr1, color=1, auc.main = FALSE)
plot(pr2, add = TRUE, color = 2)
plot(pr3, add = TRUE, color = 3)
plot(pr4, add = TRUE, color = 4)
plot(pr5, add = TRUE, color = 5)
plot(pr6, add = TRUE, color = 6)
legend("bottomright",c("l=15","l=20","l=25","l=35","l=45","l=50"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))

# STEP4B Generate a ultimate data set with all k-mer order and l info plus the acuracy

# the datasets have to be generated a priori

k_list=c(2,4,6,8,10,12)
l_list=c(15,20,25,30,35,40,45,50)
k_vector=c()
l_vector=c()
auc_vector=c()
name_vector=c()
for (k in k_list){
  for (l in l_list){
    name=paste("roc_k_",k,"_l_",l,sep="")
    name_vector <- c(name_vector, name)
    k_vector= c(k_vector,k)
    l_vector= c(l_vector,l)
    y=get(name)
    auc_vector=c(auc_vector,y[2])
  }
}
AUC=as.numeric(auc_vector)
final_plot=data.frame(name_vector,k_vector,l_vector,AUC)

# STEP 5B Ploting the data set.

ggplot(final_plot) + geom_tile(aes(k_vector, l_vector, fill = AUC)) + scale_fill_distiller(palette ="Spectral")
qplot(l_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm") + ggtitle("AUC vs l-window")
qplot(k_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm") + ggtitle("AUC vs K")
qplot(k_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm",formula=y ~ poly(x, 2)) + ggtitle("AUC vs K")
