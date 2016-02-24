##### IMPORTING LIBRARIES ######
################################
library(PRROC)
library(ggplot2)
setwd("/Users/davidmasp/GD_UPF-bioinformatica/AGB/Project/ProjectMM/Markov_Model_Project/results/")

##### FIRST LOOP TO LOAD ######
##### DATA AND PERFORM ########
##### THE CALCULATIONS ########
###############################

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

#### PLOTING roc's AUCs in terms of k-order.
plot(roc_k_2_l_45, color=1, auc.main = FALSE)
plot(roc_k_4_l_45, add = TRUE, color = 2)
plot(roc_k_6_l_45, add = TRUE, color = 3)
plot(roc_k_8_l_45, add = TRUE, color = 4)
plot(roc_k_10_l_45, add = TRUE, color = 5)
plot(roc_k_12_l_45, add = TRUE, color = 6)
legend("bottomright",c("k=2","k=4","k=6","k=8","k=10","k=12"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))

#### PLOTING presision's AUCs in terms of k-order.

plot(pr_k_2_l_45, color=1, auc.main = FALSE)
plot(pr_k_4_l_45, add = TRUE, color = 2)
plot(pr_k_6_l_45, add = TRUE, color = 3)
plot(pr_k_8_l_45, add = TRUE, color = 4)
plot(pr_k_10_l_45, add = TRUE, color = 5)
plot(pr_k_12_l_45, add = TRUE, color = 6)
legend("bottomright",c("k=2","k=4","k=6","k=8","k=10","k=12"),lwd=2,col=c("black", "red","green","blue","cyan","pink"))

#### PLOTING roc's AUC vs K in one simple chart.

K=c(2,4,6,8,10,12)
aucs=c(roc_k_2_l_45[2],roc_k_4_l_45[2],roc_k_6_l_45[2],roc_k_8_l_45[2],roc_k_10_l_45[2],roc_k_12_l_45[2])
AUC=as.numeric(aucs)
data_plot=data.frame(K,AUC)
qplot(K, AUC, data=data_plot) + geom_smooth() + ggtitle("AUC vs K")

#### PLOTING roc's AUC vs Lwindow in one simple chart.

L=c(15,20,25,35,45,50)
aucs=c(roc_k_6_l_15[2],roc_k_6_l_20[2],roc_k_6_l_25[2],roc_k_6_l_35[2],roc_k_6_l_45[2],roc_k_6_l_50[2])
AUC=as.numeric(aucs)
data_plot=data.frame(L,AUC)
qplot(L, AUC, data=data_plot) + geom_smooth() + ggtitle("AUC vs l-window")

# GETING A DATA SET TO PLOT A UNIFIED CHART.

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

#### SOME PLOTS

ggplot(final_plot) + geom_tile(aes(k_vector, l_vector, fill = AUC)) + scale_fill_distiller(palette ="Spectral")
ggplot(data=final_plot, aes(l_vector,AUC)) + geom_smooth() + 
qplot(l_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm") + ggtitle("AUC vs l-window")
qplot(k_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm") + ggtitle("AUC vs K")

qplot(k_vector,AUC,data=final_plot) + geom_point() + geom_smooth(method = "lm",formula=y ~ poly(x, 2)) + ggtitle("AUC vs K")

#### GETTING THE FINAL PLOT

ggplot(data.frame(roc_k_6_l_50$curve),aes(x=X1,y=X2,color=X3)) + geom_line(size=1.5) + labs(x="Sensitivity",y="FPR",title=format(roc_k_8_l_40$auc,digits=3), colour="Threshold") + scale_colour_gradient2(low="red", mid="orange",high="yellow")

ggplot(data.frame(pr_k_6_l_50$curve),aes(x=X1,y=X2,color=X3)) + geom_line(size=1.5) + labs(x="Recall",y="Precision",title=format(pr_k_8_l_40$auc.integral,digits=3), colour="Threshold") + scale_colour_gradient2(low="red", mid="orange",high="yellow")
