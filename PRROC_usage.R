library(PRROC);
setwd("/Users/davidmasp/GD_UPF-bioinformatica/AGB/Project/ProjectMM/Markov_Model_Project")
fg_data = read.table("GBBFOX_k_10_l_10_fg.th");
bg_data = read.table("GBBFOX_k_10_l_10_bg.th");
start=sample(0:length(fg_data$V1)-10000,1)
end=start+10000
fg=fg_data$V1[start:end]
start=sample(0:length(bg_data$V1)-10000,1)
end=start+10000
bg=bg_data$V1[start:end]

roc=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
pr = pr.curve(scores.class0 = fg,scores.class1 = bg, curve = TRUE)
#png('roc_curve.png')
roc
plot(roc)
#dev.off()
#png('precission_plot.png')
pr
plot(pr)
#dev.off()
