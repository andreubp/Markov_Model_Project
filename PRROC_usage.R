library(PRROC);
setwd("/Users/davidmasp/GD_UPF-bioinformatica/AGB/Project/ProjectMM/Markov_Model_Project")
fg_data = read.table("v2_fg.th");
bg_data = read.table("v2_bg.th");
fg=fg_data$V1
bg=bg_data$V1
roc=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
pr = pr.curve(scores.class0 = fg,scores.class1 = bg, curve = TRUE)
png('roc_curve.png')
roc
plot(roc)
dev.off()
png('precission_plot.png')
pr
plot(pr)
dev.off()
