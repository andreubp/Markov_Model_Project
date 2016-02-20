library(PRROC);
fg = rnorm(300);
bg = rnorm(500,-2);
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
