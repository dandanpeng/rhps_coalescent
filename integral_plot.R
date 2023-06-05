load("nosel_2000_100/nosel_2000_100_iter1_100_analyzed_trees.RData")

weak0204.truetraj <- mat.true.phentrajs

weak0204.err.array <- err.array
weak0204.err.std.array <- err.std.array
weak0204.qxtest_mat <- qxtest_mat

weak0204.relate.err.array <- err.array.relate
weak0204.relate.err.std.array <- err.std.array.relate
weak0204.relate.qxtest_mat <- qxtest_mat_relate

e1 <- new.env(parent = baseenv())
load("nosel_check_tsinfer/nosel_2000_100_iter1_100_analyzed_trees.RData", envir = e1)
weak0204.tsinfer.err.array <- e1$err.array.tsinfer
weak0204.tsinfer.err.std.array <- e1$err.std.array.tsinfer
weak0204.tsinfer.qxtest_mat <- e1$qxtest_mat_tsinfer

e2 <- new.env(parent = baseenv())
load("nosel_200_100/wo_rescale/nosel_200_100_iter1_94_analyzed_trees.RData", envir = e2)

weak0204.rent.err.array <- e2$err.array.rent
weak0204.rent.err.std.array <- e2$err.std.array.rent
weak0204.rent.qxtest_mat <- e2$qxtest_mat_rent

e3 <- new.env(parent = baseenv())
load("nosel_check_aw/nosel_20_100_iter1_15_analyzed_trees.RData", envir = e3)

e4 <- new.env(parent = baseenv())
load("nosel_check_aw/nosel_20_100_iter16_30_analyzed_trees.RData", envir = e4)


weak0204.aw.err.array <- abind(e3$err.array.aw, e4$err.array.aw)
weak0204.aw.std.array <- abind(e3$err.std.array.aw, e4$err.std.array.aw)

weak0204.avg.err.mat <- apply(weak0204.err.array, 1:2, mean, na.rm = TRUE) #average err of all iterations
weak0204.avg.sq.err.mat <- apply(weak0204.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.rent.avg.err.mat <- apply(weak0204.rent.err.array, 1:2, mean, na.rm = TRUE)
weak0204.rent.avg.sq.err.mat <- apply(weak0204.rent.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.relate.avg.err.mat <- apply(weak0204.relate.err.array, 1:2, mean, na.rm = TRUE)
weak0204.relate.avg.sq.err.mat <- apply(weak0204.relate.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.tsinfer.avg.err.mat <- apply(weak0204.tsinfer.err.array, 1:2, mean, na.rm = TRUE)
weak0204.tsinfer.avg.sq.err.mat <- apply(weak0204.tsinfer.err.array^2, 1:2, mean, na.rm = TRUE)

weak0204.aw.avg.err.mat <- apply(weak0204.aw.err.array, 1:2, mean, na.rm = TRUE)
weak0204.aw.avg.sq.err.mat <- apply(weak0204.aw.err.array^2, 1:2, mean, na.rm = TRUE)

pal <- c("black", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00")

xl <- c(-.25, 0)
yl.b <- c(-1,4)
yl.mse <- c(0,6)
wid <- 5
#hei <- 2.5
hei <- 4


############ Bias ####################
# plot(rev(time) - max(time), weak0204.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
#      bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
#      ylab = "Bias", col = pal[1])
#axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
plot(rev(time) - max(time), weak0204.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
     bty = "n", xlab = "Time", 
     ylab = "Bias", col = pal[1])
lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,1] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,1] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,1] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,1] , col = pal[5])
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Proportion-of-lineage")

plot(rev(time) - max(time), weak0204.avg.err.mat[,2], type = "l", xlim = xl, ylim = yl.b,
     bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
     ylab = "Bias", col = pal[1])
axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,2] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,2] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,2] , col = pal[5])
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Lineages-remaining")


plot(rev(time) - max(time), weak0204.avg.err.mat[,3], type = "l", xlim = xl, ylim = yl.b,
     bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
     ylab = "Bias", col = pal[1])
axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,3] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,3] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,3] , col = pal[5])
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Waiting-time")

############     MSE   ################
pdf("proportion-of-lineage_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
     bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
     ylab = "Mean squared error", col = pal[1])
axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))

lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,1] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,1] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,1] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,1] , col = pal[5])

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

#mtext("F", side = 3, line = .5, at = min(xl))
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Proportion-of-lineage")
#mtext("B", side = 3, line = .5, at = min(xl))
dev.off()

#weak0204.aw.avg.sq.err.mat <- apply(weak0204.aw.err.array^2 * 2.5, 1:2, mean, na.rm = TRUE)
pdf("lineages_remaining_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,2], type = "l", xlim = xl, ylim = yl.mse,
     bty = "n", xlab = "Time", xaxt = "n",
     ylab = "Mean squared error", col = pal[1])
axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))

lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,2] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,2] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,2] , col = pal[5])

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
#mtext("F", side = 3, line = .5, at = min(xl))
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Lineages-remaining")
#mtext("B", side = 3, line = .5, at = min(xl))

dev.off()


pdf("waiting-time_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,3], type = "l", xlim = xl, ylim = yl.mse,
     bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
     ylab = "Mean squared error", col = pal[1])
axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,3] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,3] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,3] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,3] , col = pal[5])

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
#mtext("F", side = 3, line = .5, at = min(xl))
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,3,5,2,1)], legend = c("tsinfer", "RELATE", "ARGWeaver", "RENT+", "True"))
title("Waiting-time")
#mtext("B", side = 3, line = .5, at = min(xl))

dev.off()


pdf("proportion-of-lineage_mse_zoom_in.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,1], type = "l", xlim = c(-.02, 0), ylim = c(0,0.8),
     bty = "n", xlab = "Time (approximate kya)", xaxt = "n",
     ylab = "Mean squared error", col = pal[1])
axis(1, at = seq(-.02, 0, by = 0.02), labels = seq(10, 0, by = -10))

lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,1] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,1] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,1] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,1] , col = pal[5])


#mtext("F", side = 3, line = .5, at = min(xl))
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,5,3,2,1)], legend = c("tsinfer", "ARGWeaver", "RELATE", "RENT+", "True"))
title("Proportion-of-lineage")
#mtext("B", side = 3, line = .5, at = min(xl))
dev.off()

#weak0204.aw.avg.sq.err.mat <- apply(weak0204.aw.err.array^2 * 2.5, 1:2, mean, na.rm = TRUE)
pdf("lineages_remaining_mse_zoomin.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,2], type = "l", xlim = c(-.02, 0), ylim = c(0,4),
     bty = "n", xlab = "Time", xaxt = "n",
     ylab = "Mean squared error", col = pal[1])
axis(1, at = seq(-.02, 0, by = 0.02), labels = seq(10, 0, by = -10))

lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,2] , col = pal[2])
lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,2] , col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,2] , col = pal[4])
lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,2] , col = pal[5])

#mtext("F", side = 3, line = .5, at = min(xl))
legend("topright", bty = "n", lty = 1, cex = 0.7, col = pal[c(4,5,3,2,1)], legend = c("tsinfer", "ARGWeaver", "RELATE", "RENT+", "True"))
title("Lineages-remaining")
#mtext("B", side = 3, line = .5, at = min(xl))

dev.off()
