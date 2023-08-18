# Dandan Peng Feb-1-2023

load("nosel_20_100_iter1_20_analyzed_trees.RData")

weak0204.truetraj <- mat.true.phentrajs

weak0204.err.array <- err.array
weak0204.err.std.array <- err.std.array
weak0204.qxtest_mat <- qxtest_mat

weak0204.rent.err.array <- err.array.rent
weak0204.rent.err.std.array <- err.std.array.rent
weak0204.rent.qxtest_mat <- qxtest_mat_rent

weak0204.relate.err.array <- err.array.relate
weak0204.relate.err.std.array <- err.std.array.relate
weak0204.relate.qxtest_mat <- qxtest_mat_relate

weak0204.tsinfer.err.array <- err.array.tsinfer
weak0204.tsinfer.err.std.array <- err.std.array.tsinfer
weak0204.tsinfer.qxtest_mat <- qxtest_mat_tsinfer

weak0204.aw.err.array <- err.array.aw
weak0204.aw.std.array <- err.std.array.aw

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

#--- multi-panel plot of bias and MSE

xl <- c(-.25, 0)
yl.b <- c(-1,4)
yl.mse <- c(0,6)
wid <- 5
#hei <- 2.5
hei <- 4

####### proportion-of-lineage ####################

pdf("proportion-of-lineage_bias.pdf", width = wid, height = hei)
#par(mfrow = c(1,2), mar = c(4,4,2,1), las = 1, mgp = c(2,.7,0), cex = .85)

plot(rev(time) - max(time), weak0204.avg.err.mat[,1], type = "l", xlim = xl, ylim = yl.b,
 bty = "n", xlab = "Time",
ylab = "Bias", col = pal[1]) #err.neut(proportioon-of-lineages)
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,1] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,1] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,1] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,1] , col = pal[5])}
#lines(rev(time) - max(time), weak0204.avg.err.mat[,2] , col = pal[2]) #err.smoothmom(lineages remaining)
#lines(rev(time) - max(time), weak0204.avg.err.mat[,3] , col = pal[3]) #err.wt_l1(waiting time)
#lines(rev(time) - max(time), weak0204.avg.err.mat[,4] , col = pal[4]) #err.sharedN.wt(sharedN waiting)
#lines(rev(time) - max(time), weak0204.avg.err.mat[,5] , col = pal[5]) #err.straight

#lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,1] , col = pal[1], lty = 2)
#lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,2] , col = pal[2], lty = 2)
#lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,3] , col = pal[3], lty = 2)
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))
title("Proportion-of-lineage")
dev.off()

########## Lineages-remaining ##############

pdf("lineage-remaining_bias.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.err.mat[,2], type = "l", xlim = xl, ylim = yl.b,
     bty = "n", xlab = "Time",
     ylab = "Bias", col = pal[1]) #err.neut(proportioon-of-lineages)
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,2] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,2] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,2] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,2] , col = pal[5])}

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))
title("Lineages-remaining")
dev.off()

########## Waiting-time ##############

pdf("waiting-time_bias.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.err.mat[,3], type = "l", xlim = xl, ylim = yl.b,
     bty = "n", xlab = "Time",
     ylab = "Bias", col = pal[1]) #err.neut(proportioon-of-lineages)
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.err.mat[,3] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.err.mat[,3] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.err.mat[,3] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.err.mat[,3] , col = pal[5])}

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))
title("Waiting-time")
dev.off()



#########################################    MSE   ############################################
pdf("proportion-of-lineage_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,1], type = "l", xlim = xl, ylim = yl.mse,
 bty = "n", xlab = "Time",
ylab = "Mean squared error", col = pal[1])
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,1] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,1] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,1] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,1] , col = pal[5])}

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
#mtext("F", side = 3, line = .5, at = min(xl))
legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))
title("Proportion-of-lineage")
#mtext("B", side = 3, line = .5, at = min(xl))

dev.off()

############ lineages-remaining  mse #########
pdf("lineages-remaining_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,2], type = "l", xlim = xl, ylim = yl.mse,
     bty = "n", xlab = "Time",
     ylab = "Mean squared error", col = pal[1])
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,2] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,2] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,2] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,2] , col = pal[5])}
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
#mtext("F", side = 3, line = .5, at = min(xl))
legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))

#mtext("B", side = 3, line = .5, at = min(xl))
title("Lineages-remaining")
dev.off()

########### waiting-time mse ########

pdf("waiting-time_mse.pdf", width = wid, height = hei)
plot(rev(time) - max(time), weak0204.avg.sq.err.mat[,3], type = "l", xlim = xl, ylim = yl.mse,
     bty = "n", xlab = "Time",
     ylab = "Mean squared error", col = pal[1])
if(!is.na(new_argv$rent)){lines(rev(time) - max(time), weak0204.rent.avg.sq.err.mat[,3] , col = pal[2])}
if(!is.na(new_argv$relate)){lines(rev(time) - max(time), weak0204.relate.avg.sq.err.mat[,3] , col = pal[3])}
if(!is.na(new_argv$tsinfer)){lines(rev(time) - max(time), weak0204.tsinfer.avg.sq.err.mat[,3] , col = pal[4])}
if(!is.na(new_argv$argweaver)){lines(rev(time) - max(time), weak0204.aw.avg.sq.err.mat[,3] , col = pal[5])}

if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}
#mtext("F", side = 3, line = .5, at = min(xl))
legend("topleft", bty = "n", lty = 1, col = pal[1:5], legend = c("true", "Rent+", "RELATE", "tsinfer", "ARGWeaver"))
title("Waiting-time")
#mtext("B", side = 3, line = .5, at = min(xl))

dev.off()


#Confidence intervals
prop.covered <- function(vec, q = 1.96, ...){
  mean(abs(vec) <= q, ...)
}

weak0204.covered.mat <- apply(weak0204.err.std.array, 1:2, prop.covered, na.rm = TRUE)
weak0204.rent.covered.mat <- apply(weak0204.rent.err.std.array, 1:2, prop.covered, na.rm = TRUE)
weak0204.relate.covered.mat <- apply(weak0204.relate.err.std.array, 1:2, prop.covered, na.rm = TRUE)
weak0204.tsinfer.covered.mat <- apply(weak0204.tsinfer.err.std.array, 1:2, prop.covered, na.rm = TRUE)


plot(rev(time) - max(time), weak0204.covered.mat[,1], type = "l", xlim = c(-.25,0), ylim = c(0,1))
lines(rev(time) - max(time), weak0204.rent.covered.mat[,1], col = pal[2])
lines(rev(time) - max(time), weak0204.relate.covered.mat[,1], col = pal[3])
lines(rev(time) - max(time), weak0204.tsinfer.covered.mat[,1], col = pal[4])
if(sel.intenses > 0){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}


plot(rev(time) - max(time), weak0204.covered.mat[,2], type = "l", xlim = c(-.25,0), ylim = c(.6,1))
lines(rev(time) - max(time), weak0204.relate.covered.mat[,2], col = "orange")

plot(rev(time) - max(time), weak0204.covered.mat[,3], type = "l", xlim = c(-.25,0), ylim = c(.6,1))
lines(rev(time) - max(time), weak0204.relate.covered.mat[,3], col = "orange")


