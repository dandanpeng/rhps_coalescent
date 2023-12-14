sample_num = 2000
scene = "recent"

#e4 <- new.env(parent = baseenv())
#load(paste(scene, "_", as.character(sample_num), "_100/scale/nosel_", as.character(sample_num), "_rent_wt_iter1_100.RData", sep = ""), envir = e4)


pal <- c("black", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#2a4ac5")
time <- c(seq(0, 4, by = 0.001))

prop.covered <- function(vec, q = 1.96, ...){
  mean(abs(vec) <= q, ...)
}

if(scene == "nosel"){
  yl.b <- c(-1, 1)
  y.pos <- 1
  if(sample_num == 20){
    yl.mse <- c(0, 1)
    yl.cv <- c(0.6, 1)
  }else if(sample_num == 200){
    yl.mse <- c(0, 2)
  }else if(sample_num == 2000){
    yl.mse <- c(0, 1.5)
    yl.cv <- c(0.4, 1)
  }
}else if(scene == "recent"){
  yl.mse = c(0, 5)
  yl.b <- c(-1, 4)
  y.pos <- 4
  yl.cv <- c(0, 1)
}

################ integral plot ###################
#integral 
e1 <- new.env(parent = baseenv())
load(paste(scene, "_2000_100/", scene, "_2000_100_iter101_1000_analyzed_trees.RData", sep = ""), envir = e1)
e2 <- new.env(parent = baseenv())
load(paste(scene, "_200_100/wo_scale/", scene, "_200_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e2)
e3 <- new.env(parent = baseenv())
load(paste(scene, "_check_tsinfer/", scene, "_2000_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e3)
e4 <- new.env(parent = baseenv())
load(paste(scene, "_check_argweaver/", scene, "_20_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e4)
e5 <- new.env(parent = baseenv())
load(paste(scene, "_check_argneedle/", scene, "_2000_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e5)
e6 <- new.env(parent = baseenv())
load(paste(scene, "_2000_100/scale/", scene, "_2000_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e6)

err.array <- abind(e6$err.array, e1$err.array)
err.std.array <- abind(e6$err.std.array, e1$err.std.array)

pdf(paste("integral_", scene, ".pdf", sep = ""), width = 15, height = 15)
par(mfrow = c(3, 3), mar = c(5, 7, 3, 1), oma = c(1, 1, 1, 1))
labels <- LETTERS[1:9]
for(i in c(1, 3, 2)){
  #plot(rev(time) - max(time), apply(e1$err.array[,i,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = c(-1, 1), 
  #     las = 1, bty = "n", xlab = "Time", ylab = "Bias", cex.lab = 1.5)
  plot(rev(time) - max(time), apply(err.array[,i,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
       las = 1, bty = "n", xlab = "", ylab = "", lwd = 1.25, cex.axis = 1.8)
  if(i == 1){
    title(ylab = "Bias", line = 4, cex.lab = 2.5)
    mtext(text = labels[1], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
    legend(x = -0.2, y = 1, bty = "n", lty = 1, lwd = 1.3, cex = 1.5, x.intersp = 1, seg.len = 1.5, ncol = 2,
           col = pal[c(1, 2, 3, 4, 5, 6)], legend = c("True", "RENT+", "RELATE", "Tsinfer", "ARGweaver", "ARG-Needle"))
  }else if(i == 2){
    mtext(text = labels[3], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }else if(i == 3){
    mtext(text = labels[2], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }
  #axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
  lines(rev(time) - max(time), apply(e2$err.array.rent[,i,], 1, mean, na.rm = TRUE), col = pal[2], lwd = 1.25)
  lines(rev(time) - max(time), apply(e6$err.array.relate[,i,], 1, mean, na.rm = TRUE), col = pal[3], lwd = 1.25)
  lines(rev(time) - max(time), apply(e3$err.array.tsinfer[,i,], 1, mean, na.rm = TRUE), col = pal[4], lwd = 1.25)
  lines(rev(time) - max(time), apply(e4$err.array.aw[,i,], 1, mean, na.rm = TRUE), col = pal[5], lwd = 1.25)
  lines(rev(time) - max(time), apply(e5$err.array.an[,i,], 1, mean, na.rm = TRUE), col = pal[6], lwd = 1.25)

  if(scene == 'recent'){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(err.array[,i,]^2, 1, mean, na.rm = TRUE),  type = "l", xlim = c(-.25, 0), ylim = yl.mse, 
       las = 1, bty = "n", xlab = "", ylab = "", lwd = 1.3, cex.axis = 1.8)
  if(i == 1){
    title(ylab = "MSE", line = 4, cex.lab = 2.5)
    mtext(text = labels[4], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }else if(i == 2){
    mtext(text = labels[6], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }else if(i == 3){
    mtext(text = labels[5], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }
  #axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
  lines(rev(time) - max(time), apply(e2$err.array.rent[,i,]^2, 1, mean, na.rm = TRUE), col = pal[2])
  lines(rev(time) - max(time), apply(e6$err.array.relate[,i,]^2, 1, mean, na.rm = TRUE), col = pal[3])
  lines(rev(time) - max(time), apply(e3$err.array.tsinfer[,i,]^2, 1, mean, na.rm = TRUE), col = pal[4])
  lines(rev(time) - max(time), apply(e4$err.array.aw[,i,]^2, 1, mean, na.rm = TRUE), col = pal[5])
  lines(rev(time) - max(time), apply(e5$err.array.an[,i,]^2, 1, mean, na.rm = TRUE), col = pal[6])
  #if(i == 1){
  #  legend(x = -0.08, y = 5, bty = "n", lty = 1, cex = 1, x.intersp = 1, seg.len = 1,
  #         col = pal[c(1, 2, 3, 4, 5, 6)], legend = c("True", "RENT+", "RELATE", "Tsinfer", "ARGweaver", "ARG-Needle"))
  #}
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

for(estimator in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(err.std.array, 1:2, prop.covered, na.rm = TRUE)[,estimator], bty = "n", type = "l", xlim = c(-.25, 0), ylim = yl.cv,
       las = 1, xlab = "",  ylab = "", lwd = 1.3, cex.axis = 1.8)
  #axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
  if(estimator == 1){
    title(ylab = "Coverage", line = 4, cex.lab = 2.5)
    mtext(text = labels[7], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }else if(estimator == 2){
    mtext(text = labels[9], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }else if(estimator == 3){
    mtext(text = labels[8], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
  }

  lines(rev(time) - max(time), apply(e2$err.std.array.rent, 1:2, prop.covered, na.rm = TRUE)[,estimator], col = pal[2])
  lines(rev(time) - max(time), apply(e6$err.std.array.relate, 1:2, prop.covered, na.rm = TRUE)[,estimator], col = pal[3])
  lines(rev(time) - max(time), apply(e3$err.std.array.tsinfer, 1:2, prop.covered, na.rm = TRUE)[,estimator], col = pal[4])
  lines(rev(time) - max(time), apply(e4$err.std.array.aw, 1:2, prop.covered, na.rm = TRUE)[,estimator], col = pal[5])
  lines(rev(time) - max(time), apply(e5$err.std.array.an, 1:2, prop.covered, na.rm = TRUE)[,estimator], col = pal[6])
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

mtext("Time (coalescent units, past to left)", side = 1, outer = TRUE, cex = 2)

dev.off()

##########################  Bias #####################################
e1 <- new.env(parent = baseenv())
load(paste(scene, "_", as.character(sample_num), "_100/wo_scale/", scene, "_", as.character(sample_num), "_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e1)

e2 <- new.env(parent = baseenv())
load(paste(scene, "_check_tsinfer/", scene, "_", as.character(sample_num), "_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e2)

if(sample_num == 20){
  e3 <- new.env(parent = baseenv())
  load(paste(scene, "_check_argweaver/", scene, "_20_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e3)
  
  if(scene == "nosel"){
    e5 <- new.env(parent = baseenv())
    load(paste(scene, "_", as.character(sample_num), "_100/time_coarse_lineages_remaining/0-4by0.005.RData", sep = ""), envir = e5)
    
    e6 <- new.env(parent = baseenv())
    load(paste(scene, "_check_argweaver/interval0.05nosel_check_argweaver.RData", sep = ""), envir = e6)
  }
}

if(sample_num == 2000){
  e7 <- new.env(parent = baseenv())
  load(paste(scene, "_check_argneedle/", scene, "_", as.character(sample_num), "_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e7)
}

pdf(paste(scene, "_", as.character(sample_num), "_100/", scene, "_", as.character(sample_num), "_integral.pdf", sep = ""), width = 15, height = 15)
par(mfrow = c(3, 3), mar = c(5, 7, 3, 1), oma = c(1, 1, 1, 1))
plot(rev(time) - max(time), apply(e1$err.array[,1,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
     las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
title(ylab = "Bias", line = 4, cex.lab = 2)
if(sample_num != 2000){
  lines(rev(time) - max(time), apply(e1$err.array.rent[,1,], 1, mean, na.rm = TRUE), col = pal[2])
}
lines(rev(time) - max(time), apply(e1$err.array.relate[,1,], 1, mean, na.rm = TRUE), col = pal[3])
lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,1,], 1, mean, na.rm = TRUE), col = pal[4])
if(sample_num == 20){
  lines(rev(time) - max(time), apply(e3$err.array.aw[,1,], 1, mean, na.rm = TRUE), col = pal[5])
  legend(x = -0.08, y = y.pos, bty = "n", lty = 1, cex = 1.5, x.intersp = 1, seg.len = 1.5, ncol = 2,
         col = pal[c(1, 2, 3, 4, 5)], legend = c("True", "RENT+", "RELATE", "Tsinfer", "ARGweaver"))
}else if(sample_num == 200){
  legend(x = -0.08, y = y.pos, bty = "n", lty = 1, cex = 1.3, x.intersp = 1, seg.len = 1,
         col = pal[c(1, 2, 3, 4)], legend = c("True", "RENT+", "RELATE", "Tsinfer"))
}else if(sample_num == 2000){
  lines(rev(time) - max(time), apply(e7$err.array.an[,1,], 1, mean, na.rm = TRUE), col = pal[6])
  legend(x = -0.2, y = y.pos, bty = "n", lty = 1, cex = 1.3, x.intersp = 1, seg.len = 1, ncol = 2,
           col = pal[c(1, 6, 3, 4)], legend = c("True", "ARG-Needle", "RELATE", "Tsinfer"))
}
if(scene == "recent"){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

plot(rev(time) - max(time), apply(e1$err.array[,3,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
     las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
lines(rev(time) - max(time), apply(e1$err.array.rent[,3,], 1, mean, na.rm = TRUE), col = pal[2])
lines(rev(time) - max(time), apply(e1$err.array.relate[,3,], 1, mean, na.rm = TRUE), col = pal[3])
lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,3,], 1, mean, na.rm = TRUE), col = pal[4])
if(sample_num == 20){
  lines(rev(time) - max(time), apply(e3$err.array.aw[,3,], 1, mean, na.rm = TRUE), col = pal[5])
}
if(sample_num == 2000){
  lines(rev(time) - max(time), apply(e7$err.array.an[,3,], 1, mean, na.rm = TRUE), col = pal[6])
}
if(scene == "recent"){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}


if(sample_num == 20 & scene == "nosel"){
  plot(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$ms_err, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
       las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$rent_err, 1, mean, na.rm = TRUE), col = pal[2])
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$relate_err, 1, mean, na.rm = TRUE), col = pal[3])
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$tsinfer_err, 1, mean, na.rm = TRUE), col = pal[4])
  lines(rev(e6$times.coarse) - max(e6$times.coarse), apply(e6$aw_err, 1, mean, na.rm = TRUE), col = pal[5])
}else{
  plot(rev(time) - max(time), apply(e1$err.array[,2,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
       las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
  if(sample_num != 2000){
    lines(rev(time) - max(time), apply(e1$err.array.rent[,2,], 1, mean, na.rm = TRUE), col = pal[2])
  }
  lines(rev(time) - max(time), apply(e1$err.array.relate[,2,], 1, mean, na.rm = TRUE), col = pal[3])
  lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,2,], 1, mean, na.rm = TRUE), col = pal[4])
  if(sample_num == 20){
    lines(rev(time) - max(time), apply(e3$err.array.aw[,2,], 1, mean, na.rm = TRUE), col = pal[5])
  }
  if(sample_num == 2000){
    lines(rev(time) - max(time), apply(e7$err.array.an[,2,], 1, mean, na.rm = TRUE), col = pal[6])
  }
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}


##########################  MSE #####################################
plot(rev(time) - max(time), apply(e1$err.array[,1,]^2, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.mse,
     las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
title(ylab = "MSE", line = 4, cex.lab = 2)
if(sample_num != 2000){
  lines(rev(time) - max(time), apply(e1$err.array.rent[,1,]^2, 1, mean, na.rm = TRUE), col = pal[2])
}
lines(rev(time) - max(time), apply(e1$err.array.relate[,1,]^2, 1, mean, na.rm = TRUE), col = pal[3])
lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,1,]^2, 1, mean, na.rm = TRUE), col = pal[4])
if(sample_num == 20){
  lines(rev(time) - max(time), apply(e3$err.array.aw[,1,]^2, 1, mean, na.rm = TRUE), col = pal[5])
}
if(sample_num == 2000){
  lines(rev(time) - max(time), apply(e7$err.array.an[,1,]^2, 1, mean, na.rm = TRUE), col = pal[6])
}
if(scene == "recent"){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

#legend(x = -0.1, y = 1, bty = "n", lty = 1, cex = 1.3, x.intersp = 1, seg.len = 1,
#       col = pal[c(1, 2, 3, 4, 5)], legend = c("True", "RENT+", "RELATE", "Tsinfer"))


plot(rev(time) - max(time), apply(e1$err.array[,3,]^2, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.mse, 
     las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
if(sample_num != 2000){
  lines(rev(time) - max(time), apply(e1$err.array.rent[,3,]^2, 1, mean, na.rm = TRUE), col = pal[2])
}
lines(rev(time) - max(time), apply(e1$err.array.relate[,3,]^2, 1, mean, na.rm = TRUE), col = pal[3])
lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,3,]^2, 1, mean, na.rm = TRUE), col = pal[4])
if(sample_num == 20){
  lines(rev(time) - max(time), apply(e3$err.array.aw[,3,]^2, 1, mean, na.rm = TRUE), col = pal[5])
}
if(sample_num == 2000){
  lines(rev(time) - max(time), apply(e7$err.array.an[,3,]^2, 1, mean, na.rm = TRUE), col = pal[6])
}
if(scene == "recent"){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}


if(sample_num == 20 & scene == "nosel"){
  plot(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$ms_err^2, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.mse, 
       las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$rent_err^2, 1, mean, na.rm = TRUE), col = pal[2])
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$relate_err^2, 1, mean, na.rm = TRUE), col = pal[3])
  lines(rev(e5$times.coarse) - max(e5$times.coarse), apply(e5$tsinfer_err^2, 1, mean, na.rm = TRUE), col = pal[4])
  lines(rev(e6$times.coarse) - max(e6$times.coarse), apply(e6$aw_err^2, 1, mean, na.rm = TRUE), col = pal[5])
}else{
  plot(rev(time) - max(time), apply(e1$err.array[,2,]^2, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.mse, 
       las = 1, bty = "n", xlab = "", ylab = "", cex.axis = 1.8)
  if(sample_num != 2000){
    lines(rev(time) - max(time), apply(e1$err.array.rent[,2,]^2, 1, mean, na.rm = TRUE), col = pal[2])
  }
  lines(rev(time) - max(time), apply(e1$err.array.relate[,2,]^2, 1, mean, na.rm = TRUE), col = pal[3])
  lines(rev(time) - max(time), apply(e2$err.array.tsinfer[,2,]^2, 1, mean, na.rm = TRUE), col = pal[4])
  if(sample_num == 20){
    lines(rev(time) - max(time), apply(e3$err.array.aw[,2,]^2, 1, mean, na.rm = TRUE), col = pal[5])
  }
  if(sample_num == 2000){
    lines(rev(time) - max(time), apply(e7$err.array.an[,2,]^2, 1, mean, na.rm = TRUE), col = pal[6])
  }
#if(sample_num == 20){
#  lines(rev(time) - max(time), apply(e3$err.array.aw[,3,]^2, 1, mean, na.rm = TRUE), col = pal[5])
}
if(scene == "recent"){
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}


##########################  Confidence Interval #####################################

for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.std.array, 1:2, prop.covered, na.rm = TRUE)[,i], bty = "n", type = "l", xlim = c(-.25, 0), ylim = yl.cv,
       xlab = "", ylab = "", cex.axis = 1.8)
  #axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
  if(i == 1){title(ylab = "Coverage", line = 4, cex.lab = 2)}
  if(sample_num != 2000){
    #if(sample_num == 20 & estimator == 3 & scene == "nosel"){
    #  lines(rev(time) - max(time), apply(e4$err.wt_l1.std.rent, 1, prop.covered, na.rm = TRUE), col = pal[2])
    #}else{
    lines(rev(time) - max(time), apply(e1$err.std.array.rent, 1:2, prop.covered, na.rm = TRUE)[,i], col = pal[2])
    #}
  }
  lines(rev(time) - max(time), apply(e1$err.std.array.relate, 1:2, prop.covered, na.rm = TRUE)[,i], col = pal[3])
  lines(rev(time) - max(time), apply(e2$err.std.array.tsinfer, 1:2, prop.covered, na.rm = TRUE)[,i], col = pal[4])
  if(sample_num == 20){
    lines(rev(time) - max(time), apply(e3$err.std.array.aw, 1:2, prop.covered, na.rm = TRUE)[,i], col = pal[5])
  }
  if(sample_num == 2000){
    lines(rev(time) - max(time), apply(e7$err.std.array.an, 1:2, prop.covered, na.rm = TRUE)[,i], col = pal[6])
  }
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

mtext("Time (coalescent units, past to left)", side = 1, outer = TRUE, cex = 2)

dev.off()


############### RENT+ scaled vs. original ########################
e1 <- new.env(parent = baseenv())
load(paste(scene, "_", as.character(sample_num), "_100/scale/", scene, "_", as.character(sample_num), "_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e1)
e2 <- new.env(parent = baseenv())
load(paste(scene, "_", as.character(sample_num), "_100/wo_scale/", scene, "_", as.character(sample_num), "_100_iter1_100_analyzed_trees.RData", sep = ""), envir = e2)

pdf(paste("RENT_", scene, "_", as.character(sample_num), "_compare.pdf", sep = ""), width = 15, height = 15)
par(mfrow = c(3, 3), mar = c(4, 5, 3, 1), oma = c(1, 1, 1, 1))
bias_y_lab <- c("Bias", NA, NA)
for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.array.rent[,i,], 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.b, 
       las = 1, bty = "n", xlab = NA, ylab = bias_y_lab[i], cex.lab = 1.8, cex.axis = 1.8, col = pal[2])
  lines(rev(time) - max(time), apply(e2$err.array.rent[,i,], 1, mean, na.rm = TRUE), col = "#6796FF")
  if(i == 1){
    legend(x = -0.08, y = y.pos, bty = "n", lty = 1, cex = 2, x.intersp = 1, seg.len = 1,
           col = c(pal[2], "#6796FF"), legend = c("scaled", "original"))
  }
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

mse_y_lab <- c("MSE", NA, NA)
for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.array.rent[,i,]^2, 1, mean, na.rm = TRUE), type = "l", xlim = c(-.25, 0), ylim = yl.mse,
       las = 1, bty = "n", xlab = NA, ylab = mse_y_lab[i], cex.lab = 1.8,cex.axis = 1.8, col = pal[2])
  lines(rev(time) - max(time), apply(e2$err.array.rent[,i,]^2, 1, mean, na.rm = TRUE), col = "#6796FF")
  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}



ci_y_lab <- c("Coverage", NA, NA)
for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.std.array.rent, 1:2, prop.covered, na.rm = TRUE)[,i], bty = "n", type = "l", xlim = c(-.25, 0), ylim = yl.cv,
       las = 1, xlab = NA,  ylab = ci_y_lab[i], cex.lab = 1.8, cex.axis = 1.8, col = pal[2])
  #axis(1, at = seq(-0.25, 0, by = 0.05), labels = seq(120, 0, by = -24))
  lines(rev(time) - max(time), apply(e2$err.std.array.rent, 1:2, prop.covered, na.rm = TRUE)[,i], col = "#6796FF")

  if(scene == "recent"){
    lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
    lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
  }
}

mtext("Time (coalescent units, past to left)", side = 1, outer = TRUE, cex = 2)

dev.off()

#################### Compare metrics based on different sample sizes #######################
e1 <- new.env(parent = baseenv())
e2 <- new.env(parent = baseenv())
e3 <- new.env(parent = baseenv())

labels <- LETTERS[1:6]
pdf("SampleSize_Compare.pdf", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 1))

scene <- "nosel"
load(paste(scene, "_20_100/wo_scale/", scene, "_20_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e1)
load(paste(scene, "_200_100/wo_scale/", scene, "_200_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e2)
load(paste(scene, "_2000_100/scale/", scene, "_2000_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e3)

for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.array[,i,]^2, 1, mean, na.rm = TRUE),  type = "l", xlim = c(-.25, 0), 
       ylim = c(0, 1), las = 1, bty = "n", xlab = "", ylab = "", lwd = 1.3, cex.axis = 1.8, col = "#E69F00")
  if(i == 1){
    mtext(text = labels[1], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)
    legend(-0.1, 1, bty = "n", lty = 1, cex = 1.5, lwd = 1.3,legend = c("20 samples", "200 samples", "2000 samples"), 
           col = c("#E69F00", "#56B4E9", "#CC79A7"))
  }
  if(i == 3){mtext(text = labels[2], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)}
  lines(rev(time) - max(time), apply(e2$err.array[,i,]^2, 1, mean, na.rm = TRUE), col = "#56B4E9")
  
  if(i == 2){mtext(text = labels[3], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)}
  lines(rev(time) - max(time), apply(e3$err.array[,i,]^2, 1, mean, na.rm = TRUE), col = "#CC79A7")
}

scene <- "recent"
load(paste(scene, "_20_100/wo_scale/", scene, "_20_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e1)
load(paste(scene, "_200_100/wo_scale/", scene, "_200_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e2)
load(paste(scene, "_2000_100/scale/", scene, "_2000_100_iter1_100_analyzed_trees.RData", sep = ''), envir = e3)

for(i in c(1, 3, 2)){
  plot(rev(time) - max(time), apply(e1$err.array[,i,]^2, 1, mean, na.rm = TRUE),  type = "l", xlim = c(-.25, 0), 
       ylim = c(0, 4), las = 1, bty = "n", xlab = "", ylab = "", lwd = 1.3, cex.axis = 1.8, col = "#E69F00")
  if(i == 1){mtext(text = labels[4], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)}
  
  lines(rev(time) - max(time), apply(e2$err.array[,i,]^2, 1, mean, na.rm = TRUE), col = "#56B4E9")
  if(i == 3){mtext(text = labels[5], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)}
  lines(rev(time) - max(time), apply(e3$err.array[,i,]^2, 1, mean, na.rm = TRUE), col = "#CC79A7")
  if(i == 2){mtext(text = labels[6], side = 3, line = 0.8, at = par("usr")[1], adj = 0, cex = 1.3)}
  
  lines(-c(.02, .02), c(-100, 100), col = "gray", lty = 2)
  lines(-c(.04, .04), c(-100, 100), col = "gray", lty = 2)
}

mtext("MSE", side = 2, outer = TRUE, cex = 2, las = 0)
mtext("Time (coalescent units, past to left)", side = 1, outer = TRUE, cex = 2)

dev.off()
