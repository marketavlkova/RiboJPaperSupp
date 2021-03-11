#!/usr/local/bin/Rscript
### Analysis of amplification efficiency

### get directory info
root.path <- getwd()
work.path <- paste0(root.path, "/qPCR/201009/")
setwd(work.path)
### load data table
data <- read.csv('201009MV_rundata.csv', header = T)
### convert data table into a matrix for easy plotting
col <- 1
row <- 1
mat <- matrix(, ncol = 4, nrow = length(data[, 1]) / 3)
for (r in 1:length(data[, 1])) {
  mat[row, col] <- data[r, 3]
  if (row %% 32 == 0) {
    col <- col + 1
    row <- 1
  } else {
    row <- row + 1
  }
}
### add column with template concentration
mat[, 4] <- c(rep(c(22.45, 2.245, 0.2245, 0.02245,
                    12.6, 1.26, 0.126, 0.0126), 2),
              rep(c(19.75, 1.975, 0.1975, 0.01975,
                    16.55, 1.655, 0.1655, 0.01655), 2))
### reorder matrix for plotting
fin <- rbind(mat[1:4,], mat[9:12,],
            mat[5:8,], mat[13:16,],
            mat[17:20,], mat[25:28,],
            mat[21:24,], mat[29:32,])
### save info for main in plotting
m <- c(expression('pA'^'-39'*'A'^'+69'*'RJ-'),
      expression('pA'^'-39'*'A'^'+69'*'RJ+'),
      expression('pA'^'-39'*'RJ-'),
      expression('pA'^'-39'*'RJ+'))
setwd(root.path)
### plotting
pdf(file = "FigureS2.pdf", width = 8, height = 8)
par(mfcol = c(2,2))

n <- 1
for (j in 0:((length(fin[, 1]) - 4) / 4)) {
  i <- j * 4 + 1
  y <- c(fin[i:(i + 3), 1], fin[i:(i + 3), 2], fin[i:(i + 3), 3])
  x <- c(log10(fin[i:(i + 3), 4]), log10(fin[i:(i + 3), 4]), log10(fin[i:(i + 3), 4]))
  fit <- lm(y~x)
  r2fit <- summary(fit)$r.squared
  slp <- coefficients(fit)[[2]]
  stder <- sqrt(diag(vcov(fit)))
  eff <- (-1 + 10^(-1 / slp)) * 100

  if ((j + 1) %% 2 != 0) {
    plot(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 1], cex = 1.5,
          pch = 16, col = rgb(0, 0, 0, 0.3), main = m[n],
          xlab = expression(paste("Template concentration (log10 ng/", mu, "l)")), ylab = "Ct values",
          xlim = c(-2, 1.5), ylim = c(15, 35))
    points(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 2], pch = 16, col = rgb(0, 0, 0, 0.3), cex = 1.5)
    points(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 3], pch = 16, col = rgb(0, 0, 0, 0.3), cex = 1.5)
    abline(fit, col = 1, lwd=0.5)
    mtext(paste0("Efficiency: ", round(eff, digits = 2), "%"), side = 3, line = -2, at = 0, adj = 0, cex = 0.75)
    mtext(paste0("r2: ", round(r2fit, digits = 3)), side = 3, line = -3, at = 0, adj = 0, cex = 0.75)
    mtext(paste0("Slope: ", round(slp, digits = 2)), side = 3, line = -4, at = 0, adj = 0, cex = 0.75)
    mtext(paste0("Std. Error: ", round(stder[2], digits = 3)), side = 3, line = -5, at = 0, adj = 0, cex = 0.75)
    if (n == 1) {
      legend("topleft", legend = c("F1", "F2"),
              lty = 1, col = c("black", "red"),
              title = "Forward primer")
    }
    n <- n + 1
  } else {
    points(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 1], cex = 1.5,
          pch = 16, col = rgb(1, 0, 0, 0.3))
    points(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 2], pch = 16, col = rgb(1, 0, 0, 0.3), cex = 1.5)
    points(x = log10(fin[i:(i + 3), 4]), y = fin[i:(i + 3), 3], pch = 16, col = rgb(1, 0, 0, 0.3), cex = 1.5)
    abline(fit, col = 2, lwd=0.5)
    mtext(paste0("Efficiency: ", round(eff, digits = 2), "%"),
          side = 3, line = -11, at = -2, adj = 0, cex = 0.75, col = "red")
    mtext(paste0("r2: ", round(r2fit, digits = 3)),
          side = 3, line = -12, at = -2, adj = 0, cex = 0.75, col = "red")
    mtext(paste0("Slope: ", round(slp, digits = 2)),
          side = 3, line = -13, at = -2, adj = 0, cex = 0.75, col = "red")
    mtext(paste0("Std. Error: ", round(stder[2], digits = 3)),
          side = 3, line = -14, at = -2, adj = 0, cex = 0.75, col = "red")
  }
}

dev.off()
