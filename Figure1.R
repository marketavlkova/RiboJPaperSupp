#!/usr/local/bin/Rscript
### Analasis of fluorescence data

library("flowCore") ### handling FCS files
### incorporate gating function
if (!exists("foo", mode="function")) source("FCgating.R")
### incorporate file reordering function
if (!exists("foo", mode="function")) source("ReFiles.R")

### save directory information
root.path <- getwd()
full.path <- c(paste0(root.path, "/FC/PromoterA_RiboJ-/"), paste0(root.path, "/FC/PromoterB_RiboJ-/"),
              paste0(root.path, "/FC/PromoterA_RiboJ+/"), paste0(root.path, "/FC/PromoterB_RiboJ+/"))
### save names of strains for legend
strains <- c("SC312", "SC392", "MG1655")
nstr <- length(strains)
### define carbon sources used
CS <- c("glucose", "galactose", "arabinose", "ribose", "lactose")
nCS <- length(CS)
### define variables for the loop below
vals <- vector()
names <- vector()
cs <- 1
rep <- 1
rj <- 0
### loop through all files in all directories and
### save max kernel density for log10 fluorescence
for (dir in 1:4) {
  ### distinguish between the 2 promoter types
  if (dir %% 2 != 0) {
    prom <- "A"
    if (dir > 2) {
      rj <- rj + 1
    }
  } else {
    prom <- "B"
  }
  setwd(full.path[dir])
  ### save file names
  Files <- Sys.glob("*.fcs")
  ### reorder files so that they are processed by strain
  InFiles <- reorder(Files, nstr)
  str <- 1
  for (fcs in InFiles) {
    ### set vector for technical replicates
    if (rep == 1) {
      set <- vector()
    }
    ### load data file
    data <- read.FCS(fcs, alter.names = T)
    ### extract gfp values from gated cells
    cells <- FCgating(data)
    gfp.log <- log10(exprs(cells$FITC.H))
    gfp.log[is.infinite(gfp.log)] <- 0
    ### get max kernel density value
    denspl <- density(gfp.log)
    peak <- which.max(denspl$y)
    gfp.peak <- denspl$x[peak]
    ### get unique name for the file processed
    name <- paste(strains[str], prom, rj, CS[cs], rep, sep = ".")
    if (rep == 3) {
      ### asign max kernel density value of 3rd tech. rep.
      set <- c(set, gfp.peak)
      ### calculate mean and stdev of the replicates
      vals <- c(vals, mean(set, na.rm = T), sd(set, na.rm = T))
      ### get names for the means and stdevs
      n <- unlist(strsplit(name, split = ".", fixed = T))[1:4]
      names <- c(names, paste(n[1], n[2], n[3], n[4], "mean", sep = "."),
                        paste(n[1], n[2], n[3], n[4], "sd", sep = "."))
      ### update values monitoring progress
      if (cs == nCS) {
        str <- str + 1
        cs <- 1
      } else {
        cs <- cs + 1
      }
      rep <- 1
    } else {
      ### assign max kernel density value of first 2 tech. reps.
      set <- c(set, gfp.peak)
      rep <- rep + 1
    }
  }
}
### reorder mean and stdev values according to the unique names
names(vals) <- names
vals <- vals[order(names)]
### create a matrices for easy plotting
vals <- matrix(vals, nrow = 12, ncol = length(vals) / 12, byrow = T)
a0 <- rbind(vals[5,], vals[9,], vals[1,])
a1 <- rbind(vals[6,], vals[10,], vals[2,])
b0 <- rbind(vals[7,], vals[11,], vals[3,])
b1 <- rbind(vals[8,], vals[12,], vals[4,])
am0 <- cbind(a0[, 5], a0[, 3], a0[, 7])
as0 <- cbind(a0[, 6], a0[, 4], a0[, 8])
bm0 <- cbind(b0[, 5], b0[, 3], b0[, 7])
bs0 <- cbind(b0[, 6], b0[, 4], b0[, 8])
am1 <- cbind(a1[, 5], a1[, 3], a1[, 7])
as1 <- cbind(a1[, 6], a1[, 4], a1[, 8])
bm1 <- cbind(b1[, 5], b1[, 3], b1[, 7])
bs1 <- cbind(b1[, 6], b1[, 4], b1[, 8])
### plotting
setwd(root.path)
pdf(file = "Figure1.pdf", width = 8, height = 6)
par(mfcol = c(1, 2),
    las = 1, xpd = NA)
  ### without RiboJ
  plot(c(1:3), am0[1,], xaxt = "n", pch = 16,
        ylim = c(2, 5), xlim = c(0.75, 3.25), type = "b",
        xlab = "", col = "magenta",
        ylab = "Modal population expression level (log10 a.u.)")
  arrows(c(1:3), am0[1,] - as0[1,],
        c(1:3), am0[1,] + as0[1,],
        length = 0.05, angle = 90, code = 3, col = "magenta")
  points(c(1:3), bm0[1,], pch = 16, type = "b", lty = 3, col = "magenta")
  arrows(c(1:3), bm0[1,] - bs0[1,],
        c(1:3), bm0[1,] + bs0[1,],
        length = 0.05, angle = 90, code = 3, col = "magenta")
  points(c(1:3), am0[2,], pch = 16, type = "b", col = "green")
  arrows(c(1:3), am0[2,] - as0[2,],
        c(1:3), am0[2,] + as0[2,],
        length = 0.05, angle = 90, code = 3, col = "green")
  points(c(1:3), bm0[2,], pch = 16, type = "b", lty = 3, col = "green")
  arrows(c(1:3), bm0[2,] - bs0[2,],
        c(1:3), bm0[2,] + bs0[2,],
        length = 0.05, angle = 90, code = 3, col = "green")
  points(c(1:3), am0[3,], pch = 16, type = "b", col = "saddlebrown")
  arrows(c(1:3), am0[3,] - as0[3,],
        c(1:3), am0[3,] + as0[3,],
        length = 0.05, angle = 90, code = 3, col = "saddlebrown")
  points(c(1:3), bm0[3,], pch = 16, type = "b", lty = 3, col = "saddlebrown")
  arrows(c(1:3), bm0[3,] - bs0[3,],
        c(1:3), bm0[3,] + bs0[3,],
        length = 0.05, angle = 90, code = 3, col = "saddlebrown")
  axis(side = 1, at = 1:3, labels = c("glucose", "galactose", "lactose"))
  text(0.5, 5.25, "A)", cex = 1.5)
  ### with RiboJ
  plot(c(1:3), am1[1,], xaxt = "n", pch = 16,
        ylim = c(2, 5), xlim = c(0.75, 3.25), type = "b",
        xlab = "", col = "magenta",
        ylab = "Modal population expression level (log10 a.u.)")
  arrows(c(1:3), am1[1,] - as1[1,],
        c(1:3), am1[1,] + as1[1,],
        length = 0.05, angle = 90, code = 3, col = "magenta")
  points(c(1:3), bm1[1,], pch = 16, type = "b", lty = 3, col = "magenta")
  arrows(c(1:3), bm1[1,] - bs1[1,],
        c(1:3), bm1[1,] + bs1[1,],
        length = 0.05, angle = 90, code = 3, col = "magenta")
  points(c(1:3), am1[2,], pch = 16, type = "b", col = "green")
  arrows(c(1:3), am1[2,] - as1[2,],
        c(1:3), am1[2,] + as1[2,],
        length = 0.05, angle = 90, code = 3, col = "green")
  points(c(1:3), bm1[2,], pch = 16, type = "b", lty = 3, col = "green")
  arrows(c(1:3), bm1[2,] - bs1[2,],
        c(1:3), bm1[2,] + bs1[2,],
        length = 0.05, angle = 90, code = 3, col = "green")
  points(c(1:3), am1[3,], pch = 16, type = "b", col = "saddlebrown")
  arrows(c(1:3), am1[3,] - as1[3,],
        c(1:3), am1[3,] + as1[3,],
        length = 0.05, angle = 90, code = 3, col = "saddlebrown")
  points(c(1:3), bm1[3,], pch = 16, type = "b", lty = 3, col = "saddlebrown")
  arrows(c(1:3), bm1[3,] - bs1[3,],
        c(1:3), bm1[3,] + bs1[3,],
        length = 0.05, angle = 90, code = 3, col = "saddlebrown")
  axis(side = 1, at = 1:3, labels = c("glucose", "galactose", "lactose"))
  text(0.5, 5.25, "B)", cex = 1.5)
  legend("topleft", legend = strains,
          pch = 16, col = c("magenta", "green", "saddlebrown"),
          title = "Strains")
dev.off()
