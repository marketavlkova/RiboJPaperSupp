#!/usr/local/bin/Rscript
### Analysis of RiboJ's cutting efficiency

### get directory info
root.path <- getwd()
data.path <- c(paste0(root.path, "/qPCR/201207/"),
              paste0(root.path, "/qPCR/201208/"),
              paste0(root.path, "/qPCR/201209/"),
              paste0(root.path, "/qPCR/201213/"),
              paste0(root.path, "/qPCR/210128/"))
### create list of all data tables
### first 3 directories
l <- list()
for (dir in 1:3) {
  setwd(data.path[dir])
  data <- read.csv(dir(pattern='^2')[1])
  col <- 1
  row <- 1
  str <- vector()
  mat <- matrix(, ncol = 3, nrow = length(data[, 1]) / 3)
  for (r in 1:length(data[, 1])) {
    mat[row, col] <- data[r, 3]
    if (col == 3) {
      str <- c(str, paste0(data[r, 1], "-", data[r, 2]))
      col <- 1
      row <- row + 1
    } else {
      col <- col + 1
    }
  }
  rownames(mat) <- str
  colnames(mat) <- c(1:3)
  l[[dir]] <- mat
}
### 4th directory (2 data files)
setwd(data.path[4])
Files <- Sys.glob("*.csv")
for (f in 1:2) {
  data <- read.csv(Files[f], header = T)
  col <- 1
  row <- 1
  str <- vector()
  mat <- matrix(, ncol = 3, nrow = length(data[, 1]) / 3)
  for (r in 1:length(data[, 1])) {
    mat[row, col] <- data[r, 3]
    if (col == 3) {
      str <- c(str, paste0(data[r, 1], "-", data[r, 2]))
      col <- 1
      row <- row + 1
    } else {
      col <- col + 1
    }
  }
  rownames(mat) <- str
  colnames(mat) <- c(1:3)
  l[[(f + 3)]] <- mat
}
### 5th directory
setwd(data.path[5])
data <- read.csv(dir(pattern='^2')[1])
col <- 1
row <- 1
str <- vector()
mat <- matrix(, ncol = 3, nrow = length(data[, 1]) / 3)
for (r in 1:length(data[, 1])) {
  mat[row, col] <- data[r, 3]
  if (col == 3) {
    spl <- data[r, 1]
    name <- unlist(strsplit(spl, split = "-", fixed = T))
    str <- c(str, paste0(name[1], "-", name[2], ".", name[3], "-", data[r, 2]))
    col <- 1
    row <- row + 1
  } else {
    col <- col + 1
  }
}
rownames(mat) <- str
colnames(mat) <- c(1:3)
l[[6]] <- mat
### create a single matrix from all data files
fin.mat <- matrix(, nrow = 63, ncol = 8)
rownames(fin.mat) <- c(rep("K12.1", 9), rep("K12.2", 9),
                        rep("G9", 9), rep("A5", 9),
                        rep("DH5alpha", 9), rep("BW25113", 9),
                        rep("BL21", 9))
colnames(fin.mat) <- c("pD9ua-lacZ", "pD9mv-lacZ",
                      "pD9ua-GFP", "pD9mv-GFP",
                      "pD9m2ua-lacZ", "pD9m2mv-lacZ",
                      "pD9m2ua-GFP", "pD9m2mv-GFP")
for (ls in 1:6) {
  for (r in 1:length(rownames(l[[ls]]))) {
    spl <- rownames(l[[ls]])[r]
    name <- unlist(strsplit(spl, split = "-", fixed = T))
    r.name <- name[2]
    c.name <- paste0(name[1], "-", name[3])
    vals <- l[[ls]][spl,]
    v <- 1
    for (i in 1:8) {
      for (j in 1:63) {
        if ((colnames(fin.mat)[i] == c.name) && (rownames(fin.mat)[j] == r.name) && is.na(fin.mat[j, i])) {
          fin.mat[j, i] <- vals[v]
          v <- v + 1
        }
      }
    }
  }
}
### bootstrapping
n = 1
s <- c("MG1655 D1", "MG1655 D2", "SC392", "SC312",
      "DH5alpha", "BW25113", "BL21 Star")
ss <- c("MG1655 D1", "", "SC392", "", "DH5alpha", "", "BL21 Star")
### matrix for fold-changes
all.boot <- matrix(, nrow = 0, ncol = 3)
### matrix for efficiencies in percentages
all.boot.per <- matrix(, nrow = 0, ncol = 3)
for (str in s) {
  boot <- matrix(, nrow = 10000, ncol = 3)
  boot.per <- matrix(, nrow = 10000, ncol = 3)
  all = fin.mat[n:(n + 8), 1:8]
  for (bt in 1:10000) {
    d9uL = mean(sample(all[, 1], size = 6, replace = T), na.rm = T)
    d9mL = mean(sample(all[, 2], size = 6, replace = T), na.rm = T)
    d9uG = mean(sample(all[, 3], size = 6, replace = T), na.rm = T)
    d9mG = mean(sample(all[, 4], size = 6, replace = T), na.rm = T)
    boot[bt, 1] <- log10(1.95766 ^ (d9uL - d9mL) / 1.95766 ^ (d9uG - d9mG))
    boot.per[bt, 1] <- 100 - 100 * (10 ^ as.numeric(boot[bt, 1]))
    m2uL = mean(sample(all[, 5], size = 6, replace = T), na.rm = T)
    m2mL = mean(sample(all[, 6], size = 6, replace = T), na.rm = T)
    m2uG = mean(sample(all[, 7], size = 6, replace = T), na.rm = T)
    m2mG = mean(sample(all[, 8], size = 6, replace = T), na.rm = T)
    boot[bt, 2] <- log10(1.95766 ^ (m2uL - m2mL) / 1.95766 ^ (m2uG - m2mG))
    boot.per[bt, 2] <- 100 - 100 * (10 ^ as.numeric(boot[bt, 2]))
    boot[bt, 3] <- str
    boot.per[bt, 3] <- str
  }
  all.boot <- rbind(all.boot, boot)
  all.boot.per <- rbind(all.boot.per, boot.per)
  n <- n + 9
}
### adjust for plotting
all.boot <- data.frame(all.boot)
all.boot <- mapply(all.boot[, 1:2], FUN = as.numeric)
all.boot.per <- data.frame(all.boot.per)
all.boot.per <- mapply(all.boot.per[, 1:2], FUN = as.numeric)
### plotting
setwd(root.path)
pdf(file = "Figure3.pdf", width = 9, height = 8)
par(fig = c(0, 1, 0, 1), las = 1)
    boxplot(all.boot.per[1:10000, 1], all.boot.per[1:10000, 2],
            all.boot.per[10001:20000, 1], all.boot.per[10001:20000, 2],
            all.boot.per[20001:30000, 1], all.boot.per[20001:30000, 2],
            all.boot.per[30001:40000, 1], all.boot.per[30001:40000, 2],
            all.boot.per[40001:50000, 1], all.boot.per[40001:50000, 2],
            all.boot.per[50001:60000, 1], all.boot.per[50001:60000, 2],
            all.boot.per[60001:70000, 1], all.boot.per[60001:70000, 2],
          at = c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20), outline = F,
          ylim = c(88, 100),
          xaxt = "n", col = c("orange", "cyan"), xlab = "", ylab = "Autocatalytic efficiency of RiboJ (%)")
    axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5), labels = s)
    legend("bottomright", legend = c("69A", "69C"),
            pch = 22, pt.bg = c("orange", "cyan"),
            title = "Promoter")
par(fig = c(0.175, 0.8, 0.06, 0.575), new = TRUE)
  boxplot(all.boot[1:10000, 1], all.boot[1:10000, 2],
          all.boot[10001:20000, 1], all.boot[10001:20000, 2],
          all.boot[20001:30000, 1], all.boot[20001:30000, 2],
          all.boot[30001:40000, 1], all.boot[30001:40000, 2],
          all.boot[40001:50000, 1], all.boot[40001:50000, 2],
          all.boot[50001:60000, 1], all.boot[50001:60000, 2],
          all.boot[60001:70000, 1], all.boot[60001:70000, 2],
        at = c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20), outline = F,
        xaxt = "n", col = c("orange", "cyan"), xlab = "", ylab = "", cex.axis = 0.75)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5), labels = ss, cex.axis = 0.7)
  title(ylab = "Fold change (log 10)", line = 2, cex.lab = 0.75)
dev.off()
