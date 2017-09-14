#read in bp quality distributions
qd1 <- read.table("R1_qd.txt")
qd2 <- read.table("R2_qd.txt")
qd3 <- read.table("R3_qd.txt")
qd4 <- read.table("R4_qd.txt")

#read in read quality frequency distributions
rf1 <- read.table("R1_rf.txt")
rf2 <- read.table("R2_rf.txt")
rf3 <- read.table("R3_rf.txt")
rf4 <- read.table("R4_rf.txt")

#plot bp quality distributions
plot(qd1$V1, qd1$V2, main = "Avg Qscore per BP position R1", xlab = "BP position", ylab = "Quality score", pch = 19, col = "steel blue", type = "b")
plot(qd2$V1, qd2$V2, main = "Avg Qscore per BP position R2", xlab = "BP position", ylab = "Quality score", pch = 19, col = "steel blue", type = "b")
plot(qd3$V1, qd3$V2, main = "Avg Qscore per BP position R3", xlab = "BP position", ylab = "Quality score", pch = 19, col = "steel blue", type = "b")
plot(qd4$V1, qd4$V2, main = "Avg Qscore per BP position R4", xlab = "BP position", ylab = "Quality score", pch = 19, col = "steel blue", type = "b")

#plot read quality frequency distribution
plot(rf1$V1, log10(rf1$V2), main = "Average Quality Score Read Frequency R1", xlab = "Avg Qscore / Read", ylab = "log10(Read Frequency)", pch = 19, col = "forest green")
plot(rf2$V1, log10(rf2$V2), main = "Average Quality Score Read Frequency R2", xlab = "Avg Qscore / Read", ylab = "log10(Read Frequency)", pch = 19, col = "forest green")
plot(rf3$V1, log10(rf3$V2), main = "Average Quality Score Read Frequency R3", xlab = "Avg Qscore / Read", ylab = "log10(Read Frequency)", pch = 19, col = "forest green")
plot(rf4$V1, log10(rf4$V2), main = "Average Quality Score Read Frequency R4", xlab = "Avg Qscore / Read", ylab = "log10(Read Frequency)", pch = 19, col = "forest green")

#read in index counts
ih_out <- read.table("ih_output_table.txt")

#plots index counts
colors <- ifelse(substring(ih_out$V1, 1,8) == substring(ih_out$V1, 10), "forest green", "steel blue")
plot(log(ih_out$V2), col = colors, pch = 19, main = "Index pairs counts", ylab = "log10(counts)", type = "h")
legend("topright", inset = 0.1, lty = 1, lwd = 5, legend = c("swapped", "not swapped"), col = c("steel blue", "forest green"))






