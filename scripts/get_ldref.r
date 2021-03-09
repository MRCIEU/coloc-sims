# Get LD matrix for simulations

# Download LD reference panel

library(simulateGP)

bfile <- "/Users/gh13047/data/ld_files/ldmat/AFR"
data(ldetect)
ldetect <- subset(ldetect, chr=="chr22" & pop=="AFR")

map <- generate_ldobj("/Users/gh13047/data/ld_files/ldmat/AFR_1kg_chr22", bfile, ldetect)

