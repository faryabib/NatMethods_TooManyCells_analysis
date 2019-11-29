library(splatter)
library(Matrix)

args = commandArgs(TRUE)
seed = strtoi(args[1])
n = strtoi(args[2])
output = args[3]

common = n
rare = (1000 - n) / 2

commonFreq = n / 1000
rareFreq = rare / 1000

# Load matrix.
params = newSplatParams(nGenes = 1000, groupCells = c(common, rare, rare), seed = seed)
sim = splatSimulate(params, method = "groups")
mat = sim@assayData$counts
barcodes = colnames(sim@assayData$counts)
genes = rownames(sim@assayData$counts)

# Save labels.
dir.create(output, recursive = TRUE)

df = data.frame(item = sim$Cell, label = sim$Group)
write.table(df, file = paste0(output, "/labels.csv"), sep = ",",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save matrix.

# expression matrix
writeMM(Matrix(mat, sparse = TRUE), paste0(output, "/matrix.mtx"))

# data frame of genes
outRowDf = data.frame(x = genes, y = genes)
write.table(outRowDf, file = paste0(output, "/genes.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# data frame of cell barcodes
write.table(barcodes, file = paste0(output, "/barcodes.tsv"), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
