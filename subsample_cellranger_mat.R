args = commandArgs(TRUE)
seed = args[1]
n = args[2]
input = args[3]
output = args[4]
library(cellrangerRkit)

set.seed(seed)

# Load matrix.
matFile = paste0(input, "/matrix.mtx")
geneFile = paste0(input, "/genes.tsv")
barcodeFile = paste0(input, "/barcodes.tsv")

print("Loading files:")
print(matFile)
print(geneFile)
print(barcodeFile)

gbm = load_cellranger_matrix_from_files( mat_fn = matFile
                                       , gene_fn = geneFile
                                       , barcode_fn = barcodeFile
                                       )

# Subset matrix.
gbmSubset = gbm[,sample(1:ncol(gbm), n, replace = FALSE)]

# Save matrix.
dir.create(output, recursive = TRUE)

# expression matrix
writeMM(exprs(gbmSubset), paste0(output, "/matrix.mtx"))

# data frame of genes
write.table(fData(gbmSubset), file = paste0(output, "/genes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# data frame of cell barcodes
write.table(pData(gbmSubset), file = paste0(output, "/barcodes.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
