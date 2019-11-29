source("~/code/single_cell/R/seurat.R", chdir = TRUE)

library(wordspace)
library(dendextend)
library(data.tree)
library(jsonlite)

args = commandArgs(TRUE)
output = args[1]
inputs = args[2:length(args)]

print("Loading files:")
print(inputs)

# Load matrix
mat = do.call("cbind", lapply(inputs, loadMatrix))

# Get clustering
hc = hclust(dist.matrix(mat, as.dist = TRUE, byrow = FALSE))

# Get tree
tree = as.dendrogram(hc)

# Get nicely formatted tree from dendrogram.
tree = as.Node(tree)

# Convert to JSON
res = toJSON(as.list(tree, mode = "explicit", unname = TRUE))

# Save matrix.
dir.create(dirname(output), recursive = TRUE)
cat(res, file=output)
