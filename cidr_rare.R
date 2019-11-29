suppressMessages(library(cidr))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(cowplot))
suppressMessages(library(Rtsne))
source("~/code/single_cell/R/seurat.R", chdir = TRUE)

set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./cidr_img/", showWarnings = FALSE)

args = commandArgs(TRUE)
output = args[1]
inputs = args[2:length(args)]

loadCells = function(x) {
  df = read.table(file = paste0(x, "/barcodes.tsv"), col.names = c("item"),
                  header = FALSE)
  df$label = tail(strsplit(x, "_")[[1]], n = 1)

  return(df)
}

tsne = function(isLabel, df) {
  p = ggplot(df, aes(x = TSNE.1, y = TSNE.2, color = factor(group))) +
    geom_point() +
    xlab("TNSE 1") +
    ylab("TNSE 2") +
    theme_cowplot(font_size = 10) +
    theme(aspect.ratio = 1)

  if(isLabel) {
    p = p + scale_color_manual(guide = guide_legend(title = ""), values=c("#e41a1c", "#377eb8", "#999999"))
    } else {
      p = p + scale_color_discrete(guide = guide_legend(title = ""))
    }
    return(p)
}

cells = do.call("rbind", lapply(inputs, loadCells))

rawNamedMat = do.call("cbind", lapply(inputs, loadNamedMatrix))
rawMat = do.call("cbind", lapply(inputs, loadMatrix))

mat = scDataConstructor(as.matrix(rawNamedMat))
mat = determineDropoutCandidates(mat)
mat = wThreshold(mat)
mat = scDissim(mat)
mat = scPCA(mat)
mat = nPC(mat)
mat = try(scCluster(mat), TRUE)

if(typeof(mat) == "character") {
  sink() # No other messages
  cat("label,cluster,item",file=stdout(),sep="\n")
} else {

# Dimensionality reduction where applicable

tsneMat = Rtsne(as.matrix(t(rawNamedMat)))
tsne1 = tsneMat$Y[,1]
tsne2 = tsneMat$Y[,2]

df = data.frame(item = colnames(rawNamedMat), cluster = mat@clusters, TSNE.1 =
                                                                        tsne1,
                TSNE.2 = tsne2)
df = merge(cells, df, by = "item")

df$group = df$cluster
clusterP = tsne(FALSE, df)
df$group = df$label
labelP = tsne(TRUE, df)

suppressMessages(ggsave(clusterP,
                        file = paste0(c("./cidr_img/", output,
                                        "_clusters.pdf"), collapse = ""),
                        useDingbats = FALSE))
suppressMessages(ggsave(labelP,
                        file = paste0(c("./cidr_img/", output,
                                        "_labels.pdf"), collapse = ""),
                        useDingbats = FALSE))

sink() # No other messages

outDf = df[,c("label", "cluster", "item")]
write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)

}
