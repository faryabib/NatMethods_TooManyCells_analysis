set.seed(0) # Set seed for entire stochastic process.

sink(file="/dev/null") # No other messages
dir.create("./seurat_img/", showWarnings = FALSE)

suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(cellrangerRkit))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(cowplot))
suppressMessages(source("~/code/single_cell/R/seurat_tabula_musa.R", chdir = TRUE))

args = commandArgs(TRUE)
output = args[1]
inputs = args[2:length(args)]

loadCells = function(x) {
    df = read.table(file = paste0(x, "/barcodes.tsv"), col.names = c("item"),
                    header = FALSE)
    df$label = tail(strsplit(x, "_")[[1]], n = 1)

    return(df)
}

cells = do.call("rbind", lapply(inputs, loadCells))

mat = suppressMessages(do.call("cbind", lapply(inputs, loadMatrix)))

sMat = suppressMessages(seuratPipeline(mat, cells))

p = TSNEPlot(object = sMat, do.return = TRUE)
p = p + theme_cowplot(font_size = 10) + theme(aspect.ratio=1)
ggsave(paste0("./seurat_img/", output, ".pdf"), plot = p, useDingbats = FALSE)

print(labels)

pL = TSNEPlot(object = sMat, group.by = "label", colors.use = c("#e41a1c",
                                                                "#377eb8",
                                                                "#999999"
                                                                ),
              do.return = TRUE)
pL = pL + theme_cowplot(font_size = 12) + theme(aspect.ratio=1)
ggsave(paste0("./seurat_img/label_", output, ".pdf"), plot = pL, useDingbats =
                                                                   FALSE)

df = GetClusters(sMat)
names(df) = c("item", "cluster")
df = merge(cells, df, by = "item")
outDf = df[,c("label", "cluster", "item")]

sink()

write.csv(outDf, stdout(), row.names = FALSE, quote = FALSE)
