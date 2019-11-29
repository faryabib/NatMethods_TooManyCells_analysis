suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(cellrangerRkit))
suppressMessages(library(plyr))

# Load a 10x matrix from a directory.
loadMatrix = function(x) {
  mat = load_cellranger_matrix_from_files(mat_fn = paste0(x, "/matrix.mtx")
                                        , gene_fn = paste0(x, "/genes.tsv")
                                        , barcode_fn = paste0(x, "/barcodes.tsv")
                                          )

  return(exprs(mat))
}

# Load a 10x matrix from a directory with row and column names.
loadNamedMatrix = function(x) {
  gbm = load_cellranger_matrix_from_files(mat_fn = paste0(x, "/matrix.mtx")
                                        , gene_fn = paste0(x, "/genes.tsv")
                                        , barcode_fn = paste0(x, "/barcodes.tsv")
                                          )

  mat = exprs(gbm)
  colnames(mat) = pData(gbm)$barcode
  rownames(mat) = fData(gbm)$id

  return(mat)
}

# Execute the Seurat pipeline using guide from Tabula Musa online vignettes
# (less strict filtering) on a matrix with cell barcodes and labels.
seuratPipeline = function(mat, cells) {

  sMat = CreateSeuratObject( raw.data = mat
                          , min.cells = 3
                          , min.genes = 200
                          , project = "Random"
                            )

  labelDf = data.frame(label = cells$label)
  rownames(labelDf) = cells$item

  sMat = AddMetaData(object = sMat, metadata = labelDf)

  sMat = FilterCells(object = sMat, subset.names = c("nGene", "nUMI"),
                        low.thresholds = c(500, 1000))

  sMat = NormalizeData(object = sMat, scale.factor = 1e6)

  sMat = ScaleData(object = sMat)

  sMat = FindVariableGenes(object = sMat, x.high.cutoff = Inf, y.cutoff = 0.5, x.low.cutoff = 0.1)

  sMat = RunPCA(object = sMat, do.print = TRUE)

  sMat = ProjectPCA(object = sMat, do.print = FALSE)

  sMat = FindClusters(object = sMat, reduction.type = "pca", dims.use = 1:10,
                      resolution = 1, print.output = 0, save.SNN = TRUE) # From vignette

  sMat = RunTSNE(object = sMat, dims.use = 1:10, perplexity=30)

  return(sMat)

}
