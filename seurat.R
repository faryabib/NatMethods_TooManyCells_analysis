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

# Execute the Seurat pipeline using guide from online vignettes (more strict
# filterings) on a matrix with cell barcodes and labels.
seuratPipeline = function(mat, cells) {

  sMat = CreateSeuratObject( raw.data = mat
                          , min.cells = 3
                          , min.genes = 200
                          , project = "Random"
                            )

  labelDf = data.frame(label = cells$label)
  rownames(labelDf) = cells$item

  sMat = AddMetaData(object = sMat, metadata = labelDf)

  mito.genes = grep(pattern = "^MT-", x = rownames(x = sMat@data), value = TRUE)
  percent.mito = Matrix::colSums(sMat@raw.data[mito.genes,
                                              ])/Matrix::colSums(sMat@raw.data)
  sMat = AddMetaData(object = sMat, metadata = percent.mito, col.name =
                                                              "percent.mito")

  pdf("./seurat_img/qc.pdf", useDingbats=FALSE)
  VlnPlot(object = sMat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  dev.off()

  sMat = FilterCells(object = sMat, subset.names = c("nGene", "percent.mito"),
                        low.thresholds = c(200, -Inf), high.thresholds = c(Inf,
                                                                            0.05))

  sMat = NormalizeData( object = sMat
                    , normalization.method = "LogNormalize"
                    , scale.factor = 10000
                      )

  sMat = FindVariableGenes(object = sMat, mean.function = ExpMean,
                                                    dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                                                    x.high.cutoff = 3, y.cutoff =
                                                                        0.5)

  sMat = ScaleData(object = sMat, vars.to.regress = c("nUMI", "percent.mito"))

  sMat = RunPCA(object = sMat, pc.genes = sMat@var.genes, do.print = TRUE,
                pcs.print = 1:5,
                genes.print = 5)

  sMat = ProjectPCA(object = sMat, do.print = FALSE)

  sMat = JackStraw(object = sMat, num.replicate = 100, display.progress = FALSE)

  sMat = FindClusters(object = sMat, reduction.type = "pca", dims.use = 1:10,
                      resolution = 0.6, print.output = 0, save.SNN = TRUE)  # From vignette

  sMat = RunTSNE(object = sMat, dims.use = 1:10)

  return(sMat)

}
