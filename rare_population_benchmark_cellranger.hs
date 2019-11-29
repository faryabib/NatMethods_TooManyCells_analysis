#!/usr/bin/env stack
{- stack
  script
  --resolver lts-10.0
  --package async
  --package bytestring
  --package cassava
  --package containers
  --package diversity
  --package foldl
  --package inline-r
  --package system-filepath
  --package text
  --package text-show
  --package turtle
  --package vector
-}

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}

import Data.Maybe (fromJust, isJust)
import Data.Monoid ((<>))
import Math.Diversity.Diversity (diversity)
import TextShow (showt)
import Turtle
import qualified Control.Foldl as Fold
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Csv as CSV
import qualified Data.List as List
import qualified Data.Map as Map
import qualified Data.Text as T
import qualified Data.Vector as V
import qualified Filesystem.Path.CurrentOS as FP

import qualified Foreign.R as R
import Foreign.R (SEXP, SEXPTYPE)
import Language.R.Instance as R
import Language.R.QQ

cluster :: T.Text -> IO ()
cluster outputLabel = sh $ do
    let path = "./cellranger_matrices/common_" <> outputLabel <> ".h5"
        idString = outputLabel <> "_cellranger"

    stdout
        . inproc "cellranger"
            ["reanalyze"
            , "--id", idString
            , "--matrix", path
            ]
        $ mempty

-- | Get a cellranger matrix from an input matrix path.
makeMat :: String -> R s (R.SomeSEXP s)
makeMat path = do
    let matFile = path <> "/matrix.mtx"
        geneFile = path <> "/genes.tsv"
        barcodeFile = path <> "/barcodes.tsv"

    [r| library(cellrangerRkit)

        gbm = load_cellranger_matrix_from_files( mat_fn = matFile_hs
                                               , gene_fn = geneFile_hs
                                               , barcode_fn = barcodeFile_hs
                                               )
    |]

-- | Aggregate cellranger matrices from a list of matrix paths and save.
makeTotalMat :: String -> [FP.FilePath] -> R s ()
makeTotalMat outputLabel
  files = do
    let filesString = fmap (T.unpack . format fp) $ files
        outputFile = "./cellranger_matrices/common_" <> outputLabel <> ".h5"

    [r| library(cellrangerRkit)

        dir.create(file.path("./cellranger_matrices"), showWarnings = FALSE)
        unlink(outputFile_hs, recursive=TRUE)

        matList = lapply(filesString_hs, makeMat_hs)
        totalMat = concatenate_gene_bc_matrices(matList)
        save_cellranger_matrix_h5(totalMat, outputFile_hs, "mm10")
    |]

    return ()

-- | Plot a tsne.
tsne :: R s (R.SomeSEXP s) -> R s (R.SomeSEXP s)
tsne df = [r| suppressMessages(library(ggplot2))
              suppressMessages(library(cowplot))
              p = ggplot(df_hs, aes(x = TSNE.1, y = TSNE.2, color = group)) +
                      geom_point() +
                      xlab("TNSE 1") +
                      ylab("TNSE 2") +
                      scale_color_discrete(guide = guide_legend(title = "")) +
                      theme(aspect.ratio = 1)
          |]

-- | Plot TSNE.
getTsnePlot :: String -> IO ()
getTsnePlot outputLabel = R.runRegion $ do
    let tsnePath = "./"
                <> outputLabel
                <> "_cellranger/outs/analysis/tsne/2_components/projection.csv"
        labelPath = "../labels/labels_b_monocytes_cd4.csv" :: String
        clusterPath = "./" <> outputLabel <> "_cellranger_clustering.csv"
        outputLabelPath = "./"
                       <> outputLabel
                       <> "_cellranger/outs/analysis/tsne/2_components/"
                       <> outputLabel
                       <> "_labels.pdf"
        outputClusterPath = "./"
                         <> outputLabel
                         <> "_cellranger/outs/analysis/tsne/2_components/"
                         <> outputLabel
                         <> "_clusters.pdf"

    [r| suppressMessages(library(ggplot2))
        suppressMessages(library(cowplot))

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

        tsneDf = read.csv(tsnePath_hs)
        labelDf = read.csv(labelPath_hs)
        clusterDf = read.csv(clusterPath_hs)

        tsneDf$item = tsneDf$Barcode
        mainDf = Reduce(function(x, y) merge(x, y, by = "item"), list(tsneDf, clusterDf))

        mainDf$group = mainDf$cluster
        clusterP = tsne(FALSE, mainDf)
        mainDf$group = mainDf$label
        labelP = tsne(TRUE, mainDf)

        suppressMessages(ggsave(clusterP, file = outputClusterPath_hs, useDingbats = FALSE))
        suppressMessages(ggsave(labelP, file = outputLabelPath_hs, useDingbats = FALSE))
    |]

    return ()

-- | Get the benchmark output.
getBenchmark :: T.Text -> IO ()
getBenchmark outputLabel = sh $ do
    let path = fromText
             $ "./"
            <> outputLabel
            <> "_cellranger/outs/analysis/clustering/graphclust/clusters.csv"
        labelPath = "../labels/labels_b_monocytes_cd4.csv"
        outputPath =
            fromText $ "./" <> outputLabel <> "_cellranger_clustering.csv"

    tmp <- using $ mktempfile "/tmp" "cellrangerOutput.csv"

    output tmp
        . sed ("Barcode,Cluster" *> pure "item,cluster")
        . input
        $ path

    output outputPath
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", labelPath, format fp tmp]
      $ mempty

main :: IO ()
main = R.withEmbeddedR R.defaultConfig . sh $ do
    let totalNum       = 1000 :: Int -- Max sample size.
        percent        = 0.005 -- Incremental percent for a rare population.
        maxPercent     = 0.05 -- Maximum percent for a rare population.
        runs = 10 :: Int -- Number of runs
        seeds = [0..runs - 1]
        inc            = round $ fromIntegral totalNum * percent -- Incremental number.
        maxNum         = round $ fromIntegral totalNum * maxPercent -- Max number.
        outputRoot     = "/home/gw/research/10x_t_cells/data/b_monocytes_cd4/random"
        commonPop = ["b"]
        rarePop = ["monocytes", "cd4"]

    sh $ do -- Common population.
        commonPop <- select commonPop
        rareNum <- select [inc,(inc * 2)..maxNum]
        seed <- select seeds

        let commonNum = totalNum - (2 * rareNum)
            outputCommon = outputRoot <> "_seed_" <> showt seed <> "/" <> showt commonNum <> "_" <> commonPop
            outputRare = fmap (\x -> outputRoot <> "_seed_" <> showt seed <> "/" <> showt rareNum <> "_" <> x) rarePop
            files = fmap fromText $ outputCommon : outputRare
            outputLabel = showt commonNum <> "_seed_" <> showt seed

        liftIO $ R.runRegion $ makeTotalMat (T.unpack outputLabel) files
        liftIO $ cluster outputLabel
        liftIO $ getBenchmark outputLabel
        liftIO $ getTsnePlot (T.unpack outputLabel)
