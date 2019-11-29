#!/usr/bin/env stack
{- stack
  script
  --resolver lts-10.0
  --package async
  --package system-filepath
  --package safe
  --package text
  --package text-show
  --package turtle
  --package foldl
-}

{-# LANGUAGE OverloadedStrings #-}

import Data.Monoid ((<>))
import Safe
import System.Environment (getArgs)
import TextShow (showt)
import Turtle
import qualified Control.Foldl as Fold
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP

data Clustering = Too
                | TooBothNorm
                | Seurat
                | Phenograph
                | Monocle
                | RaceID
                | BackSPIN
                | Cidr
                | DimRed
                deriving (Read, Show)

clusterToo :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterToo outputLabel labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = outputLabel <> "_too-many-cells_clustering_pruned.csv"
        args = ["make-tree"]
            <> csvCommands
            -- <> ["--prior", outputLabel <> "_too-many-cells_clustering"]
            <> [ "-l", labelPath
               , "-o", outputLabel <> "_too-many-cells_clustering_prune"
               , "--draw-leaf", "DrawItem DrawLabel"
               , "--min-distance-search", "1"
               , "--smart-cutoff", "5"
               , "--draw-max-leaf-node-size", "300"
               , "--draw-mark", "MarkModularity"
               , "--draw-collection", "PieRing"
               , "--draw-colors", "[\"#999999\", \"#e41a1c\", \"#377eb8\"]"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . inproc "csvcut" ["-c", "label,cluster,item"]
      . inproc "csvjoin" ["-c", "item", labelPath, format fp tmp]
      $ mempty

clusterSeurat :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterSeurat outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_seurat_clustering.csv"
        args = "/home/gw/code/single_cell/R/seurat_rare.R"
             : (outputLabel <> "_seurat_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterMonocle :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterMonocle outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_monocle_clustering.csv"
        args = "/home/gw/code/single_cell/R/monocle_rare.R"
             : (outputLabel <> "_monocle_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterRaceID :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterRaceID outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_raceid_clustering.csv"
        args = "/home/gw/code/single_cell/R/raceid_rare.R"
             : (outputLabel <> "_raceid_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterCidr :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterCidr outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_cidr_clustering.csv"
        args = "/home/gw/code/single_cell/R/cidr_rare.R"
             : (outputLabel <> "_cidr_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterPhenograph :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterPhenograph outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_phenograph_clustering.csv"
        args = "/home/gw/code/single_cell/python/phenograph_rare.py"
             : (outputLabel <> "_phenograph_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "/home/gw/venv3/bin/python3" args
        $ mempty

clusterBackSPIN :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterBackSPIN outputLabel _ files = sh $ do
    let outputFile = outputLabel <> "_backspin_clustering.csv"
        args = "/home/gw/code/single_cell/python/backspin_rare.py"
             : (outputLabel <> "_backspin_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . sed (",NA" *> pure "")
        . inproc "/home/gw/venv3/bin/python3" args
        $ mempty

clusterDimRed :: T.Text -> T.Text -> [T.Text] -> IO ()
clusterDimRed outputLabel labelPath files = sh $ do
    mktree "reduce_dims"
    let outputFile = "reduce_dims/" <> outputLabel <> "_dimensionality_reduction.csv"
        args = "/home/gw/code/single_cell/python/reduce_dims.py"
             : ("reduce_dims/" <> outputLabel)
             : labelPath
             : files

    liftIO . print $ args

    stdout
        . sed (",NA" *> pure "")
        . inproc "/home/gw/venv3/bin/python3" args
        $ mempty

main :: IO ()
main = sh $ do
    args <- liftIO getArgs

    let totalNum       = 1000 -- Max sample size.
        percent        = 0.005 -- Incremental percent for a rare population.
        maxPercent     = 0.05 -- Maximum percent for a rare population.
        runs = 10 :: Int -- Number of runs
        seeds = [0..runs - 1]
        inc            = round $ fromIntegral totalNum * percent -- Incremental number.
        maxNum         = round $ fromIntegral totalNum * maxPercent -- Max number.
        outputRoot     = "/home/gw/research/10x_t_cells/data/b_monocytes_cd4/random"
        commonPop = ["b"]
        rarePop = ["monocytes", "cd4"]
        labelPath = "../labels/labels_b_monocytes_cd4.csv"
        clusterType    = maybe Too read . headMay $ args
        cluster Too         = clusterToo
        cluster TooBothNorm = clusterTooBothNorm
        cluster Seurat      = clusterSeurat
        cluster Phenograph  = clusterPhenograph
        cluster Monocle     = clusterMonocle
        cluster RaceID      = clusterRaceID
        cluster BackSPIN    = clusterBackSPIN
        cluster Cidr        = clusterCidr
        cluster DimRed      = clusterDimRed

    sh $ do -- Common population.
        commonPop <- select commonPop
        rareNum <- select [inc,(inc * 2)..maxNum]
        seed <- select seeds

        let commonNum = totalNum - (2 * rareNum) :: Int
            outputCommon = outputRoot <> "_seed_" <> showt seed <> "/" <> showt commonNum <> "_" <> commonPop
            outputRare = fmap (\x -> outputRoot <> "_seed_" <> showt seed <> "/" <> showt rareNum <> "_" <> x) rarePop
            outputLabel = showt commonNum <> "_seed_" <> showt seed

        liftIO
            . (cluster clusterType) outputLabel labelPath
            $ outputCommon : outputRare
