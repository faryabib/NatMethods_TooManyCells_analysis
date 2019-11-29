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

import Data.Function (on)
import Data.Monoid ((<>))
import Control.Monad (mapM_)
import Safe
import System.Environment (getArgs)
import TextShow (showt)
import Turtle
import qualified Control.Foldl as Fold
import qualified Data.List as List
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP

data Clustering = Too
                | TooNorm T.Text
                | TooPrior Int
                | TooBothNorm
                | Seurat
                | Phenograph
                | Monocle
                | RaceID
                | BackSPIN
                | Cidr
                deriving (Eq, Ord, Read, Show)

-- | Join clustering results with labels.
assignLabels :: FP.FilePath -> Shell Line -> Shell Line
assignLabels labelPath = inproc "csvcut" ["-c", "label,cluster,item"]
                       . inproc "csvjoin" ["-c", "item", format fp labelPath, "-"]

clusterToo :: T.Text -> [T.Text] -> IO ()
clusterToo labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-l", labelPath
               , "-o", "too-many-cells_clustering"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . assignLabels (FP.fromText labelPath)
      . input
      $ tmp

clusterTooBothNorm :: T.Text -> [T.Text] -> IO ()
clusterTooBothNorm labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells_both_norm_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-l", labelPath
               , "-o", "too-many-cells_both_norm_clustering"
               , "--normalization", "BothNorm"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooBothNormOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . assignLabels (FP.fromText labelPath)
      . input
      $ tmp

clusterTooPrior :: Clustering -> T.Text -> [T.Text] -> IO ()
clusterTooPrior (TooPrior cut) labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells_cut" <> showt cut <> "_clustering.csv"
        args = ["make-tree"]
            -- <> csvCommands
            <> ["--prior", "too-many-cells_clustering"]
            <> [ "-l", labelPath
               , "-o", "too-many-cells_cut" <> showt cut <> "_clustering"
               ]

    liftIO . print $ csvCommands

    tmp <- using $ mktempfile "/tmp" "tooPriorOutput.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . assignLabels (FP.fromText labelPath)
      . input
      $ tmp
clusterTooPrior _ _ _ = error "clusterTooPrior needs TooPrior type."

clusterTooNorm :: Clustering -> T.Text -> [T.Text] -> IO ()
clusterTooNorm (TooNorm normType) labelPath files = sh $ do
    let csvCommands = concat $ zipWith (\x y -> [x, y]) (repeat "--matrix-path") files
        outputFile = "too-many-cells_" <> normType <> "_clustering.csv"
        args = ["make-tree"]
            <> csvCommands
            <> [ "-l", labelPath
               , "-o", "too-many-cells_" <> normType <> "_clustering"
               , "--normalization", normType
               ]

    liftIO . print $ csvCommands

    tmp <- using . mktempfile "/tmp" $ "too" <> normType <> "Output.csv"

    output tmp
        . sed ("cell,cluster,path" *> pure "item,cluster,path")
        . inproc "too-many-cells" args
        $ mempty

    output (FP.fromText outputFile)
      . assignLabels (FP.fromText labelPath)
      . input
      $ tmp
clusterTooNorm _ _ _ = error "clusterTooNorm needs TooNorm type."

clusterSeurat :: T.Text -> [T.Text] -> IO ()
clusterSeurat labelPath files = sh $ do
    let outputFile = "seurat_clustering.csv"
        args = "/home/gw/code/single_cell/R/seurat_rare.R"
             : ("seurat_clustering")
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterMonocle :: T.Text -> [T.Text] -> IO ()
clusterMonocle labelPath files = sh $ do
    let outputFile = "monocle_clustering.csv"
        args = "/home/gw/code/single_cell/R/monocle_rare.R"
             : "monocle_clustering"
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterRaceID :: T.Text -> [T.Text] -> IO ()
clusterRaceID labelPath files = sh $ do
    let outputFile = "raceid_clustering.csv"
        args = "/home/gw/code/single_cell/R/raceid_rare.R"
             : "raceid_clustering"
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterCidr :: T.Text -> [T.Text] -> IO ()
clusterCidr labelPath files = sh $ do
    let outputFile = "cidr_clustering.csv"
        args = "/home/gw/code/single_cell/R/cidr_rare.R"
             : "cidr_clustering"
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "Rscript" args
        $ mempty

clusterPhenograph :: T.Text -> [T.Text] -> IO ()
clusterPhenograph labelPath files = sh $ do
    let outputFile = "phenograph_clustering.csv"
        args = "/home/gw/code/single_cell/python/phenograph_rare.py"
             : "phenograph_clustering"
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "/home/gw/venv3/bin/python3" args
        $ mempty

clusterBackSPIN :: T.Text -> [T.Text] -> IO ()
clusterBackSPIN labelPath files = sh $ do
    let outputFile = "backspin_clustering.csv"
        args = "/home/gw/code/single_cell/python/backspin_rare.py"
             : "backspin_clustering"
             : files

    liftIO . print $ args

    output (FP.fromText outputFile)
        . assignLabels (FP.fromText labelPath)
        . sed (",NA" *> pure "")
        . inproc "/home/gw/venv3/bin/python3" args
        $ mempty

main :: IO ()
main = do
    args <- liftIO getArgs

    let inFile = "../data/datasets/total_b_monocytes_cd4"
        labelPath = "../labels/labels_b_monocytes_cd4.csv" :: T.Text
        clusterType    = maybe Too read . headMay $ args
        cluster Too         = clusterToo
        cluster TooBothNorm = clusterTooBothNorm
        cluster (TooPrior normType) = clusterTooPrior (TooPrior normType)
        cluster (TooNorm normType) = clusterTooNorm (TooNorm normType)
        cluster Seurat      = clusterSeurat
        cluster Phenograph  = clusterPhenograph
        cluster Monocle     = clusterMonocle
        cluster RaceID      = clusterRaceID
        cluster BackSPIN    = clusterBackSPIN
        cluster Cidr        = clusterCidr
        algorithms = [ Too
                     , TooBothNorm
                     , TooNorm "UQNorm"
                     , TooNorm "TotalMedNorm"
                     , TooNorm "LogCPMNorm"
                     , Seurat
                     , Phenograph
                     , Monocle
                     , Cidr
                     , RaceID
                     , BackSPIN
                     ]

    mapM_ (\ a -> (cluster a) labelPath [inFile]) $ algorithms
