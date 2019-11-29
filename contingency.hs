#!/usr/bin/env stack
{- stack
  script
  --resolver lts-10.0
  --package async
  --package lens
  --package safe
  --package system-filepath
  --package text
  --package text-show
  --package turtle
  --package foldl
-}

{-# LANGUAGE OverloadedStrings #-}

import Data.Maybe (catMaybes, isJust)
import Data.Monoid ((<>))
import TextShow (showt)
import Turtle
import Safe
import qualified Control.Foldl as Fold
import qualified Control.Lens as L
import qualified Data.Foldable as F
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP

newtype Label = Label T.Text
newtype TextN = TextN T.Text
newtype Seed = Seed T.Text

data Clustering = Too
                | TooBothNorm
                | Cellranger
                | Seurat
                | Phenograph
                | Monocle
                | RaceID
                | BackSPIN
                | Cidr
                deriving (Read, Show)

-- | Convert the too-many-cells output to the contingency table input format.
tooToContingency :: Clustering -> Shell Line -> Shell Line
tooToContingency Too = sed (ends (("," <> star (notChar ',')) *> pure ""))
tooToContingency Cellranger = tooToContingency Too
tooToContingency Seurat = tooToContingency Too
tooToContingency Monocle = tooToContingency Too
tooToContingency RaceID = tooToContingency Too
tooToContingency Cidr = tooToContingency Too
tooToContingency BackSPIN = tooToContingency Too
tooToContingency Phenograph = tooToContingency Too

-- | Add labels to the contingency table.
labelContingency :: Label -> TextN -> Seed -> WithHeader Line -> Shell (Maybe Line)
labelContingency (Label label) _ _ (Header h) = return Nothing
labelContingency (Label label) (TextN n) (Seed seed) (Row _ row) =
    select . fmap Just . textToLines $ lineToText row <> "," <> label <> "," <> n <> "," <> seed

-- | Get the contingency table from the too-many-cells output.
contingency :: FP.FilePath -> Maybe (Shell Line)
contingency file = do
    n <- headMay . T.splitOn "_" . format fp . FP.basename $ file
    seed <- flip atMay 2 . T.splitOn "_" . format fp . FP.basename $ file

    let (label, clustering)
            | (\x -> T.isInfixOf "too-many-cells" x && T.isInfixOf "both_norm" x)
                . format fp
                $ file = (Label "too-many-cells_both_norm_maxmed", Too)
            | T.isInfixOf "too-many-cells" . format fp $ file = (Label "too-many-cells", Too)
            | T.isInfixOf "seurat" . format fp $ file = (Label "seurat", Seurat)
            | T.isInfixOf "phenograph" . format fp $ file = (Label "phenograph", Phenograph)
            | T.isInfixOf "monocle" . format fp $ file = (Label "monocle", Monocle)
            | T.isInfixOf "cellranger" . format fp $ file = (Label "cellranger", Cellranger)
            | T.isInfixOf "raceid" . format fp $ file = (Label "raceid", RaceID)
            | T.isInfixOf "backspin" . format fp $ file = (Label "backspin", BackSPIN)
            | T.isInfixOf "cidr" . format fp $ file = (Label "cidr", Cidr)
            | otherwise = error "Unknown clustering algorithm."

    Just
        . join
        . fmap (maybe (select $ textToLines "") return)
        . mfilter isJust
        . join
        . fmap (labelContingency label (TextN n) (Seed seed))
        . header
        . inproc "clustering-contingency" []
        . tooToContingency clustering
        . input
        $ file

main :: IO ()
main = sh $ do
    csvs <- fold (find (ends $ "_clustering.csv") ".") Fold.list

    let outputFinal = "contingency_total.csv"

    output outputFinal
        . F.foldl' ((<|>)) mempty
        . (:) (pure "label1,label2,count,percent,label,n,seed")
        . catMaybes
        . fmap contingency
        $ csvs
