#!/usr/bin/env stack
{- stack
  script
  --resolver lts-10.0
  --package async
  --package system-filepath
  --package text
  --package text-show
  --package turtle
-}

{-# LANGUAGE OverloadedStrings #-}

import TextShow (showt)
import Turtle
import qualified Data.Text as T
import qualified Filesystem.Path.CurrentOS as FP

-- | Subsample a cellranger matrix.
writePop :: Int -> Int -> T.Text -> T.Text -> IO ()
writePop seed n input output = sh $ do
    procs
        "Rscript"
        [ "/home/gw/code/single_cell/R/subsample_cellranger_mat.R"
        , showt seed
        , showt n
        , input
        , output
        ]
        mempty

main :: IO ()
main = sh $ do
    let totalNum = 1000 -- Max sample size.
        percent = 0.005 -- Incremental percent for a rare population.
        maxPercent = 0.05 -- Maximum percent for a rare population.
        runs = 10 -- Number of runs
        seeds = [0..runs - 1]
        inc = round $ fromIntegral totalNum * percent -- Incremental number.
        maxNum = round $ fromIntegral totalNum * maxPercent -- Max number.
        inputRoot = "/home/gw/research/10x_t_cells/data/datasets/"
        outputRoot = "/home/gw/research/10x_t_cells/data/b_monocytes_cd4/random"
        commonPop = ["b"]
        rarePop = ["monocytes", "cd4"]

    sh $ do -- Common population.
        dataSet <- select commonPop
        rareNum <- select . fmap (* 2) $ [inc,(inc * 2)..maxNum]
        seed <- select seeds

        let commonNum = totalNum - rareNum
            input = inputRoot <> dataSet
            output = outputRoot <> "_seed_" <> showt seed <> "/" <> showt commonNum <> "_" <> dataSet

        liftIO $ writePop seed commonNum input output

    sh $ do -- Rare population.

        dataSet <- select rarePop
        rareNum <- select [inc,(inc * 2)..maxNum]
        seed <- select seeds

        let input = inputRoot <> dataSet
            output = outputRoot <> "_seed_" <> showt seed <> "/" <> showt rareNum <> "_" <> dataSet

        liftIO $ writePop seed rareNum input output
