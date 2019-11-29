#!/usr/bin/env stack
{- stack
  runghc
  --resolver lts-10.0
  --package async
  --package system-filepath
  --package safe
  --package text
  --package text-show
  --package turtle
  --package foldl
  --package diversity
  --package lens
  --package vector
  --package cassava
  --package containers
  --package bytestring
  --package nlp-scores
  --package exact-combinatorics
-}

{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE ScopedTypeVariables #-}

import Data.Maybe (catMaybes, mapMaybe, fromMaybe)
import Data.Monoid ((<>))
import Safe
import System.Environment (getArgs)
import TextShow (showt)
import Turtle
import Math.Diversity.Diversity (diversity)
import Math.Combinatorics.Exact.Binomial (choose)
import qualified Control.Lens as L
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.Csv as CSV
import qualified Data.Foldable as F
import qualified Data.List as List
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified Data.Set as Set
import qualified Data.Text as T
import qualified Data.Text.IO as T
import qualified Data.Text.Read as T
import qualified Data.Vector as V
import qualified Filesystem.Path.CurrentOS as FP
import qualified NLP.Scores as NLP

newtype Cluster = Cluster { unCluster :: Int } deriving (Eq, Ord, Show)
newtype Size = Size { unSize :: Int } deriving (Eq, Ord)
newtype Diversity = Diversity { unDiversity :: Double } deriving (Eq, Ord)
newtype Label = Label { unLabel :: T.Text } deriving (Eq, Ord, Show)
newtype ClusterMap = ClusterMap { unClusterMap :: Map.Map Cluster (Seq.Seq Label) }
newtype InfoMap = InfoMap { unInfoMap :: Map.Map Cluster Info }
newtype Entropy = Entropy { unEntropy :: Double } deriving (Show)
newtype Purity = Purity { unPurity :: Double } deriving (Show)
newtype NMI = NMI { unNMI :: Double } deriving (Show)

data Info = Info { _diversity :: !Diversity
                 , _size :: !Main.Size
                 , _label :: !Label
                 , _maxDiv :: !Diversity
                 }

-- | Parse the clustering file names.
parseFileName :: T.Text -> Label
parseFileName = Label . fst . T.breakOn "_clustering.csv" . T.drop 2

-- | Read a file to group labels per cluster.
contentsToClusterMap :: B.ByteString-> ClusterMap
contentsToClusterMap = ClusterMap
                     . Map.fromListWith (Seq.><)
                     . fmap (\ !x -> ( Cluster . either error fst . T.decimal $ (Map.!) x "cluster"
                                     , Seq.singleton . Label $ (Map.!) x "label"
                                     )
                                )
                     . V.toList
                     . either error snd
                     . (\x -> CSV.decodeByName x :: Either String (CSV.Header, V.Vector (Map.Map T.Text T.Text)))

-- | Get the diversity of each cluster.
clusterMapToInfoMap :: Label -> ClusterMap -> InfoMap
clusterMapToInfoMap l clusterMap =
  InfoMap
    . fmap (\xs -> Info { _diversity = Diversity . diversity 1 . F.toList $ xs
                        , _size = Size $ Seq.length xs
                        , _label = l
                        , _maxDiv = maxDiv
                        }
           )
    . unClusterMap
    $ clusterMap
  where
    maxDiv = Diversity
           . diversity 0
           . F.toList
           . List.foldl' (Seq.><) Seq.empty
           . Map.elems
           . unClusterMap
           $ clusterMap

-- | Output diversity map.
outputInfoMap :: InfoMap -> Shell Line
outputInfoMap = select
                   . catMaybes
                   . fmap (\((Cluster !c), info) -> textToLine . T.intercalate ","
                            $ [ showt c
                              , showt . unDiversity . _diversity $ info
                              , showt . unDiversity . _maxDiv $ info
                              , showt . unSize . _size $ info
                              , unLabel . _label $ info
                              ]
                           )
                   . Map.toAscList
                   . unInfoMap

-- | Get the entropy of a clustering.
clusteringToEntropy :: ClusterMap -> Entropy
clusteringToEntropy (ClusterMap clusterMap) =
  Entropy
    . negate
    . sum
    . fmap (\c -> sum $ fmap (\l -> (mj c / m) * info l c) labs)
    $ clusters
  where
    clusters = Map.keys clusterMap
    labs = List.nub . F.toList . mconcat . Map.elems $ clusterMap
    info x c
      | p x c == 0 = 0
      | otherwise = (p x c) * (logBase 2 $ p x c)
    m :: Double
    m = fromIntegral . sum . fmap Seq.length . Map.elems $ clusterMap
    mj :: Cluster -> Double
    mj x = fromMaybe 0 $ fromIntegral . Seq.length <$> Map.lookup x clusterMap
    p :: Label -> Cluster -> Double
    p x c = fromMaybe 0 $ (xj x c /) . fromIntegral . Seq.length <$> Map.lookup c clusterMap
    xj :: Label -> Cluster -> Double
    xj x c = fromMaybe 0 $ fromIntegral . Seq.length . mfilter (== x) <$> Map.lookup c clusterMap

-- | Get the purity of a clustering.
clusteringToPurity :: ClusterMap -> Purity
clusteringToPurity (ClusterMap clusterMap) =
  Purity
    . sum
    . fmap (\ c -> (mj c / m) * (maximum $ fmap (\l -> p l c) labs))
    $ clusters
  where
    clusters = Map.keys clusterMap
    labs = List.nub . F.toList . mconcat . Map.elems $ clusterMap
    m :: Double
    m = fromIntegral . sum . fmap Seq.length . Map.elems $ clusterMap
    mj :: Cluster -> Double
    mj x = fromMaybe 0 $ fromIntegral . Seq.length <$> Map.lookup x clusterMap
    p :: Label -> Cluster -> Double
    p x c = fromMaybe 0 $ (xj x c /) . fromIntegral . Seq.length <$> Map.lookup c clusterMap
    xj :: Label -> Cluster -> Double
    xj x c = fromMaybe 0 $ fromIntegral . Seq.length . mfilter (== x) <$> Map.lookup c clusterMap

-- | Get the NMI of a clustering.
clusteringToNMI :: ClusterMap -> NMI
clusteringToNMI (ClusterMap clusterMap) = NMI $ mi / (min hl hc)
  where
    mi = Map.foldl' (+) 0 $ Map.mapWithKey clusterVal $ clusterMap
    convertLnToLog2 = (/ (logBase (exp 1) 2))
    hc = getEntropy clusters
    hl = getEntropy labels
    n = List.genericLength clusters
    all :: ([Cluster], [Label])
    all@(clusters, labels) = unzip
                   . concatMap (\(!c, xs) -> zip (repeat c) . F.toList $ xs)
                   . Map.toList $ clusterMap
    getEntropy :: (Eq a, Ord a) => [a] -> Double
    getEntropy = convertLnToLog2 . logBase (exp 1) . diversity 1
    lsm :: Map.Map Label Double
    lsm = getFreq labels
    csm :: Map.Map Cluster Double
    csm = getFreq clusters
    nLab = fromIntegral $ Map.size lsm
    nClust = fromIntegral $ Map.size csm
    clusterVal :: Cluster -> Seq.Seq Label -> Double
    clusterVal c xs = Map.foldl' (+) 0
                    . Map.mapWithKey (\l v -> (v / n) * logBase 2 ((n * v) / ((look l lsm) * (look c csm))))
                    . getFreq
                    . F.toList
                    $ xs
    getFreq :: Ord k => [k] -> Map.Map k Double
    getFreq = Map.fromListWith (+) . flip zip (repeat 1)
    look :: Ord k => k -> Map.Map k a -> a
    look x = Map.findWithDefault (error "Not found") x

-- | Check if file is valid.
isValidFile :: Bool -> [String] -> FP.FilePath -> Bool
isValidFile invert xs file
  | null xs = True
  | invert && (not . any (\x -> T.isInfixOf (T.pack x) $ format fp file) $ xs) = True
  | not invert && (any (\x -> T.isInfixOf (T.pack x) $ format fp file) xs) = True
  | otherwise = False

main :: IO ()
main = (>>) (putStrLn "label,type,value") $ sh $ do

  file <- find (has "clustering.csv") "."
  size <- liftIO . fmap kilobytes $ du file
  guard (size > 0)

  contents <- liftIO . B.readFile . T.unpack . format fp $ file

  let label = parseFileName . format fp $ file
      clusterMap = contentsToClusterMap contents
      divMap = clusterMapToInfoMap label clusterMap
      entropy = clusteringToEntropy clusterMap
      purity = clusteringToPurity clusterMap
      nmiMin = clusteringToNMI clusterMap
      outputPath = FP.fromText
                 . (<> "diversity.csv")
                 . fst
                 . T.breakOn "clustering"
                 . format fp
                 $ file

  liftIO $ T.putStrLn $ T.intercalate "," [unLabel label, "Entropy", showt $ unEntropy entropy]
  liftIO $  T.putStrLn $ T.intercalate "," [unLabel label, "Purity", showt $ unPurity purity]
  liftIO $  T.putStrLn $ T.intercalate "," [unLabel label, "NMI", showt $ unNMI nmiMin]
