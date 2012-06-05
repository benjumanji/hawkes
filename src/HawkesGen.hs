{-# LANGUAGE BangPatterns #-}

module Main
  where

import qualified Data.Vector.Unboxed as U

import Control.Arrow ((&&&))
import Control.Monad.Loops (unfoldrM)
import Data.List (foldl')
import System.Random.Mersenne 
  ( MTGen
  , newMTGen
  , random
  )

main :: IO () 
main = do gen <- newMTGen Nothing
          evs <- generateHawkes gen hc 1000
          print evs
  where hc = HC 1.2 0.6 0.8 

-- | This is the context that a user can pass around in order to
--   generate hawkes processes. It is for simple univariate processes
data HawkesContext = HC { hcLambda :: {-# UNPACK #-} !Double 
                        , hcAlpha  :: {-# UNPACK #-} !Double 
                        , hcBeta   :: {-# UNPACK #-} !Double 
                        } deriving (Show, Eq)

-- | This is an internal structure for the unfold
data HawkesSeed = HS { hsGen     :: MTGen
                     , hsHc      :: HawkesContext
                     , hsHorizon :: Double 
                     , hsEvents  :: [Double]
                     , hsLStar   :: Double
                     } 
                       
untilM :: Monad m => (a ->  m (Maybe a)) -> a -> m a
untilM f x = f x >>= maybe (return x) (\x -> untilM f x)
             
generateHawkes :: MTGen         -- ^ Generator for even timings 
               -> HawkesContext -- ^ Context for the process
               -> Double        -- ^ The event time horizon
               -> IO [Double]   -- ^ GeneratedProcess
generateHawkes g hc t = 
  do let lstar = hcLambda hc
     ev <- (\x -> -(1 / lstar * log x)) `fmap` random g
     let hs = HS g hc t [ev] lstar
     if ev <= t then hsEvents `fmap` untilM generateHawkes' hs
                else return []


generateHawkes' :: HawkesSeed -> IO (Maybe HawkesSeed)
generateHawkes' hs@(HS g hc t evs lstar) =
  do ev' <- eventloop g hc t ev evs lstar
     return $! update hs `fmap` ev'
  where ev = head evs
        lstar = intensity hc ev (tail evs) + (hcAlpha hc)
        update hs x = hs { hsEvents = x:evs }


eventloop :: MTGen         -- ^ PRNG for event times
          -> HawkesContext -- ^ Context
          -> Double        -- ^ Event horizon
          -> Double        -- ^ Last generated event
          -> [Double]      -- ^ Event times
          -> Double        -- ^ Max intensity
          -> IO (Maybe Double)  -- ^ Newly generated event
eventloop g hc t ev evs lstar = do
  ev' <- (\x -> ev - (1 / lstar) * log x) `fmap` random g 
  if ev' >= t 
    then return Nothing
    else do d <- random g
            let i = intensity hc ev' evs
            if d <= i / lstar
              then return $ Just ev'
              else eventloop g hc t ev' evs i

intensity :: HawkesContext -- ^ Process Context
          -> Double        -- ^ Event time        
          -> [Double]      -- ^ Events
          -> Double        -- ^ Intensity
intensity (HC l a b) ev evs = foldl' binop l evs
  where binop z x = z + a * (exp $ -b * (ev - x))

