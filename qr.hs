{-# LANGUAGE BangPatterns #-}


import Control.Monad
import Control.Monad.ST
import Data.List (sortOn)
import Data.Time.Clock
import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Devel
import Prelude hiding ((<>))

type MD = Matrix Double

type VD = Vector Double

sortBySameOrder :: VD -> MD -> (VD, MD)
sortBySameOrder eigenvalues eigenvectors = 
    let indices = [0 .. size eigenvalues - 1]
        sortedIndices = map snd $ sortOn fst $ zip (toList eigenvalues) indices
        sortedEigenvalues = fromList [eigenvalues `atIndex` i | i <- sortedIndices]
        sortedEigenvectors = eigenvectors ¿ sortedIndices
    in (sortedEigenvalues, sortedEigenvectors)

wilkinson_shift :: Matrix Double -> (Double)
wilkinson_shift !mat =
  let n = rows mat
   in if n == 1
        then mat `atIndex` (0, 0)
        else
          let subMat = subMatrix (n - 2, n - 2) (2, 2) mat
              a = subMat `atIndex` (0, 0)
              b = subMat `atIndex` (0, 1)
              c = subMat `atIndex` (1, 0)
              d = subMat `atIndex` (1, 1)
              trace = a + d
              det = a * d - b * c
              mul1 = trace / 2 + sqrt ((trace / 2) ** 2 - det)
              mul2 = trace / 2 - sqrt ((trace / 2) ** 2 - det)
           in if abs (d - mul1) < abs (d - mul2)
                then mul1
                else mul2




sign :: Double -> Double
sign x
  | x > 0 = 1.0
  | x < 0 = -1.0
  | otherwise = 0.0

vecNorm :: VD -> Double
vecNorm v = sqrt $ sumElements $ cmap (^ 2) v


vec_dot :: VD -> Double
vec_dot !x = c where
  vsize = size x
  a = at' x 0
  b = if vsize > 1 then (at' x 1) else 0.0
  c=if vsize > 1
    then  (a * a + b * b)
    else  (a * a)

vec_norm :: VD -> Double
vec_norm x = n where
  vsize = size x
  a = at' x 0
  b = if vsize > 1 then (at' x 1) else 0.0
  n=if vsize > 1
     then  (sqrt (a * a + b * b))
     else  (sqrt (a * a))

reflection_vector :: VD -> (VD, Double)
reflection_vector x = (vec, c)
  where
    x0=at' x  0
    v0=x0+(sign x0)*vec_norm x
    v=vjoin[vector [v0],subVector 1 ((size x) -1) x]
    eps=1e-10
    vec=v
    c= 2.0/(vec_dot(v)+eps)






qr_factorization_householder :: MD -> (MD, MD)
qr_factorization_householder !a = (q, r)
  where
    !n = rows a
    (q, r) = runST $ do
      r' <- thawMatrix a :: ST s (STMatrix s Double)
      q' <- thawMatrix (ident n) :: ST s (STMatrix s Double)
      forM_ [0 .. n - 2] $ \j -> do
        temp <- unsafeFreezeMatrix r'
        let subR = subMatrix (j, j) (n - j, 1) temp
        let (v, c) = reflection_vector $ flatten subR
        let v'=subVector 0 2 v
        let subR = subMatrix (j, j) (2, n - j) temp
        let newR = subR-(scale c (v' `outer` (v'<# subR)))
        setMatrix r' j j newR

        temp2 <- unsafeFreezeMatrix q'
        let subQ = subMatrix (0, j) (n, 2) temp2
        let newQ = subQ-(scale c ((subQ#>v') `outer` v'))
        setMatrix q' 0 j newQ

      changedQ <- unsafeFreezeMatrix q'
      changedR <- unsafeFreezeMatrix r'
      return (changedQ, changedR)

-----------Matrix A         eigenvectors  epsilon         return  MatA       eigenvectors
runLoop :: Matrix Double -> Matrix Double ->Double-> (Matrix Double, Matrix Double)
runLoop a eigenvectors epsilon = go a eigenvectors where
  go !mat !vectors |(((minElement $ cmap abs $ takeDiag $ subMatrix (0, 1) ((((rows mat) - 1)), ((rows mat) - 1)) mat)) <= epsilon) =(mat,vectors)|otherwise= 
    go newMat newVec where 
       u = wilkinson_shift mat
       ui = scale u (ident (rows mat))
       (q, r) = qr_factorization_householder (mat - ui)
       newMat = ((r Numeric.LinearAlgebra.<> q) + ui)
       newVec = vectors Numeric.LinearAlgebra.<> q


myEigenRecursive :: Matrix Double -> Double -> (Vector Double, Matrix Double)
myEigenRecursive a !epsilon =
    if rows a == 1
      then (Numeric.LinearAlgebra.fromList [a `atIndex` (0, 0)], ident 1)
      else
                      (eigenvalues, eigenvectors Numeric.LinearAlgebra.<> (v1 ||| v2))
  where
    n = rows a
    q = ident n

    (finalA, eigenvectors) = runLoop a q epsilon

    diagArrPosition = minIndex $ cmap abs $ takeDiag $ subMatrix (0, 1) (n - 1, n - 1) finalA

    upperMat = subMatrix (0, 0) (diagArrPosition + 1, diagArrPosition + 1) finalA
    lowMat = subMatrix (diagArrPosition + 1, diagArrPosition + 1) (n - diagArrPosition - 1, n - diagArrPosition - 1) finalA

    (eigenvaluesUpper, eigenvectorUpper) = myEigenRecursive upperMat epsilon
    (eigenvaluesLower, eigenvectorLower) = myEigenRecursive lowMat epsilon

    eigenvalues = vjoin [eigenvaluesUpper, eigenvaluesLower]

    v1 = eigenvectorUpper === konst 0.0 (rows eigenvectorLower, cols eigenvectorUpper)
    v2 = konst 0.0 (rows eigenvectorUpper, cols eigenvectorLower) === eigenvectorLower



createV :: VD -> VD
createV x = runST $ do
  v <- thawVector x :: ST s (STVector s Double)
  let first = at' x 0
  let res = first + (sign first) * (vecNorm x)
  writeVector v 0 res
  v' <- unsafeFreezeVector v
  let vNorm = vecNorm v'
  let n = if vNorm == 0 then 1.0 else vNorm
  return $ scale (1 / n) v'

applyMask :: [[Double]] -> Double -> [[Double]]
applyMask h epsilon =
  let mask i j = abs (i - j) <= 1
      applyMaskToElement i j x = if mask i j then x else if abs x < epsilon then 0 else x
   in [[applyMaskToElement i j x | (j, x) <- zip [0 ..] row] | (i, row) <- zip [0 ..] h]

hessen :: MD -> Double -> (MD, MD)
hessen a epsilon = (h, q)
  where
    n = rows a
    (h, q) = runST $ do
      h' <- thawMatrix a :: ST s (STMatrix s Double)
      q' <- thawMatrix (ident n) :: ST s (STMatrix s Double)
      forM_ [0 .. (n - 3)] $ \i -> do
        temp <- unsafeFreezeMatrix h'
        let b = subMatrix (i + 1, i) (n - i - 1, 1) temp
        let v = createV $ flatten b
        c <- unsafeFreezeMatrix h'
        let d = subMatrix (i + 1, i) (n - i - 1, n - i) c
        let vProduct = v <# d
        let out = d - 2.0 * (v `outer` vProduct)
        setMatrix h' (i + 1) i out
        e <- unsafeFreezeMatrix h'
        let m = subMatrix (0, i + 1) (n, n - i - 1) e
        let v2 = m #> v
        let out2 = m - 2.0 * (v2 `outer` v)
        setMatrix h' 0 (i + 1) out2
        l <- unsafeFreezeMatrix q'
        let o = subMatrix (0, i + 1) (n, n - i - 1) l
        let v3 = o #> v
        let out3 = o - 2.0 * (v3 `outer` v)
        setMatrix q' 0 (i + 1) out3
      f <- unsafeFreezeMatrix h'
      let aList = toLists f
      let newH = fromLists $ applyMask aList epsilon
      g <- unsafeFreezeMatrix q'
      return (newH, g)
      
main = do
  w<-loadMatrix "inv_matrix(800 x 800).txt"
  tt1 <- getCurrentTime
  let (!g, !k) = hessen w 1e-6
  let (a, j) = myEigenRecursive g 1e-6
  let vectors=k<>j
  disp 2 (subMatrix (0,0) (1,1) vectors)
  tt2 <- getCurrentTime
  print "my eig recursive:"
  print (diffUTCTime tt2 tt1)
  let (sorted_values, sorted_vectors)= sortBySameOrder a vectors
  tt1 <- getCurrentTime
  let (eigVal,eigVectors)=eig w
  print (subMatrix (0,0) (1,1) eigVectors)  
  tt2 <- getCurrentTime
  print "Hmatrix eig:"
  print (diffUTCTime tt2 tt1)
  let realEigenvectors = cmap realPart eigVectors
  let realEigenvalues = cmap realPart eigVal
  let (sortEigVal,sortEigVectors)=sortBySameOrder realEigenvalues realEigenvectors
  let absSortedVectors = cmap abs sorted_vectors
  let absSortedEigenvectors = cmap abs sortEigVectors
  let absSortedValues = cmap abs sorted_values
  let absSortedEigenvalues = cmap abs sortEigVal

-- Compute the differences
  let vectorDiff = absSortedVectors - absSortedEigenvectors
  let valueDiff = absSortedValues - absSortedEigenvalues

-- Compute the norms of the differences
  let vectorNormDifference = norm_2 vectorDiff
  let valueNormDifference = norm_2 valueDiff
  saveMatrix "vectorDiff_Haskell.txt" "%.8f" vectorDiff
  saveMatrix "valueDiff_Haskell.txt" "%.8f" (asColumn valueDiff)
-- Print the results
  putStrLn $ "\nEigenvector norm difference : \t" ++ show vectorNormDifference
  putStrLn $ "Eigenvalue norm difference  : \t" ++ show valueNormDifference


