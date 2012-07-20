# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Andreas Maunz
# Year: 2012


# Determine class maximally associated to a pattern


# Create level vector from bb table by filling a vector with class values
#
# Params: data frame with pattern and one 'level' col (bb)
#         levels to be recognized

getLevelVecFromBB <- function (x, levels) {
  sp <- c()
  for (l in unique(levels)) { 
    addLevel <- x[levels==l]
    addLevel <- addLevel[complete.cases(addLevel)]
    sp <- c(sp,rep(l,sum(addLevel)))
  }
  sp
}

# Compare local (pattern) prop table to global prop table and determine maximum element's name
# Break ties in favor of the dominant global class
# In case of local and global ties, return one of the dominant global classes with uniform probability
# 
# Params: Local prop table (propTabLocal)
#         Global prop table (propTabGlobal)

getMaxClass <- function (propTabLocal, propTabGlobal) {
  TabDiff <- propTabLocal - propTabGlobal
  MaxElemsLocal <- ( round(TabDiff,10) == round(max(TabDiff),10) )
  MaxElemsLocal <- names(MaxElemsLocal[MaxElemsLocal == TRUE])
  if (length(MaxElemsLocal) == 1)
    return(MaxElemsLocal)
  else
    propTabGlobal <- propTabGlobal[MaxElemsLocal]
    MaxElemsGlobal <- ( round(propTabGlobal,10) == round(max(propTabGlobal),10) )
    MaxElemsGlobal <- names(MaxElemsGlobal[MaxElemsGlobal == TRUE])
    if (length(MaxElemsGlobal) == 1)
      return(MaxElemsGlobal)
    else
      return(MaxElemsGlobal[trunc(runif(1)*length(MaxElemsGlobal)+1)])
}


# Squared error (with Yates correction)
#
# Params: Value to test (x)
#         Expected value (spx)

squaredErr <- function(x, spx) {
  (x-spx-0.5)^2 / spx
}


#' Computes the cumulative sum in terms of logarithmic in- and output 
#' Useful to avoid numerical underflow when summing products of probabilities 
#' When using this function, one can sum sums of log probabilities 
#' See also: http://goo.gl/aJopi 
#' @param logx a vector of log numbers (need not be probabilities) 
#' @return the log of the sum of the exponentiated input 
#' @examples { 
#'   x=c(1,2,3) 
#'   exp(logsum(log(x))) 
#'   # 6 
#' } 

logsum <- function(logx) {
  mypi=max(logx)
  mysum=0
  for (i in 1:length(logx)) {
    mysum = mysum + exp(logx[i]-mypi)
  }
  mypi + log(mysum)
}


#' Computes the weighted mean in terms of log
#' @param logx a vector of log numbers (need not be probabilities) 
#' @param logw a vector of log weights (need not be probabilities) 
#' @return the log of the weighted mean of the exponentiated input 
#' @examples { 
#'   weighted.mean(c(1,2),c(2,1)) # 1.33333333
#'   exp(wlogmean(log(c(1,2)),log(c(2,1)))) # 1.33333333
#' } 

wlogmean <- function (logx,logw) {
    names(logx)=NULL
    names(logw)=NULL
    logsum(logx+logw)-logsum(logw)
}

