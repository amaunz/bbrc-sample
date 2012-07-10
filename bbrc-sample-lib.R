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

