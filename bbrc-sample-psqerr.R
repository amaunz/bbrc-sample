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


# Probabilistic squared error
# Integrates over squaredError(x) * p(x), with X~Poisson(lambda)
# Integral approximation: rectangular in limits from, to
# 
# Params: Integral bounds (from, to)
#         Poisson parameter (lambda)
#         Expected value of x (spx)

# Rectangular integration with step width 1
rectIntegrate = function(f, from, to, ...) { 
  n=seq(floor(from), ceiling(to))
  sum(do.call(f, c(list(x=n, ...))))
}

# Inner integrand; weighted squared error values
# Params: point to estimate (x)
wSquaredErr = function(x, lambda, spx) {
  dpois(x,lambda)*squaredErr(x,spx)
}

# Squared error (Yates correction)
squaredErr = function(x, spx) {
  (x-spx-0.5)^2 / spx
}
