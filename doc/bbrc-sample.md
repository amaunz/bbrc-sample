1.1

Out-Of-Bag Predictive Graph Mining

Andreas Maunz\
Institute for Physics, Hermann-Herder-Str. 3, 79104 Freiburg, Germany

Abstract
========

Class-correlated subgraph descriptors are calculated repeatedly on
bootstrap samples of a database of molecular graphs, where each database
graph is associated with one from a finite set of target classes. The
association between descriptors and the target class are aggregated over
the *out-of-bag* instances. Accordingly, the process follows the
out-of-bag estimation approach, however, the application to graph mining
is a novelty. In fact, class-correlated graph mining is an unstable,
supervised process and should therefore be seen in close context to
model learning. The implementation presented is completely parallelized:
several bootstrap samples are processed at a time and matched onto the
out-of-bag instances. A probabilistic chi-square test is performed on
the most frequently sampled descriptors (patterns) using a poisson
maximum likelihood estimate. Patterns surviving the test are mapped back
on the database molecules, using a parallelized molecular fragment
matching method. Inside each sample, the efficient BBRC algorithm is
employed for mining the graphs.

Introduction
============

A deficiency of any data mining approach is that of selection bias: no
training set is a faithful representation of true state of affairs. For
the case of target-class correlated graph mining, this means that there
is no way of inferring the exact relationships between chemical
fragments and chemical reactivity. Moreover, class-correlated graph
mining is an unstable process, where slight changes in the training data
lead to quite different results. It has been shown that estimates
obtained from out-of-bag instances (the latter obtained via bootstrap
sampling) can drastically improve on estimates obtained from a single
model, built using the whole training data, in the case of unstable
prediction models \cite{breiman96oob}.

In contrast to traditional out-of-bag-estimation, however, the focus is
not on estimating the quality of certain model characteristics (such as
node class probability), but instead the frequencies of class-correlated
subgraphs on the target classes. The latter can be seen as a sort of
model statistics, too, in that they are the outputs of supervised
selection processes. However, no bagged predictor is built, since graph
mining is about descriptor generation, and thus is a preprocessing step
to model building.

Methods
=======

Bootstrapping Graph Databases
-----------------------------

Let $G$ a set of graphs. A *graph database* consists of $G$ and a
function $g: G \rightarrow H$, $H \subset \mathbb{N}_0$, mapping graphs
to a finite set of *target classes*. Define the class sizes as
$h_i=\vert\{G_j \in G \; \vert\; g(G_i)=i\}\vert$, $\forall i \in H$.

Consider a pattern generating process $F: (G,g) \rightarrow (X,k)$,
$k=\left(k_1,\ldots,k_{\vert H\vert}\right)$. $X$ are called patterns,
$k_i: X \rightarrow \mathbb{N}_0$ are referred to as *pattern support
functions* on the target classes.

Run non-parametric bootstrapping, by drawing $n$ samples with
replacement from $G$. Ensure that each sample comprises exactly $h_i$
graphs $G$ with $g(G)=i$, drawn with uniform probability $1/h_i$ inside
each class $i \in H$ (stratification).

The result of running $F$ on the bootstrap samples is a pattern set
$X= X^{(1)}\cup\ldots\cup,X^{(n)}$, and sets of function vectors
$k^{(1)},\ldots,k^{(n)}$.

Frequencies
-----------

Define

-   the frequency of support value $j$ on class $i$ as

    $$w_{i,j}:=\sum_{l=1}^n \delta_{k_i^{(l)}(x),\, j}\, , \; \forall i \in H, j \in \mathbb{N}_0$$

-   the frequency of support value $j$ on class $i$, given the sum of
    support value frequencies equals $h$, as

    $$w_{i,j,h}:=\sum_{l=1}^n \delta_{k_i^{(l)}(x),\, j} * \delta_{\sum_{i=1}^{\vert H \vert}k_i^{(l)}(x),\, h}\, , \; \forall i \in H, \; j,h \in \mathbb{N}_0$$

-   the cumulative frequency of support value $j$ over classes as

    $$w_j(x):=\sum_{i \in H}w_{i,j}(x)$$

The remainder of this section drops dependency on $x$ for better
readability.

Probability Distributions
-------------------------

Consider the finite set of support values $J \subset \mathbb{N}_0$ with
$w_j > 0, \; j \in J$, and the categorical distribution for support
value $j \in J$ with parameter $W_J:=\sum_{j \in J} w_j$:

$$\begin{aligned*}
  p(j|W_J) &= w_j/W_J\end{aligned*}$$

Consider the Poisson distribution for support value $k$ in class $i$,
given the sum of support values equals $j$.

$$\begin{aligned*}
  p(k\vert\lambda_{i,j}) &= \operatorname{Pois}\left(\lambda_{i,j}\right), \text{where}\\
  \lambda_{i,j}&=\sum_{l \in \mathbb{N}} l*w_{i,l,j} / \sum_{l \in \mathbb{N}} w_{i,l,j}\end{aligned*}$$

Parameter $\lambda_{i,j}$ is the sample mean, which is in turn the
Maximum Likelihood Estimate of the Poisson distribution.

Significance Testing
--------------------

[ss:significanceTesting] Consider the probabilistic version of the
$\chi^2$ distribution test, defined as

$$\begin{aligned*}
  \chi^2 = \sum_{i \in H}\; \left( \sum_{j \in J}p(j|W_j) \int_0^{\infty}dk\; p(k|\lambda_{i,j}) \frac{(k-E_j(k_i)-0.5)^2}{E_j(k_i)} \right)\end{aligned*}$$

where

$$\begin{aligned*}
  E_j(k_i) = \frac{j * h_i}{\vert G\vert}\end{aligned*}$$

the expected support value on class $i$ when the sum of support values
is $j$. Note that it is derived from $j$, the current support value, and
the target class size. Then,
$p(\chi^2)=\operatorname{Chi-Square}\left(\vert H\vert-1\right)$, i.e.
the degrees of freedom equals the number of target classes minus one.

Algorithm
=========

Implementation
--------------

Patterns mined by BBRC are canonical, therefore a hash table can be used
to gather results. Here, “canonical” means that two isomorphic subgraphs
will be always represented by the same pattern string (so called *SMARTS
strings*, which describe molecular fragments). The algorithm using $n$
bootstrap samples is shown in Algorithm [alg:bbrc-sample].

Calculate BBRC descriptors on bootstrap samples of a database
[alg:bbrc-sample]

[1] $n, minSamplingSupport, minFrequencyPerSample$ $hash \gets \{\}$
$i:=1 \to n$ Done in parallel
$res \gets BBRC(drawSample(n), minFrequencyPerSample)$
$insert(hash,res)$ $ans \gets \emptyset$ $pattern \in keys(hash)$
$length(hash[pattern]) \geq minSamplingSupport$
$(p=SignificanceTest(hash[pattern]))>1-\alpha$
$ans\gets ans \cup (pattern,p)$ $ans$

Line 1 creates an initially empty hash table to gather results from BBRC
mining in line 3. Importantly, the result $res$ consists of patterns and
support values *per class*, i.e. each pattern is used as key in the
hash, where the values stored are the class-specific support values.
Importantly, these values correspond to the *out-of-bag* instances, not
the *in-bag* instances (bootstrap samples). The necessary step of
matching the patterns mined from the *in-bag* instances (line 3) on the
*out-of-bag* instances is not shown. It is understood to happen inside
the BBRC step and implemented using a parallelized molecular fragment
matching method. On termination of the loop in line 5, each hash entry
has at most $n$ support values per class.

Post-processing the results is very fast and incurs negligible overhead
compared to the graph mining step. It consists of removing patterns that
(line 8) were not generated often enough by the sampling process (have
not enough entries in the hash table), or which (line 9) do not
significantly deviate from the overall distribution of classes, as
assessed by the probabilistic $\chi^2$ test described in section
[ss:significanceTesting].

For better readability, the listing does not show how the results in
$ans$ are used further. They are processed as follows: the caller of
Algorithm [alg:bbrc-sample] matches the patterns back onto the graphs of
the original database ($G$), which yields an instantiation matrix, with
compounds in the rows and patterns (molecular fragments) in the columns.
This matching is again done in parallel, and matrix entries can either
be of type binary (occurrence vs no occurrence) or frequency (how many
times the pattern occurs).
