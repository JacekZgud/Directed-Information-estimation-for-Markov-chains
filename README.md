
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Markov.dir.info

<!-- badges: start -->
<!-- badges: end -->

The goal of this collection is to simulate the specific class of Markov
chains and calculate relevant information theory-based measures. The
code is oriented around `MarkovProcess` S4 class object.

## Installation

You can install the development version of `Markov.dir.info` like so:

``` r
install_github("JacekZgud/Directed-Information-estimation-for-Markov-chains")
```

### Comparison to Ay Polani information flow using an example from their article.

Script for it may be found in `mutual_info_comparison.R` followed by
some other visualisations.

## Example

This is a basic example which shows you how to use the library.

### Initializing `MarkovProcess` class:

``` r
library(Markov.dir.info)

# parametrise the code
node_dim = 3
nodes = 3
work_names = tail(LETTERS, nodes)

# define parent structure
ParentStructure = matrix(nrow = nodes, ncol = nodes, data = 0)
rownames(ParentStructure) = colnames(ParentStructure) = work_names

# define the parent structure of nodes
diag(ParentStructure) = 1
ParentStructure[2, 3] = 1
ParentStructure[3, 2] = 1
ParentStructure[1, 2] = 1
ParentStructure[2, 1] = 1

#------------------------------------------------------------------------
# initialize Markov Process class 
process = MarkovProcess(node_dim, nodes, ParentStructure, work_names)
process@trans_prob
#> $X
#> Key: <X, Y>
#>        X     Y     prob_0     prob_1    prob_2
#>    <int> <int>      <num>      <num>     <num>
#> 1:     0     0 0.71940299 0.13233831 0.1482587
#> 2:     0     1 0.35377918 0.35931935 0.2869015
#> 3:     0     2 0.21518987 0.45316456 0.3316456
#> 4:     1     0 0.25820763 0.59716060 0.1446318
#> 5:     1     1 0.05003127 0.36710444 0.5828643
#> 6:     1     2 0.01873327 0.72524532 0.2560214
#> 7:     2     0 0.07771261 0.60923754 0.3130499
#> 8:     2     1 0.47762803 0.06253369 0.4598383
#> 9:     2     2 0.17537860 0.33707865 0.4875427
#> 
#> $Y
#> Key: <X, Y, Z>
#>         X     Y     Z    prob_0     prob_1     prob_2
#>     <int> <int> <int>     <num>      <num>      <num>
#>  1:     0     0     0 0.3441721 0.40520260 0.25062531
#>  2:     0     0     1 0.5580601 0.20969945 0.23224044
#>  3:     0     0     2 0.5227098 0.03849115 0.43879908
#>  4:     0     1     0 0.1682641 0.76570820 0.06602769
#>  5:     0     1     1 0.4931657 0.34390377 0.16293056
#>  6:     0     1     2 0.4786730 0.11690363 0.40442338
#>  7:     0     2     0 0.4095745 0.28073286 0.30969267
#>  8:     0     2     1 0.5769459 0.18419489 0.23885918
#>  9:     0     2     2 0.2304688 0.59765625 0.17187500
#> 10:     1     0     0 0.5318471 0.18152866 0.28662420
#> 11:     1     0     1 0.2322159 0.38511856 0.38266558
#> 12:     1     0     2 0.1947917 0.50416667 0.30104167
#> 13:     1     1     0 0.3346810 0.41554358 0.24977538
#> 14:     1     1     1 0.3647830 0.18567104 0.44954591
#> 15:     1     1     2 0.0996420 0.48031026 0.42004773
#> 16:     1     2     0 0.7203166 0.23482850 0.04485488
#> 17:     1     2     1 0.5450501 0.40044494 0.05450501
#> 18:     1     2     2 0.3258741 0.13006993 0.54405594
#> 19:     2     0     0 0.1649151 0.30719483 0.52789006
#> 20:     2     0     1 0.3087432 0.28734062 0.40391621
#> 21:     2     0     2 0.3871560 0.43486239 0.17798165
#> 22:     2     1     0 0.5068871 0.15371901 0.33939394
#> 23:     2     1     1 0.3655914 0.61146953 0.02293907
#> 24:     2     1     2 0.4235214 0.35486064 0.22161795
#> 25:     2     2     0 0.5743415 0.28632116 0.13933730
#> 26:     2     2     1 0.2109153 0.32916430 0.45992041
#> 27:     2     2     2 0.2500000 0.41766827 0.33233173
#>         X     Y     Z    prob_0     prob_1     prob_2
#> 
#> $Z
#> Key: <Y, Z>
#>        Y     Z     prob_0     prob_1    prob_2
#>    <int> <int>      <num>      <num>     <num>
#> 1:     0     0 0.23817568 0.37795608 0.3838682
#> 2:     0     1 0.41187050 0.48021583 0.1079137
#> 3:     0     2 0.49596478 0.48275862 0.0212766
#> 4:     1     0 0.60285285 0.03153153 0.3656156
#> 5:     1     1 0.28823888 0.30773918 0.4040219
#> 6:     1     2 0.44585091 0.32536334 0.2287857
#> 7:     2     0 0.06445312 0.66601562 0.2695312
#> 8:     2     1 0.48950652 0.25411231 0.2563812
#> 9:     2     2 0.28968714 0.04634994 0.6639629

# example of quick manual conditional probability changes:

process2 = MarkovProcess(node_dim, nodes, ParentStructure, work_names)
process2@trans_prob$X$prob_0 = c(0.05, 0.05, 0.90, 0.90, 0.05, 0.05, 0.90, 0.90, 0.1)
process2@trans_prob$X$prob_1 = c(0.90, 0.90, 0.05, 0.05, 0.90, 0.90, 0.05, 0.05, 0.1)
process2@trans_prob$X$prob_2 = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.8) 
process2@trans_prob$X
#> Key: <X, Y>
#>        X     Y prob_0 prob_1 prob_2
#>    <int> <int>  <num>  <num>  <num>
#> 1:     0     0   0.05   0.90   0.05
#> 2:     0     1   0.05   0.90   0.05
#> 3:     0     2   0.90   0.05   0.05
#> 4:     1     0   0.90   0.05   0.05
#> 5:     1     1   0.05   0.90   0.05
#> 6:     1     2   0.05   0.90   0.05
#> 7:     2     0   0.90   0.05   0.05
#> 8:     2     1   0.90   0.05   0.05
#> 9:     2     2   0.10   0.10   0.80
```

### Calculating transition matrix and stationary probability

``` r
#------------------------------------------------------------------------
# transition matrix

process = stationary.probability(process)

process = trans.matrix(process)
```

### Calculate transfer entropy

``` r

#------------------------------------------------------------------------
#markov process simulation

m = 10 

process = simulate(process, m)
process@simulation
#>       X Y Z
#>  [1,] 0 0 1
#>  [2,] 2 0 1
#>  [3,] 2 1 1
#>  [4,] 0 1 0
#>  [5,] 1 1 0
#>  [6,] 1 2 0
#>  [7,] 1 2 2
#>  [8,] 1 0 0
#>  [9,] 1 0 1
#> [10,] 1 2 1

#-------------------------------------------------------------------------
# simulate marginalized markov process
n_2 = 10 

process = simulate.marginalized(process, c('Y', 'Z'), n_2)

process@marg_sim
#> $sim_target
#>       Y Z
#>  [1,] 2 1
#>  [2,] 0 1
#>  [3,] 2 1
#>  [4,] 2 0
#>  [5,] 0 2
#>  [6,] 2 0
#>  [7,] 1 2
#>  [8,] 1 0
#>  [9,] 0 2
#> [10,] 0 1
#> 
#> $ft
#>      X         1         2         3         4          5         6         7
#> [1,] 0 0.3357679 0.1366435 0.2320264 0.1610724 0.09883333 0.2903232 0.1129166
#> [2,] 1 0.2849749 0.5300234 0.5616110 0.4383024 0.54991371 0.5256752 0.5597933
#> [3,] 2 0.3792571 0.3333331 0.2063626 0.4006252 0.35125296 0.1840016 0.3272901
#>              8         9        10
#> [1,] 0.1848112 0.3675508 0.3944260
#> [2,] 0.2780165 0.1573009 0.3830546
#> [3,] 0.5371723 0.4751483 0.2225193
#> 
#> $mt
#>   Y Z          1          2          3          4           5           6
#> 1 0 0 0.18187205 0.21001625 0.12448763 0.23667534 0.039431348 0.146195259
#> 2 0 1 0.13272261 0.10902350 0.14514496 0.12286275 0.407457267 0.142302486
#> 3 0 2 0.09332656 0.10999692 0.03261685 0.12395974 0.164894730 0.006271690
#> 4 1 0 0.11422872 0.14724425 0.13532255 0.16425859 0.016941612 0.215149061
#> 5 1 1 0.08335934 0.07643734 0.15777782 0.08526981 0.175063325 0.209420240
#> 6 1 2 0.05861579 0.07711982 0.03545569 0.08603115 0.070846742 0.009229767
#> 7 2 0 0.14975014 0.13224602 0.15206033 0.08857259 0.008080164 0.134620463
#> 8 2 1 0.10928138 0.06865147 0.17729304 0.04597975 0.083495032 0.131035895
#> 9 2 2 0.07684339 0.06926443 0.03984113 0.04639028 0.033789778 0.005775138
#>             7          8          9          10
#> 1 0.038880788 0.11076873 0.23898873 0.201718698
#> 2 0.401768142 0.08083439 0.01250003 0.196347490
#> 3 0.162592386 0.05684032 0.14494086 0.008653613
#> 4 0.016605080 0.17754550 0.20473668 0.148827845
#> 5 0.171585826 0.12956527 0.01070852 0.144864974
#> 6 0.069439425 0.09110642 0.12416783 0.006384626
#> 7 0.008967257 0.15753668 0.15912744 0.145418240
#> 8 0.092661657 0.11496368 0.00832298 0.141546157
#> 9 0.037499439 0.08083901 0.09650693 0.006238356
#> 
#> $target_name
#> [1] "Y" "Z"
#> 
#> $sim_length
#> [1] 10
#calculate transfer entropy from V/target -----> target
process = trans_entropy(process, c('Y', 'Z'), sim.length = n_2)

process@trans_entropy_table
#>   target origin transfer.entropy
#> 1   Y, Z      X       -0.2476936
```
