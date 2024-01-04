PCA for the people
================
Karl Rohe
2024-01-03

The current approach to PCA (principal components analysis) is *matrix
first*. This note begins to explore an alternative path, one that is
*model first*. Essentially, any discussion of matrices is an
intermediate step in current workflows. So, this is code that makes all
of that transparent. In the end, the hope is that this makes PCA legible
for folks that have not yet learned linear algebra (just like linear
models). In this way, the hope of this project is that PCA becomes more
legible and usable for a wider audience.

I am currently very inspired by this approach because (despite the fact
that I love matrices and linear algebra) I find that this *model first*
way of thinking is so much easier.

``` r
# be sure these packages are installed:
library(tidyverse)
library(irlba)
library(Matrix)
library(gdim)
# library(softImpute)

# here is my code:
source("pca_4_the_people.R")
```

To illustrate the code, consider the popular data example `nycflights13`
which contains a row for every flight departing from the 3 main New York
City airports during 2013 (LGA, JFK, and EWR). It includes things like
date, destination, and information about delays.

``` r
library(nycflights13)
flights
```

    ## # A tibble: 336,776 × 19
    ##     year month   day dep_time sched_dep_time dep_delay arr_time sched_arr_time
    ##    <int> <int> <int>    <int>          <int>     <dbl>    <int>          <int>
    ##  1  2013     1     1      517            515         2      830            819
    ##  2  2013     1     1      533            529         4      850            830
    ##  3  2013     1     1      542            540         2      923            850
    ##  4  2013     1     1      544            545        -1     1004           1022
    ##  5  2013     1     1      554            600        -6      812            837
    ##  6  2013     1     1      554            558        -4      740            728
    ##  7  2013     1     1      555            600        -5      913            854
    ##  8  2013     1     1      557            600        -3      709            723
    ##  9  2013     1     1      557            600        -3      838            846
    ## 10  2013     1     1      558            600        -2      753            745
    ## # ℹ 336,766 more rows
    ## # ℹ 11 more variables: arr_delay <dbl>, carrier <chr>, flight <int>,
    ## #   tailnum <chr>, origin <chr>, dest <chr>, air_time <dbl>, distance <dbl>,
    ## #   hour <dbl>, minute <dbl>, time_hour <dttm>

There are so many matrices “inside” of this data, but you don’t think of
them when you see this data. Many applications are like this. The data
does not look like matrix. Instead, it looks like a spreadsheet or a SQL
database or a tibble. That is how the users think about their data. *And
underneath it*, there are so many possible matrices.

## A formula to make the matrix transparent

For example, in the flights example, you could make a matrix where every
row is a day of the year and every column is a destination airport. Each
element of the matrix is the number of flights to that destination on
that day. To do that, we propose using a formula, like in `lm`.

``` r
formula = 1 ~ (month & day)*(dest)
```

This formula and the flights data will make a sparse matrix. The
left-hand side is `1` to denote we are just going to count the
occurrences over the 336,776 rows of `flights`. If there was another
column (e.g. `size`) for which it was sensible to sum multiple values
(e.g. number of people on the flight, or number of gallons of fuel
consumed), then you could do that with `size ~ (month & day)*(dest)`. On
the right hand side, there is an “interaction” in the parlance of linear
models and the notation of `lm`. The first terms in the interaction will
index the rows of a matrix. Here, there are two columns of `flights`, in
particular `month & day`. There are 365 unique combinations of these
values. Each will be a row. The second term in the interaction does the
same thing, but those terms will become the columns.

All suggestions welcome for the formula syntax. Maybe it should also
allow `(month & day) ~ dest`? Maybe the ampersand `&` is confusing and
another symbol would be better? Other thoughts? Relevant footnote: the
syntax should also allow low-rank matrix completion, e.g. predicting
arrival delay: `arr_delay ~ (month & day)*(dest)`.

#### Technical detail here. You can skip to “Let’s do the analysis” without missing anything.

Hopefully, you rarely need to access the actual sparse matrix. But
perhaps it is instructive to see it anyways:

``` r
sparse_matrix_data = make_sparse_matrix_raw(formula, flights)
```

This creates a list of three elements. First, the sparse matrix of
counts `sparse_matrix_data$A`:

``` r
dim(sparse_matrix_data$A)
```

    ## [1] 365 105

``` r
sparse_matrix_data$A[1:5,1:8]
```

    ## 5 x 8 sparse Matrix of class "dgCMatrix"
    ##                            
    ## [1,] 20 31 3 40 47 39 14 39
    ## [2,] 20 31 3 51 45 41 14 43
    ## [3,] 19 32 3 49 46 40 17 41
    ## [4,] 20 32 3 48 46 40 17 42
    ## [5,] 13 33 3 35 26 38 11 39

Second and third, tibbles that give basic information about the row and
column “universes”…

``` r
sparse_matrix_data$row_universe
```

    ## # A tibble: 365 × 3
    ##    month   day row_id
    ##    <int> <int>  <int>
    ##  1     1     1      1
    ##  2     1     2      2
    ##  3     1     3      3
    ##  4     1     4      4
    ##  5     1     5      5
    ##  6     1     6      6
    ##  7     1     7      7
    ##  8     1     8      8
    ##  9     1     9      9
    ## 10     1    10     10
    ## # ℹ 355 more rows

``` r
sparse_matrix_data$column_universe
```

    ## # A tibble: 105 × 2
    ##    dest  col_id
    ##    <chr>  <int>
    ##  1 IAH        1
    ##  2 MIA        2
    ##  3 BQN        3
    ##  4 ATL        4
    ##  5 ORD        5
    ##  6 FLL        6
    ##  7 IAD        7
    ##  8 MCO        8
    ##  9 PBI        9
    ## 10 TPA       10
    ## # ℹ 95 more rows

The universes give the “row names” and “column names,” but in a tidy
format that allows for things like multiple variables `(month & day)`
indexing the rows. Right now, the code recomputes this object for
various different functions. For bigger spreadsheets, that can be slow.
Future extensions should allow the argument `sparse_matrix_data` for
more familiar users.

# Let’s do the analysis

### Examining the matrix sparsity

The code inside `pca_4_the_people.R` contains a few helper functions.
First, if lots of rows or columns have very few non-zero elements, this
can cause “localization issues”. The matrix needs to be “dense enough”
for the PCA to find good stuff. So, `diagnose` prints some simple
diagnostics and plots the “degree distribution” for the rows and
columns. Here, “degree” is the number of non-zero elements in that row
or column.

``` r
# inspect "degree distributions" with this funciton:
diagnose(formula, flights)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    ## # A tibble: 6 × 3
    ##   measurement      dest `month & day`
    ##   <chr>           <dbl>         <dbl>
    ## 1 number_of_items   105           365
    ## 2 average_degree    297            86
    ## 3 median_degree     365            86
    ## 4 percent_le_1        2             0
    ## 5 percent_le_2        2             0
    ## 6 percent_le_3        2             0

For example, if either average degree was less than 10, then I might be
worried. `percent_le_x` gives the percent of rows/columns that have
row/col sums less than or equal to `x`. If these values are large, then
the matrix is very sparse and you might have trouble. Issues with
sparsity will likely manifest in localization; something that will be
evaluated in functions below.

### Picking k with cross-validated eigenvalues

When doing a PCA, you need to pick the number of dimensions that you
want. The way that we do this in my lab is with cross-validated
eigenvalues. It gives you a Z-score and a p-value. [Here is the arxiv
paper](https://arxiv.org/abs/2108.03336). Alex Hayes made it a proper R
package on CRAN
[gdim](https://cran.r-project.org/web/packages/gdim/gdim.pdf). For this
example, it picks `k=4`.

``` r
cv_eigs = pick_dim(formula, flights, dimMax = 10,num_bootstraps = 5) 
plot(cv_eigs)
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
cv_eigs
```

    ## Estimated graph dimension:    4
    ## 
    ## Number of bootstraps:         10
    ## Edge splitting probabaility:  0.1
    ## Significance level:       0.05
    ## 
    ##  ------------ Summary of Tests ------------
    ##   k           z        pvals         padj
    ##   1 166.7382015 0.000000e+00 0.000000e+00
    ##   2  13.0456829 3.362584e-39 3.362584e-39
    ##   3   8.2594938 7.314574e-17 7.314574e-17
    ##   4   3.9947133 3.238629e-05 3.238629e-05
    ##   5   0.4845076 3.140128e-01 3.140128e-01
    ##   6  -3.0294648 9.987751e-01 9.987751e-01
    ##   7  -5.3966457 1.000000e+00 1.000000e+00
    ##   8  -6.5058484 1.000000e+00 1.000000e+00
    ##   9  -7.2124035 1.000000e+00 1.000000e+00
    ##  10  -8.7714105 1.000000e+00 1.000000e+00

Notice that the top-line of the printout says that the estimated graph
dimension is 4. So, we will use `k=6` and see that in this example they
become harder to interpret. This is what we would expect if it was just
noise… but also, maybe they are not just noise?

### Let’s get the people some PCA

For right now, I named the function `pca_count`, because if there are
multiple rows of the data that have the same values for `(month & day)`
and also for `dest`, then the value inside the matrix is a *sum* of the
values on the left hand side of the formula. Right now, that is just
`1`. So, it *counts* how many times that entry appears. Some people call
this a cross-tab or a contingency table. When we do PCA to this matrix
of counts, some folks call that Correspondence Analysis.

``` r
embeddings = pca_count(formula, tib = flights, k = 6)
# In some settings, the verb "sum" is a more sensible than "count"... pca_sum is the identical function
# embeddings = pca_sum(formula, tib = flights, k = 6) 
```

In particular, the code takes a square root of every count. Then,
computes the normalized and regularized Laplacian `L` (using the number
of non-zero entries as the degree). Then, computes the leading `k`
singular vectors. This is all done with sparse linear algebra via the
packages `Matrix` and `irlba`.

``` r
names(embeddings)
```

    ## [1] "row_features"    "column_features" "middle_B"        "settings"

The `row_features` and `column_features` are the PC’s and loadings (I
don’t prefer those old terms). `middle_B` gives the singular values.
`settings` contains some details that are handy in later functions.

``` r
sample_n(embeddings$row_features, size = 3)
```

    ## # A tibble: 3 × 11
    ##   month   day row_id degree weighted_degree pc_1_rows pc_2_rows pc_3_rows
    ##   <int> <int>  <int>  <int>           <dbl>     <dbl>     <dbl>     <dbl>
    ## 1     5    15    227     86             967     1.03      0.231    -0.888
    ## 2    11     9     71     80             715     0.890    -0.112     1.61 
    ## 3    12    22    114     89             895     0.970    -0.706     1.59 
    ## # ℹ 3 more variables: pc_4_rows <dbl>, pc_5_rows <dbl>, pc_6_rows <dbl>

``` r
sample_n(embeddings$column_features, size=3)
```

    ## # A tibble: 3 × 10
    ##   dest  col_id degree weighted_degree pc_1_columns pc_2_columns pc_3_columns
    ##   <chr>  <int>  <int>           <dbl>        <dbl>        <dbl>        <dbl>
    ## 1 MSP       16    365            7185        1.50         0.528       -0.278
    ## 2 HNL       47    365             707        0.470       -0.267        0.455
    ## 3 SFO       12    365           13331        2.05         1.51         0.822
    ## # ℹ 3 more variables: pc_4_columns <dbl>, pc_5_columns <dbl>,
    ## #   pc_6_columns <dbl>

Notice that these features are in a wide and tidy form, making it easy
to `lubridate::make_date` (for `row_features`) and left-join with
airports (to get latitude and longitude) for `column_features`.

##### Diagnostic plots

You can `plot(embeddings)`. It makes these five plots, each described
after all plots are displayed.

``` r
plot(embeddings) 
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## Press [Enter] to continue to the next plot...

![](README_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

    ## Press [Enter] to continue to the next plot...

    ## `geom_smooth()` using formula = 'y ~ s(x, bs = "cs")'

![](README_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

    ## Press [Enter] to continue to the next plot...

![](README_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

    ## Press [Enter] to continue to the next plot...

![](README_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

These are the five plots:

1)  Screeplot: The top `k` singular values of `L`.  
2)  Better screeplot: its singular values `2:k` (because the first one
    is usually dominant and difficult to see an elbow past it).  
3)  A “localization plot” which is very similar (maybe exact?) to [this
    stuff](https://github.com/karlrohe/LocalizationDiagnostic); for each
    row (and column) compute its degree and its leverage score. Take the
    log of both. Fit a linear model `log(leverage)~log(degree)` and plot
    the residuals against `log(degree)`. If there is localization, I
    suspect that there will be a big curl on the right side.
4)  Pairs plot of `row_features`. This is the plot emphasized in the
    varimax paper. In these example plots below, we do not see very
    clear radial streaks.
5)  A pairs plot for `column_features`. In both pairs plots, if there
    are more than 1000 points, then the code samples 1000 points with
    probability proportional to their leverage scores. It will plot up
    to `k=10` dimensions. If `k` is larger, then it plots the first 5
    and the last 5.

### Interpreting the features with ggplot (for the people)

It is fun and easy to use these `row_features` and `column_features`.

First, let’s do the rows (i.e. dates). To interpret pc’s, it is best to
plot it with contextual information and/or in the *native space*. For
dates, the native space is a time series or a sequence. Let’s plot it
there. I give my interpretation after the plots.

``` r
embeddings$row_features %>% 
  mutate(dd = make_date(day = day, month=month, year = 2013)) %>% 
  select(dd, contains("pc_")) %>% 
  pivot_longer(contains("pc_")) %>% 
  ggplot(aes(x = dd, y = value)) + geom_line() + facet_wrap(~name, scales= "free") + geom_smooth()
```

    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

I always think of the first pc as the “mean”. What we see is that
flights are more or less constant throughout the year (see y-axis). I
presume that the oscillations are for the weekends; on average
(i.e. across destinations), there are more flights during the work week
and fewer flights on the weekends. The second pc gives a seasonal effect
(fewer flights in winter, more in summer); remember after `pc_1`, some
destinations will have negative values of this (i.e. more in the winter,
fewer in the summer). The third pc is positive on weekend destinations
(more flights on the weekends and fewer during the weekdays relative to
`pc_1`). The last three are harder to interpret. My uninformed guess is
that it is some artifact of airline decisions. If you have a guess, I’d
love to hear it.

Now, let’s do the columns (i.e. destinations). The “native space” for
destinations is a map. Let’s plot it there. Be sure you have `maps`
installed.

``` r
airports %>% sample_n(size = 3)
```

    ## # A tibble: 3 × 8
    ##   faa   name                   lat    lon   alt    tz dst   tzone            
    ##   <chr> <chr>                <dbl>  <dbl> <dbl> <dbl> <chr> <chr>            
    ## 1 ITH   Ithaca Tompkins Rgnl  42.5  -76.5  1099    -5 A     America/New_York 
    ## 2 Y51   Municipal Airport     43.6  -90.9  1292    -6 A     America/Chicago  
    ## 3 KGK   Koliganek Airport     59.7 -157.    269    -9 A     America/Anchorage

``` r
# first, get the lat and lon for the airports:
airport_dat = embeddings$column_features %>% 
  left_join(airports %>% select(dest=faa, lat,lon)) %>% 
  select(lat, lon, contains("_col")) %>% 
  pivot_longer(contains("pc_")) %>% 
  drop_na()
```

    ## Joining with `by = join_by(dest)`

``` r
library(maps)
```

    ## 
    ## Attaching package: 'maps'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

``` r
usa_map <- map_data("state")
p <- ggplot() + 
  geom_polygon(data = usa_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  coord_fixed(1.3, xlim = c(-125, -65), ylim = c(25, 50)) # i'm only keeping lower 48 states, dropping Anchorage and Honolulu.


p + geom_point(data = airport_dat, aes(x = lon, y = lat, size = abs(value), color = value)) +
  facet_wrap(~ name)  +
  scale_color_gradient2(low = "red", high = "blue", mid = "white")
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Here `pc_1` likely aligns with larger and smaller airports. `pc_2` is
negative on Denver and Florida and positive in Maine. Looking back at
the time series plots, I interpret this to mean that people go to Denver
(skiing) and Florida (beach) in the winter and Maine (coastline) in the
summer. `pc_3` picks up the “fly-over zone”… looking back at the time
series, folks only travel here the work week; the blue areas are more
weekend (vacation) destinations and the red areas are the fly-over. The
other pc’s are difficult for me to interpret (my guess is that they are
weird artifacts of airline things). We do see that the last three are
heavily localized on a few airports, looking back at the pairs plots you
can see this localization. Given that, my sense is that they are not so
interesting, but if I needed to make sense of them, I would print out
their most extreme elements and dig into those airports. Making this
function is a todo item.
