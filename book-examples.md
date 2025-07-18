# xtsf3g Illustration and Replication

In this article, the functionalityIn this article, the functionality of the commands **xtsf3gpss1**,
**xtsf3gpss2**, **xtsf3gpss3**, and **xtsf3gkss** are showcased. The
subset of banking is used:

``` stata
use ../../banks00_07, clear
```

# Define specifications

Here are the specification and formula for first derivatives of the cost
function with respect to input prices and outputs to check monotonicity
assumptions and compute returns to scale.

``` stata
global spetech "lny1 lny2 lnw1 trend c.half#(c.lny1#c.lny1 c.lny2#c.lny2 c.lnw1#c.lnw1) c.lny1#(c.lny2 c.lnw1) c.lny2#c.lnw1 c.trend#(c.lny1 c.lny2 c.lnw1 c.trend#c.half)"
global year_c = 2001
global itermax = 1000
```

## Model 1

`xtsf3gpss1` fits the PSS Type I estimator of the stochastic frontier
model for panel data, where some regressors are allowed to be correlated
with effects.

It allows using factor variables (see `fvvarlist`). Unbalanced panels
are supported.

The estimation goes in two stages. First, `xtsf3gpss1` will go over the
grid from `gr0` to `gr1` with a step/increment `gri` to find a bandwidth
that results in the smallest MSE. In each of these grid points,
`xtsf3gpss1` will do `reps` bootstrap replications. This grid search is
extremely time-consuming. If bandwidth has been found previously for
this specification (and only this specification), there is an option to
specify this bandwidth using the `bandwidth`. In the second step, the
PSS Type I estimator is obtained using the optimal bandwidth. The second
step is very fast.

``` stata
timer clear 1
timer on 1
xtsf3gpss1 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)
timer off 1
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M1

. timer clear 1

. timer on 1

. xtsf3gpss1 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)

Description of the panel data:------------------------------------------------

      id:  1155, 2040, ..., 3217331                          n =        500
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         2       2       4         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     57.80   57.80 |  111111
       34      6.80   64.60 |  11111.
       33      6.60   71.20 |  1111..
       27      5.40   76.60 |  111...
       24      4.80   81.40 |  11....
       15      3.00   84.40 |  11.111
       14      2.80   87.20 |  ..1111
       11      2.20   89.40 |  .11111
        5      1.00   90.40 |  .1111.
       48      9.60  100.00 | (other patterns)
 ---------------------------+---------
      500    100.00         |  XXXXXX
  
Calculating optimal bandwidth for the PSS (some regressors are correlated with effects) estimator.
Please be patient!
It may take long time to compute estimates if the data size is large.
  
Choosing the bandwidth which has the smallest MSE for the PSS estimator.
Going over the grid, which contains 9 grid points
(in each grid point, 9 bootstrap replications are used)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
.........  
  
Optimal bandwidth (for the within panel data estimator) is 0.1000

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.8214
 Adj R-squared    = 0.7763
 AIC              = -3.7229
 BIC              = -3.6885
 Root MSE         = 0.0937
-----------------------------

PSS Type I estimator: some regressors are correlated with effects
Park, Sickles, and Simar (1998), Journal of Econometrics, 84(2):273–301

 Cost Stochastic Frontier
----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |     0.2596     0.0388     6.68   0.000       0.1834      0.3357
                  lny2 |    -1.9301     0.0873   -22.10   0.000      -2.1013     -1.7590
                  lnw1 |    -0.0077     0.0470    -0.16   0.870      -0.0998      0.0845
                 trend |    -0.5000     0.0324   -15.45   0.000      -0.5634     -0.4366
                       |
  c.half#c.lny1#c.lny1 |     0.0423     0.0014    30.82   0.000       0.0396      0.0449
                       |
  c.half#c.lny2#c.lny2 |     0.2722     0.0068    39.93   0.000       0.2588      0.2855
                       |
  c.half#c.lnw1#c.lnw1 |    -0.0226     0.0042    -5.42   0.000      -0.0308     -0.0145
                       |
         c.lny1#c.lny2 |    -0.0487     0.0027   -17.76   0.000      -0.0540     -0.0433
                       |
         c.lny1#c.lnw1 |    -0.0010     0.0017    -0.57   0.572      -0.0044      0.0024
                       |
         c.lny2#c.lnw1 |     0.0155     0.0035     4.40   0.000       0.0086      0.0225
                       |
        c.trend#c.lny1 |     0.0023     0.0013     1.78   0.075      -0.0002      0.0049
                       |
        c.trend#c.lny2 |     0.0103     0.0025     4.20   0.000       0.0055      0.0151
                       |
        c.trend#c.lnw1 |    -0.0105     0.0016    -6.50   0.000      -0.0137     -0.0073
                       |
c.trend#c.trend#c.half |     0.0768     0.0024    32.53   0.000       0.0722      0.0814
----------------------------------------------------------------------------------------

. timer off 1

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -3.7229442

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -3.6885239

. eststo M1
```

## Model 2

`xtsf3gpss2` PSS Type II estimator of the stochastic frontier model for
panel data, where error follow AR(1).

It allows using factor variables (see `fvvarlist`). Unbalanced panels
are supported.

The estimation is performed in two stages. First, `xtsf3gpss2` will go
over the grid from `gr0` to `gr1` with a step/increment `gri` to find a
bandwidth that results in the smallest MSE. In each of these grid
points, `xtsf3gpss2` will do `reps` bootstrap replications. This grid
search is extremely time-consuming. If bandwidth has been found
previously for this specification (and only this specification), there
is an option to specify this bandwidth using the `bandwidth.` In the
second step, the PSS Type II estimator is obtained using the optimal
bandwidth. The second step is very fast.

``` stata
timer clear 2
timer on 2
xtsf3gpss2 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)
timer off 2
estadd scalar AIC = e(aicW)
estadd scalar BIC = e(bicW)
eststo M2

. timer clear 2

. timer on 2

. xtsf3gpss2 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)

Description of the panel data:------------------------------------------------

      id:  1155, 2040, ..., 3217331                          n =        500
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         2       2       4         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     57.80   57.80 |  111111
       34      6.80   64.60 |  11111.
       33      6.60   71.20 |  1111..
       27      5.40   76.60 |  111...
       24      4.80   81.40 |  11....
       15      3.00   84.40 |  11.111
       14      2.80   87.20 |  ..1111
       11      2.20   89.40 |  .11111
        5      1.00   90.40 |  .1111.
       48      9.60  100.00 | (other patterns)
 ---------------------------+---------
      500    100.00         |  XXXXXX

IDs with 2 or fewer observations have been excluded from estimation


      id:  1155, 2040, ..., 3217331                          n =        471
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         3       3       5         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     61.36   61.36 |  111111
       34      7.22   68.58 |  11111.
       33      7.01   75.58 |  1111..
       27      5.73   81.32 |  111...
       15      3.18   84.50 |  11.111
       14      2.97   87.47 |  ..1111
       11      2.34   89.81 |  .11111
        5      1.06   90.87 |  .1111.
        5      1.06   91.93 |  1.1111
       38      8.07  100.00 | (other patterns)
 ---------------------------+---------
      471    100.00         |  XXXXXX
  
  Calculating optimal bandwidth for the PSS (AR(1) error) estimator
  Please be patient!
  
  Going over the grid, which contains 9 grid points
  (in each grid point, 9 bootstrap replications are used)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
.........  
  
  Optimal bandwidth (for the  within  estimator) is .6
  Optimal bandwidth (for the    GLS   estimator) is .4

Sample:----------------------
 Number of obs    = 2488
 Number of groups = 471
Diagnostics:-----------------
Within:----------------------
 R-squared        = 0.8108
 Adj R-squared    = 0.7651
 AIC              = -2.2908
 BIC              = -2.2557
 Root MSE         = 0.1918
GLS:-------------------------
 R-squared        = 0.8163
 Adj R-squared    = 0.7719
 AIC              = -2.3204
 BIC              = -2.2853
 Root MSE         = 0.1890
-----------------------------

PSS Type 2 estimator: AR(1) error
Park, Sickles, and Simar (2003), Journal of Econometrics, 117(2):279–309

 Cost Stochastic Frontier
----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
Within                 |
                  lny1 |    -0.1337     0.1166    -1.15   0.251      -0.3622      0.0947
                  lny2 |    -2.9294     0.3208    -9.13   0.000      -3.5582     -2.3005
                  lnw1 |    -0.0834     0.1486    -0.56   0.575      -0.3747      0.2079
                 trend |    -0.5785     0.0263   -21.96   0.000      -0.6302     -0.5269
                       |
  c.half#c.lny1#c.lny1 |     0.0514     0.0039    13.32   0.000       0.0438      0.0589
                       |
  c.half#c.lny2#c.lny2 |     0.3360     0.0260    12.94   0.000       0.2851      0.3868
                       |
  c.half#c.lnw1#c.lnw1 |    -0.0270     0.0122    -2.22   0.027      -0.0508     -0.0031
                       |
         c.lny1#c.lny2 |    -0.0230     0.0092    -2.50   0.012      -0.0410     -0.0050
                       |
         c.lny1#c.lnw1 |     0.0061     0.0045     1.36   0.175      -0.0027      0.0150
                       |
         c.lny2#c.lnw1 |     0.0170     0.0113     1.51   0.132      -0.0051      0.0391
                       |
        c.trend#c.lny1 |     0.0030     0.0011     2.63   0.009       0.0008      0.0053
                       |
        c.trend#c.lny2 |     0.0107     0.0020     5.34   0.000       0.0068      0.0146
                       |
        c.trend#c.lnw1 |    -0.0099     0.0018    -5.57   0.000      -0.0134     -0.0064
                       |
c.trend#c.trend#c.half |     0.0913     0.0013    68.17   0.000       0.0886      0.0939
-----------------------+----------------------------------------------------------------
GLS                    |
                  lny1 |    -0.0566     0.0523    -1.08   0.279      -0.1592      0.0460
                  lny2 |    -2.8499     0.1577   -18.08   0.000      -3.1589     -2.5409
                  lnw1 |    -0.0521     0.0715    -0.73   0.466      -0.1921      0.0880
                 trend |    -0.5407     0.0121   -44.66   0.000      -0.5644     -0.5169
                       |
  c.half#c.lny1#c.lny1 |     0.0538     0.0017    31.19   0.000       0.0505      0.0572
                       |
  c.half#c.lny2#c.lny2 |     0.3374     0.0129    26.15   0.000       0.3121      0.3627
                       |
  c.half#c.lnw1#c.lnw1 |    -0.0139     0.0060    -2.33   0.020      -0.0257     -0.0022
                       |
         c.lny1#c.lny2 |    -0.0299     0.0041    -7.21   0.000      -0.0380     -0.0217
                       |
         c.lny1#c.lnw1 |    -0.0016     0.0020    -0.79   0.427      -0.0056      0.0024
                       |
         c.lny2#c.lnw1 |     0.0169     0.0055     3.10   0.002       0.0062      0.0276
                       |
        c.trend#c.lny1 |     0.0041     0.0005     7.79   0.000       0.0031      0.0051
                       |
        c.trend#c.lny2 |     0.0090     0.0009     9.63   0.000       0.0072      0.0108
                       |
        c.trend#c.lnw1 |    -0.0096     0.0008   -11.97   0.000      -0.0112     -0.0080
                       |
c.trend#c.trend#c.half |     0.0843     0.0006   143.53   0.000       0.0831      0.0854
----------------------------------------------------------------------------------------

. timer off 2

. estadd scalar AIC = e(aicW)

added scalar:
                e(AIC) =  -2.2908232

. estadd scalar BIC = e(bicW)

added scalar:
                e(BIC) =  -2.2557394

. estadd scalar shat = e(shatW)

added scalar:
               e(shat) =  .191812

. estadd scalar RSS = e(RSSW)

added scalar:
                e(RSS) =  24.312897

. eststo M2
```

## Model 3

`xtsf3gpss3` PSS Type II estimator of the stochastic frontier model for
panel data, where error follow AR(1).

It allows using factor variables (see `fvvarlist`). Unbalanced panels
are supported.

The estimation is performed in two stages. First, `xtsf3gpss3` will go
over the grid from `gr0` to `gr1` with a step/increment `gri` to find a
bandwidth that results in the smallest MSE. In each of these grid points
`reps` reps bootstrap replications will be performed. This grid search
is extremely time-consuming. If bandwidth has been found previously for
this specification (and only this specification), there is an option to
specify this bandwidth using the `bw`. In the second step, the PSS Type
III estimator is obtained using the optimal bandwidth. The second step
is very fast.

``` stata
timer clear 3
timer on 3
xtsf3gpss3 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)
timer off 3
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M3

. timer clear 3

. timer on 3

. xtsf3gpss3 lnc $spetech if year > $year_c, cost gr0(0.1) gr1(0.9) gri(0.1) reps(9)

Description of the panel data:------------------------------------------------

      id:  1155, 2040, ..., 3217331                          n =        500
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         2       2       4         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     57.80   57.80 |  111111
       34      6.80   64.60 |  11111.
       33      6.60   71.20 |  1111..
       27      5.40   76.60 |  111...
       24      4.80   81.40 |  11....
       15      3.00   84.40 |  11.111
       14      2.80   87.20 |  ..1111
       11      2.20   89.40 |  .11111
        5      1.00   90.40 |  .1111.
       48      9.60  100.00 | (other patterns)
 ---------------------------+---------
      500    100.00         |  XXXXXX

IDs with 4 or fewer observations have been excluded from estimation


      id:  1155, 2040, ..., 3153297                          n =        359
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         5       5       6         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     80.50   80.50 |  111111
       34      9.47   89.97 |  11111.
       15      4.18   94.15 |  11.111
       11      3.06   97.21 |  .11111
        5      1.39   98.61 |  1.1111
        4      1.11   99.72 |  1111.1
        1      0.28  100.00 |  111.11
 ---------------------------+---------
      359    100.00         |  XXXXXX
  
Calculating optimal bandwidth for the PSS (dynamic panel data model) estimator.
Please be patient!
It may take long time to compute estimates if the data size is large.
  
Choosing the bandwidth which has the smallest MSE for the PSS estimator.
Going over the grid, which contains 9 grid points
(in each grid point, 9 bootstrap replications are used)
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
.........  
  
Optimal bandwidth (for the dynamic panel data estimator) is 0.5000

Sample:----------------------
 Number of obs    = 2084
 Number of groups = 359
Diagnostics:-----------------
 R-squared        = 0.9066
 Adj R-squared    = 0.8864
 AIC              = -3.1583
 BIC              = -3.1177
 Root MSE         = 0.1242
-----------------------------

PSS Type 3 estimator: dynamic panel data model
Park, Sickles, and Simar (2007), Journal of Econometrics, 136(1):281–301

 Cost Stochastic Frontier
----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |     0.4531     0.1484     3.05   0.002       0.1622      0.7440
                  lny2 |    -2.8233     0.5082    -5.56   0.000      -3.8193     -1.8272
                  lnw1 |     0.0246     0.2276     0.11   0.914      -0.4214      0.4706
                 trend |    -0.5087     0.0384   -13.23   0.000      -0.5840     -0.4334
                       |
  c.half#c.lny1#c.lny1 |     0.0088     0.0047     1.87   0.062      -0.0004      0.0179
                       |
  c.half#c.lny2#c.lny2 |     0.3408     0.0416     8.18   0.000       0.2592      0.4224
                       |
  c.half#c.lnw1#c.lnw1 |    -0.0605     0.0177    -3.42   0.001      -0.0951     -0.0259
                       |
         c.lny1#c.lny2 |    -0.0500     0.0119    -4.20   0.000      -0.0734     -0.0267
                       |
         c.lny1#c.lnw1 |     0.0134     0.0061     2.21   0.027       0.0015      0.0252
                       |
         c.lny2#c.lnw1 |     0.0005     0.0173     0.03   0.976      -0.0334      0.0345
                       |
        c.trend#c.lny1 |     0.0028     0.0016     1.70   0.088      -0.0004      0.0060
                       |
        c.trend#c.lny2 |    -0.0010     0.0030    -0.34   0.734      -0.0068      0.0048
                       |
        c.trend#c.lnw1 |     0.0042     0.0026     1.65   0.098      -0.0008      0.0093
                       |
c.trend#c.trend#c.half |     0.0866     0.0017    49.98   0.000       0.0832      0.0900
                       |
                   lnc |
                   L1. |     0.5126          .        .       .            .           .
----------------------------------------------------------------------------------------

. timer off 3

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -3.1583365

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -3.1177268

. eststo M3
```

## Model 4

`xtsf3gkss` fits the KSS estimator of the stochastic frontier model for
panel data, where arbitrary temporal heterogeneity is allowed.

``` stata
timer clear 4
timer on 4
xtsf3gkss lnc $spetech if year > $year_c, cost gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)
timer off 4
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M4

. timer clear 4

. timer on 4

. xtsf3gkss lnc $spetech if year > $year_c, cost gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)

Description of the panel data:------------------------------------------------

      id:  1155, 2040, ..., 3217331                          n =        500
    year:  2002, 2003, ..., 2007                             T =          6
           Delta(year) = 1 unit
           Span(year)  = 6 periods
           (id*year uniquely identifies each observation)

Distribution of T_i:   min      5%     25%       50%       75%     95%     max
                         2       2       4         6         6       6       6

     Freq.  Percent    Cum. |  Pattern
 ---------------------------+---------
      289     57.80   57.80 |  111111
       34      6.80   64.60 |  11111.
       33      6.60   71.20 |  1111..
       27      5.40   76.60 |  111...
       24      4.80   81.40 |  11....
       15      3.00   84.40 |  11.111
       14      2.80   87.20 |  ..1111
       11      2.20   89.40 |  .11111
        5      1.00   90.40 |  .1111.
       48      9.60  100.00 | (other patterns)
 ---------------------------+---------
      500    100.00         |  XXXXXX
  
Select L such that Cl is less than 2.33
Please be patient!
It may take long time to compute estimates if the data size is large.
  
Going over the grid, which contains 5 grid points
----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 
.
Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.5480
 Adj R-squared    = 0.4338
 AIC              = -1.8907
 BIC              = -1.8563
 Root MSE         = 0.2343
-----------------------------

Kneip-Sickles-Song Estimator
Kneip, Sickles, and Song (2012), Econometric Theory, 28(3):590–628

 Cost Stochastic Frontier
----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |    -0.0813     0.1752    -0.46   0.643      -0.4247      0.2621
                  lny2 |    -2.7225     0.6251    -4.36   0.000      -3.9476     -1.4974
                  lnw1 |    -0.2761     0.2035    -1.36   0.175      -0.6749      0.1227
                 trend |    -0.4539     0.0857    -5.30   0.000      -0.6218     -0.2859
                       |
  c.half#c.lny1#c.lny1 |     0.0121     0.0054     2.23   0.026       0.0014      0.0228
                       |
  c.half#c.lny2#c.lny2 |     0.2557     0.0519     4.92   0.000       0.1539      0.3575
                       |
  c.half#c.lnw1#c.lnw1 |    -0.0396     0.0145    -2.72   0.007      -0.0681     -0.0111
                       |
         c.lny1#c.lny2 |     0.0052     0.0150     0.34   0.731      -0.0243      0.0346
                       |
         c.lny1#c.lnw1 |     0.0203     0.0049     4.13   0.000       0.0107      0.0299
                       |
         c.lny2#c.lnw1 |     0.0258     0.0175     1.47   0.141      -0.0085      0.0600
                       |
        c.trend#c.lny1 |    -0.0115     0.0033    -3.49   0.000      -0.0180     -0.0051
                       |
        c.trend#c.lny2 |     0.0358     0.0064     5.60   0.000       0.0233      0.0483
                       |
        c.trend#c.lnw1 |     0.0008     0.0046     0.16   0.869      -0.0082      0.0097
                       |
c.trend#c.trend#c.half |     0.0400     0.0042     9.44   0.000       0.0317      0.0482
----------------------------------------------------------------------------------------

. timer off 4

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -1.8907164

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -1.8562961

. eststo M4
```

## Results of Models 1-4

Use *estout* command for this:

``` stata
estout M1 M2 M3 M4, ///
 cells(b(star fmt(%9.4f)) se(par)) ///
 stats(AIC BIC shat RSS ll N sumTi, ///
 labels("AIC" "BIC" "\$\hat\sigma\$" "RSS" "log-likelihood" "\$N\$" "\$\sum T_{i}\$") ///
 fmt(%9.4f %9.4f %9.4f %9.2f %9.2f %9.0f %9.0f)) ///
 starlevels(* 0.10 ** 0.05 *** 0.01) ///
 varlabels(_cons Constant ) ///
 substitute("_ " "Frontier") ///
 legend label collabels(none) mlabels(none) replace

. estout M1 M2 M3 M4, ///
>  cells(b(star fmt(%9.4f)) se(par)) ///
>  stats(AIC BIC shat RSS N sumTi, ///
>  labels("AIC" "BIC" "\$\hat\sigma\$" "RSS" "\$N\$" "\$\sum T_{i}\$") ///
>  fmt(%9.4f %9.4f %9.4f %9.2f %9.0f %9.0f)) ///
>  starlevels(* 0.10 ** 0.05 *** 0.01) ///
>  varlabels(_cons Constant ) ///
>  substitute("_ " "Frontier") ///
>  legend label collabels(none) mlabels(none) replace

------------------------------------------------------------------------------------
main                                                                                
lny1                       0.2596***      -0.1337          0.4531***      -0.0813   
                         (0.0388)        (0.1166)        (0.1484)        (0.1752)   
lny2                      -1.9301***      -2.9294***      -2.8233***      -2.7225***
                         (0.0873)        (0.3208)        (0.5082)        (0.6251)   
lnw1                      -0.0077         -0.0834          0.0246         -0.2761   
                         (0.0470)        (0.1486)        (0.2276)        (0.2035)   
trend                     -0.5000***      -0.5785***      -0.5087***      -0.4539***
                         (0.0324)        (0.0263)        (0.0384)        (0.0857)   
half # lny1 # lny1         0.0423***       0.0514***       0.0088*         0.0121** 
                         (0.0014)        (0.0039)        (0.0047)        (0.0054)   
half # lny2 # lny2         0.2722***       0.3360***       0.3408***       0.2557***
                         (0.0068)        (0.0260)        (0.0416)        (0.0519)   
half # lnw1 # lnw1        -0.0226***      -0.0270**       -0.0605***      -0.0396***
                         (0.0042)        (0.0122)        (0.0177)        (0.0145)   
lny1 # lny2               -0.0487***      -0.0230**       -0.0500***       0.0052   
                         (0.0027)        (0.0092)        (0.0119)        (0.0150)   
lny1 # lnw1               -0.0010          0.0061          0.0134**        0.0203***
                         (0.0017)        (0.0045)        (0.0061)        (0.0049)   
lny2 # lnw1                0.0155***       0.0170          0.0005          0.0258   
                         (0.0035)        (0.0113)        (0.0173)        (0.0175)   
trend # lny1               0.0023*         0.0030***       0.0028*        -0.0115***
                         (0.0013)        (0.0011)        (0.0016)        (0.0033)   
trend # lny2               0.0103***       0.0107***      -0.0010          0.0358***
                         (0.0025)        (0.0020)        (0.0030)        (0.0064)   
trend # lnw1              -0.0105***      -0.0099***       0.0042*         0.0008   
                         (0.0016)        (0.0018)        (0.0026)        (0.0046)   
trend # trend # half       0.0768***       0.0913***       0.0866***       0.0400***
                         (0.0024)        (0.0013)        (0.0017)        (0.0042)   
L.lnc                                                      0.5126                   
                                                              (.)                   
------------------------------------------------------------------------------------
GLS                                                                                 
lny1                                      -0.0566                                   
                                         (0.0523)                                   
lny2                                      -2.8499***                                
                                         (0.1577)                                   
lnw1                                      -0.0521                                   
                                         (0.0715)                                   
trend                                     -0.5407***                                
                                         (0.0121)                                   
half # lny1 # lny1                         0.0538***                                
                                         (0.0017)                                   
half # lny2 # lny2                         0.3374***                                
                                         (0.0129)                                   
half # lnw1 # lnw1                        -0.0139**                                 
                                         (0.0060)                                   
lny1 # lny2                               -0.0299***                                
                                         (0.0041)                                   
lny1 # lnw1                               -0.0016                                   
                                         (0.0020)                                   
lny2 # lnw1                                0.0169***                                
                                         (0.0055)                                   
trend # lny1                               0.0041***                                
                                         (0.0005)                                   
trend # lny2                               0.0090***                                
                                         (0.0009)                                   
trend # lnw1                              -0.0096***                                
                                         (0.0008)                                   
trend # trend # half                       0.0843***                                
                                         (0.0006)                                   
------------------------------------------------------------------------------------
AIC                       -3.7229         -2.2908         -3.1583         -1.8907   
BIC                       -3.6885         -2.2557         -3.1177         -1.8563   
$\hat\sigma$               0.0937          0.1918          0.1242          0.2343   
RSS                     512068.91           24.31       456358.91        1.26e+06   
$N$                           500             471             359             500   
$\sum T_{i}$                 2546            2488            2084            2546   
------------------------------------------------------------------------------------
* p<0.10, ** p<0.05, *** p<0.01
```

## Execution Speed, Models 1-4

It took so many seconds:

``` stata
timer list 1
timer list 2
timer list 3
timer list 4

. timer list 1
   1:     22.09 /        1 =      22.0870

. timer list 2
   2:      4.75 /        1 =       4.7540

. timer list 3
   3:      1.79 /        1 =       1.7900

. timer list 4
   4:     21.82 /        1 =      21.8250
``

That is it for now.
