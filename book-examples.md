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
```

## Execution Speed, Models 1-4

It took so many seconds:

``` stata
timer list 1
timer list 2
timer list 3
timer list 4
```

That is it for now.
