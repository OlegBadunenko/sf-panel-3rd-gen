{smcl}
{right:version 1.3  31Aug2020}
{cmd:help xtsf3gpss1}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 20 22 2}{...} {phang} {bf:xtsf3gpss1} {hline 2} PSS Type I estimator of the stochastic frontier model for panel data, where some regressors are allowed to be correlated with effects. Unbalanced panels are supported{p_end} {p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 10 17 2} {cmd:xtsf3gpss1} {it:{help varlist:depvar}} {it:{help varlist:indepvars}} {ifin}, [{cmd:}{it:{help xtsf3gpss1##options:options}}] 

{synoptset 17 tabbed}{...}
{marker Specification}{...}
{synopthdr:Specification}
{synoptline}
{syntab :Model}
{synopt :{it:{help varname:depvars}}}left-hand-side variable{p_end}
{synopt :{it:{help varname:indepvars}}}right-hand-side variables. {it:indepvars} may contain factor variables; see {help fvvarlist}{p_end}

{synoptset 17 tabbed}{...}
{synopthdr :options}
{synoptline}
{syntab :Cost frontier}
{synopt :{opt cost}}fit cost frontier model; default is production frontier
model{p_end}

{syntab:Bandwidth}
{synopt :{opt gr0(#)}}grid begin{p_end}
{synopt :{opt gr1(#)}}grid end{p_end}
{synopt :{opt gri(#)}}grid step/increment{p_end}
{synopt:{opt reps(#)}}perform # bootstrap replications to find optimal bandwidth; default is reps(399){p_end}
{synopt:{opt bw(#)}}specify the bandwidth to use{p_end}

{syntab :Reporting}
{synopt :{opt lev:el(#)}}set confidence level; default as set by set level{p_end}
{synopt :{opt {ul on}nolog{ul off}}}suppress display of a log{p_end}
{synopt :{help ereturn##display_options :{it:display_options}}}further options for displaying output{p_end}

{syntab:Miscellaneous}
{synopt:{opt nodots}}prevent the display of replication dots {p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd} 
{cmd:xtsf3gpss1} fits the PSS Type I estimator of the stochastic frontier model for panel data, where some regressors are allowed to be correlated with effects. 

{pstd}
It allows using factor variables (see {help fvvarlist}). Unbalanced panels are supported. 

{pstd}The estimation goes in two stages.  First, {cmd:xtsf3gpss1} will go over the grid from {cmd:gr0} to {cmd:gr1} with a step/increment {cmd:gri} to find a bandwidth that results in the smallest MSE.  In each of these grid points, {cmd:reps} bootstrap replications will be performed.  This grid search is extremely time-consuming. If bandwidth has been found previously for this specification (and only this specification), there is an option to specify this bandwidth using the {cmd:bandwidth}. In the second step, the PSS Type I estimator is obtained using the optimal bandwidth.  The second step is very fast.

{pstd}See {help xtsf3gpss1_postestimation:xtsf3gpss1 postestimation} for features available after estimation.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang} {opt cost} specifies that frontier fit a cost frontier model

{dlgtab:Bandwidth}

{phang} {opt gr0(#)}, {opt gr1(#)} and {opt gri(#)} specify the lower limit, upper limit, and increment of the grid for search, respectively. The default values are 0.1, 2, and 0.1, which implies that {cmd:xtsf3gpss1} will go over 0.1, 0.2, 0.3,..., 2.0.

{phang} {opt bandwidth(#)} specifies bandwidth that has been found using {cmd:xtsf3gpss1} previously for this same specification.  Do not use this option for a new specification.

{dlgtab:Reporting}

{phang} {opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang} {opt nolog} prevents {cmd:xtsf3gpss1} from displaying any output on the screen (akin to using {cmd:quietly}).

{title:Example}

{pstd}Load data{p_end}

{phang2}{cmd:. webuse xtfrontier1, clear}{p_end}

{pstd}PSS1 estimator when optimal bandwidth is chosen on the [0.1 0.6] grid with interval of 0.5. In each grid point, 9 bootstrap replications are used. In reality large number such as for example 499 should be used.{p_end}

{phang2}{cmd:. xtsf3gpss1 lnwidgets lnmachines lnworkers, gr0(0.1) gr1(0.6) gri(0.5) reps(9)}{p_end}

{pstd}Using factor variables for interaction{p_end}

{phang2}{cmd:. xtsf3gpss1 lnwidgets c.lnmachines##c.lnworkers, gr0(0.1) gr1(0.6) gri(0.05) reps(9)}{p_end}

{pstd}List the stored results{p_end}

{phang2}{cmd:. ereturn list}{p_end}

{pstd}Using if{p_end}

{phang2}{cmd:. quietly summarize lnwidgets}{p_end}

{phang2}{cmd:. xtsf3gpss1 lnwidgets c.lnmachines##c.lnworkers if lnwidgets < r(mean), gr0(0.1) gr1(0.6) gri(0.05) reps(9)}{p_end}

{pstd}Suppose you noted that the optimal bandwidth{p_end}

{phang2}{cmd:. display `e(bandwidth)'}{p_end}

{pstd}which is 0.1 (could be different in your case as too few `reps` are used); you did something else and now wish to look at the estimation results again:{p_end}

{phang2}{cmd:. xtsf3gpss1 lnwidgets c.lnmachines##c.lnworkers if lnwidgets < r(mean), b(0.1)}{p_end}

{pstd}Display the 99% confidence interval and format coefficients, standard errors, and confidence limits in the coefficient table to have 4 places after the decimal point
{p_end}

{phang2}{cmd:. xtsf3gpss1 lnwidgets c.lnmachines##c.lnworkers, gr0(0.1) gr1(0.6) gri(0.05) reps(9) level(99) cformat(%9.4f)}{p_end}

{title:Saved results}

{pstd}
{cmd:xtsf3gpss1} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sumTi})}the sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(bandwidth)}}optimal bandwidth{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}R-squared adjusted{p_end}
{synopt:{cmd:e(aic)}}AIC{p_end}
{synopt:{cmd:e(bic)}}BIC{p_end}
{synopt:{cmd:e(cp})}Mallows's Cp = RSS/shat^2-NT+2*(k+1), where k is the number of regressors{p_end}
{synopt:{cmd:e(shat})}standard error of the regression{p_end}
{synopt:{cmd:e(RSS})}RSS{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtsf3gpss1}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(predict)}}program used to implement {opt predict}{p_end}
{synopt:{cmd:e(properties)}}{opt b V}{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}vector of estimated coefficients{p_end}
{synopt:{cmd:e(V)}}estimated variance-covariance matrix{p_end}
{synopt:{cmd:e(residuals)}}residuals{p_end}
{synopt:{cmd:e(alpha)}}Individual effects, Nx1{p_end}
{synopt:{cmd:e(alpha_p)}}Individual effects, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(eff)}}Measures of technical efficiency, Nx1{p_end}
{synopt:{cmd:e(eff_p)}}Measures of technical efficiency, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}

{synoptset 20 tabbed}{...}{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Park, Byeong U., Sickles, Robin C. and Simar, Léopold 1998, "Stochastic panel frontiers: A semiparametric approach," {it:Journal of Econometrics}, 84(2), 273-301.

{marker codes}{...}
{title:Codes, Programming Kit}

{pstd} {cmd:xtsf3gpss1} is an implementation of the MATLAB code {p_end}

{pmore2}Written by Wonho Song{p_end}
{pmore2}Last Update: September 1, 2017{p_end}
{pmore2}E-mail: whsong@cau.ac.kr, whsong73@hotmail.com{p_end}

{title:Author}

{psee} Oleg Badunenko{p_end}{psee} Brunel University London{p_end}{psee}E-mail: oleg.badunenko@brunel.ac.uk {p_end}

{title:Disclaimer}
 
{pstd} This software is provided "as is" without warranty of any kind, either expressed or implied. The entire risk as to the quality and 
performance of the program is with you. Should the program prove defective, you assume the cost of all necessary servicing, repair or 
correction. In no event will the copyright holders or their employers, or any other party who may modify and/or redistribute this software, 
be liable to you for damages, including any general, special, incidental or consequential damages arising out of the use or inability to 
use the program.{p_end}

{title:Also see}

{p 7 14 2}Help: {help xtsf1g}, {help xtsf2g}, {help xtsf2gbi}, {help xtsf3gpss2}, {helpb xtsf3gpss3}, {help xtsf3gkss} (if installed){p_end}

