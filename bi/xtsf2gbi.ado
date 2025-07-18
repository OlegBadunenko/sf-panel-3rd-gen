*! version 1.0.0  11Apr2020
*! version 1.1.0  12Apr2020
*! version 1.2.0  15Apr2020
*! version 1.3.0  16Apr2020
*! version 1.3.1  12Aug2020
*! version 1.4.0  18Sep2020

/*
Origin:

MATLAB code

% Originally written by Joe J. Qian (jhqian@sjtu.edu.cn)
% Modified by Liu Junrong (liujunrong_529@hotmail.com)
% Updated Sep 2015 by Wonho Song (whsong73@hotmail.com)

Translated by Oleg Badunenko
oleg.badunenko@brunel.ac.uk

*/


// if(c(MP)){
// 	set processors 1
// }


capture program drop xtsf2gbi
program define xtsf2gbi, eclass
version 11

  if !replay() {
	syntax varlist(numeric fv min=2) [if] [in]                  ///
		[, Distribution(string) ITERate(integer 102) TRACElevel(string)] ///
    [COST LEVEL(string) NOLOG ] ///
    [noCI] [noPValues] [noOMITted] [noEMPTYcells] ///
    [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] ///
    [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] ///
    [SFORMAT(passthru)] [nolstretch]
	marksample touse
	// display 100

	// handle the lists
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'
	
	_rmcoll `indepvars' if `touse', expand `constant'
	local indepvars `r(varlist)'
	
	tempname b V n1 nt1 dist R2 R2adj aic bic Ebie alpha rez fitted ///
  shat1 RSS1 mypanelvar mytimevar nrep convstatus Cp1
          
	
	//display 101
  if "`distribution'" == ""{
    di as error "distribution() must be specified"
		exit 198
  }
  if "`distribution'" != "texponential" & "`distribution'" != "thnormal"  & /// 
   "`distribution'" != "dthnormal"{
    di as error "distribution() is not appropriate"
		exit 198
  }
  else if "`distribution'" == "texponential" {
    local dist 0
  }
  else if "`distribution'" == "thnormal" {
    local dist 1
  }
  else if "`distribution'" == "dthnormal" {
    local dist 2
  }
  //display "|`dist'|"
  //display "|`iterate'|"
  
  // handle production/cost function
	if "`cost'" == "" { 
		local myprod = -1 
		local function "Production" 
	}
	else {
		local myprod = 1
		local function "Cost"
	}
  
  if `iterate' < 0  {
		display as error " iterate() must be greater than 0"
		exit 198
	}
   else {
    local iter = `iterate'
  }
  //display "|`iter'|"
   
  if "`tracelevel'" == ""{
    local tracelevel "value"
  }
  if "`tracelevel'" != "none" & "`tracelevel'" != "value"  & ///
   "`tracelevel'" != "tolerance" & "`tracelevel'" != "step" & ///
   "`tracelevel'" != "paramdiffs" & "`tracelevel'" != "params" & ///
   "`tracelevel'" != "gradient" & "`tracelevel'" != "hessian" {
    di as error "tracelevel() can be specified as follows"
    di as input "none, value, tolerance, step, paramdiffs, params, gradient, or hessian"
    di as text "see help for " as result "mf_optimize##i_tracelevel" as text " for details"
		exit 198
  }

	//display 102
  display
	display as result "Description of the panel data:" as input "{hline 48}
	xtdescribe if `touse'
  quietly xtset
  local mypanelvar `r(panelvar)'
  local mytimevar `r(timevar)'
  //display "|`mytimevar'|"
	
	//display 104

	mata: bie_robyty("`depvar'", "`indepvars'", "`touse'", "`mypanelvar'",     ///
      "`mytimevar'", "`cost'", "`b'", "`V'", "`n1'", "`nt1'", "`dist'",      ///
      "`R2'", "`R2adj'",  "`aic'", "`bic'",   "`Ebie'",  "`alpha'",          ///
      "`iter'", "`tracelevel'", "`convstatus'",                              ///
      "`rez'", "`fitted'", "`shat1'", "`RSS1'", "`Cp1'")
	
	//display 105
  
  
  //display "`indepvars'"
	
	matrix colnames `b'        = `indepvars'
  //matrix list `b'
	//display 106
  matrix rownames `V'        = `indepvars'
	matrix colnames `V'        = `indepvars'
  //matrix list `V'
	
	//display 107
	//matrix colnames `alpha'   = "alpha"
	//display 108

	matrix colnames `rez'     = "residuals"
	//display 113
	matrix colnames `Ebie'   = "Efficiency"
	//display 115
	ereturn post `b' `V', esample(`touse') buildfvinfo
	//display 118
	ereturn scalar r2         = `R2'
	//display 119
	ereturn scalar r2_a       = `R2adj'
	//display 122
	ereturn scalar aic        = `aic'
	//display 123
	ereturn scalar bic        = `bic'
  //display 123
	ereturn scalar N           = `n1'
	//display 127
	ereturn scalar NT          = `nt1'
  //display `shat1'
  ereturn scalar shat        = `shat1'
  //display 128
  //display `RSS1'
  ereturn scalar RSS         = `RSS1'
  ereturn scalar cp          = `Cp1'
  ereturn scalar converged   = `convstatus'
	//display 129
	ereturn matrix eff         = `Ebie'
	//display 132
	ereturn matrix alpha       = `alpha'
	//display 135
	ereturn matrix residuals   = `rez'
	//display 137
	ereturn matrix xb          = `fitted'
	//di 10
	ereturn local predict "xtsf2gbi_p"
	ereturn local cmd   "xtsf2gbi"
  ereturn local cmdline "`0'"

  }
	if replay() {
//     display "replay here"
    syntax, [LEVel(real `c(level)')] [noCI] [noPValues] [noOMITted] [noEMPTYcells] [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] [SFORMAT(passthru)] [nolstretch]
  }
  if "`nolog'" == "" {
    display
    display as result "Sample:" as input "{hline 22}
    display as input " Number of obs    " as text "= " as result `nt1'
    display as input " Number of groups " as text "= " as result `n1'
    display as result "Diagnostics:" as input "{hline 17}
    display as input " R-squared        " as text "= " as result  %5.4f `R2'
    display as input " Adj R-squared    " as text "= " as result  %5.4f `R2adj'
    display as input " AIC              " as text "= " as result  %5.4f `aic'
    display as input " BIC              " as text "= " as result  %5.4f `bic'
    display as input " Root MSE         " as text "= " as result  %5.4f `shat1'
    display as input "{hline 29}"
    
    //di as text "{hline 79}"
    display as input "Stochastic frontier with bounded inefficiency"
    display as input "Pavlos Almanidis, Junhui Qian, and Robin C. Sickles (2014) in Handbook "
    display as input "'Festschrift in Honor of Peter Schmidt: Econometric Methods and Applications'"
display
    display as input " `function'" as text " Stochastic Frontier"
    ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
  }
	
end


mata mata clear

mata:

void bie_robyty( string scalar depvar,                                       ///
  string scalar indepvars, string scalar touse,  string scalar mypanelvar,   ///
  string scalar mytimevar, string scalar mycost,                             ///
  string scalar bsname,    string scalar vsname,                             ///
  string scalar n1name,    string scalar nt1name,  string scalar dist0,      ///
  string scalar R2name,    string scalar R2adjname,                          ///
  string scalar aicname,   string scalar bicname,                            ///
  string scalar Eff_name,  string scalar aname,                              ///
  string scalar iter0,     string scalar tracelevel,                         ///
  string scalar convstatus,                                                  ///
  string scalar rezname,   string scalar xbname,                             ///
  string scalar shatname,  string scalar RSSname, string scalar Cpname)
{
	//200

	real scalar p, n, nt,  Ti_max, dist
	
	real vector y, Ti0, ids0
	
	real matrix x, ids
	
	//201
	
	y         = st_data(., depvar, touse)
	x         = st_data(., indepvars, touse)
	
	//202
	Ti0       = st_data(., mytimevar, touse)
	Ti        = Ti0 :- min(Ti0) :+ 1

	if (mycost == ""){
    cost = 0
  }
  else {
    cost = 1
  } 	
	//203

	ids0      = st_data(., mypanelvar, touse)
	ids       = panelsetup(ids0, 1)
	ids       = ids, ids[,2] - ids[,1] :+ 1
	
	//205

	p         = cols(x)
	n         = rows(ids)
	nt        = rows(y)

  Ti_max    = max(ids[,3])

  dist      = strtoreal(dist0)
  iter      = strtoreal(iter0)
  //"iter"
  //iter

  // 206

  // initial values

  S  = optimize_init()
  if (dist == 0) {
    theta0 = initipanelexpt(y, x, ids, Ti, Ti_max, n, nt, p)
    optimize_init_evaluator(S, &ipanelexpt())
  }
  else if (dist == 1) {
    theta0 = initipanelhalftgam(y, x, ids, Ti, Ti_max, n, nt, p)
    optimize_init_evaluator(S, &ipanelhalftgam())
  }
  else if (dist == 2) {
    theta0 = initipaneldoublytgam(y, x, ids, Ti, Ti_max, n, nt, p)
    optimize_init_evaluator(S, &ipaneldoublytgam())
  }
  
  optimize_init_valueid(S, "log-likelihood")
  optimize_init_iterid(S, "iter")
  optimize_init_evaluatortype(S, "gf1")
	optimize_init_technique(S, "nr")
  optimize_init_tracelevel(S, tracelevel)
  optimize_init_conv_maxiter(S, iter)
  //optimize_init_singularHmethod("hybrid")
	optimize_init_argument(S, 1, y)
	optimize_init_argument(S, 2, x)
	optimize_init_argument(S, 3, ids)
	optimize_init_argument(S, 4, Ti)
	optimize_init_argument(S, 5, Ti_max)
	optimize_init_argument(S, 6, n)
	optimize_init_argument(S, 7, nt)
	optimize_init_argument(S, 8, p)
	optimize_init_params(S, theta0)
   
  if (tracelevel != "none"){
    printf("\n\n{input:Log-likelihood maximization using} {result:optimize()} \n\n")
  }
  //211
	bh        = optimize(S)
  //212
  bhx       = bh[1..p]
  //213
  xb        = x * bhx'
  //214
  res0      = y - xb
  //215
  vh = optimize_result_V_oim(S)
  //216
  bcov      = vh[|1,1 \ p, p|]
  //217
  if (dist == 0) {
    //2171
    sigma_u= bh[p+1]
    //2172
    sigma_v= bh[p+2]
    //2173
    B       = bh[(p+3)..(p+Ti_max+2)]'
    Bse     = sqrt(diagonal( vh[|(p+3),(p+3) \ (p+Ti_max+2),(p+Ti_max+2)|] ))
    //2174
    BB        = J(nt, 1, .)
    sample_all= 1::nt
    //2175
    for (t = 1; t <= Ti_max; t++) {
      sample_i0 = Ti :== t
      sample_i  = sample_all :* sample_i0
      BB[selectindex(sample_i)] = J(sum(sample_i0), 1, B[t])
    }
    //2176
    ustar = -res0 :- sigma_v^2/sigma_u
    //2177
    u     = ustar :+ sigma_v * 
     (normalden(-ustar/sigma_v) - normalden((BB:-ustar)/sigma_v)) :/ 
     (normal((BB-ustar)/sigma_v) - normal(-ustar/sigma_v))
     //2178
  }
  else if (dist == 1) {
    //2171
    sigma  = bh[p+1]
    //2172
    gamma  = bh[p+2]
    //2173
    B      = bh[(p+3)..(p+Ti_max+2)]'
    Bse     = sqrt(diagonal( vh[|(p+3),(p+3) \ (p+Ti_max+2),(p+Ti_max+2)|] ))
    //2174
    lamda  = sqrt(gamma/(1-gamma))
    //2175
    sigma_u= sigma/sqrt(1+1/lamda^2)
    //2176
    sigma_v= sqrt(sigma^2-sigma_u^2)
    //2177
    BB     = J(nt, 1, .)
    sample_all= 1::nt
    //2178
    for (t = 1; t <= Ti_max; t++) {
      sample_i0 = Ti :== t
      sample_i  = sample_all :* sample_i0
      BB[selectindex(sample_i)] = J(sum(sample_i0), 1, B[t])
    }
    ustar  = -res0*sigma_u^2/sigma^2
    sigmastar = sigma_u*sigma_v/sigma
    u      = ustar + sigmastar *
     (normalden(-ustar/sigmastar)-normalden((BB-ustar)/sigmastar)) :/
     (normal((BB-ustar)/sigmastar) - normal(-ustar/sigmastar))
  }
  else if (dist == 2) {
    //2171
    sigma  = bh[p+1]
    //2172
    gamma  = bh[p+2]
    mu     = bh[p+3]
    //2173
    B      = bh[(p+4)..(p+Ti_max+3)]'
    Bse     = sqrt(diagonal( vh[|(p+4),(p+4) \ (p+Ti_max+3),(p+Ti_max+3)|] ))
    //2174
    lamda  = sqrt(gamma/(1-gamma))
    //2175
    sigma_u= sigma/sqrt(1+1/lamda^2)
    //2176
    sigma_v= sqrt(sigma^2-sigma_u^2)
    //2177
    BB     = J(nt, 1, .)
    sample_all= 1::nt
    //2178
    for (t = 1; t <= Ti_max; t++) {
      sample_t0 = Ti :== t
      sample_t  = sample_all :* sample_t0
      BB[selectindex(sample_t)] = J(sum(sample_t0), 1, B[t])
    }
    ustar  = (mu*sigma_v^2:-res0*sigma_u^2)/sigma^2
    sigmastar = sigma_u*sigma_v/sigma
    u      = ustar + sigmastar *
     (normalden(-ustar/sigmastar)-normalden((BB-ustar)/sigmastar)) :/
     (normal((BB-ustar)/sigmastar) - normal(-ustar/sigmastar))
  }
   
  // Calculation of Efficiency
  //218
  u         = offinf(u, ids, n, nt)
  //219
  utmin     = J(Ti_max, 1, .)
  for (t = 1; t <= Ti_max; t++) {
    sample_t0 = Ti :== t
    sample_t  = sample_all :* sample_t0
    utmin[t]  = min( u[selectindex(sample_t)] )
  }
  //220
  TE        = J(nt, 1, .)
  for (i = 1; i <= n; i++) {
    ind_ti = panelsubmatrix(Ti,  i, ids)
    TE[|ids[i,1] \ ids[i,2]|] = 
     exp(- panelsubmatrix(u,  i, ids) + utmin[selectindex(ind_ti)] )
  }
  //221
  BTE       = exp(-(B-utmin))
  //222
  BTEci     = exp(-(B+1.96*Bse-utmin)), exp(-(B-1.96*Bse-utmin))
  //223
   
  xb        = x * bhx'
  res0      = y - xb
  //res       = y - x * bhx'
  ek        = res0 :- mean(res0)

  // Diagnostics 
  ey        = y :- mean(y)
  R2        = 1 - cross(ek,ek) / cross(ey,ey)
  R2a       = 1 - (1-R2) * (nt-1) / (nt-p-n)
  
  // Score
  shat2     = variance(res0)
  shat      = sqrt(shat2)
  //998
  RSS       = cross(res0, res0)
  Cp        = RSS/shat2 - nt + 2*(p+1)
  //999
  aic       = log((nt-1)/nt*shat2)+1+2*(p+1)/nt
  //9991
	bic       = log((nt-1)/nt*shat2)+1+(p+1)*log(nt)/nt
  
  //265

	//xb        = x * bh
	//bh        = beta'
	//210
	st_matrix(bsname,        bhx)
	//211
  //bcov
	st_matrix(vsname,        bcov)
	//216
	st_numscalar(n1name,     n)
	//217
	st_numscalar(nt1name,    nt)
	//218
	st_numscalar(R2name,     R2)
	//219
	st_numscalar(R2adjname,  R2a)
	//222
	st_numscalar(aicname,    aic)
	//223
	st_numscalar(bicname,    bic)
	//227
	st_matrix(aname,         u)
	//231
	st_matrix(rezname,       res0)
	//233
	st_matrix(Eff_name,      TE)
	//235
	st_matrix(xbname,        xb)
	//236
	st_numscalar(shatname,   shat)
  //237
  st_numscalar(RSSname,    RSS)
  //238
  st_numscalar(convstatus, optimize_result_converged(S))
}


end

// begin functions for bie

/*

  OLS, very limited
  
  bh, e, sse, s2, bvar, se, y1, R2, R2bar  all must be set before executing it.

*/
  
capture mata mata drop bieols()
mata:

void bieols(real vector y, real matrix x, real scalar n, real scalar p,    ///
            bh, e, sse, s2, bvar, se, y1, R2, R2bar)
{
  real matrix xxi

  xxi       = cross(x,x)
  bh        = invsym( xxi ) * cross(x,y)
  e         = y - x*bh
  sse       = sum(e:^2)
  s2        = sse/(n - p)
  bvar      = xxi * s2
  se        = sqrt(diagonal(bvar))
  R2        = 1 - sse / sum(y:^2)
  y1        = y :- mean(y)
  R2bar     = 1 - (sse/(n-p)) / (sum(y1:^2)/(n-1))
}  
end

/*

  within, very limited
  
  bh, e, sse, s2, bvar, se, y1, R2, R2bar  all must be set before executing it.

*/

capture mata mata drop biepanelfe()
mata:

void biepanelfe(real vector y, real matrix x,  real matrix ids,             ///
                real scalar n, real scalar nt, real scalar p,              ///
                bh, e, sse, s2, bvar, se, y1, R2, R2bar)
{
  real colvector y_demean
  real matrix x_demean

  y_demean   = J(nt,1,.)
  x_demean   = J(nt,p,.)
  for (i = 1; i <= n; i++) {
    y_demean[| ids[i,1], .\ ids[i,2], .|] = 
     y[| ids[i,1], .\ ids[i,2], .|] - 
     J(ids[i,3], 1, mean( panelsubmatrix(y,  i, ids) ))
    x_demean[| ids[i,1], .\ ids[i,2], .|] = 
     x[| ids[i,1], .\ ids[i,2], .|] - 
     J(ids[i,3], 1, mean( panelsubmatrix(x,  i, ids) ))
  }
  bh = e = sse = s2 = bvar = se = R2 = y1 = R2bar = .
  bieols(y_demean, x_demean, nt-n, p, bh, e, sse, s2, bvar, se, y1, R2, R2bar)
}  
end

mata
bh1 = e1 = sse1 = s21 = bvar1 = se1 = R21 = y11 = R2bar1 = .
biepanelfe(y, x, ids, n, nt, p, bh1, e1, sse1, s21, bvar1, se1, y11, R21, R2bar1)
end

/*

  This is very approximate if unbalanced so use only for initial values in ll
  
  bh, e, sse, s2, bvar, se, y1, R2, R2bar, sigmau, sigmav  all must be set 
  before executing it.

*/

capture mata mata drop biepanelre()
mata:

void biepanelre(real vector y, real matrix x,  real matrix ids,             ///
                real scalar n, real scalar nt, real scalar p,              ///
                bh, e, sse, s2, bvar, se, y1, R2, R2bar, sigmau, sigmav)
{
  real matrix Xand1, x_q_demean
  real colvector y_q_demean, u, v

  // pooled regression
  Xand1     = x, J(nt,1,1)
  bh = e = sse = s2_po = bvar = se = R2 = y1 = R2bar = .
  bieols(y, Xand1, nt, cols(Xand1), bh, e, sse, s2_po, bvar, se, y1, R2, R2bar)
  
  // fixed effect model
  s2_fe     = .
  biepanelfe(y, x, ids, n, nt, p, bh, e, sse, s2_fe, bvar, se, y1, R2, R2bar)

  // Greene's
  thetai    = 1 :- sqrt(s2_fe :/ (s2_fe :+ (s2_po - s2_fe) * ids[,3]))
  
  y_q_demean = J(nt,1,.)
	x_q_demean = J(nt,p,.)
  for (i = 1; i <= n; i++) {
    y_q_demean[| ids[i,1], .\ ids[i,2], .|] = 
     (y[| ids[i,1], .\ ids[i,2], .|] - 
     thetai[i] * J(ids[i,3], 1, mean( panelsubmatrix(y,  i, ids) )) ) / sqrt(s2_fe)
    x_q_demean[| ids[i,1], .\ ids[i,2], .|] = 
     (x[| ids[i,1], .\ ids[i,2], .|] - 
     thetai[i] * J(ids[i,3], 1, mean( panelsubmatrix(x,  i, ids) )) ) / sqrt(s2_fe)
  }

  bh = e = sse = s2 = bvar = se = R2 = y1 = R2bar = .
  bieols(y_q_demean, x_q_demean, nt, p, bh, e, sse, s2, bvar, se, y1, R2, R2bar)

  e         = y - x * bh
  u         = J(n, 1, .)
  v         = J(nt, 1, .)
  for (i = 1; i <= n; i++) {
    ei     = panelsubmatrix(e,  i, ids)
    u[i]   = mean( ei )
    v[| ids[i,1], .\ ids[i,2], .|] = ei - J(ids[i,3], 1, u[i])
  }

  sigmau    = sqrt( variance ( u ) )
  sigmav    = sqrt( variance ( v ) )
  sse       = sum(e:^2)
  s2        = sse/(nt - p)
  R2        = 1 - sse / (sum(y:^2))
}
end

/*

  get initial values for three set-ups
  
  in all three cases, the inputs are the same;
  the outputs are slightly different

*/

capture mata mata drop initipanelexpt()
mata:

real rowvector initipanelexpt(real vector y,   real matrix x, ///
                              real matrix ids, real colvector Ti0, ///
                              real scalar Ti_max, ///
                              real scalar n, real scalar nt, real scalar p)
{
  real colvector e, B, tymch, sample_all, sample_i

  // run GLS
  bh1 = e1 = sse1 = s21 = bvar1 = se1 = R21 = y11 = R2bar1 = sigmau1 = sigmav1 = .
  biepanelre(y, x, ids, n, nt, p, bh1, e1, sse1, s21, bvar1, se1, 
   y11, R21, R2bar1, sigmau1, sigmav1)
    
  e         = y - x * bh1
  //Ti_max    = max(Ti0)
   
  B         = J(1, Ti_max, .)
  sample_all = 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i = sample_all :* (Ti0 :== t)
    et = e[selectindex(sample_i)]
    B[t]   = max (-et )
  }
  //tymch = bh1 \ log(sigmau1) \ log(sigmav1) \ log(B)
  tymch = bh1', sigmau1,  sigmav1, B
  return(tymch)
}
end

mata
theta = initipanelexpt(y, x, ids, Ti0, Ti_max, n, nt, p)
theta
end

capture mata mata drop initipanelhalftgam()
mata:

real rowvector initipanelhalftgam(real vector y,   real matrix x, ///
                              real matrix ids, real colvector Ti0, ///
                              real scalar Ti_max, ///
                              real scalar n, real scalar nt, real scalar p)
{
  real colvector e, B, tymch, sample_all, sample_i

  // run GLS
  bh1 = e1 = sse1 = s21 = bvar1 = se1 = R21 = y11 = R2bar1 = sigmau1 = sigmav1 = .
  biepanelre(y, x, ids, n, nt, p, bh1, e1, sse1, s21, bvar1, se1, 
   y11, R21, R2bar1, sigmau1, sigmav1)
    
  e         = y - x * bh1
  //Ti_max    = max(Ti0)
   
  B         = J(1, Ti_max, .)
  sample_all = 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i = sample_all :* (Ti0 :== t)
    et = e[selectindex(sample_i)]
    B[t]   = max (-et )
  }
  sigma1    = sqrt(sigmau1^2 + sigmav1^2);
  gamma1    = sigmau1^2 / sigma1^2;
  tymch     = bh1', sigma1, gamma1, B
  return(tymch)
}
end

capture mata mata drop initipaneldoublytgam()
mata:

real rowvector initipaneldoublytgam(real vector y,   real matrix x, ///
                              real matrix ids, real colvector Ti0, ///
                              real scalar Ti_max, ///
                              real scalar n, real scalar nt, real scalar p)
{
  real colvector e, B, tymch, sample_all, sample_i

  // run GLS
  bh1 = e1 = sse1 = s21 = bvar1 = se1 = R21 = y11 = R2bar1 = sigmau1 = sigmav1 = .
  biepanelre(y, x, ids, n, nt, p, bh1, e1, sse1, s21, bvar1, se1, 
   y11, R21, R2bar1, sigmau1, sigmav1)
    
  e         = y - x * bh1
  mu1       = mean(e)
   
  B         = J(1, Ti_max, .)
  sample_all = 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i = sample_all :* (Ti0 :== t)
    et = e[selectindex(sample_i)]
    B[t]   = max (-et )
  }
  sigma1    = sqrt(sigmau1^2 + sigmav1^2);
  gamma1    = sigmau1^2 / sigma1^2;
  tymch     = bh1', sigma1, gamma1, mu1, B
  return(tymch)
}
end

/*

  log-likelihood and gradient functions (it-specific) for three set-ups
  
  required for Mata's optimize

*/


capture mata mata drop ipanelexpt()
mata:
void ipanelexpt(real scalar todo, real vector theta,                         ///
	real vector y, real matrix x,                                             ///
	real matrix ids, real colvector Ti0, real scalar Ti_max,                  ///
  real scalar n, real scalar nt, real scalar p,                             ///
	ll, grad, Hess)
{
  real scalar sigu, sigv
	
  real vector  alpha, B, e, BB, sample_all, sample_i0, sample_i, A1, A2,    ///
   phiA1, PhiA1, phiA2, PhiA2, dA1A2, EXPmBBoverSigu, tymch1, g_sigu,       ///
   g_sigv, EXPmBoverSigu, tymch2, me, mA1, mA2, mphiA1, mPhiA1, mPhiA2,     ///
   dmA1A2, tymch3, Ti0i
    
  real matrix g_a, g_B
  //1
	//Ti_max    = max(Ti0)
  //2
  alpha     = theta[1 .. p]
  // 3
  sigu      = theta[p+1]
  //4
  sigv      = theta[p+2]
  //5
  B         = theta[(p+3) .. (p+Ti_max+2)]
  //B         = B :* (B :>= 0) + 1e-8 * (B :< 0)
  //6
  e = y - x * alpha'
  //7
  // BB is B[t] written into respecful time periods in each panel
  BB        = J(nt, 1, .)
  sample_all= 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i0 = Ti0 :== t
    sample_i  = sample_all :* sample_i0
    BB[selectindex(sample_i)] = J(sum(sample_i0), 1, B[t])
  }
  //9
  A1        = (BB+e)/sigv :+ sigv/sigu
  //10
  A2        = e/sigv :+ sigv/sigu
  //11
  phiA1     = normalden(A1)
  PhiA1     = normal(A1)
  phiA2     = normalden(A2)
  PhiA2     = normal(A2)
  dA1A2     = PhiA1 - PhiA2
  //12
  dA1A2 = dA1A2 :* (dA1A2 :> 0) + 1e-8 * (dA1A2 :<= 0)
  //13
  EXPmBBoverSigu = exp(-BB/sigu)
  //14
  ll        = -log(sigu) :- log(1:-EXPmBBoverSigu) :+ (1/2) * (sigv/sigu)^2 :+
   e/sigu + log(dA1A2)
  //15
	if (todo>=1) {
    //101
		tymch1    = (-phiA1+phiA2) :/ dA1A2
    //102
    g_a       = -x/sigu :+ tymch1 :* x/sigv
    //103
    g_sigu    = -1/sigu :+ BB/sigu^2 :* EXPmBBoverSigu :/ (1:-EXPmBBoverSigu) -
     1/sigu^2*e :- sigv^2/sigu^3 :+ tymch1*sigv/sigu^2
    //104
    g_sigv    = sigv/sigu^2 :+ (phiA1 :* (-(BB+e)/sigv^2 :+ 1/sigu) - 
     phiA2 :* (-e/sigv^2 :+ 1/sigu)) :/ dA1A2
     //105
    EXPmBoverSigu = exp(-B/sigu)
    tymch2    = EXPmBoverSigu :/ (1:-EXPmBoverSigu)
    //107
    g_B = J(nt, Ti_max, 0)
    for (i = 1; i <= n; i++){
      me     = panelsubmatrix(e,  i, ids)
      //1071
      mA1    = (panelsubmatrix(BB,  i, ids) + me) / sigv :+ sigv/sigu
      //1072
      mA2    = me/sigv :+ sigv/sigu
      //1073
      mphiA1 = normalden(mA1)
      //1074
      mPhiA1 = normal(mA1)
      //1075
      mPhiA2 = normal(mA2)
      //1076
      dmA1A2 = mPhiA1 - mPhiA2
      //1077
      dmA1A2 = dmA1A2 :* (dmA1A2 :> 0) + 1e-8 * (dmA1A2 :<= 0)
      //1078
      tymch3 = mphiA1/sigv :/ dmA1A2
      //1079
      Ti0i   = panelsubmatrix(Ti0, i, ids)
      g_B[ids[i,1]::ids[i,2], Ti0i ] =
       -1/sigu * J(ids[i,3], 1, 1/ids[i,3] * tymch2[selectindex(Ti0i)]) + 
       J(ids[i,3], 1, 1/ids[i,3]*tymch3')
       //10710
    }
    //108
    grad = g_a, g_sigu, g_sigv, g_B
    //109
	}
}
end

capture mata mata drop ipanelhalftgam()
mata:
void ipanelhalftgam(real scalar todo, real vector theta,                     ///
	real vector y, real matrix x,                                             ///
	real matrix ids, real colvector Ti0, real scalar Ti_max,                  ///
  real scalar n, real scalar nt, real scalar p,                             ///
	ll, grad, Hess)
{
  real scalar sigma, gamma, lamda, sigma_u, lamdag
	
  real vector alpha, B, e, esq, BB, sample_i0, sample_i, A1, A3, A4, PhiA1, ///
   PhiA3, PhiA4, phiA1, phiA3, phiA4, dA3A4, dA1A2, d12, d34, Dlam, Dsig,   ///
   tymch2, g_sig, g_lam, me, BBi, mA3, mA4, mphiA3, mPhiA3, mPhiA4,         ///
   dmPhiA3mPhiA4, md34, tymch3
    
  real matrix g_a, g_B
  //1
	//Ti_max    = max(Ti0)
  //2
  alpha     = theta[1 .. p]
  // 3
  sigma     = theta[p+1]
  //4
  gamma     = theta[p+2]
  //5
  B         = theta[(p+3) .. (p+Ti_max+2)]
  //B         = B :* (B :>= 0) + 1e-8 * (B :< 0)
  //6
  e         = y - x * alpha'
  esq       = e:^2
  
  lamda     = sqrt(gamma/(1-gamma))
  sigma_u   = sigma/sqrt(1+1/lamda^2)
  //7
  // BB is B[t] written into respecful time periods in each panel
  BB        = J(nt, 1, .)
  sample_all= 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i0 = Ti0 :== t
    sample_i  = sample_all :* sample_i0
    BB[selectindex(sample_i)] = J(sum(sample_i0), 1, B[t])
  }
  //9
  A1        = BB / sigma_u
  A3        = ( (BB+e)*lamda + BB/lamda ) / sigma
  A4        = e * lamda / sigma
  PhiA1     = normal(A1)
  PhiA3     = normal(A3)
  PhiA4     = normal(A4)

  phiA1     = normalden(A1)
  phiA3     = normalden(A3)
  phiA4     = normalden(A4)
  
  dA3A4     = PhiA3 - PhiA4
  dA3A4     = dA3A4 :* (dA3A4 :> 0) + 1e-8 * (dA3A4 :<= 0)
  dA1A2     = PhiA1 :- 1/2
  dA1A2     = dA1A2 :* (dA1A2 :> 0) + 1e-8 * (dA1A2 :<= 0)
  
  ll        =  -esq / (2*sigma^2) + log(dA3A4) - log(dA1A2) :- log(sigma) :-
   log(2*pi())/2

  //15
	if (todo>=1) {
    d12    = 1 :/ dA1A2
    d34    = 1 :/ dA3A4
    Dlam   = sigma / (sqrt(1+1/lamda^2)*lamda)^3
    Dsig   = 1 / sqrt(1+1/lamda^2)
    //101
		tymch1 = d12 :* phiA1 / sigma_u
    tymch2 = tymch1 :* BB / sigma_u
    //102
    g_a    =  (e/sigma - d34:*(phiA3 - phiA4)*lamda ) :* x / sigma
    //103      
    g_sig  = esq/sigma^3 - d34 :* (phiA3:*A3 - phiA4:*e*lamda/sigma)/sigma +
     tymch2*Dsig :- 1/sigma
    //104
    g_lam  = d34 :* (phiA3:*((BB+e)/sigma - BB/(lamda^2*sigma)) - 
     phiA4:*e/sigma) + tymch2*Dlam
    //107
    g_B = J(nt, Ti_max, 0)
    for (i = 1; i <= n; i++){
      me  = panelsubmatrix(e,  i, ids)
      BBi = panelsubmatrix(BB,  i, ids)
      mA3 = ( (BBi + me) *lamda + BBi/lamda ) / sigma
      mA4 = me*lamda / sigma
      mphiA3 = normalden(mA3)
      mPhiA3 = normal(mA3)
      //mphiA4 = normalden(mA4)
      mPhiA4 = normal(mA4)
      dmPhiA3mPhiA4 = mPhiA3 - mPhiA4
      dmPhiA3mPhiA4 = dmPhiA3mPhiA4 :* (dmPhiA3mPhiA4 :> 0) + 
       1e-8 * (dmPhiA3mPhiA4 :<= 0)
      md34 = 1 :/ dmPhiA3mPhiA4
      tymch3 = md34 :* mphiA3 * (lamda/sigma + 1/(lamda*sigma)) - 
       panelsubmatrix(tymch1,  i, ids)
      g_B[ids[i,1]::ids[i,2], panelsubmatrix(Ti0, i, ids) ] = 
       J(ids[i,3], 1, 1/ids[i,3]*tymch3')
    }
    lamdag = 1/(2*lamda) * (1/(1-gamma)^2)
    g_gam  = g_lam * lamdag
    //108
    grad   = g_a, g_sig, g_gam, g_B
    //colsum(grad)
    //109
	}
}
end

capture mata mata drop ipaneldoublytgam()
mata:
void ipaneldoublytgam(real scalar todo, real vector theta,                   ///
	real vector y, real matrix x,                                             ///
	real matrix ids, real colvector Ti0, real scalar Ti_max,                  ///
  real scalar n, real scalar nt, real scalar p,                             ///
	ll, grad, Hess)
{
  real scalar sigma, gamma, lamda, mu, sigma_u, lamdag
   
  real vector alpha, B, e, BB, sample_all, sample_i0, sample_i, A1, A2, A3, ///
   A4, PhiA1, PhiA2, PhiA3, PhiA4, phiA1, phiA2, phiA3, phiA4,dA3A4, dA1A2, ///
   emsq, d12, d34, Dlam, Dsig, tymch2, tymch3, g_sig, g_lam, g_mu, me, BBi, ///
   mA3, mA4, mphiA3, mPhiA3, mPhiA4, dmPhiA3mPhiA4, md34, tymch4, g_gam
	
  real matrix g_a, g_B
  //2000
  //Ti0'
  //1
	//Ti_max    = max(Ti0)
  //2
  alpha     = theta[1 .. p]

  // 3
  sigma     = theta[p+1]
  //4
  gamma     = theta[p+2]
  mu        = theta[p+3]
  //5
  B         = theta[(p+4) .. (p+Ti_max+3)]
  //B         = B :* (B :>= 0) + 1e-8 * (B :< 0)
  //6
  e         = y - x * alpha'
  
  lamda     = sqrt(gamma/(1-gamma))
  sigma_u   = sigma/sqrt(1+1/lamda^2)
  //7
  // BB is B[t] written into respecful time periods in each panel
  BB        = J(nt, 1, .)
  sample_all= 1::nt
  for (t = 1; t <= Ti_max; t++) {
    sample_i0 = Ti0 :== t
    sample_i  = sample_all :* sample_i0
    BB[selectindex(sample_i)] = J(sum(sample_i0), 1, B[t])
  }
  //9
  A1        = (BB :- mu) / sigma_u
  A2        = -mu / sigma_u
  A3        = ( (BB + e)*lamda +  (BB :- mu)/lamda ) / sigma
  A4        = (e*lamda :- mu/lamda) / sigma
  PhiA1     = normal(A1)
  PhiA2     = normal(A2)
  PhiA3     = normal(A3)
  PhiA4     = normal(A4)

  phiA1     = normalden(A1)
  phiA2     = normalden(A2)
  phiA3     = normalden(A3)
  phiA4     = normalden(A4)
  
  dA3A4     = PhiA3 - PhiA4
  dA3A4     = dA3A4 :* (dA3A4 :> 0) + 1e-8 * (dA3A4 :<= 0)
  dA1A2     = PhiA1 :- PhiA2
  dA1A2     = dA1A2 :* (dA1A2 :> 0) + 1e-8 * (dA1A2 :<= 0)
  
  emsq      = (e :+ mu):^2
  ll        = -emsq / (2*sigma^2) + log(dA3A4) - log(dA1A2) :- 
   log(sigma) :- log(2*pi())/2
  
  //15
	if (todo>=1) {
    d12    = 1 :/ dA1A2
    d34    = 1 :/ dA3A4
    Dlam   = sigma / (sqrt(1+1/lamda^2)*lamda)^3
    Dsig   = 1 / sqrt(1+1/lamda^2)
    //101
		tymch1 = d12 :* phiA1/sigma_u
    tymch2 = tymch1 :* (BB:-mu) / sigma_u
    tymch3 = d12 :* phiA2*mu/sigma_u^2
    tymch4 = tymch2 + tymch3
    //102
    g_a    =  ((e :+ mu)/sigma - d34 :* (phiA3 - phiA4)*lamda ) :* x / sigma
    //103      
    g_sig  = emsq / sigma^3 - d34 :* (phiA3:*A3/sigma - phiA4:*A4/sigma) +
     tymch4 * Dsig :- 1/sigma
    //104
    g_lam  = d34 :* ( phiA3:*((BB:+e)/sigma - (BB:-mu)/(lamda^2*sigma)) - 
     phiA4:*(e/sigma :+ mu/(sigma*lamda^2)) ) + tymch4* Dlam 
    g_mu  = -((e:+mu)/sigma^2 + d34:*(phiA3 - phiA4)/(sigma*lamda)) + 
     d12:*(phiA1 :- phiA2)/sigma_u
    //107
    g_B = J(nt, Ti_max, 0)
    for (i = 1; i <= n; i++){
      me  = panelsubmatrix(e,  i, ids)
      BBi = panelsubmatrix(BB,  i, ids)
      //1071
      mA3 = ( (BBi + me) *lamda + (BBi :- mu)/lamda ) / sigma
      //1072
      mA4 = (me*lamda :- mu/lamda) / sigma
      //1073
      mphiA3 = normalden(mA3)
      //1074
      mPhiA3 = normal(mA3)
      //1075
      mPhiA4 = normal(mA4)
      //1076
      dmPhiA3mPhiA4 = mPhiA3 - mPhiA4
      //1077
      dmPhiA3mPhiA4 = dmPhiA3mPhiA4 :* (dmPhiA3mPhiA4 :> 0) + 
       1e-8 * (dmPhiA3mPhiA4 :<= 0)
      //1078
      md34 = 1 :/ dmPhiA3mPhiA4
      tymch4 = md34 :* mphiA3 * (lamda/sigma + 1/(lamda*sigma)) - 
       panelsubmatrix(tymch1,  i, ids)
      //1079
      g_B[ids[i,1]::ids[i,2], panelsubmatrix(Ti0, i, ids) ] = 
       J(ids[i,3], 1, 1/ids[i,3]*tymch4')
       ///10710
    }
    lamdag = 1/(2*lamda) * (1/(1-gamma)^2)
    g_gam  = g_lam * lamdag
    //108
    grad   = g_a, g_sig, g_gam, g_mu, g_B
    //colsum(grad)
    //109
	}
}
end

capture mata mata drop offinf()
mata:

real colvector offinf(real vector q, real matrix ids,                        ///
 real scalar n, real scalar nt)
{
  real colvector tymch
  
  tymch     = J(nt, 1, .)
  for (i = 1; i <= n; i++){
    qi     = panelsubmatrix(q,  i, ids)
    medqi  = median( qi )
    qi     = qi :* (qi :!= .) + medqi * (qi :== .)
    tymch[|ids[i,1] \ ids[i,2]|] = qi
  }
  return(tymch)
}
end

capture mata mata drop median()
mata:

real scalar median(real vector q)
{
  real vector qs, qc
  real scalar n
  
  n  = length(q)
  qc = q
  if( cols(q) != 1 ) qc = q'
  qs = sort(qc, 1)
  // if number of elements are even
  if ( mod(n,2) == 0 ) {
    tymch = (qs[n/2] + qs[n/2+1]) / 2
  }
  // if number of elements are odd
  else {
    tymch = qs[n/2+1]
  }
  return(tymch)
}
end

// end functions for bie

