*! version 1.0.0  17Oct2019
*! version 1.1.0  27Feb2020
*! version 1.2.0  03Mar2020
*! version 1.2.1  05Mar2020
*! version 1.3.0  07Mar2020
*! version 1.3.1  10Mar2020
*! version 1.3.2  7Apr2020
*! version 1.3.3  9Apr2020
*! version 1.3.4  16Apr2020
*! version 1.3.5  5Aug2020
*! version 1.4.0  31Aug2020

/*
Origin:

MATLAB code

% PSS1 Estimator
% Written by Park, Sickles, and Simar
% Updated by Wonho Song, May 2014
% E-mail: whsong@cau.ac.kr, whsong73@hotmail.com

Implemented by Oleg Badunenko
oleg.badunenko@brunel.ac.uk, obadunenko@gmail.com

*/


// if(c(MP)){
// 	set processors 1
// }

capture program drop xtsf3gpss1
program define xtsf3gpss1, eclass
version 11

  if !replay() {
    
	syntax varlist(numeric fv min=2) [if] [in]                  ///
		[, NOLOG GR0(real 0.1) GR1(real 2) GRI(real 0.1)      ///
		REPS(integer 399) noDOTS BW(string) COST ///
    LEVEL(string) NOLOG ] ///
    [noCI] [noPValues] [noOMITted] [noEMPTYcells] ///
    [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] ///
    [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] ///
    [SFORMAT(passthru)] [nolstretch]
	
  marksample touse
	//display 0

	// handle the lists
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'
	
	_rmcoll `indepvars' if `touse', expand `constant'
	local indepvars `r(varlist)'
	
	tempname b V sd R2 R2adj PSS1 EPSS1 EPSS1_p              ///
	         mytrace sgrid alpha1 alpha1_p rez bandwidth1 aic bic shat1 RSS1 ///
            mypanelvar Cp1
	
	// display 1
	if "`dots'" == "nodots" {
		local mytrace 0
	}
	else {
		local mytrace 1
	}
	// sort data
	
	// data balanced?
	
	//display 2
	display
  display as result "Description of the panel data:" as input "{hline 48}
	xtdescribe if `touse'
	local n  = r(N)
  local nt = r(sum)
  quietly xtset
  local mypanelvar `r(panelvar)'

  // handle production/cost function
	if "`cost'" == "" { 
		local myprod = -1 
		local function "Production" 
	}
	else {
		local myprod = 1
		local function "Cost"
	}
	
  local bandwidth `bw'
	if "`bandwidth'" == ""{
		local bandwidth -999
	}
   
  // handle level
  if ("`level'" == "") {
    local mylevel 95
  }
  else 	if `level' < 10 |  `level' > 99.99 {
		display "{p 1 1 7}{error}level() must be between 10 and 99.99 inclusive{p_end}"
		exit 198
	}
  else {
    local mylevel `level'
  }
   
   //xtset
	
	mata: pss1_work("`depvar'", "`indepvars'", "`touse'", "`mypanelvar'", ///
      "`cost'", "`b'", "`V'",  "`n1'", "`nt1'",              ///
			"`sd'", "`R2'", "`R2adj'", "`PSS1'", "`EPSS1'", "`EPSS1_p'",             ///
			"`reps'", "`gr0'", "`gr1'", "`gri'",  ///
			"`bandwidth'", "`bandwidth1'", "`sgrid'", "`alpha1'", "`alpha1_p'", ///
			"`rez'", "`aic'", "`bic'", "`shat1'", "`RSS1'", "`Cp1'", "`mytrace'")

	matrix colnames `b'        = `indepvars'
	matrix rownames `V'        = `indepvars'
	matrix colnames `V'        = `indepvars'
	if "`bandwidth'" == "-999"{
		matrix colnames `sgrid'= "bandwidth" "MSE"
	}
	//di 9
	matrix colnames `alpha1'   = "alpha"
	matrix colnames `alpha1_p' = "alpha"
	matrix colnames `rez'      = "residuals"
	matrix colnames `EPSS1'    = "Efficiency"
	matrix colnames `EPSS1_p'  = "Efficiency"
	ereturn post `b' `V', esample(`touse') buildfvinfo depname("`depvar'")
	ereturn scalar N           = `n'
	ereturn scalar sumTi       = `nt'
	ereturn scalar bandwidth   = `bandwidth1'
	ereturn scalar r2          = `R2'
	ereturn scalar r2_a        = `R2adj'
	ereturn scalar aic         = `aic'
	ereturn scalar bic         = `bic'
  ereturn scalar cp          = `Cp1'
  ereturn scalar shat        = `shat1'
  ereturn scalar RSS         = `RSS1'
	ereturn matrix eff_p       = `EPSS1_p'
	ereturn matrix eff         = `EPSS1'
	if "`bandwidth'" == "-999"{
		ereturn matrix bw_grid = `sgrid'
	}
	ereturn matrix alpha_p     = `alpha1_p'
	ereturn matrix alpha       = `alpha1'
	ereturn matrix residuals   = `rez'
	ereturn local predict "xtsf3gpss1_p"
	ereturn local cmd   "xtsf3gpss1"
  ereturn local cmdline "`0'"
  
  }
	if replay() {
//     display "replay here"
    syntax, [LEVel(real `c(level)')] [noCI] [noPValues] [noOMITted] [noEMPTYcells] [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] [SFORMAT(passthru)] [nolstretch]
  }
  if "`nolog'" == "" {
    if "`cformat'" == ""{
      local cformat = "cformat(%9.4f)"
    }
    display
    display as result "Sample:" as input "{hline 22}
    display as input " Number of obs    " as text "= " as result `nt'
    display as input " Number of groups " as text "= " as result `n'
    display as result "Diagnostics:" as input "{hline 17}
    display as input " R-squared        " as text "= " as result  %5.4f `R2'
    display as input " Adj R-squared    " as text "= " as result  %5.4f `R2adj'
    display as input " AIC              " as text "= " as result  %5.4f `aic'
    display as input " BIC              " as text "= " as result  %5.4f `bic'
    display as input " Root MSE         " as text "= " as result  %5.4f `shat1'
    display as input "{hline 29}"
    display
    //di as text "{hline 79}"
    display as input "PSS Type I estimator: some regressors are correlated with effects"
    display as text "Park, Sickles, and Simar (1998), Journal of Econometrics, 84(2):273–301"
    display
    display as input " `function'" as text " Stochastic Frontier"
    ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
  }
	
end

mata mata clear

//capture mata mata drop pss1_work()
mata:
// mata drop pss1_work()
void pss1_work( string scalar depvar,                        ///
	string scalar indepvars, string scalar touse,  string scalar mypanelvar,  ///
  string scalar mycost, ///
	string scalar bs,   string scalar vs,                    ///
	string scalar n1,   string scalar nt1,               ///
	string scalar sd,   string scalar R2, ///
	string scalar R2adj,string scalar PSS1, ///
	string scalar Eff, string scalar Eff_p, //
	string scalar NB,   string scalar gr0,                   ///
	string scalar gr1,  string scalar gri,                   ///
	string scalar bw0,  string scalar bw1,                   ///
	string scalar sg1,  string scalar a1,                    ///
	string scalar a1_p, string scalar rez1,                  ///
	string scalar aic,  string scalar bic,                  ///
  string scalar shatname,  string scalar RSSname, ///
	string scalar Cpname, string scalar mytrace)
{
	
	//200

	real vector y, bh, se1
	real matrix x, vh, v1, x_sd_work, xp
	real scalar tr, p1, min_i, w, n, t, nt, p
	
	bw   = strtoreal(bw0)

	y    = st_data(., depvar, touse)
	//rows(y)
	x    = st_data(., indepvars, touse)
	//rows(x)
	p    = cols(x)
	
	// panel variable
	//stata("xtset")
   //stata("return list")
	//stata("local mypanelvar `r(panelvar)'")
	//mypanelvar = st_local("mypanelvar1")
   //"mypanelvar"
   //mypanelvar
   //mypanelvar = st_local("`r(panelvar)'")
   //"mypanelvar"
   //mypanelvar
	ids0  = st_data(., mypanelvar, touse)
	ids   = panelsetup(ids0, 1)
	ids   = ids, ids[,2] - ids[,1] :+ 1
	//ids
	//st_view(myid = ., ., mypanelvar)
	//mypanelvar[1..10,]
	
	//32
	//t
	//n    = strtoreal(n)
	n    = rows(ids)
	//33
	//n
	//nt   = strtoreal(nt)
	nt   = rows(y)
	//t    = strtoreal(t)
	t    = round(nt/n)
	//34
	//nt
	gr0  = strtoreal(gr0)
	//35
	//gr0
	gr1  = strtoreal(gr1)
	//36
	//gr1
	gri  = strtoreal(gri)
	//37
	//gri
	NB   = strtoreal(NB)
	//38
	//NB
	tr   = strtoreal(mytrace)
	//p1   = 4
	//bw1  = 0.1
  if (mycost == ""){
    cost = 0
  }
  else {
    cost = 1
  }
	
	// positions of time-constant variables
	
	// get the variance in each column
	//391
	x_sd_work = J(n,p,.)
	for (i=1; i<=rows(ids); i++) {
		xp = panelsubmatrix(x,i,ids)
		for(j=1;j<=p;j++){
			x_sd_work[i,j] = variance(xp[,j])
		}
	}
	//392
	//x_sd_work
	x_sd_mean = mean(x_sd_work)
	//393
	// is this x_j time-constant or time-varying?
	x_time_constant = x_sd_mean :== 0
	x_time_varying = x_sd_mean :!= 0
	
	//select(x_sd_work, x_time_varying)
	
	// number of time-constant regressors
	p1 = sum(x_time_varying)
	//p1
	//394
	// do this only if p1 < p
	
	if(p1 < p){
		tmp3 = ""
		invtokens(tmp3)
		tmp3 = "Note:", strofreal(p-p1) 
		tmp3 = tmp3, "regressors are time-constant"
    printf ("{text}Note: {result}%f {input}regressor(s) is(are) time-constant\n", p-p1)
		//tmp3 = invtokens(tmp3)
		//invtokens(tmp3)
		
		// get the indicies on which time-constant or time-varying
		x_time_constant_which = J(1,p-p1,.)
		x_time_varying_which = J(1,p1, .)
		jj = 1
		jjj = 1
		for(j=1;j<=p;j++){
			if(x_time_constant[,j] == 1){
				x_time_constant_which[,jj] = j
				jj++
			} else {
				x_time_varying_which[,jjj] = j
				jjj++
			}
		}
		//x_time_constant_which
		//x_time_varying_which
		// make indices with the order
		x_indicies = x_time_varying_which, x_time_constant_which
		x_indicies = 1..p\ x_indicies
		//x_indicies
		x = x[, x_indicies[2,]]
	}
	
	//395
	
	// Finding bandwidth for PSS estimators using bootstrap
	// if it is not given
	if (bw == -999){
		//20
		// this is the longest
		// the result is matrix with 2 columns.
		// col 1: bandwidth
		// col 2: MSE
		sgrid = bootbpss1(y, x, ids, n, t, nt, p, p1, NB, gr0, gri, gr1, tr)
		//21
		// this is the optimal bandwidth which gives the smallest MSE
		min_i = J(0,0,.); w  = J(0,0,.)
		//111
		minindex(sgrid[.,2], 1, min_i, w)
		""
		""
		//tmp3 = "Optimal bandwidth is", strofreal(sgrid[min_i,1])
		//tmp3 = invtokens(tmp3)
		// this is just dispaying in MATA
		//invtokens(tmp3)
    printf ("{text}Optimal bandwidth (for the within panel data estimator) is {result}%5.4f\n", sgrid[min_i,1])
		//112
		//st_numscalar(bw, sgrid[min_i,1])
		s1   = sgrid[min_i,1]
		//sgrid
	} else {
		""
    printf ("{text}Bandwidth bw ({input}%f{text}) is provided\n", bw)
    printf ("{error}Make sure it was chosen to minimize MSE for the PSS estimator\n")
		//tmp3 = "Bandwidth", strofreal(bw) 
		//tmp3 = tmp3, "is provided.  Make sure it is found using optimality criterion"
		//tmp3 = invtokens(tmp3)
		// this is just dispaying in MATA
		//invtokens(tmp3)
		s1   = bw
	}
	
	//s1
	// this gives beta_hat and vcov in the last two arguments, bh, vh
	pss1onlyB(y, x, ids, n, t, nt, p, p1, s1, bh, vh)
		
	//22
	//bh
	//23
	//vh
	//vh = J(cols(x), cols(x), 1)
	//24
	//se1
	//pss1post(y, x, bh, n, t, nt, PSS1, EPSS1, R2p1, alphapss1, residual)
	pss1post(cost, y, x, ids, n, p, t, nt, bh, PSS1, EPSS1, EPSS1_p,        ///
		R2p1, alphapss1, alphapss1_p, residual, residual0, shat, RSS, Cp,     ///
      aic_pss1, bic_pss1)
	bh = bh'
	//variance(residual)
	// AIC BIC
	//R2p1
	//rows(EPSS1)
	//81
	//PSS1
	//82
	//EPSS1
	//83
	//R2p1
	//84
	//alphapss1
	//85
	//residual[1:10]
	if(p1 < p){		
		// recover b from b_t and vcov from vcov_t

		b_t_t = bh
		b_t_t[x_indicies[2,]] = bh[x_indicies[1,]]
		bh = b_t_t
		
		vcov_t_t = vh
		vcov_t_t[x_indicies[2,],x_indicies[2,]] = vh[x_indicies[1,],x_indicies[1,]]
		vh = vcov_t_t
		// done recovering b from b_t and vcov from vcov_t
	}
	st_matrix(bs, bh)
	st_matrix(vs, vh)
	//86
	//s1
	st_numscalar(bw1, s1);
	//1
	//n
	//st_numscalar(ns, n)
	//2
	//nt
	//st_numscalar(nt1, nt)
	//3
	st_numscalar(R2, R2p1[1])
	//4
	st_numscalar(R2adj, R2p1[2])
	//5
	st_numscalar(aic, aic_pss1)
	//6
	st_numscalar(bic, bic_pss1)
	//87
	if (bw == -999){
		st_matrix(sg1, sgrid)
	}
	//88
	st_matrix(a1, alphapss1)
	//89
	st_matrix(a1_p, alphapss1_p)
	//891
	st_matrix(rez1, residual0)
	//892
	st_matrix(Eff, EPSS1)
	st_matrix(Eff_p, EPSS1_p)
   //236
	st_numscalar(shatname,   shat)
   //237
  st_numscalar(RSSname,    RSS)
  st_numscalar(Cpname,    Cp)
}

//end

//capture mata mata drop pss1onlyB()
//mata
// mata drop pss1onlyB()
void pss1onlyB(real vector y,  real matrix x,        ///
	real matrix ids, real scalar n,        ///
	real scalar t,  real scalar nt,
	real scalar p, real scalar p1,         ///
	real scalar bw, b_pss1, v_pss1)
{
	real scalar v, stilde2, sbarm, sb2, theta,            ///
	            xp, fI0, I0
	
	real vector ybar, ybar_p, utilde, sbar, ssbar, ystar, bht,       ///
	            nxbbar, zbbar, kernel, ker, fpdf,            ///
					cat, kerp, kerprime, wpdw, correc1,          ///
					correc2, correc3, y_p_work
					
	real matrix xbar, xbar_p, xxbar_p, xtilde, ytilde, ttilde,                ///
	            z, zbar, nx, nxbar, xqv,                     ///
					tmp1, tmp2, x2hat, xx, xxbar, xstar,         ///
					zstar, phatz, nxbtil, nzbtil,                ///
					SBXZ, SBX, SWX, SWZ, SWXZ,                   ///
					Ihat, Ihatinv, x_p_work
					
	// not used
	// real matrix ttildecv, nxtilde, ztilde, tglscov, SWINV, SW
					
	//p         = cols(x)
	//nt        = rows(y)
	//n         = rows(ids)
	//t         = round(nt/n)
	//nt        = rows(x)
	//t         = nt/n
	v         = p - p1
	//41
	//ybar      = mean(reshape(y,n))'
	//xbar      = reshape(mean(reshape(x,n*p)),p)	
	ybar      = J(n, 1, .)
	xbar      = J(n, p, .)	
	//for (i=1; i<=n; i++) {
	//	yp = panelsubmatrix(y,i,ids)
	//	xp = panelsubmatrix(x,i,ids)
	//	length_p = ids[i,2] - ids[i,1] + 1
	//	ybar[i,.] = mean(yp)
	//	xbar[i,.] = mean(xp)
	//}
	//42
	//xtilde    = x - xbar # J(t,1,1)
	//ytilde    = y - ybar # J(t,1,1)
	xbar_p    = J(nt, p, .)
	ybar_p    = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		//1234
		x_p_work = panelsubmatrix(x, i, ids)
		y_p_work = panelsubmatrix(y, i, ids)
		//length_p = ids[i,2] - ids[i,1] + 1
		ybar[i,.] = mean(y_p_work)
		xbar[i,.] = mean(x_p_work)
		ybar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, ybar[i,.])
		xbar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, xbar[i,.])
	}

	ytilde    = y - ybar_p
	xtilde    = x - xbar_p
	//43
	ttilde    = invsym(quadcross(xtilde, xtilde))*quadcross(xtilde, ytilde)
	//ttilde
	//44
	utilde    = ytilde - xtilde*ttilde

	stilde2   = quadcross(utilde, utilde) / (n*(t-1))
	//ttildecv  = stilde2*invsym(xtilde'*xtilde)
	sbar      = ybar - xbar*ttilde;

	sbarm     = mean(sbar)
	ssbar     = sbar :- sbarm
	sb2       = ssbar'*ssbar/n
	theta     = (stilde2/(t*sb2))^(.5)
	//431
	if(p1 == p){
		z      = J(nt,0,.)
		zbar   = J(n,0,.)
		zbar_p = J(nt,0,.)
	}
	else {
		z      = x[.,p1+1..p]
		zbar   = xbar[.,p1+1..p]
		zbar_p = J(nt,v,.)
		for (i=1; i<=n; i++) {
			//	z_p_work = panelsubmatrix(zbar,i,ids)
			//length_p = ids[i,2] - ids[i,1] + 1
			zbar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, zbar[i,.])
		}
	}
	//432
	nx        = x[.,1..p1]

	nxbar     = xbar[.,1..p1]
	nxbar_p   = xbar_p[.,1..p1]
	//45
	//nxbar

	//nxtilde   = xtilde[.,1..p1]

	// ztilde    = z - zbar # J(t,1,1)

	// Hausman - Taylor

	//xqv       = nxbar # J(t,1,1), xtilde
	xqv       = nxbar_p, xtilde
	//112
	//mean(xtilde)
	// tmp1 is where there is discrepancy with the MATLAB code
	tmp1      = invsym(quadcross(xqv, xqv))
	//113
	//mean(tmp1)
	tmp2      = quadcross(xqv, z)
	//114
	//mean(tmp2)
	x2hat     = xqv*tmp1*tmp2
	//115
	//mean(x2hat)

	xx        = nx, x2hat, J(nt,1,1)
	//116
	//mean(xx)

	//xp        = cols(xx)
	//46
	//xp

	//xxbar     = reshape(mean(reshape(xx,n*xp)),xp)
	
	xxbar_p = J(nt,cols(xx),.)
	for (i=1; i<=n; i++) {
		x_p_work = panelsubmatrix(xx,i,ids)
		//length_p = ids[i,2] - ids[i,1] + 1
		xxbar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, mean(x_p_work))
	}
	
	//xstar     = xx - xxbar # J(t,1,1) + (theta*xxbar) # J(t,1,1)
	//ystar     = y - ybar # J(t,1,1) + (theta*ybar) # J(t,1,1)
	xstar     = xx - (1-theta)*xxbar_p		
	ystar     = y  - (1-theta)*ybar_p
	//115
	//xstar[1..10,]
	//116
	//mean(ystar)

	bht       = invsym(quadcross(xstar, xstar)) * quadcross(xstar, ystar)
	// tglscov   = stilde2*invsym(xstar'*xstar)
	//117
	//invsym(xstar'*xstar)
	// Get Hausman Parameter estimates
	
	bht       = bht[1..p]

	nxbbar    = mean(nxbar)
	zbbar     = mean(zbar)

	zstar     = J(n,1,1), zbar

	phatz     = zstar * invsym(quadcross(zstar, zstar)) * quadcross(zstar, nxbar)


	//nxbtil    = nx - nxbar_p - (nxbar_p - nxbbar # J(nt,1,1))
	nxbtil    = nx - nxbar_p - (nxbar_p :- nxbbar)
	//nzbtil    = z - zbar # J(t,1,1) - (zbar # J(t,1,1) -  zbbar # J(nt,1,1))
	//nzbtil    = z - zbar_p - (zbar_p - zbbar # J(nt,1,1))
	nzbtil    = z - zbar_p - (zbar_p :- zbbar)

	SBXZ      = quadcross(nxbar-phatz, nxbar-phatz)/(n) 

	SBX       = SBXZ
	SWX       = quadcross(nxbtil, nxbtil) / (n) 
	SWZ       = quadcross(nzbtil, nzbtil) / (n) 
	SWXZ      = quadcross(nzbtil, nxbtil) / (n) 
	// SW        = (SWX\ SWXZ), (SWXZ'\ SWZ)
	// SWINV     = invsym(SW)

	kernel    = estpdfk(utilde,utilde,bw,1)
	ker       = estpdfk(utilde,utilde,bw,2)/bw
	fpdf      = ker:/ kernel

	fI0       = mean(fpdf:^2)
	cat       = sbar, zbar

	// Binwidth of the distribution of the effects and regressors 
	
	//70
	//cat
	//70
	//bw
	kernel    = estpdfk(cat,cat,bw,1)
	//71
	//kernel
	ker       = estpdfk(sbar,sbar,bw,2)/bw
	//72
	//ker
	kerp      = estpdfk(zbar,zbar,bw,1)
	//73
	//kerp
	kerprime  = ker:*kerp
	//74
	//kerprime


	wpdw      = kerprime :/ kernel
	I0        = mean(wpdw:^2)

	Ihat      = ((SWX*fI0 + I0*SBX), fI0*SWXZ')\ (fI0*SWXZ, J(v,v,0))
	
	if(p1 == p){
		//768
	}
	else {
		Ihat[p1+1..p, p1+1..p] = SWZ*fI0
	}

	Ihatinv   = invsym(Ihat)

	correc1   = (nxbtil)'*fpdf
	//61
	//correc1
	correc2   = correc1 + nxbbar'*sum(wpdw)
	//62	
	//correc2
	correc3   = (correc2 - nxbar'*wpdw) \ ((nzbtil)'*fpdf)
	//63
	//correc3
	//64
	b_pss1    = bht + Ihatinv*correc3/n
	//65
	//b_pss1
	//se_pss1   = sqrt(diagonal(Ihatinv/n))
	v_pss1    = Ihatinv/n
	//66
	
	//return
	
	// return(b_pss1)
}

//end

//mata
//mata drop pss1post()
void pss1post(real scalar mycost, real vector y, real matrix x,              ///
	real matrix ids, real scalar n, real scalar p,                             ///
	real scalar t,   real scalar nt,                                           ///
  real vector b_pss1,  PSS1, EPSS1, EPSS1_p, R2p1,                           ///
	alphapss1, alphapss1_p, residual, residual0, shat, RSS, Cp, aic, bic)
{
	real vector ep1, ey
	
	//91
	residual0 = y - x*b_pss1
   //residual0[1..10]
	//92
	residual = residual0 :- mean(residual0)
	//93

	ep1      = residual

	// Calculation of R2 
	ey       = y :- mean(y)
	//94
	R2p1     = 1 - (ep1'*ep1)/(ey'*ey)
	//95
	//R2p1
	//95
	R2p1     = R2p1, 1-(1-R2p1)*(nt-1)/(nt-p-n)
	//96
	//R2p1
	//PSS1     = mean(reshape(residual,n))'
	//97
	//alphapss1= PSS1
  //PSS1     = PSS1 :-mean(PSS1)
	
	
	alphapss1= J(n, 1, .)
	for (i=1; i<=n; i++) {
		alphapss1[i,.] = mean( panelsubmatrix(residual,i,ids) )
	}
	//97
	PSS1     = alphapss1 :- mean(alphapss1)
	//PSS1
	//alphapss1
	
	//98
  if (mycost == 1){
    EPSS1    = exp(min(PSS1) :- PSS1)
  }
  else {
    EPSS1    = exp(PSS1 :- max(PSS1))
  }
	
	//99
	
	
	EPSS1_p     = J(nt, 1, .)
	PSS1_p      = J(nt, 1, .)
	alphapss1_p = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		//length_p = ids[i,2] - ids[i,1] + 1
		PSS1_p[ids[i,1]::ids[i,2],.]      = J(ids[i,3], 1, PSS1[i,])
		EPSS1_p[ids[i,1]::ids[i,2],.]     = J(ids[i,3], 1, EPSS1[i,])
		alphapss1_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, alphapss1[i,])
	}
	//995
	//residual = y - x*b_pss1 - PSS1 # J(t,1,1)
	residual = y - x*b_pss1 - PSS1_p
   //997
  shat2    = variance(residual)
  shat     = sqrt(shat2)
   //998
  RSS      = cross(residual, residual)
  Cp       = RSS/shat2 - nt + 2*(p+1)
   //999
  aic      = log((nt-1)/nt*shat2)+1+2*(p+1)/nt
   //9991
	bic      = log((nt-1)/nt*shat2)+1+(p+1)*log(nt)/nt

	//100
	//sqrt(variance(residual))
	//return
}

//mata
//mata drop bpss1()
real scalar bpss1(real vector y,  real matrix x,        ///
	real matrix ids, real scalar n,        ///
	real scalar t,  real scalar nt,
	real scalar p, real scalar p1,          ///
	real scalar bw,  real scalar NB)
{
	real scalar crit, crit1, NB1, sighat
	real vector b_pss1, eps, yb
	
	//nt        = rows(x)
	//t         = nt/n
	
	// ESTIMATION
	
	//b_pss1    = pss1onlyB(y, x, p1, n, t, nt, bw)
	//51
	//pss1onlyB(y, x, p1, n, t, nt, bw, b_pss1, vh)
	pss1onlyB(y, x, ids, n, t, nt, p, p1, bw, b_pss1, vh)
	//52
	//pss1post(y, x, b_pss1, n, t, nt, PSS1, EPSS1, R2p1, alphapss1_p, rez)
	pss1post(1, y, x, ids, n, p, t, nt, b_pss1, PSS1, EPSS1, EPSS1_p, ///
		R2p1, alphapss1, alphapss1_p, rez, rez_o, shat, RSS, Cp, aic, bic)
     
	crit      = 0
	NB1       = 0

	sighat    = sqrt(variance(rez))
	//53
	//sighat
	//54
	//bw
	// BOOTSTRAP DATA GENERATION
	
	for(indb = 1; indb<=NB; indb++){
		eps    = rnormal(nt, 1, 0, sighat)
		//yb     = x*b_pss1 + alphapss1 # J(t,1,1) + eps
		yb     = x*b_pss1 + alphapss1_p + eps
		
		//55
		//mean(yb)
		// RERUN OF PROGRAM
		
		//b_pss1b = pss1onlyB(ystar, x, p1, n, t, nt, bw)
		//pss1onlyB(yb, x, p1, n, t, nt, bw, b_pss1b, vh)
		pss1onlyB(yb, x, ids, n, t, nt, p, p1, bw, b_pss1b, vh)
	
		if( missing(b_pss1b) ) NB1++
		//56
		//NB1
	
		// END OF BOOTSTRAP
		//114
		//b_pss1, b_pss1b, b_pss1b - b_pss1, (b_pss1b - b_pss1):^2
		//116
		//sum( (b_pss1b - b_pss1):^2 )
		crit   = crit + sum( (b_pss1b - b_pss1):^2 )	
	}
	
	// end of bootloop
	
	crit1     = crit/(NB-NB1)
	//117
	//crit1
	return(crit1)
	
}

//end

// y,x,n,NB,gr_pss,,p1,EST,datatype

//mata
//mata drop bootbpss1()
real matrix bootbpss1(real vector y,  real matrix x,        ///
	real matrix ids, real scalar n,        ///
	real scalar t,  real scalar nt,
	real scalar p, real scalar p1,          ///
	real scalar NB,  real scalar gr0,         ///
	real scalar gri, real scalar gr1, real scalar tr)
{
						//real scalar zp,  real scalar k1,          ///
                  //real scalar j1,  real scalar l1,          ///
	real scalar sl, min_i, w
	// real vector b_pss1, eps, ystar,
	
	// nt        = rows(x)
	// t         = nt/n
	// s1        = 0
	//101
	sl        = ceil( (gr1 - gr0) / gri + 1 )
	//102
	//sl
	sgrid     = J(sl, 2, .)
	//103

  ""
  printf("{text:Calculating optimal bandwidth for the PSS (some regressors are correlated with effects) estimator.}\n")
  printf ("{error:Please be patient!}\n")
  printf ("{input:It may take long time to compute estimates if the data size is large.}\n")
  ""
  printf ("{text:Choosing the bandwidth which has the smallest MSE for the PSS estimator.}\n")
  printf ("{result}Going over the grid, which contains {input}%f {result}grid points\n", sl)
  printf ("{result}(in each grid point, {input}%f {result}bootstrap replications are used)\n", NB)
	if (tr == 1) stata("_dots 0")
	for(j=1; j<=sl; j++) {
		if (tr == 1) {
			st_numscalar("q982769726w", j)
			stata("_dots q982769726w 0")
			//st_local("iter1", strofreal(j))
			//stata("_dots `=scalar(iter1)' 0")
		}
		sgrid[j,1] = gr0 + (j-1)*gri
		//117
		//j
		//sgrid[j,1]
		//sgrid[j,2] = bpss1(y,x,p1,n,t,nt,sgrid[j,1],NB)
		sgrid[j,2] = bpss1(y, x, ids, n, t, nt, p, p1, sgrid[j,1], NB)
	}
	//110
	//sgrid
	return(sgrid)
}

//mata drop reshape()
real matrix reshape(numeric matrix z, numeric scalar np)
{
	return(rowshape(vec(z),np)')
}

//mata drop kerneln()
real matrix kerneln(numeric matrix z)
{
	real scalar k
	real vector kk
	k = cols(z)
	if (k == 1) {
		kk = (1/sqrt((2*pi())^k))*exp(-.5*((z:*z)))
	}
	else {
		kk = (1/sqrt((2*pi())^k))*exp(-.5*colsum((z:*z)')')
	}
	return(kk)
}

//mata drop kernpri()
real matrix kernpri(numeric matrix z)
{
	real vector kk
	kk = -z:*kerneln(z)
	return(kk)
}

//mata drop colprod()
real matrix colprod(numeric matrix z)
{
	real vector kk
	kk = z[1,]
	for(q = 2; q<=rows(z); q++) {
		kk = kk :* z[q,]
	}
	return(kk)
}

//mata drop kerneln4()
real matrix kerneln4(numeric matrix z)
{
	real scalar k
	real vector kk
	k  = cols(z)
	kk = colprod(((1.5:-.5*(z:*z))*(1/sqrt((2*pi())^k)):*exp(-.5*(z:*z)))')
	return(kk)
}

//mata drop kern4pri()
real matrix kern4pri(numeric matrix z)
{
	real vector kk	
	kk = (1.5 :- .5*(z:*z)) :* (kernpri(z)) - (z:*kerneln(z))
	return(kk)
}

//mata drop kernelu()
real matrix kernelu(numeric matrix z)
{
	real vector kk
	kk=colprod((abs(z):<(.5))')
	return(kk)
}

//mata drop estpdfk()
real colvector estpdfk(numeric matrix x, numeric matrix x0, ///
numeric scalar h, numeric scalar c)
{
	real scalar n0
	real colvector pdf
	real matrix w	
	n0  = rows(x0)
	pdf = J(n0,1,0)
	for(i = 1; i<=n0; i++) {
		psi = ( x :- x0[i,.] # J(n0,1,1) )/h
		//psi = ( x :- x0[i,.] )/h
		if (c == 1){
			w = kerneln(psi);
		}
		else if (c == 2){
			w = kernpri(psi);
		}
		pdf[i] = mean(w)/h
	}
	return(pdf)
}

end

