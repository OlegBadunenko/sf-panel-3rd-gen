*! version 1.0.0  19Mar2020
*! version 1.1.0  20Mar2020
*! version 1.2.0  22Mar2020
*! version 1.3.0  24Mar2020
*! version 1.4.0  25Mar2020
*! version 1.4.1  27Mar2020
*! version 1.4.2  30Mar2020
*! version 1.5.0  7Ap2020
*! version 1.5.1  9Apr2020
*! version 1.5.2  16Apr2020
*! version 1.5.3  12Aug2020
*! version 1.6.0  31Aug2020

/*
Origin:

MATLAB code

% PSS2 Estimator
% Written by Park, Sickles, and Simar
% Updated by Wonho Song, May 2014
% E-mail: whsong@cau.ac.kr, whsong73@hotmail.com

Implemented by Oleg Badunenko
oleg.badunenko@brunel.ac.uk, obadunenko@gmail.com

*/


// if(c(MP)){
// 	set processors 1
// }


capture program drop xtsf3gpss2
program define xtsf3gpss2, eclass
version 11

  if !replay() {
    
	syntax varlist(numeric fv min=2) [if] [in]                  ///
		[, GR0(real 0.1) GR1(real 0.2) GRI(real 0.1)      ///
		REPS(integer 699) noDOTS BWWithin(string) BWGls(string) COST ///
    LEVEL(string) NOLOG ] ///
    [noCI] [noPValues] [noOMITted] [noEMPTYcells] ///
    [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] ///
    [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] ///
    [SFORMAT(passthru)] [nolstretch]
    
	marksample touse

	// handle the lists
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'
	
	_rmcoll `indepvars' if `touse', expand `constant'
	local indepvars `r(varlist)'
	
	tempname bW VW bG VG n1 nt1 R2W R2Wadj R2G R2Gadj                        ///
	         PSS2W EPSS2W EPSS2W_p PSS2G EPSS2G EPSS2G_p                     ///
			 alphaW alphaW_p alphaG alphaG_p                                 ///
			 bwW bwG sgrid aicW bicW aicG bicG rezW rezG mytrace             ///
			 b2 V12 V2 fittedW fittedG shatW RSSW CpW shatG RSSG CpG mypanelvar
			 
	if "`dots'" == "nodots" {
		local mytrace 0
	}
	else {
		local mytrace 1
	}
	
	if "`bwwithin'" == "" & "`bwgls'" == ""{
		//display 1031
		local bwwithin -999
		local bwgls -999
		local bw_provided 0
	} 
	else if "`bwwithin'" != "" & "`bwgls'" != ""{
		// all is great
		//display 1032
		local bw_provided 1
	}
	else {
		display as error "Either both bww and bwg are provided or none"
		exit
	}
	
	display
  display as result "Description of the panel data:" as input "{hline 48}
	xtdescribe if `touse'
  quietly xtset
  local mypanelvar `r(panelvar)'
  local mytimevar `r(timevar)'
  
  // handle production/cost function
	if "`cost'" == "" { 
		local myprod = -1 
		local function "Production" 
	}
	else {
		local myprod = 1
		local function "Cost"
	}
	
	mata: pss2_work("`depvar'", "`indepvars'", "`touse'",               ///
    "`mypanelvar'", "`mytimevar'", "`cost'",                           ///
		"`bW'", "`VW'", "`bG'", "`VG'",  "`n1'", "`nt1'",                   ///
		"`R2W'", "`R2Wadj'", "`R2G'", "`R2Gadj'",                          ///
		"`EPSS2W'", "`EPSS2G'",  "`EPSS2W_p'", "`EPSS2G_p'",               ///
		"`alphaW'", "`alphaW_p'","`alphaG'", "`alphaG_p'",                  ///
		"`reps'", "`gr0'", "`gr1'", "`gri'", "`sgrid'",                     ///
		"`bwwithin'", "`bwW'", "`bwgls'", "`bwG'",                           ///
		"`aicW'", "`bicW'","`aicG'", "`bicG'",                                ///
    "`shatW'", "`RSSW'", "`CpW'","`shatG'", "`RSSG'", "`CpG'",              ///
		"`rezW'", "`rezG'", "`fittedW'", "`fittedG'", "`mytrace'")    ///

	/*
		This is in case two estimators are displayed:
	*/
	
	local indepvarsW = ""
	local indepvarsG = ""
	foreach lname of local indepvars{
		local indepvarsW = "`indepvarsW' Within:`lname'"
		local indepvarsG = "`indepvarsG' GLS:`lname'"
	}
	//display "`indepvarsW'"
	//display "`indepvarsG'"
	
	matrix `b2' = `bW' , `bG'
	matrix colnames `b2'        = `indepvarsW' `indepvarsG'
	//matrix list `b2'
	matrix `V12' = J(rowsof(`VW'), rowsof(`VW'), 0)
	matrix `V2' = `VW', `V12' \ `V12', `VG'
	matrix rownames `V2'        = `indepvarsW' `indepvarsG'
	matrix colnames `V2'        = `indepvarsW' `indepvarsG'
	
	//display 106
	
	if "`bw_provided'" == "0"{
		matrix colnames `sgrid'= "bandwidth" "MSE_W" "MSE_G"
	}
	matrix colnames `alphaW'   = "alpha"
	matrix colnames `alphaW_p' = "alpha"
	matrix colnames `alphaG'   = "alpha"
	matrix colnames `alphaG_p' = "alpha"
	matrix colnames `rezW'     = "residuals"
	matrix colnames `rezG'     = "residuals"
	matrix colnames `EPSS2W'   = "Efficiency"
	matrix colnames `EPSS2W_p' = "Efficiency"
	matrix colnames `EPSS2G'   = "Efficiency"
	matrix colnames `EPSS2G_p' = "Efficiency"
	ereturn post `b2' `V2', esample(`touse') buildfvinfo depname("`depvar'")
	ereturn scalar N           = `n1'
	ereturn scalar sumTi       = `nt1'
	ereturn scalar bandwidthW  = `bwW'
	ereturn scalar bandwidthG  = `bwG'
	ereturn scalar r2W         = `R2W'
	ereturn scalar r2G         = `R2G'
	ereturn scalar r2W_a       = `R2Wadj'
	ereturn scalar r2G_a       = `R2Gadj'
	ereturn scalar aicW        = `aicW'
	ereturn scalar bicW        = `bicW'
	ereturn scalar aicG        = `aicG'
	ereturn scalar bicG        = `bicG'
  ereturn scalar cpW         = `CpW'
  ereturn scalar cpG         = `CpG'
  ereturn scalar shatW       = `shatW'
  ereturn scalar shatG       = `shatG'
  ereturn scalar RSSW        = `RSSW'
  ereturn scalar RSSG        = `RSSG'
	ereturn matrix effG_p      = `EPSS2G_p'
	ereturn matrix effW_p      = `EPSS2W_p'
	ereturn matrix effG        = `EPSS2G'
	ereturn matrix effW        = `EPSS2W'
	if "`bw_provided'" == "0"{
		ereturn matrix bw_grid = `sgrid'
	}
	ereturn matrix alphaG_p    = `alphaG_p'
	ereturn matrix alphaW_p    = `alphaW_p'
	ereturn matrix alphaG      = `alphaG'
	ereturn matrix alphaW      = `alphaW'
	ereturn matrix residualsG  = `rezG'
	ereturn matrix residualsW  = `rezW'
	ereturn matrix xbG         = `fittedG'
	ereturn matrix xbW         = `fittedW'
	ereturn local predict "xtsf3gpss2_p"
	ereturn local cmd   "xtsf3gpss2"
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
    display as input " Number of obs    " as text "= " as result `nt1'
    display as input " Number of groups " as text "= " as result `n1'
    display as result "Diagnostics:" as input "{hline 17}
    display as result "Within:" as input "{hline 22}
    display as input " R-squared        " as text "= " as result  %5.4f `R2W'
    display as input " Adj R-squared    " as text "= " as result  %5.4f `R2Wadj'
    display as input " AIC              " as text "= " as result  %5.4f `aicW'
    display as input " BIC              " as text "= " as result  %5.4f `bicW'
    display as input " Root MSE         " as text "= " as result  %5.4f `shatW'
    display as result "GLS:" as input "{hline 25}
    display as input " R-squared        " as text "= " as result  %5.4f `R2G'
    display as input " Adj R-squared    " as text "= " as result  %5.4f `R2Gadj'
    display as input " AIC              " as text "= " as result  %5.4f `aicG'
    display as input " BIC              " as text "= " as result  %5.4f `bicG'
    display as input " Root MSE         " as text "= " as result  %5.4f `shatG'  
    display as input "{hline 29}"
    display
    //di as text "{hline 79}"
    display as input "PSS Type 2 estimator: AR(1) error"
    display as input "Park, Sickles, and Simar (2003), Journal of Econometrics, 117(2):279–309"
    display
    display as input " `function'" as text " Stochastic Frontier"
    ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
  }

end



mata mata clear

mata:

void pss2_work( string scalar depvar,                                        ///
    string scalar indepvars, string scalar touse,                            ///
    string scalar mypanelvar, string scalar mytimevar,                       ///
    string scalar mycost,                                                      ///
    string scalar bsWname,   string scalar vsWname,                          ///
    string scalar bsGname,   string scalar vsGname,                          ///
    string scalar n1name,    string scalar nt1name ,                         ///
    string scalar R2Wname,   string scalar R2Wadjname,                       ///
    string scalar R2Gname,   string scalar R2Gadjname,                       ///
    string scalar EffW_name,string scalar EffG_name,                         ///
    string scalar EffW_pname,string scalar EffG_pname,                       ///
    string scalar aWname,    string scalar aW_pname,                         ///
    string scalar aGname,    string scalar aG_pname,                         ///
    string scalar NB,    string scalar gr00,                                 ///
    string scalar gr10,  string scalar gri0, string scalar sg1name,          ///
    string scalar bwW0,  string scalar bwWname,                              ///
    string scalar bwG0,  string scalar bwGname,                              ///
    string scalar aicWname,  string scalar bicWname,                         ///
    string scalar aicGname,  string scalar bicGname,                         ///
    string scalar shatWname, string scalar RSSWname,                         ///
    string scalar CpWname,                                                   ///
    string scalar shatGname, string scalar RSSGname,                         ///
    string scalar CpGname,                                                   ///
    string scalar rezWname,  string scalar rezGname,                         ///
    string scalar xbWname,   string scalar xbGname,                          ///
    string scalar mytrace)
{
	//200

	real vector y, Ti0, Ti, ids0
	real matrix x, ids, ids1a, ids1b, ids2a, ids2b
	real scalar bwW, bwG, n, nt, p
	
	//201
	
	bwW       = strtoreal(bwW0)
	bwG       = strtoreal(bwG0)

	y         = st_data(., depvar, touse)
	x         = st_data(., indepvars, touse)
	
	//202
	
	// panel variable
	//stata("quietly xtset")
	//stata("local mypanelvar `r(panelvar)'")
	//mypanelvar= st_local("mypanelvar")
	
	//203

	ids0      = st_data(., mypanelvar, touse)
	ids       = panelsetup(ids0, 1)
	ids       = ids, ids[,2] - ids[,1] :+ 1
	ids1a     = ids[,1] :+ 1, ids[,2], ids[,3] :- 1
	ids1b     = ids[,1], ids[,2] :- 1, ids[,3] :- 1
	ids2a     = ids[,1] :+ 2, ids[,2], ids[,3] :- 2
	ids2b     = ids[,1], ids[,2] :- 2, ids[,3] :- 2
	
	//204
	
	//stata("local mytimevar `r(timevar)'")
	//mytimevar = st_local("mytimevar")
	Ti0       = st_data(., mytimevar, touse)
	Ti        = Ti0 :- min(Ti0) :+ 1
	
	//205

	p         = cols(x)
	n         = rows(ids)
	nt        = rows(y)
	//t         = round(nt/n)
	gr0       = strtoreal(gr00)
	gr1       = strtoreal(gr10)
	gri       = strtoreal(gri0)
	NB        = strtoreal(NB)
	tr        = strtoreal(mytrace)
  if (mycost == ""){
    cost = 0
  }
  else {
    cost = 1
  }

	//206
	
	// Finding bandwidth for PSS estimators using bootstrap
	// if it is not given
	if (bwW == -999){
		//2061
		// this is the longest
		// the result is matrix with 3 columns.
		// col 1: bandwidth
		// col 2: MSE_W
		// col 3: MSE_G
		sgrid = bootbpss2( y, x, ids, ids1a, ids1b, ids2a, ids2b, Ti,        ///
                           n, nt, p, NB, gr0, gri, gr1, tr )
						   
		//2062
		// this is the optimal bandwidth which gives the smallest MSE
		s1W   = sort(sgrid, 2)[1,1]
		s1G   = sort(sgrid, 3)[1,1]
		
		//2063

		""
		""
		tmp3  = "Optimal bandwidth (for the  within  estimator) is", strofreal(s1W)
		tmp3  = invtokens(tmp3)
		// this is just dispaying in MATA
		invtokens(tmp3)
		tmp3  = "Optimal bandwidth (for the    GLS   estimator) is", strofreal(s1G)
		tmp3  = invtokens(tmp3)
		// this is just dispaying in MATA
		invtokens(tmp3)
		//sgrid
	} else {
		//2064
		""
		tmp3  = "Bandwidths bww (", strofreal(bwW), ") and bwg (", strofreal(bwG)
		tmp3  = tmp3, ") are provided."
		tmp3  = invtokens(tmp3)
		// this is just dispaying in MATA
		invtokens(tmp3)
// 		tmp3  = "Make sure they were found using optimality criterion"
// 		invtokens(tmp3)
    printf ("{error} Make sure they were chosen to minimize MSE for the PSS estimators\n")

		s1W   = bwW
		s1G   = bwG
	}
	//207
		
	// ESTIMATION
	pss2onlyB(y, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,               ///
              s1W, s1G, bhW, vhW, bhG, vhG )
			  
	//208
	// Getting rhat and shat2
	pss2post( cost, y, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,           ///
              bhW, bhG,                                                      ///
              EPSS2W, EPSS2W_p, R2p2W, EPSS2G, EPSS2G_p, R2p2G,              ///
              alphapss2W, alphapss2W_p, alphapss2G, alphapss2G_p,            ///
              aic_pss2W, bic_psswW, aic_pss2G, bic_psswG,                    ///
			  residualW, residualG, rhat, shat, shatW, RSSW, CpW, shatG, RSSG, CpG )
	//209
	xbW  = x*bhW
	xbG  = x*bhG
	bhW  = bhW'
	bhG  = bhG'
	//210
	st_matrix(bsWname,       bhW)
	//211
	st_matrix(vsWname,       vhW)
	//212
	st_matrix(bsGname,       bhG)
	//213
	st_matrix(vsGname,       vhG)
	//214
	st_numscalar(bwWname,    s1W);
	//215
	st_numscalar(bwGname,    s1G);
	//216
	st_numscalar(n1name,     n)
	//217
	st_numscalar(nt1name,    nt)
	//218
	st_numscalar(R2Wname,    R2p2W[1])
	//219
	st_numscalar(R2Wadjname, R2p2W[2])
	//220
	st_numscalar(R2Gname,    R2p2G[1])
	//221
	st_numscalar(R2Gadjname, R2p2G[2])
	//222
	st_numscalar(aicWname,   aic_pss2W)
	//223
	st_numscalar(bicWname,   bic_psswW)
	//224
	st_numscalar(aicGname,   aic_pss2G)
	//225
	st_numscalar(bicGname,   bic_psswG)
	//226
	if (bwW == -999){
		st_matrix(sg1name,   sgrid)
	}
	//227
	st_matrix(aWname,        alphapss2W)
	//228
	st_matrix(aW_pname,      alphapss2W_p)
	//229
	st_matrix(aGname,        alphapss2G)
	//230
	st_matrix(aG_pname,      alphapss2G_p)
	//231
	st_matrix(rezWname,      residualW)
	//232
	st_matrix(rezGname,      residualG)
	//233
	st_matrix(EffW_pname,    EPSS2W_p)
	//234
	st_matrix(EffG_pname,    EPSS2G_p)
	st_matrix(EffW_name,     EPSS2W)
	//234
	st_matrix(EffG_name,     EPSS2G)
	//235
	st_matrix(xbWname,       xbW)
	//236
	st_matrix(xbGname,       xbG)
   //236
	st_numscalar(shatWname,   shatW)
   //237
  st_numscalar(RSSWname,    RSSW)
  st_numscalar(CpWname,     CpW)
   //236
	st_numscalar(shatGname,   shatG)
   //237
  st_numscalar(RSSGname,    RSSG)
  st_numscalar(CpGname,     CpG)
}
end

mata 
//mata clear
//mata drop pss2onlyB()
void pss2onlyB( real vector y,     real matrix x,                     ///
    real matrix ids,   real matrix ids1a, real matrix ids1b,                 ///
    real matrix ids2a, real matrix ids2b,                                    ///
    real scalar n,     real scalar nt,    real scalar p,                     ///
    real scalar bwW,   real scalar bwG,                                      ///
    b_pss2W, v_pss2W,  b_pss2G, v_pss2G )
{
	real scalar rtilde, rtilde_g_1, trho, stilde2, I_f
	
	real vector ybar, ybar_p, y_p_work, ytilde, btilde, restilde,            ///
	            C_i0, C_i1, C_i2, ct0, ebig, wtilde, wtilde_p,               ///
                ztilde, tymch1, ystar1, btildegls, xtilbar, ztilde1,         ///
                term1, temp1, temp2, wpdw, term3, correc
	
	real matrix xbar, xbar_p, x_p_work, xtilde, xstar1, xtilcent, x1,        ///
	            zxrho, temp, K, Kprime, S1, S2, Ihat
				
	numeric scalar sq_1_rt2
	numeric vector ystar, ebig_sqrt
	numeric matrix xstar
	
	// the WITHIN estimator
	ybar      = J(n, 1, .)
	xbar      = J(n, p, .)	
	ybar_p    = J(nt, 1, .)
	xbar_p    = J(nt, p, .)
	for (i=1; i<=n; i++) {
		x_p_work = panelsubmatrix(x, i, ids)
		y_p_work = panelsubmatrix(y, i, ids)
		ybar[i,.] = mean(y_p_work)
		xbar[i,.] = mean(x_p_work)
		ybar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, ybar[i,.])
		xbar_p[ids[i,1]::ids[i,2],.] = J(ids[i,3], 1, xbar[i,.])
	}
	
	ytilde    = y - ybar_p
	xtilde    = x - xbar_p

	//btilde    = invsym(quadcross(xtilde, xtilde))*quadcross(xtilde, ytilde)
	btilde    = qrsolve(cross(xtilde,xtilde), cross(xtilde,ytilde))
	
	// the OLS residuals
	restilde  = y - x * btilde
	
	// the OLS WITHIN residuals
	//uw        = ytilde - xtilde*btilde
	
	// \tilde \rho = rtilde
	C_i0      = J(n, 1, .)
	C_i1      = J(n, 1, .)
	C_i2      = J(n, 1, .)
	for (i=1; i<=n; i++) {
		C_i0[i] = cross(panelsubmatrix(restilde, i, ids),
		                panelsubmatrix(restilde, i, ids)) / ids[i,3]
		C_i1[i] = cross(panelsubmatrix(restilde, i, ids1a),
		                panelsubmatrix(restilde, i, ids1b)) / ids1a[i,3]
		C_i2[i] = cross(panelsubmatrix(restilde, i, ids2a),
		                panelsubmatrix(restilde, i, ids2b)) / ids2a[i,3]
	}
	rtilde    = sum(C_i1-C_i2) / sum(C_i0-C_i1)
	// this could be outside [-1,1]l then need to deal with complex numbers
	rtilde_g_1= rtilde > 1 | rtilde < -1
	//rtilde_comp = iscomplex(...)
	
	// get the weights c_t and e_t
	// trho      = (1-rtilde^2) + (t-1) * (1-rtilde)^2
	ebig      = J(nt, 1, .)
	// calculate \tilde w_{i}, \tilde x_{i}, \tilde z_{it}and \tilde\sigma^2
	wtilde    = J(n, 1, .)
	xtilde    = J(n, p, .)
	ytilde    = J(n, 1, .)
	wtilde_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		
		trho           = (1-rtilde^2) + (ids[i,3]-1) * (1-rtilde)^2
		ct0            = J(ids[i,3], 1, ((1-rtilde)^2) / trho)
		ct0[1]         = (1-rtilde)/trho
		ct0[ids[i,3]]  = (1-rtilde)/trho
		ebig[ ids[i,1]::ids[i,2]] = ct0 * trho * (1+rtilde) / ((ids[i,3]-1) * (1-rtilde))		

		wtilde[i]      = cross(ct0, panelsubmatrix(restilde, i, ids))
		xtilde[i,]     = cross(ct0, panelsubmatrix(x, i, ids))
		ytilde[i]      = cross(ct0, panelsubmatrix(y, i, ids))
		wtilde_p[ ids[i,1]::ids[i,2]] = J(ids[i,3], 1, wtilde[i])
	}
	ztilde    = restilde - wtilde_p
	//eltype(ebig[i])
	// it can be complex
	ebig_sqrt = ebig
	if ( rtilde_g_1 ) ebig_sqrt = C(ebig)
	tymch1    = sqrt(ebig_sqrt) :* ztilde
	stilde2   = tymch1' * tymch1 / n
	//stilde2
	//iscomplex(stilde2)
	//Re(stilde2)
	//Im(stilde2)
	if( Im(stilde2) == 0 ) stilde2 = Re(stilde2)
	//iscomplex(stilde2)
	stilde2 = Re(stilde2)
	
	xstar1    = x[ids[,1],] - xtilde
	ystar1    = y[ids[,1]]  - ytilde
	if ( rtilde_g_1 ) {
		sq_1_rt2  = sqrt( C(1-rtilde^2) )
	} else {
		sq_1_rt2  = sqrt(1-rtilde^2)
	}
	xstar1    = xstar1 * sq_1_rt2
	ystar1    = ystar1 * sq_1_rt2
	
	xstar     = J(nt, p, .)
	ystar     = J(nt, 1, .)
	if ( rtilde_g_1 ) {
		xstar     = C( J(nt, p, .) )
		ystar     = C( J(nt, 1, .) )
	}
	for (i=1; i<=n; i++) {
		xstar[(ids[i,1]+1)::ids[i,2],] =
			panelsubmatrix(x, i, ids1a) - 
			panelsubmatrix(x, i, ids1b) * rtilde - 
			J(ids1a[i,3],1,xtilde[i,] * (1-rtilde))
		xstar[ids[i,1],] = xstar1[i,]
		ystar[(ids[i,1]+1)::ids[i,2]] =
			panelsubmatrix(y, i, ids1a) - 
			panelsubmatrix(y, i, ids1b) * rtilde - 
			J(ids1a[i,3],1,ytilde[i] * (1-rtilde))
		ystar[ids[i,1]] = ystar1[i]
	}
	// GLS within estimator
	btildegls = cholsolve(xstar'*xstar, xstar'*ystar)
	if ( rtilde_g_1 ) {
		btildegls = Re(btildegls)
	}
	
	// \tilde x_{\cdot} = xtilbar
	xtilbar   = mean(xtilde)
	xtilcent  = xtilde - J(n, 1, xtilbar)

	//calculate I_f, \Sigma_1, \Sigma_2, \hat I and  \hat \beta
	
	ztilde1   = ztilde[ids[,1]]
	x1        = x[ids[,1],]
	term1     = ((1-rtilde^2)/stilde2)*(ztilde1'*x1)

	zxrho     = J(n, p, .)
	for (i=1; i<=n; i++) {
		temp1 = panelsubmatrix(ztilde, i, ids1a) - 
		        panelsubmatrix(ztilde, i, ids1b) * rtilde
		temp2 = panelsubmatrix(x, i, ids1a) - 
		        panelsubmatrix(x, i, ids1b) * rtilde
		zxrho[i,] = cross(temp1, temp2)
	}
	term2     = colsum(zxrho) / stilde2
	// term2
	
	// density estimator
	// with current bandwith s chosen above by the s loop
		
	//temp      = J(1,n,1) # exp(wtilde/s)
	temp      = J(1,n,1) # exp(wtilde/bwW)
	temp      = temp :/ temp'
	// temp

	// temp_{ij}=exp((\tilde w_i -\tilde w_j)/s)

	//K         = (temp :/ ((1 :+ temp) :^ 2)) / s
	//Kprime    = (temp :* (1 :- temp) :/ ((1 :+ temp) :^3 )) / (s^2)
	K         = (temp :/ ((1 :+ temp) :^ 2)) / bwW
	Kprime    = (temp :* (1 :- temp) :/ ((1 :+ temp) :^3 )) / (bwW^2)

	// K is  a (n x n) matrix, so is Kprime (derivatives)

	// sum over j
	K         = colsum(K')
	// now K is a (1 x n) vector: \hat f(w_i),i=1,...,n

	Kprime    = colsum(Kprime')
	// now Kprime is a (1 x n) vector: \hat f^{{1)}(w_i),i=1,...,n

	wpdw      = Kprime :/ K
	// wpdw is a (1 x n) vector of {\hat f{{1)}\over \hat f}

	term3     = -(wpdw * xtilcent)

	I_f       = mean(wpdw':^2)

	// Sigma_1 et Sigma_2

	S1        = (xstar'*xstar)/n
	S2        = (xtilcent'*xtilcent)/n
	if ( rtilde_g_1 ){
		S1    = Re(S1)
		S2    = Re(S2)
	}

	// calculate \hat I and \hat \beta

	Ihat      = S1/stilde2 + S2*I_f
	correc    = invsym(Ihat) * (term1+term2+term3)'/n
	//bhat      = btilde + correc
	//bhat      = Re(bhat)
	b_pss2W   = btilde + correc
	//b_pss2W

	// Variance of beta gls
	//Vbwgls    = stilde2 * invsym(S1) / n
	// Variance of beta hat, efficient estimator
	//Vbhat     = invsym(Ihat) / n
	v_pss2W   = invsym(Ihat) / n
	
	// GLS =====================================================================
	
	// redo the same thing with the GLS within as first step
	// save OLS values of \tilde beta, \tilde rho and \tilde sigma^2
	
	//btilde1   = btilde
	//rtilde1   = rtilde
	//stilde21  = stilde2

	btilde    = btildegls
	// NOW the GLS estimator !!!!!!

	restilde  = y - x * btilde
	// NOW the GLS residuals !!!!!!!
	
	C_i0      = J(n, 1, .)
	C_i1      = J(n, 1, .)
	C_i2      = J(n, 1, .)
	for (i=1; i<=n; i++) {
		C_i0[i] = cross(panelsubmatrix(restilde, i, ids),
		                panelsubmatrix(restilde, i, ids)) / ids[i,3]
		C_i1[i] = cross(panelsubmatrix(restilde, i, ids1a),
		                panelsubmatrix(restilde, i, ids1b)) / ids1a[i,3]
		C_i2[i] = cross(panelsubmatrix(restilde, i, ids2a),
		                panelsubmatrix(restilde, i, ids2b)) / ids2a[i,3]
	}
	rtilde = sum(C_i1-C_i2) / sum(C_i0-C_i1)
	// this could be outside [-1,1]l then need to deal with complex numbers
	rtilde_g_1= rtilde > 1 | rtilde < -1
	//rtilde_comp = iscomplex(...)
	
	// get the weights c_t and e_t
	// trho      = (1-rtilde^2) + (t-1) * (1-rtilde)^2
	ebig      = J(nt, 1, .)
	// \tilde w_{i}, \tilde x_{i}, \tilde z_{it}and \tilde\sigma^2
	wtilde    = J(n, 1, .)
	xtilde    = J(n, p, .)
	ytilde    = J(n, 1, .)
	wtilde_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		
		trho           = (1-rtilde^2) + (ids[i,3]-1) * (1-rtilde)^2
		ct0            = J(ids[i,3], 1, ((1-rtilde)^2) / trho)
		ct0[1]         = (1-rtilde)/trho
		ct0[ids[i,3]]  = (1-rtilde)/trho
		et0            = ct0 * trho * (1+rtilde) / ((ids[i,3]-1) * (1-rtilde))
		ebig[ids[i,1]::ids[i,2]] = et0		

		wtilde[i]      = cross(ct0, panelsubmatrix(restilde, i, ids))
		xtilde[i,]     = cross(ct0, panelsubmatrix(x, i, ids))
		ytilde[i]      = cross(ct0, panelsubmatrix(y, i, ids))
		wtilde_p[ids[i,1]::ids[i,2]] = J(ids[i,3], 1, wtilde[i])
	}
	ztilde    = restilde - wtilde_p
	//eltype(ebig)
	// it can be complex
	ebig_sqrt = ebig
	if ( rtilde_g_1 ) ebig_sqrt = C(ebig)
	tymch1    = sqrt(ebig_sqrt) :* ztilde
	stilde2   = tymch1' * tymch1 / n
	//stilde2
	//iscomplex(stilde2)
	//Re(stilde2)
	//Im(stilde2)
	if ( rtilde_g_1 ) stilde2 = Re(stilde2)
	// iscomplex(stilde2)

	xstar1    = x[ids[,1],] - xtilde
	ystar1    = y[ids[,1]]  - ytilde
	sq_1_rt2  = sqrt(1 - rtilde^2)
	if ( rtilde_g_1 ) {
		sq_1_rt2  = sqrt( C(1-rtilde^2) )
	} else {
		sq_1_rt2  = sqrt(1-rtilde^2)
	}
	xstar1    = xstar1 * sq_1_rt2
	ystar1    = ystar1 * sq_1_rt2

	xstar     = J(nt, p, .)
	ystar     = J(nt, 1, .)
	if ( rtilde_g_1 ) {
		xstar     = C( J(nt, p, .) )
		ystar     = C( J(nt, 1, .) )
	}
	for (i=1; i<=n; i++) {
		xstar[(ids[i,1]+1)::ids[i,2],] =
			panelsubmatrix(x, i, ids1a) - 
			panelsubmatrix(x, i, ids1b) * rtilde - 
			J(ids1a[i,3],1,xtilde[i,] * (1-rtilde))
		xstar[ids[i,1],] = xstar1[i,]
		ystar[(ids[i,1]+1)::ids[i,2]] =
			panelsubmatrix(y, i, ids1a) - 
			panelsubmatrix(y, i, ids1b) * rtilde - 
			J(ids1a[i,3],1,ytilde[i] * (1-rtilde))
		ystar[ids[i,1]] = ystar1[i]
	}
	
	// calculate \tilde x_{\cdot} = xtilbar

	xtilbar   = mean(xtilde)
	xtilcent  = xtilde - J(n, 1, xtilbar)

	// calculate I_f, \Sigma_1, \Sigma_2, \hat I and  \hat \beta
	
	ztilde1   = ztilde[ids[,1]]
	x1        = x[ids[,1],]
	term1     = ((1-rtilde^2)/stilde2)*(ztilde1'*x1)

	zxrho     = J(n, p, .)
	for (i=1; i<=n; i++) {
		temp1 = panelsubmatrix(ztilde, i, ids1a) - 
				panelsubmatrix(ztilde, i, ids1b) * rtilde
		temp2 = panelsubmatrix(x, i, ids1a) - 
				panelsubmatrix(x, i, ids1b) * rtilde
		zxrho[i,] = cross(temp1, temp2)
	}
	term2     = colsum(zxrho) / stilde2
	
	// density estimator
	// with current bandwith chosen above
	
	//temp      = J(1,n,1) # exp(wtilde / s)
	temp      = J(1,n,1) # exp(wtilde / bwG)
	temp      = temp :/ temp'
	// temp

	// temp_{ij}=exp((\tilde w_i -\tilde w_j)/s)

	//K         = (temp :/ ((1 :+ temp) :^ 2)) / s
	//Kprime    = (temp :* (1 :- temp) :/ ((1 :+ temp) :^3 )) / (s^2)
	K         = (temp :/ ((1 :+ temp) :^ 2)) / bwG
	Kprime    = (temp :* (1 :- temp) :/ ((1 :+ temp) :^3 )) / (bwG^2)
	
	// K is  a (n x n) matrix, so is Kprime (derivatives)

	// sum over j
	K         = colsum(K')
	// now K is a (1 x n) vector: \hat f(w_i),i=1,...,n

	Kprime    = colsum(Kprime')
	// now Kprime is a (1 x n) vector: \hat f^{{1)}(w_i),i=1,...,n

	wpdw      = Kprime :/ K
	// wpdw is a (1 x n) vector of {\hat f{{1)}\over \hat f}

	term3     = -(wpdw * xtilcent)

	I_f       = mean(wpdw':^2)

	// Sigma_1 et Sigma_2

	S1        = (xstar'*xstar)/n
	S2        = (xtilcent'*xtilcent)/n
	if ( rtilde_g_1 ){
		S1    = Re(S1)
		S2    = Re(S2)
	}

	// calculate \hat I and \hat \beta

	Ihat      = S1/stilde2 + S2*I_f
	correc    = invsym(Ihat) * (term1+term2+term3)'/n	
	//bglshat   = btilde + correc
	b_pss2G   = btilde + correc

	//bglshat   = Re(bglshat)
	
	// Variance of beta gls hat, efficient estimator
	//Vbglshat     = invsym(Ihat)/n
	v_pss2G     = invsym(Ihat)/n
}
end

mata 
//mata clear
//mata drop pss2post()
void pss2post(real scalar mycost, real vector y, real matrix x,              ///
    real matrix ids,   real matrix ids1a, real matrix ids1b,                 ///
    real matrix ids2a, real matrix ids2b,                                    ///
    real scalar n,     real scalar nt,    real scalar p,                     ///
    real vector b_pss2W, real vector b_pss2G,                                ///
    EPSS2W, EPSS2W_p, R2W, EPSS2G, EPSS2G_p, R2G,                            ///
    alphapss2W, alphapss2W_p, alphapss2G, alphapss2G_p,                      ///
    aicW, bicW, aicG, bicG, residualW, residualG,                            ///
    rhat, shat, shatW, RSSW, CpW, shatG, RSSG, CpG )
{
	real scalar ey2, R2p2w, R2p2w_adj, trho, shatW2, shatG2
	
	real vector ep2w, reshat, ey, C_i0, C_i1, C_i2, ebig,                    ///
                wtilde, xtilde, wtilde_p, ct0, ztilde, tymch1, ep2g	
	
	// Within ==================================================================
	ep2w      = y - x * b_pss2W
	residualW = ep2w

	reshat    = ep2w
	ep2w      = ep2w :- mean(ep2w)

	// Calculation of R2 
	ey        = y :- mean(y)
	ey2       = cross(ey, ey)
	R2p2w     = 1 - cross(ep2w, ep2w) / ey2
	R2p2w_adj = 1 - (1 - R2p2w) * (nt - 1) / (nt-p-n)
	R2W       = R2p2w, R2p2w_adj

	// wtilde1   = wtilde
	
	// redo the calculations with \hat \beta 
	// for better estimating \hat rho, \hat sigma and \hat alpha

	// get \hat rho= rhat
	// the within residuals are given by uw

	// reshat -- the efficient residuals
	
	C_i0      = J(n, 1, .)
	C_i1      = J(n, 1, .)
	C_i2      = J(n, 1, .)
	for (i=1; i<=n; i++) {
		C_i0[i] = cross(panelsubmatrix(reshat, i, ids),
		                panelsubmatrix(reshat, i, ids)) / ids[i,3]
		C_i1[i] = cross(panelsubmatrix(reshat, i, ids1a),
		                panelsubmatrix(reshat, i, ids1b)) / ids1a[i,3]
		C_i2[i] = cross(panelsubmatrix(reshat, i, ids2a),
		                panelsubmatrix(reshat, i, ids2b)) / ids2a[i,3]
	}
	rhat      = sum(C_i1-C_i2) / sum(C_i0-C_i1)
	//rhat
	
	// revised on Sep. 2003
	if (rhat >=  1) rhat =  0.99    
	if (rhat <= -1) rhat = -0.99
	
	// Get the NEW weights c_t and e_t
	// trho      = (1-rhat^2) + (t-1) * (1-rhat)^2
	//ebig      = J(nt, 1, .)
	// recalculate \tilde w_{i}, \tilde x_{i}, \tilde z_{it} and \hat \sigma^2
	wtilde    = J(n, 1, .)
	//xtilde    = J(n, p, .)
	//ytilde    = J(n, 1, .)
	//wtilde_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		// Get the NEW weights c_t and e_t
		trho           = (1-rhat^2) + (ids[i,3]-1) * (1-rhat)^2
		ct0            = J(ids[i,3], 1, ((1-rhat)^2) / trho)
		ct0[1]         = (1-rhat)/trho
		ct0[ids[i,3]]  = (1-rhat)/trho
		//ebig[ids[i,1]::ids[i,2]] = ct0 * trho * (1+rhat) / ((ids[i,3]-1) * (1-rhat))
		// recalculate \tilde w_{i}, \tilde x_{i}, \tilde z_{it} and \hat \sigma^2
		wtilde[i]      = cross(ct0, panelsubmatrix(reshat, i, ids))
		//xtilde[i,]     = cross(ct0, panelsubmatrix(x, i, ids))
		//wtilde_p[ids[i,1]::ids[i,2]] = J(ids[i,3], 1, wtilde[i])
	}
	//ztilde    = reshat - wtilde_p
	//eltype(ebig[i])
	// it cannot be complex
	//ebig_sqrt = ebig
	//tymch1    = sqrt(ebig) :* ztilde
	//shat2   = tymch1' * tymch1 / n
	//shat2

	// INDIVIDUAL EFFECTS PSS2W for WITHIN

	alphapss2W= wtilde

	PSS2W     = wtilde
	PSS2W     = PSS2W :- mean(PSS2W)

  if (mycost == 1){
    EPSS2W    = exp(min(PSS2W) :- PSS2W)
  }
  else {
    EPSS2W    = exp(PSS2W :- max(PSS2W))
  }
  
	EPSS2W_p  = J(nt, 1, .)
	alphapss2W_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		EPSS2W_p[ ids[i,1]::ids[i,2]]     = J(ids[i,3], 1, EPSS2W[i])
		alphapss2W_p[ ids[i,1]::ids[i,2]] = J(ids[i,3], 1, alphapss2W[i])
	}
   
  shatW2    = variance(residualW)
  residualW = residualW - alphapss2W_p
  aicW      = log((nt-1)/nt*shatW2)+1+2*(p+1)/nt;
  bicW      = log((nt-1)/nt*shatW2)+1+(p+1)*log(nt)/nt;
  shatW     = sqrt(shatW2)
  RSSW      = cross(residualW, residualW)
  CpW       = RSSW/shatW2 - nt + 2*(p+1)
	
	// GLS =====================================================================
	
	//ep2g      = y - x * bglshat
	ep2g      = y - x * b_pss2G
	
	residualG = ep2g
	
	reshat    = ep2g
	ep2g      = ep2g :- mean(ep2g)

	//Calculation of R2 
	//ey        = y :- mean(y)
	R2p2g     = 1 - cross(ep2g, ep2g) / ey2
	R2p2g_adj = 1 - (1 - R2p2g) * (nt - 1) / (nt-p-n)
	R2G       = R2p2g, R2p2g_adj

	//wtilde1   = wtilde

	// redo the calculations with \hat \beta 
	// for better estimating \hat rho, \hat sigma and \hat alpha
	
	// get \hat rho= rhat
	// the within residuals are given by uw

	// reshat -- the efficient residuals
	
	C_i0      = J(n, 1, .)
	C_i1      = J(n, 1, .)
	C_i2      = J(n, 1, .)
	for (i=1; i<=n; i++) {
		C_i0[i] = cross(panelsubmatrix(reshat, i, ids),
		                panelsubmatrix(reshat, i, ids)) / ids[i,3]
		C_i1[i] = cross(panelsubmatrix(reshat, i, ids1a),
		                panelsubmatrix(reshat, i, ids1b)) / ids1a[i,3]
		C_i2[i] = cross(panelsubmatrix(reshat, i, ids2a),
		                panelsubmatrix(reshat, i, ids2b)) / ids2a[i,3]
	}
	rhat      = sum(C_i1-C_i2) / sum(C_i0-C_i1)
	//rhat
	
	// revised on Sep. 2003
	if (rhat >=  1) rhat =  0.99    
	if (rhat <= -1) rhat = -0.99
	
	// Get the NEW weights c_t and e_t
	// trho      = (1-rhat^2) + (t-1) * (1-rhat)^2
	ebig      = J(nt, 1, .)
	// recalculate \tilde w_{i}, \tilde x_{i}, \tilde z_{it} and \hat \sigma^2
	wtilde    = J(n, 1, .)
	xtilde    = J(n, p, .)
	//ytilde    = J(n, 1, .)
	wtilde_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		
		trho           = (1-rhat^2) + (ids[i,3]-1) * (1-rhat)^2
		ct0            = J(ids[i,3], 1, ((1-rhat)^2) / trho)
		ct0[1]         = (1-rhat)/trho
		ct0[ids[i,3]]  = (1-rhat)/trho
		et0            = ct0 * trho * (1+rhat) / ((ids[i,3]-1) * (1-rhat))
		ebig[ids[i,1]::ids[i,2]] = et0		
		
		wtilde[i]      = cross(ct0, panelsubmatrix(reshat, i, ids))
		xtilde[i,]     = cross(ct0, panelsubmatrix(x, i, ids))
		wtilde_p[ids[i,1]::ids[i,2]] = J(ids[i,3], 1, wtilde[i])
	}
	ztilde    = reshat - wtilde_p
	//eltype(ebig[i])
	// it cannot be complex
	//ebig_sqrt = ebig
	tymch1    = sqrt(ebig) :* ztilde
	shat      = sqrt(tymch1' * tymch1 / n)
	//shat2

	//  INDIVIDUAL EFFECTS PSS2G for GLS

	alphapss2G= wtilde

	PSS2G     = wtilde
	PSS2G     = PSS2G :- mean(PSS2G)

  if (mycost == 1){
    EPSS2G    = exp(min(PSS2G) :- PSS2G)
  }
  else {
    EPSS2G    = exp(PSS2G :- max(PSS2G))
  }
	
	EPSS2G_p  = J(nt, 1, .)
	alphapss2G_p  = J(nt, 1, .)
	for (i=1; i<=n; i++) {
		EPSS2G_p[ ids[i,1]::ids[i,2]]     = J(ids[i,3], 1, EPSS2G[i])
		alphapss2G_p[ ids[i,1]::ids[i,2]] = J(ids[i,3], 1, alphapss2G[i])
	}
  
  //RSS       = cross(residual, residual)
  shatG2    = variance(residualG)
  residualG = residualG - alphapss2G_p
  aicG      = log((nt-1)/nt*shatG2)+1+2*(p+1)/nt
  bicG      = log((nt-1)/nt*shatG2)+1+(p+1)*log(nt)/nt
  shatG     = sqrt(shatG2)
  RSSG      = cross(residualG, residualG)
  CpG       = RSSG/shatG2 - nt + 2*(p+1)

}
end

mata 
//mata clear
//mata drop bpss2()
real rowvector bpss2(     real vector y,     real matrix x,                     ///
    real matrix ids,   real matrix ids1a, real matrix ids1b,                 ///
    real matrix ids2a, real matrix ids2b, real vector Ti, real scalar t_max, ///
    real scalar n,     real scalar nt,    real scalar p,                     ///
    real scalar bwW,   real scalar bwG,   real scalar NB )
{
	//400
	
	// ESTIMATION
	pss2onlyB(y, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,               ///
              bwW, bwG, beta_W, v_W, beta_G, v_G )
	//401
	// Getting rhat and shat
	pss2post( 1, y, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,               ///
              beta_W, beta_G,                                                ///
              EPSS2W, EPSS2W_p, R2p2W, EPSS2G, EPSS2G_p, R2p2G,              ///
              alphapss2W, alphapss2W_p, alphapss2G, alphapss2G_p,            ///
              aicW, bicW, aicG, bicG, residualW, residualG,                  ///
              rhat, shat, shatW, RSSW, CpW, shatG, RSSG, CpG )
	//402		  
	crit1     = 0
	crit2     = 0
	NB1       = 0
	
	sigeps1b  = shat / (sqrt(1-rhat^2))
	//t_max     = max(ids[,3])
	
	//403
	
	// BOOTSTRAP DATA GENERATION
	
	for(indb = 1; indb <= NB; indb++){
		// time period 1
		//4031
		epsil     = J(t_max, n, .)
		temp      = rnormal(1, n, 0, sigeps1b)
		epsil[1,] = temp
		// time period 2 to t_max
		for(jt = 2; jt <= t_max; jt++){
			temp       = rhat*temp + rnormal(1, n, 0, shat)
			epsil[jt,] = temp
		}
		//4032
		// every row is a time period
		// every column is a firm
		// use Ti to select relevant elements of epsil[,i]
		epsilon = J(nt, 1, .)
		for (i=1; i <= n; i++) {
			epsilon[ ids[i,1]..ids[i,2] ] = epsil[panelsubmatrix(Ti, i, ids), i]
		}
		//4033

		// RERUN OF PROGRAM
		// With Within as initial
		yb     = x*beta_W + alphapss2W_p + epsilon
		pss2onlyB(yb, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,          ///
              bwW, bwG, beta_Wb, v_W, beta_G, v_G )
		// With GLS as initial
		yb     = x*beta_G + alphapss2G_p + epsilon
		pss2onlyB(yb, x, ids, ids1a, ids1b, ids2a, ids2b, n, nt, p,          ///
              bwW, bwG, beta_W, v_W, beta_Gb, v_G )
		// Check	  
		if( missing(beta_Wb) | missing(beta_Gb) ) NB1++
		// Add
		crit1   = crit1 + sum( (beta_W - beta_Wb):^2 )
		crit2   = crit2 + sum( (beta_G - beta_Gb):^2 )
		
	
	}
	
	//404
	
	// end of bootloop
	
	crit  = crit1, crit2	
	crit  = crit / (NB-NB1)
	
	return( crit )
}
end

mata 
//mata clear
//mata drop bootbpss2()
real matrix bootbpss2( real vector y,     real matrix x,                     ///
    real matrix ids,   real matrix ids1a, real matrix ids1b,                 ///
    real matrix ids2a, real matrix ids2b, real vector Ti,                    ///
    real scalar n,     real scalar nt,    real scalar p,                     ///
    real scalar NB,    real scalar gr0,   real scalar gri, real scalar gr1,  ///
    real scalar tr )
{
	t_max = max(Ti)
	sl        = ceil( (gr1 - gr0) / gri + 1 )
	//102
	//sl
	sgrid     = J(sl, 3, .)
	//103
	//sgrid
	//stata("scalar list `iter'")
	//st_local("iter1", strofreal(sl))
	//st_local("myNB", strofreal(NB))
	//"can it do this?"
	""
	"Calculating optimal bandwidth for the PSS (AR(1) error) estimator"
	"Please be patient!"
	""
	tmp3 = "Going over the grid, which contains", strofreal(sl)
	tmp3 = invtokens(tmp3)
	tmp3 = tmp3, "grid points"
	invtokens(tmp3)
	tmp3 = "(in each grid point,", strofreal(NB)
	tmp3 = invtokens(tmp3)
	tmp3 = tmp3, "bootstrap replications are used)"
	invtokens(tmp3)
	if (tr == 1) stata("_dots 0")
	for(j = 1; j <= sl; j++) {
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
		sgrid[j,2..3] = bpss2(y, x, ids, ids1a, ids1b, ids2a, ids2b,         ///
		                      Ti, t_max, n, nt, p,                           ///
                              sgrid[j,1], sgrid[j,1], NB )
	}
	//110
	//sgrid
	return(sgrid)
}
end
