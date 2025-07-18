*! version 1.0.0  27Mar2020
*! version 1.1.0  28Mar2020
*! version 1.2.0  30Mar2020
*! version 1.3.0  7Apr2020
*! version 1.3.1  8Apr2020
*! version 1.3.2  9Apr2020
*! version 1.3.3  16Apr2020
*! version 1.3.4  12Aug2020
*! version 1.4.0  31Aug2020
*! version 1.4.1  18Jul2025

/*
Origin:

MATLAB code

% PSS3 Estimator
% Written by Wonho Song, October 2002
% Updated August 2004
% Updated May 2014
% E-mail: whsong@cau.ac.kr, whsong73@hotmail.com

Translated by Oleg Badunenko
oleg.badunenko@brunel.ac.uk, obadunenko@gmail.com

*/


// if(c(MP)){
// 	set processors 1
// }


capture program drop xtsf3gpss3
program define xtsf3gpss3, eclass
version 11

  if !replay() {
	syntax varlist(numeric fv min=2) [if] [in]                  ///
		[, GR0(real 0.1) GR1(real 0.2) GRI(real 0.1)      ///
		REPS(integer 699) noDOTS BW(string) COST ///
    LEVEL(string) NOLOG ] ///
    [noCI] [noPValues] [noOMITted] [noEMPTYcells] ///
    [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] ///
    [fvwrap(passthru)] [fvwrapon(passthru)] ///
		[CFORMAT(passthru)] [PFORMAT(passthru)] ///
    [SFORMAT(passthru)] [nolstretch]
    
	marksample touse
	//display 100

	// handle the lists
	gettoken depvar indepvars : varlist
	_fv_check_depvar `depvar'
	
	_rmcoll `indepvars' if `touse', expand `constant'
	local indepvars `r(varlist)'
	
	tempname b V n1 nt1 R2 R2adj Epss3 Epss3_p alpha alpha_p bw_opt sgrid    ///
			 aic bic rez mytrace fitted shat1 RSS1 Cp1 mypanelvar ///
			 Ti
			 
	if "`dots'" == "nodots" {
		local mytrace 0
	}
	else {
		local mytrace 1
	}
	
	if "`bw'" == ""{
		//display 1031
		local bw -999
		local bw_provided 0
	} 
	else {
		//display 1032
		local bw_provided 1
	}
	
	display
  display as result "Description of the panel data:" as input "{hline 48}
	  quietly xtset
  local mypanelvar `r(panelvar)'
  local mytimevar `r(timevar)'
	xtdescribe if `touse'
// 	di 11
	if r(min) <= 4 {
		display
		display as error "IDs with 4 or fewer observations have been excluded from estimation"
		display
		// 		generate sample = `touse'
		quietly bysort `mypanelvar' (`touse'): egen `Ti' = count(`mypanelvar') if `touse'

		* Exclude IDs with 2 or fewer observations
		quietly replace `touse' = 0 if `Ti' <= 4

// 		* Optional: display a message
		xtdescribe if `touse'
	}
// 	di 12
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
		
	//display 104
   
   local AB 1
   local cn 0

	mata: pss3_robyty("`depvar'", "`indepvars'", "`touse'", "`mypanelvar'", ///
    "`cost'", "`AB'", "`cn'",                                           ///
    "`b'", "`V'", "`n1'", "`nt1'",                                       ///
    "`R2'", "`R2adj'",  "`aic'", "`bic'",                                ///
    "`Epss3'", "`Epss3_p'", "`alpha'", "`alpha_p'",                      ///
    "`reps'", "`gr0'", "`gr1'", "`gri'", "`sgrid'",                      ///
    "`bw'", "`bw_opt'", "`rez'", "`fitted'", "`shat1'", "`RSS1'",        ///
    "`Cp1'", "`mytrace'")

	/*
		This is in case two estimators are displayed:
	*/
	
	//display 105

  //display "`indepvars'"
		
	matrix colnames `b'        = `indepvars' "L.`depvar'"
	matrix rownames `V'        = `indepvars' "L.`depvar'"
	matrix colnames `V'        = `indepvars' "L.`depvar'"
	
	if "`bw_provided'" == "0"{
		matrix colnames `sgrid'= "bandwidth" "MSE"
	}
	matrix colnames `alpha'   = "alpha"
	matrix colnames `alpha_p' = "alpha"
	matrix colnames `rez'     = "residuals"
	matrix colnames `Epss3_p' = "Efficiency"
	matrix colnames `Epss3'   = "Efficiency"
	ereturn post `b' `V', esample(`touse') buildfvinfo depname("`depvar'")
	ereturn scalar N           = `n1'
	ereturn scalar sumTi          = `nt1'
	ereturn scalar bandwidth  = `bw_opt'
	ereturn scalar r2         = `R2'
	ereturn scalar r2_a       = `R2adj'
	ereturn scalar aic        = `aic'
	ereturn scalar bic        = `bic'
  ereturn scalar cp          = `Cp1'
  ereturn scalar shat        = `shat1'
  ereturn scalar RSS         = `RSS1'
	ereturn matrix eff_p       = `Epss3_p'
	ereturn matrix eff         = `Epss3'
	if "`bw_provided'" == "0"{
		ereturn matrix bw_grid = `sgrid'
	}
	ereturn matrix alpha_p     = `alpha_p'
	ereturn matrix alpha       = `alpha'
	ereturn matrix residuals   = `rez'
	ereturn matrix xb          = `fitted'
	ereturn local predict "xtsf3gpss3_p"
	ereturn local cmd   "xtsf3gpss3"
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
    display as input " R-squared        " as text "= " as result  %5.4f `R2'
    display as input " Adj R-squared    " as text "= " as result  %5.4f `R2adj'
    display as input " AIC              " as text "= " as result  %5.4f `aic'
    display as input " BIC              " as text "= " as result  %5.4f `bic'
    display as input " Root MSE         " as text "= " as result  %5.4f `shat1'
    display as input "{hline 29}"
    display
    //di as text "{hline 79}"
    display as input "PSS Type 3 estimator: dynamic panel data model"
    display as input "Park, Sickles, and Simar (2007), Journal of Econometrics, 136(1):281–301"
    display
    display as input " `function'" as text " Stochastic Frontier"
    ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
  }
  
end


mata mata clear

mata:

void pss3_robyty( string scalar depvar,                                        ///
    string scalar indepvars, string scalar touse,  string scalar mypanelvar, ///
    string scalar mycost, ///
    string scalar AB, string scalar cn,        ///
    string scalar bsname, string scalar vsname,     ///
    string scalar n1name,    string scalar nt1name,                          ///
    string scalar R2name,    string scalar R2adjname,                        ///
    string scalar aicname,   string scalar bicname,                          ///
    string scalar Eff_name,  string scalar Eff_pname,                        ///
    string scalar aname,     string scalar a_pname,                          ///
    string scalar NB,        string scalar gr00,                             ///
    string scalar gr10,      string scalar gri0, string scalar sg1name,      ///
    string scalar bw0,       string scalar bwname,                           ///
    string scalar rezname,   string scalar xbname,                           ///
    string scalar shatname,  string scalar RSSname,                          ///
    string scalar Cpname,    string scalar mytrace )
{
	//200

	real scalar bw, subboth, p, n, nt, gr0, gr1, gri, tr, s1, Ti_max
	
	real vector y, ids0, tymch1, tymch2, xb
	
	real matrix x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,  ids4a1b3,       ///
    ids4a0b4, ids4, ids1
	
	//201
	
	bw        = strtoreal(bw0)
	
	y         = st_data(., depvar, touse)
	x         = st_data(., indepvars, touse)
	
	ids0      = st_data(., mypanelvar, touse)
	ids       = panelsetup(ids0, 1)
	ids       = ids, ids[,2] - ids[,1] :+ 1
	ids1a     = ids[,1] :+ 1, ids[,2], ids[,3] :- 1
	ids1b     = ids[,1], ids[,2] :- 1, ids[,3] :- 1
	//ids2a     = ids[,1] :+ 2, ids[,2], ids[,3] :- 2
	//ids2b     = ids[,1], ids[,2] :- 2, ids[,3] :- 2

	// skip first 4 and last 0
	ids4a4b0    = ids[,1] :+ 4, ids[,2] :- 0, ids[,3] :- 4
	// skip first 3 and last 1
	ids4a3b1    = ids[,1] :+ 3, ids[,2] :- 1, ids[,3] :- 4
	// skip first 2 and last 2
	ids4a2b2    = ids[,1] :+ 2, ids[,2] :- 2, ids[,3] :- 4
	// skip first 1 and last 3
	ids4a1b3    = ids[,1] :+ 1, ids[,2] :- 3, ids[,3] :- 4
	// skip first 0 and last 3
	ids4a0b4    = ids[,1] :+ 0, ids[,2] :- 4, ids[,3] :- 4
	
// 	ids4a0b4
	
	//204
	
	//stata("local mytimevar `r(timevar)'")
	//mytimevar = st_local("mytimevar")
	//Ti0       = st_data(., mytimevar, touse)
	//Ti        = Ti0 :- min(Ti0) :+ 1
	
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
  AB        = strtoreal(AB)
  cn        = strtoreal(cn)
  Ti_max    = max(ids[,3])
  if (mycost == ""){
    cost = 0
  }
  else {
    cost = 1
  } 
  
  //2031
	// subtract from both sides
	subboth   = 4
	tymch1    = (0::n-1) * subboth
	tymch2    = (1::n)   * subboth
   
	ids4      = ids[,1] - tymch1, ids[,2] - tymch2, ids[,3] :- subboth
	//2032
	// subtract from both sides
	subboth   = 1
	tymch1    = (0::n-1) * subboth
	tymch2    = (1::n)   * subboth
	ids1      = ids[,1] - tymch1, ids[,2] - tymch2, ids[,3] :- subboth

	//206
	
	// Finding bandwidth for PSS estimators using bootstrap
	// if it is not given
	if (bw == -999){
		//2061
		// this is the longest
		// the result is matrix with 3 columns.
		// col 1: bandwidth
		// col 2: MSE_W
		// col 3: MSE_G
		sgrid = bootbpss3( y, x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,     ///
                         ids4a1b3, ids4a0b4, ids4, ids1, n, nt, p, Ti_max,   ///
                         AB, cn, NB, gr0, gri, gr1, tr )
                           
				   
		//2062
		// this is the optimal bandwidth which gives the smallest MSE
		s1   = sort(sgrid, 2)[1,1]
		
		//2063

		""
		""
		//tmp3  = "Optimal bandwidth (for the dynamic panel data estimator) is", strofreal(s1)
		//tmp3  = invtokens(tmp3)
      printf ("{text}Optimal bandwidth (for the dynamic panel data estimator) is {result}%5.4f\n", s1)
		// this is just dispaying in MATA
		//invtokens(tmp3)
		//sgrid
	} else {
		//2064
		""
      printf ("{text}Bandwidth bw ({input}%f{text}) is provided\n", bw)
      printf ("{error}Make sure it was chosen to minimize MSE for the PSS estimator\n")
		/*
      tmp3  = "Bandwidth bw (", strofreal(bw), ")"
		tmp3  = tmp3, " is provided."
		tmp3  = invtokens(tmp3)
		// this is just dispaying in MATA
		invtokens(tmp3)
		tmp3  = "Make sure it was found using optimality criterion"
		invtokens(tmp3)
      */
		s1   = bw
	}
	
   //207
	// ESTIMATION
   pss3onlyB(y, x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,                 ///
    ids4a1b3, ids4a0b4, ids4, ids1, n, nt, p, Ti_max, s1, AB, cn,            ///
    bh, vh, b_hat_tilde)

	//208
	// Post estimation
   pss3post(cost, y, x, ids, ids1, n, nt, p, bh, b_hat_tilde,                ///
    EPSS3, EPSS3_p, R2p3, alphapss3, alphapss3_p, aic, bic,                  ///
    residual, xb, shat, RSS, Cp )  
 
	//209

	//xb        = x * bh
	bh        = bh'
	//210
	st_matrix(bsname,       bh)
	//211
	st_matrix(vsname,       vh)
	//212
	st_numscalar(bwname,     s1);
	//216
	st_numscalar(n1name,     n)
	//217
	st_numscalar(nt1name,    nt)
	//218
	st_numscalar(R2name,     R2p3[1])
	//219
	st_numscalar(R2adjname,  R2p3[2])
	//222
	st_numscalar(aicname,    aic)
	//223
	st_numscalar(bicname,    bic)
	//226
	if (bw == -999){
		st_matrix(sg1name,    sgrid)
	}
	//227
	st_matrix(aname,         alphapss3)
	//228
	st_matrix(a_pname,       alphapss3_p)
	//231
	st_matrix(rezname,       residual)
	//233
	st_matrix(Eff_name,      EPSS3)
	//234
	st_matrix(Eff_pname,     EPSS3_p)
	//235
	st_matrix(xbname,        xb)
	//236
	st_numscalar(shatname,   shat)
   //237
  st_numscalar(RSSname,    RSS)
  st_numscalar(Cpname,     Cp)
}

//mata clear
//mata drop pss3onlyB()
void pss3onlyB( real vector y,     real matrix x,                            ///
    real matrix ids,      real matrix ids1b,    real matrix ids4a4b0,        ///
    real matrix ids4a3b1, real matrix ids4a2b2, real matrix ids4a1b3,        ///
    real matrix ids4a0b4, real matrix ids4,     real matrix ids1,            ///
    real scalar n,        real scalar nt,       real scalar p,               ///
    real scalar Ti_max,   real scalar bw,                                    ///
    real scalar AB,       real scalar cn,                                    ///
    b_pss3, v_pss3, thetatilde )
{
	real scalar nt4, sighat2, gammai, c_tilde, Iw_hat, ksi, I22_A, I22_B,     ///
    I22_C, I22_D, I22_E, I22_F, I22_G, I22_H, I_22
	
	real vector tymch, dropn1, dropnTi, yit, yil, vyit, vyil, ivy, dy2, dy3,  ///
    y2, ybar, vybar, vylbar, ivybar, yitwith, yilwith, vyitwith, vyilwith,   ///
    ivywith, betai, zit, zbar, zitwith, gammaij, xxiw, ziw_tilde, zzitw,     ///
    ctg,cc_with, c_with, wkp, lig_A, lig_B, zitt, yill, lib_B, lig_C,        ///
    li_gamma, I12_A, I12_B, I12_C, I_12, sli
	
	real matrix xit, vxit, xbar, vxbar, vxxit, xitwith, vxitwith, lowgam,     ///
    xiw_tilde, xitw, xitw_with, xxitw, lib_A, xitt, xit_btn, li_beta,        ///
    xiw_btn, li_star, S_wtn, S_btn, I_11, Ihat, Ihati
	
	// INITIAL IV ESTIMATION
   
   //300
	
	dropn1     = 1::nt
	dropn1[ids[,1]] = J(n, 1, 0)
   
   //301
	
	dropnTi    = 1::nt
	dropnTi[ids[,2]] = J(n, 1, 0)
   
   //302
	
	xit        = select(x, dropn1)
	yit        = select(y, dropn1)
	yil        = select(y, dropnTi)
   
   //303
	
	nt4        = nt - 4 * n
	
	vxit       = J(nt4, p, .)
	vyit       = J(nt4, 1, .)
	vyil       = J(nt4, 1, .)
	ivy        = J(nt4, 1, .)
	for (i = 1; i <= n; i++) {
		vxit[ ids4[i,1]::ids4[i,2], . ]  = panelsubmatrix(x, i, ids4a4b0)
		vyit[ ids4[i,1]::ids4[i,2], . ]  = panelsubmatrix(y, i, ids4a4b0)
		
		tymch  = panelsubmatrix(y, i, ids4a3b1)
		
		vyil[ ids4[i,1]::ids4[i,2], . ]  = tymch

		dy2    = panelsubmatrix(y, i, ids4a2b2) - panelsubmatrix(y, i, ids4a1b3)
		dy3    = panelsubmatrix(y, i, ids4a1b3) - panelsubmatrix(y, i, ids4a0b4)
		y2     = panelsubmatrix(y, i, ids4a2b2)
		//y3     = panelsubmatrix(y, i, ids4a1b3)
		
		//x1     = panelsubmatrix(x, i, ids4a3b1)
		//x2     = panelsubmatrix(x, i, ids4a2b2)
		
		//dx1    = panelsubmatrix(x, i, ids4a3b1) - panelsubmatrix(x, i, ids4a2b2)
		//dx2    = panelsubmatrix(x, i, ids4a2b2) - panelsubmatrix(x, i, ids4a1b3)
		//dx3    = panelsubmatrix(x, i, ids4a1b3) - panelsubmatrix(x, i, ids4a0b4)
		
		//iv     = dy2, dy3, y2, x1, dx1, dx2
		iv     = dy2, dy3, y2
		
		ivy[ ids4[i,1]::ids4[i,2], . ]  = iv * invsym(iv'*iv) * iv' * tymch
	}
   
   //304

	// Preliminary Definitions
	
	ybar       = J(n, 1, .)
	xbar       = J(n, p, .)
	for (i = 1; i <= n; i++) {
		ybar[ i, . ] = mean(panelsubmatrix(y, i, ids))
		xbar[ i, . ] = mean(panelsubmatrix(x, i, ids))
	}
   
   //305
	
	vxbar      = J(n, p, .)
	vybar      = J(n, 1, .)
	vylbar     = J(n, 1, .)
	ivybar     = J(n, 1, .)
	for (i = 1; i <= n; i++) {
		vxbar[i,.]   = mean( panelsubmatrix(vxit, i, ids4) )
		vybar[i,.]   = mean( panelsubmatrix(vyit, i, ids4) )
		vylbar[i,.]  = mean( panelsubmatrix(vyil, i, ids4) )
		ivybar[i,.]  = mean( panelsubmatrix(ivy,  i, ids4) )
	}
	
   //306
	
	xitwith    = xit
	yitwith    = yit
	yilwith    = yil
	vxitwith   = vxit
	vyitwith   = vyit
	vyilwith   = vyil
	ivywith    = ivy
	for (i = 1; i <= n; i++) {
		xitwith[|ids1[i,1],. \ids1[i,2],.|] = 
		 xitwith[ ids1[i,1]::ids1[i,2], . ] - J(ids1[i,3], 1, xbar[i, ])
		yitwith[|ids1[i,1],. \ids1[i,2],.|] = 
		 yitwith[ ids1[i,1]::ids1[i,2], . ] - J(ids1[i,3], 1, ybar[i, ])
		yilwith[|ids1[i,1],. \ids1[i,2],.|] = 
		 yilwith[ ids1[i,1]::ids1[i,2], . ] - J(ids1[i,3], 1, ybar[i, ])
		vxitwith[|ids4[i,1],. \ids4[i,2],.|] = 
		 vxitwith[ ids4[i,1]::ids4[i,2], . ] - J(ids4[i,3], 1, vxbar[i, ])
		vyitwith[|ids4[i,1],. \ids4[i,2],.|] = 
		 vyitwith[ ids4[i,1]::ids4[i,2], . ] - J(ids4[i,3], 1, vybar[i, ])
		vyilwith[|ids4[i,1],. \ids4[i,2],.|] = 
		 vyilwith[ ids4[i,1]::ids4[i,2], . ] - J(ids4[i,3], 1, vylbar[i, ])							 
		ivywith[|ids4[i,1],. \ids4[i,2],.|] = 
		 ivywith[ ids4[i,1]::ids4[i,2], . ] - J(ids4[i,3], 1, ivybar[i, ])
	}
   
   //307
	
	// Within OLS / Arelleno and Bond IV Estimator  estimator

	if (AB == 1){
      
      //308
		
      vxxit = vxit, ivy
		//tymch      = panelr4pss3( vyit, vxxit, ids4 )
		panelr4pss3( vyit, vxxit, ids4, thetatilde, sighat2)
      //tymch
      //309
		//p1         = cols(vxxit)
		//tymch
		//thetatilde = tymch[| 1 \ p+1 |]
		//thetatilde
		//setilde    = tymch[| p1+1 \ 2*p1 |]
		//setilde
		//sighat2    = tymch[2*(p+1)+1]
		//sighat2

		// Arelleno-Bond estimator for beta coefficient
		betai      = thetatilde[| 1 \ p |]
		// Arelleno-Bond estimator for lagged dependent variable
		gammai     = thetatilde[p+1]
		// standard error for Arelleno-Bond estimator for beta coefficient
		//se_beta    = setilde[| 1 \ p |]
		// standard error for Arelleno-Bond estimator for lagged dependent variable
		//se_gamma   = setilde[p+1]
	} 
	else {
		//%%%%  Arelleno and Bond IV using OLS %%%%%%%%%
		//%abx=[vxitwith ivywith];
		//%abx2=[vxitwith vyilwith];

		//%thetatilde=inv(abx'*abx)*abx'*vyitwith;   

		//%betai=thetatilde(1:k);
		//%gammai=thetatilde(k+1);

		ss2        = vyitwith - abx2 * thetatilde
		sighat2    = sum(ss2:^2)/(n*(nt/n-4))

		abxwith    = invsym(abx'*abx)
		setilde    = sqrt(diagonal(sighat2*abxwith))
	}
	
   //311
   
	// Definitions
	
	zit        = yit - yil * gammai - xit * betai
	zbar       = J(n,1,.)
	zitwith    = zit
	for (i = 1; i <= n; i++) {
		zbar[i,]   = mean( panelsubmatrix(zit, i, ids1) )
		zitwith[| ids1[i,1], . \ ids1[i,2], .|] = 
		 zitwith[| ids1[i,1], . \ ids1[i,2], .|] - J(ids1[i,3], 1, zbar[i,])
	}
	
   //312
   //gammai
   //Ti_max
	//gammaij = J(1, n, gammai :^ (0::Ti_max-2))
	gammaij = gammai :^ (0::Ti_max-2)
   //gammaij
	
	lowgam  = J(Ti_max-1, Ti_max-1, 0)
	
   //314
   
	//lowertriangle(A, .)
	
	for(j = 1; j <= Ti_max-1; j++){
		lowgam[|j,j \Ti_max-1,j|] = gammaij[|1,1 \Ti_max-j,1|]
	}
	
	xiw_tilde  = J(n, p, .)
	xitw       = J(nt-n, p, .)
	xitw_with  = J(nt-n, p, .)
	for (i = 1; i <= n; i++) {
		xxitw         = lowgam[|1,1 \ ids1b[i,3],ids1b[i,3]|] * panelsubmatrix(x, i, ids1b)
		xxiw          = colsum(xxitw) / ids1b[i,3]
		xiw_tilde[i,] = xxiw
		xitw[| ids1[i,1], . \ ids1[i,2],. |]      = xxitw
		xitw_with[| ids1[i,1], . \ ids1[i,2],. |] = xxitw - J(ids1[i,3], 1, xxiw)
	}
	
	ziw_tilde  = J(n, 1, .)
	for (i = 1; i <= n; i++) {
		zzitw = lowgam[|1,1 \ ids1[i,3], ids1[i,3]|] * panelsubmatrix(zit, i, ids1)
		ziw_tilde[i,] = colsum(zzitw) / ids1[i,3]
	}
	
	ctg        = colsum(lowgam')'
	c_tilde    = sum(ctg) / Ti_max

	cc_with    = ctg - J(Ti_max-1, 1, c_tilde)
	c_with     = J(nt-n, 1, .)
	for (i = 1; i <= n; i++) {
		c_with[| ids1[i,1], . \ ids1[i,2],. |] = cc_with[|1 \ ids1[i,3]|]
	}

	// log_likelihood
	
	// li_beta (sum over i)
	
	wkp        = J(n, 1, .)
	for (i = 1; i <= n; i++) {
		wkp[i] = wprime(zbar[i], zbar, bw) / wkernel(zbar[i], zbar, bw, cn)
	}
	
	lib_A      = J(n, p, .)
	lig_A      = J(n, 1, .)
	lig_B      = J(n, 1, .)
	
	for (i = 1; i <= n; i++) {
		zitt       = panelsubmatrix(zitwith, i, ids1)
		xitt       = panelsubmatrix(xit, i, ids1)
		yill       = panelsubmatrix(yil, i, ids1)
		lib_A[i,]  = colsum(zitt :* xitt / sighat2)
		lig_A[i,]  = sum(zitt :* yill /sighat2)
		lig_B[i,]  = (c_tilde/(ids[i,3]-1) * sighat2) * sum(zitt :^ 2)
	}
	
	xit_btn    = xbar - J(n, 1, mean(xbar))

	lib_B      = wkp :* xit_btn

	li_beta    = lib_A - lib_B
	
	// li_gamma (sum over i)
	
	xiw_btn    = xiw_tilde - J(n, 1, mean(xiw_tilde))
	lig_C      = wkp :* (xiw_btn * betai + ziw_tilde - c_tilde * zbar)
	
	li_gamma   = lig_A + lig_B - lig_C
	
	// li_star
	
	li_star    = li_beta, li_gamma

	// Partition of Ihat
	S_wtn      = cross(xitwith, xitwith) / n
	S_btn      = cross(xit_btn, xit_btn) / n
	Iw_hat     = mean(wkp:^2)

	ksi        = 0
	for (ki = 1; ki <= Ti_max-1; ki++) {
		for (kj = 1; kj <= Ti_max-1; kj++) {
			for (kk = 1; kk <= min(ki\kj)-1; kk++){
				ksi = ksi + gammai^(abs(ki-kj)+2*kk)
			}
		}
	}
	
	// Calculation of I11 matrix
	I_11       = S_wtn / sighat2 + Iw_hat * S_btn
	
	// Calculation of I12 matrix
	I12_A      = colsum(xitw*betai :* xitwith) / n
	ctg_long   = J(nt-n, 1, .)
	for (i = 1; i <= n; i++) {
		ctg_long[| ids1[i,1], . \ ids1[i,2],. |] = ctg[|1 \ ids1[i,3]|]
	}

	I12_B      = colsum( ctg_long :* xitwith) / n * mean(zbar)

	I12_C      = Iw_hat * mean(xiw_btn*betai :* xit_btn)

	I_12       = (I12_A + I12_B) / sighat2 + I12_C

	// Calculation of I22 matrix
	tymch      = xitw_with*betai
	I22_A      = cross(tymch, tymch) / n / sighat2
	I22_B      = 2*(colsum(c_with:*xitw_with)/n) * betai * mean(zbar) / sighat2
	tymch      = J(nt-n, 1, .)
	tymch1     = zbar :^ 2
	for (i = 1; i <= n; i++) {
		tymch[| ids1[i,1], . \ ids1[i,2],. |] = J(ids1[i,3], 1, tymch1[i])
	}
	I22_C      = sum(c_with:^2:*tymch) / n / sighat2
	I22_D      = mean((( -ids[,3]*c_tilde :+ ksi ) * sighat2 ) :/ (ids[,3]:^2) )
	tymch      = xiw_btn*betai
	I22_E      = cross(tymch, tymch) / n

	I22_F      = 0
	for (ki = 1; ki <= Ti_max-1; ki++) {
		for (kj = 1; kj <= ki-1; kj++) {
			I22_F = I22_F + gammai^(2*kj)
		}
	}
	I22_F      = I22_F * mean( -1:/ids[,3] :+ 1 )

	I22_G      = mean( sum(ctg:^2) :/ ids[,3] :+ 2*c_tilde^2:/(ids[,3]:-1) )

	//This part added on August 9, 2004
	tymch      = xiw_tilde*betai+c_tilde*zbar
	I22_H      = cross(tymch, tymch) / n / sighat2


	I_22       = I22_A + I22_B + I22_C + Iw_hat * (I22_D + I22_E) + I22_F - I22_G + I22_H

	Ihat       = I_11, I_12'\I_12, I_22
	
	// FINAL STEP
	
	sli        = colsum(li_star)'
   //334
   //Ihat
   Ihati      = invsym(Ihat)
   //335
   //invsym(Ihat)
   //336
   //Ihatiqr    = qrinv(Ihat)
   //(Ihatiqr+Ihatiqr')/2
	b_pss3     = thetatilde + (1/n) * Ihati * sli
   //337
   //thetatilde
   //338
   //b_pss3

	//b_pss3p     = b_pss3_full[|1 \ p|]
	//b_pss3l    = b_pss3_full[p+1]
	
	v_pss3     = Ihati/n
	//v_pss3p    = v_pss3[|1, 1 \ p, p|]

	//se_pss=sqrt(  diagonal(invsym(Ihat)/n)  );
	//se_pss     = sqrt(  abs(diagonal(invsym(Ihat)/n))  )

	//se_pss3    = se_pss[|1 \ p|]
	//se_pss3l   = se_pss[p+1]
	// se_pss(k+1) is the standard error of lagged dependent variable
    //from the PSS3 estimator
}
end


mata 
//mata clear
//mata drop pss3post()
void pss3post(real scalar mycost, real vector y, real matrix x,              ///
    real matrix ids,    real matrix ids1,                                    ///
    real scalar n,      real scalar nt,    real scalar p,                    ///
    real vector b_pss3, real vector b_hat_tilde,                             ///
    EPSS3, EPSS3_p, R2, alphapss3, alphapss3_p,                              ///
    aic, bic, residual, xb, shat, RSS, Cp )
{
	real scalar eyey, res_sqsum, R2p3, R2p3a
	
	real vector dropn1, dropnTi, yit, yil, rez, ep3, ey, rez1, PSS3
                
  real matrix xit
   
//    401
	
  dropn1     = 1::nt
	dropn1[ids[,1]] = J(n, 1, 0)
   
   //402
	
	dropnTi     = 1::nt
	dropnTi[ids[,2]] = J(n, 1, 0)
   
//    403
	
	xit        = select(x, dropn1)
	yit        = select(y, dropn1)
	yil        = select(y, dropnTi)
   
//    404
   
  xb0        = (xit, yil) * b_pss3
	rez0       = yit - xb0
	ep3        = rez0 :- mean(rez0)
   //rez0
   
//    405

	// Calculation of R2 
	ey         = y :- mean(y)
  eyey       = cross(ey, ey)
  res_sqsum  = cross(ep3, ep3) 
	R2p3       = 1 - res_sqsum / eyey
	R2p3a      = 1 - (1-R2p3) * (nt-1) / (nt-p-n)
  R2         = R2p3, R2p3a
   //aic        = log(res_sqsum/nt)+1+2*(p+1)/nt
   //bic        = log(res_sqsum/nt)+1+(p+1)*log(nt)/nt
   
//    411

	// INDIVIDUAL EFFECTS ALPHA

	rez1      = yit - (xit, yil) * b_hat_tilde
	alphapss3 = J(n,1,.)
	for(i = 1; i <= n; i++){
		alphapss3[i] = mean( panelsubmatrix(rez1, i, ids1) )
	}
	PSS3       = alphapss3 :- mean(alphapss3)
   
//    417

  if (mycost == 1){
    EPSS3    = exp(min(PSS3) :- PSS3)
  }
  else {
    EPSS3    = exp(PSS3 :- max(PSS3))
  }
  
//   "alphapss3"
//   alphapss3
//   "EPSS3"
//   EPSS3
   
   
  EPSS3_p    = J(nt, 1, .)
	alphapss3_p= J(nt, 1, .)
  residual   = J(nt, 1, .)
  xb         = J(nt, 1, .)
	for (i=1; i<=n; i++) {
//       4170
      //i
		EPSS3_p[| ids[i,1] \ ids[i,2] |]      = J(ids[i,3], 1, EPSS3[i])
		alphapss3_p[| ids[i,1] \ ids[i,2] |]  = J(ids[i,3], 1, alphapss3[i])
//       4171
//       alphapss3_p[| ids[i,1] \ ids[i,2] |]
      //4172
      //residual[| ids[i,1] \ ids[i,2] |]
      //4173
      //rez0[| ids1[i,1] \ ids1[i,2] |]
      //4174
      //.\ rez0[| ids1[i,1] \ ids1[i,2] |] 
      //4175
      //.\ alphapss3_p[| ids1a[i,1] \ ids[i,2] |]
      //4176
      //rez0[| ids1[i,1] \ ids1[i,2] |] - alphapss3_p[| ids1a[i,1] \ ids[i,2] |]
      //4177
      residual[| ids[i,1] \ ids[i,2] |] =  .\ rez0[| ids1[i,1] \ ids1[i,2] |] 
      //- .\ alphapss3_p[| ids1a[i,1] \ ids[i,2] |]
//       4178
      xb[| ids[i,1] \ ids[i,2] |]       = .\  xb0[| ids1[i,1] \ ids1[i,2] |]
//       4179
      //residual[ ids[i,1] ] = .
      //alphapss3_p_[| ids1[i,1] \ ids1[i,2] |] = J(ids1[i,3], 1, alphapss3[i])
	}
   
//   420
   
  //rows(residual)
  residual  = residual - alphapss3_p
   
  //421
  shat2     = variance(residual)
  shat      = sqrt(shat2)
  //422
  RSS       = cross(residual, residual)
  Cp        = RSS/shat2 - nt + 2*(p+1)

  aic       = log((nt-1)/nt*shat2)+1+2*(p+1)/nt
  bic       = log((nt-1)/nt*shat2)+1+(p+1)*log(nt)/nt
}
end

mata 
//mata clear
//mata drop bpss3()
real scalar bpss3( real vector y,     real matrix x,                         ///
    real matrix ids,      real matrix ids1b,    real matrix ids4a4b0,        ///
    real matrix ids4a3b1, real matrix ids4a2b2, real matrix ids4a1b3,        ///
    real matrix ids4a0b4, real matrix ids4,     real matrix ids1,            ///
    real scalar n,        real scalar nt,       real scalar p,               ///
    real scalar Ti_max,   real scalar bw,                                    ///
    real scalar AB,       real scalar cn,                                    ///
    real scalar NB )
{
   real scalar crit3, NB1
   
   real vector xbPlusA, epsil, yb
   
	//600
	
	// ESTIMATION
  pss3onlyB(y, x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,                 ///
    ids4a1b3, ids4a0b4, ids4, ids1, n, nt, p, Ti_max, bw, AB, cn,            ///
    bh, vh, b_hat_tilde)

	//601
	// Post estimation
  pss3post( 1, y, x, ids, ids1, n, nt, p, bh, b_hat_tilde,                     ///
    EPSS3, EPSS3_p, R2p3, alphapss3, alphapss3_p, aic, bic, 
    residual, xb, shat, RSS, Cp )  
   
   //602
   
  xbPlusA   = x * bh[|1 \ p|] +  alphapss3_p
   
   //xbPlusA[|1,1\10,1|]
   
   //603

  gam       = bh[p+1]
   
  //604
	 
	crit3     = 0
	NB1       = 0
		
	//406
	
	// BOOTSTRAP DATA GENERATION
	
	for(indb = 1; indb <= NB; indb++){
		// time period 1
		//6061
		epsilon = rnormal(nt, 1, 0, shat)
      
      yb      = J(nt, 1, .)
      
      for(i = 1; i <= n; i++){
         //epsiloni       = panelsubmatrix(epsilon,  i, ids)
         yb[ ids[i,1] ] = xbPlusA[ ids[i,1] ] + epsilon[ ids[i,1] ]
         for(t = 1; t < ids[i,3]; t++){
            yb[ ids[i,1] + t ] = 
             yb[ ids[i,1] + t - 1 ] * gam + 
             xbPlusA[ ids[i,1] + t ] + 
             epsilon[ ids[i,1] + t ]
         }
      }
      
      //yy0 = y, yb
      //yy0[|1,1\10,2|]
      
      
      //6062

		// RERUN OF PROGRAM
		pss3onlyB(yb, x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,             ///
       ids4a1b3, ids4a0b4, ids4, ids1, n, nt, p, Ti_max, bw, AB, cn,         ///
       bhb, vh, b_hat_tilde)
		// Check	  
		if( missing(bhb) ) NB1++
		// Add
      //6063
      //bh, bhb, bh - bhb, (bh - bhb):^2
		crit3   = crit3 + sum( (bh - bhb):^2 )
	}
	
	//607
	
	// end of bootloop
	
	crit3  = crit3 / (NB-NB1)
	
	return( crit3 )
}
end

mata 
//mata clear
//mata drop bootbpss3()
real matrix bootbpss3( real vector y,     real matrix x,                     ///
    real matrix ids,      real matrix ids1b,    real matrix ids4a4b0,        ///
    real matrix ids4a3b1, real matrix ids4a2b2, real matrix ids4a1b3,        ///
    real matrix ids4a0b4, real matrix ids4,     real matrix ids1,            ///
    real scalar n,        real scalar nt,       real scalar p,               ///
    real scalar Ti_max,   real scalar AB,       real scalar cn,              ///
    real scalar NB,    real scalar gr0,   real scalar gri, real scalar gr1,  ///
    real scalar tr )
{
   //501
	sl        = ceil( (gr1 - gr0) / gri + 1 )
	//502
	//sl
	sgrid     = J(sl, 2, .)
	//503
	//sgrid
	//stata("scalar list `iter'")
	//st_local("iter1", strofreal(sl))
	//st_local("myNB", strofreal(NB))
	//"can it do this?"
	//""
	//"Calculating optimal bandwidth for the PSS (dynamic panel data model) estimator"
	//"Please be patient!"
	""
   printf("{text:Calculating optimal bandwidth for the PSS (dynamic panel data model) estimator.}\n")
   printf ("{error:Please be patient!}\n")
   printf ("{input:It may take long time to compute estimates if the data size is large.}\n")
   ""
   printf ("{text:Choosing the bandwidth which has the smallest MSE for the PSS estimator.}\n")
   printf ("{result}Going over the grid, which contains {input}%f {result}grid points\n", sl)
   printf ("{result}(in each grid point, {input}%f {result}bootstrap replications are used)\n", NB)
   /*
	tmp3 = "Going over the grid, which contains", strofreal(sl)
	tmp3 = invtokens(tmp3)
	tmp3 = tmp3, "grid points"
	invtokens(tmp3)
	tmp3 = "(in each grid point,", strofreal(NB)
	tmp3 = invtokens(tmp3)
	tmp3 = tmp3, "bootstrap replications are used)"
	invtokens(tmp3)
   */
   //504
	if (tr == 1) stata("_dots 0")
	for(j = 1; j <= sl; j++) {
      //5041
		if (tr == 1) {
			st_numscalar("q982769726w", j)
			stata("_dots q982769726w 0")
			//st_local("iter1", strofreal(j))
			//stata("_dots `=scalar(iter1)' 0")
		}	
		sgrid[j,1] = gr0 + (j-1)*gri
		//5042
		//j
		//sgrid[j,1]
		//sgrid[j,2] = bpss1(y,x,p1,n,t,nt,sgrid[j,1],NB)
		sgrid[j,2] = bpss3(y, x, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,  ///
 ids4a1b3, ids4a0b4, ids4, ids1, n, nt, p, Ti_max, sgrid[j,1], AB, cn, NB )
      //5043
	}
	//505
	//sgrid
	return(sgrid)
}
end

capture mata mata drop wkernel()
mata
real scalar wkernel(numeric scalar z, numeric  colvector zbar, ///
                   numeric scalar bn, numeric scalar cn)
{
	real vector u, ww
	u          = -zbar :+ z
	ww         = mean( Kernel1(u, bn)) + cn
	return(ww)
}
end


capture mata mata drop Kernel1()
mata
real colvector Kernel1(numeric colvector u, numeric scalar bn)
{
	real vector kbn
	kbn = Kernel2(u/bn) / bn
	return(kbn)
}
end
	
capture mata mata drop Kernel2()
mata:
real colvector Kernel2(numeric colvector u)
{
	real vector ku, emu
	emu         = exp(-u)
	ku          = emu :* (emu :+ 1) :^ (-2)
	return(ku)
}
end



capture mata mata drop wprime()
mata
real scalar wprime(numeric scalar z, numeric  colvector zbar, numeric scalar bn)
{
	real vector u, wp
	u          = -zbar :+ z
	wp         = mean( Krnl1(u, bn) )
	return(wp)
}
end

capture mata mata drop Krnl1()
mata
real colvector Krnl1(numeric colvector u, numeric scalar bn)
{
	real vector kbn
	kbn = Krnl2(u/bn) / bn^2
	return(kbn)
}
end
	
capture mata mata drop Krnl2()
mata:
real colvector Krnl2(numeric colvector u)
{
	real vector ku, emu
	emu         = exp(-u)
	ku          = emu :* (emu :- 1) :/ (emu :+ 1) :^ 3
	return(ku)
}
end
	

capture mata mata drop panelr4pss3()
mata
void panelr4pss3( real vector y, real matrix x, real matrix ids,
             bg, s2g)
{
	nt         = rows(x)
	p          = cols(x)
	n          = rows(ids)
	
	//myP        = J(nt, nt, 0)
	//myZ        = J(nt, n, 0)
	ym         = J(n,1,.)
	xm         = J(n,p,.)
	y_demean   = J(nt,1,.)
	x_demean   = J(nt,p,.)
	for (i = 1; i <= n; i++) {
		//myP[| ids[i,1], ids[i,1] \ ids[i,2], ids[i,2] |] = 
		// J(ids[i,3], ids[i,3], 1 / ids[i,3])
		//myZ[| ids[i,1], i \ ids[i,2], i |] = 
		// J(ids[i,3], 1, 1 )
		ym[i, .]    = mean( panelsubmatrix(y,  i, ids) )
		xm[i, .]    = mean( panelsubmatrix(x,  i, ids) )
		y_demean[| ids[i,1], .\ ids[i,2], .|] = 
		 panelsubmatrix(y,  i, ids) - J(ids[i,3], 1, ym[i, .])
		x_demean[| ids[i,1], .\ ids[i,2], .|] = 
		 panelsubmatrix(x,  i, ids) - J(ids[i,3], 1, xm[i, .])
	}
	
	// see for which variance is 0
   toinclude0 = toinclude = 1..p
	for (j = 1; j <= p; j++) {
      if ( variance(xm[,j] ) < 1e-8 ) toinclude[j] = 0
	}
   //toinclude
   // if variance is zero for more than 1 of X's
   // allow of zero variant X's to be included
   if( sum(toinclude0 :!= toinclude) > 1 ){
      for (j = p; j >= 1; j++) {
         if (toinclude[j] == 0){
            toinclude[j] = toinclude0[j]
            break
         }
      }
   }
   else{
      toinclude = toinclude0
   }  
   //toinclude
   // between
   //stata("xtreg y x1 x2 x1_sd x2_sd t t2, be")
   xm1        = select(xm, toinclude)
   //xm1        = xm
	bb         = qrsolve(xm1' * xm1, xm1' * ym)
	//bb
	eb         = ym - xm1 * bb
	//eb
	//ssr_b_star = sum(ids[,3] :* eb:^2)
   ssr_b      = sum( eb:^2 )

	
	// within
	
	bw         = qrsolve(x_demean' * x_demean, x_demean' * y_demean)
	//bw
	ew         = y_demean - x_demean * bw
	s2e        = cross(ew, ew) / (nt - n - sum(bw :!= 0))

   SA = 0
	//for unbalanced panels
	if(SA){
      //stata("xtreg y x1 x2 x1_sd x2_sd t t2, re sa")
      //stata("ereturn list")
      //stata("display e(sigma_u)")
      ctr        = trace( invsym(x' * myP * x ) * x' * myZ * myZ' * x)
      s2u        = ( ssr_b_star - (n-p)*s2e ) / (nt-ctr)
      if(s2u < 0) s2u = 0
      //sqrt(s2u)
      //sigma_u = st_numscalar("e(sigma_u)")
      //sigma_u
      thetai     = 1 :- sqrt(s2e :/ (ids[,3] * s2u :+ s2e) )
      //14
      //thetai
      //stata("display e(thta_min)")
      //stata("display e(thta_max)")
   }
   else {
      Tbar = n / sum(1:/ids[,3])
      s2u = ssr_b / (n-p) - s2e / Tbar
      if(s2u < 0) s2u = 0
      //14
      //thetai
      //stata("display e(thta_min)")
      //stata("display e(thta_max)")
   }
   thetai     = 1 :- sqrt(s2e :/ (ids[,3] * s2u :+ s2e) )
	
	y_q_demean = J(nt,1,.)
	x_q_demean = J(nt,p,.)
	for (i = 1; i <= n; i++) {
		y_q_demean[| ids[i,1], . \ ids[i,2], . |] = 
		 panelsubmatrix(y,  i, ids) - thetai[i] * J(ids[i,3], 1, ym[i, .])
		x_q_demean[| ids[i,1], . \ ids[i,2], . |] = 
		 panelsubmatrix(x,  i, ids) - thetai[i] * J(ids[i,3], 1, xm[i, .])
	}

	//19
	//x_q_demean[|1,. \ 20, .|]
	//x_q_demean
	xx         = x_q_demean' * x_q_demean
	//20
	//xx 
	bg         = qrsolve(xx, x_q_demean' * y_q_demean)
	//21
	//bg
	eg         = y_q_demean - x_q_demean * bg
	//22
	//eg
	s2g        = cross(eg, eg) / (nt - sum(bg :!= 0) )
	//23
	//s2g
	//seg        = sqrt(s2g*(diagonal(invsym(xx))))
	//24
	//seg
	//rez        = bg \ seg \ s2g
	//return(rez)
}

	//panelr4pss3( vyit, vxxit, ids4)

end

