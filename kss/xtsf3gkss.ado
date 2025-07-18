*! version 1.0.0  3Apr2020
*! version 1.1.0  6Apr2020
*! version 1.2.0  7Apr2020
*! version 1.2.1  8Apr2020
*! version 1.2.2  9Apr2020
*! version 1.2.3  15Apr2020
*! version 1.3.0  16Apr2020
*! version 1.3.4  12Aug2020
*! version 1.4.0  18Sep2020


/*
Origin:

MATLAB code

% Kneip, Sickles, and Song Estimator
% Written by Wonho Song
% Last Update: September 1, 2017
% E-mail: whsong@cau.ac.kr, whsong73@hotmail.com

Translated by Oleg Badunenko
oleg.badunenko@brunel.ac.uk, obadunenko@gmail.com

*/


// if(c(MP)){
// 	set processors 1
// }


capture program drop xtsf3gkss
program define xtsf3gkss, eclass
version 11

  if !replay() {
	syntax varlist(numeric fv min=2) [if] [in]                  ///
		[, GR0(real 0.1) GR1(real 0.2) GRI(real 0.1)      ///
		LI(real 1) LS(real 7) GMEAN IMEAN TMEAN] [COST ///
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
	
	tempname b V n1 nt1 mygmean myimean mytmean R2 R2adj aic bic Ekss alpha   ///
            lstar rez fitted shat1 RSS1 mypanelvar mytimevar Cp1
	
	//display 101
	
	   
	if "`gmean'" == ""{
		//display 1031
		local mygmean 0
	} 
	else {
		local mygmean 1
	}
   if "`imean'" == ""{
		//display 1031
		local myimean 0
	} 
	else {
		local myimean 1
	}
	if "`tmean'" == ""{
		//display 1031
		local mytmean 0
	} 
	else {
		local mytmean 1
	}
	
	//display 102
  display
	display as result "Description of the panel data:" as input "{hline 48}
	xtdescribe if `touse'
  quietly xtset
  local mypanelvar `r(panelvar)'
  local mytimevar `r(timevar)'
  //display "|`mytimevar'|"
  
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
  

	mata: kss_robyty("`depvar'", "`indepvars'", "`touse'", "`mypanelvar'",     ///
      "`mytimevar'",  "`cost'", "`b'", "`V'", "`n1'", "`nt1'",               ///
      "`mygmean'", "`myimean'", "`mytmean'",                                 ///
      "`R2'", "`R2adj'",  "`aic'", "`bic'", "`Ekss'",  "`alpha'",            ///
      "`gr0'", "`gr1'", "`gri'", "`li'",  "`ls'", "`lstar'"   ,              ///
      "`rez'", "`fitted'", "`shat1'", "`RSS1'", "`Cp1'")

	
	// display 105

   
  //display "`indepvars'"
	
	matrix colnames `b'        = `indepvars'
  matrix rownames `V'        = `indepvars'
	matrix colnames `V'        = `indepvars'
	matrix colnames `alpha'   = "alpha"

	matrix colnames `rez'     = "residuals"
	matrix colnames `Ekss'   = "Efficiency"
	ereturn post `b' `V', esample(`touse') buildfvinfo depname("`depvar'")
	ereturn scalar N           = `n1'
	ereturn scalar sumTi          = `nt1'
	ereturn scalar L          = `lstar'
	ereturn scalar r2         = `R2'
	ereturn scalar r2_a       = `R2adj'
	ereturn scalar aic        = `aic'
	ereturn scalar bic        = `bic'
  ereturn scalar cp          = `Cp1'
  ereturn scalar shat        = `shat1'
  ereturn scalar RSS         = `RSS1'
	ereturn matrix eff         = `Ekss'
	ereturn matrix alpha       = `alpha'
	ereturn matrix residuals   = `rez'
	ereturn matrix xb          = `fitted'
	ereturn local predict "xtsf3gkss_p"
	ereturn local cmd   "xtsf3gkss"
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
    display as input "Kneip-Sickles-Song Estimator"
    display as input "Kneip, Sickles, and Song (2012), Econometric Theory, 28(3):590–628"
    display
    display as input " `function'" as text " Stochastic Frontier"
    ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
  }
	
end


mata mata clear

mata:

void kss_robyty( string scalar depvar,                                       ///
    string scalar indepvars, string scalar touse,  string scalar mypanelvar, ///
    string scalar mytimevar, string scalar mycost,                           ///
    string scalar bsname,    string scalar vsname,                           ///
    string scalar n1name,    string scalar nt1name,                          ///
    string scalar gmean0,    string scalar imean0,  string scalar tmean0,    ///
    string scalar R2name,    string scalar R2adjname,                        ///
    string scalar aicname,   string scalar bicname,                          ///
    string scalar Eff_name,  string scalar aname,                            ///
    string scalar gr00,      string scalar gr10,    string scalar gri0,      ///
    string scalar li0,       string scalar ls0,    string scalar Lname ,     ///
    string scalar rezname,   string scalar xbname,                           ///
    string scalar shatname,  string scalar RSSname,
    string scalar Cpname )
{
// 	200

	real scalar bw, subboth, p, n, nt, gr0, gr1, gri, LI, LS, Ti_max, gmean,  ///
   imean, tmean, Cp
	
	real vector y, y0, Ti0,  ids0, tymch1, tymch2, xb, res0
	
	real matrix x, x0, ids, ids1b, ids4a4b0, ids4a3b1, ids4a2b2,  ids4a1b3,   ///
    ids4a0b4, ids4, ids1
	
// 	201
	
	
	y = y0    = st_data(., depvar, touse)
	x = x0    = st_data(., indepvars, touse)
	
// 	202
	
	// panel variable
	//stata("quietly xtset")
	//stata("local mypanelvar `r(panelvar)'")
	//mypanelvar= st_local("mypanelvar")
   
   //stata("local mytimevar `r(timevar)'")
	//mytimevar = st_local("mytimevar")
	Ti0       = st_data(., mytimevar, touse)
	Ti        = Ti0 :- min(Ti0) :+ 1

	if (mycost == ""){
    cost = 0
  }
  else {
    cost = 1
  } 
  
// 	203

	ids0      = st_data(., mypanelvar, touse)
	ids       = panelsetup(ids0, 1)
	ids       = ids, ids[,2] - ids[,1] :+ 1
	
	//205

	p         = cols(x)
	n         = rows(ids)
	nt        = rows(y)
	//t         = round(nt/n)
	gr0       = strtoreal(gr00)
   if (gr0 < 0 |  gr0 > 1) {
      //_error("gr0() must be between 0 and 1 inclusive")
      //errprintf(" argument gr0() must be between 0 and 1 inclusive\n")
      printf("{error:argument} {text: gr0 = } {input:%3.2f} {err: must be between 0 and 1 inclusive;} {result:respecify}\n", gr0)
      exit(error(498))
   }
	gr1       = strtoreal(gr10)
   if (gr1 < 0 |  gr1 > 1) {
      //_error("gr0() must be between 0 and 1 inclusive")
      //errprintf(" argument gr0() must be between 0 and 1 inclusive\n")
      printf("{error:argument} {text: gr1 = } {input:%3.2f} {err: must be between 0 and 1 inclusive;} {result:respecify}\n", gr1)
      exit(error(498))
   }
   if (gr0 > gr1) {
      //_error("gr0() must be between 0 and 1 inclusive")
      //errprintf(" argument gr0() must be between 0 and 1 inclusive\n")
      printf("{error:argument} {text: gr0 = } {input:%3.2f} {err: must smaller than} {text: gr1 = } {input:%3.2f}{error:;} {result:respecify}\n", gr0, gr1)
      exit(error(498))
   }
	gri       = strtoreal(gri0)
   if (gri > gr1-gr0) {
      //_error("gr0() must be between 0 and 1 inclusive")
      //errprintf(" argument gr0() must be between 0 and 1 inclusive\n")
      printf("{error:argument} {text: gri = } {input:%3.2f} {err: must smaller than} {text: gr1 - gr0 = } {input:%3.2f}{error:;} {result:respecify}\n", gri, gr1-gr0)
      exit(error(498))
   }
   LI        = strtoreal(li0)
   LS        = strtoreal(ls0)
   Ti_max    = max(ids[,3])
   gmean     = strtoreal(gmean0)
   imean     = strtoreal(imean0)
   tmean     = strtoreal(tmean0)
   
   // elimination of global mean 
   
//    206

   if (gmean ==1){
      y      = y :- mean(y)
      x      = x - J(nt, 1, mean(x))
   }

   // elimination of individual effects 
   
//    207

   if (imean == 1){
      for (i = 1; i <= n; i++) {
         y[| ids[i,1] ,. \ ids[i,2],. |]  = y[| ids[i,1] ,. \ ids[i,2],. |] - 
          J(ids[i,3], 1, mean( panelsubmatrix(y, i, ids) ))
         x[| ids[i,1] ,. \ ids[i,2],. |]  = x[| ids[i,1] ,. \ ids[i,2],. |] - 
          J(ids[i,3], 1, mean( panelsubmatrix(x, i, ids) ))
      }
   }

   // elimination of time-fixed effects
   
//    208
//    Ti_max
   //ids
   
   ytbar     = J(Ti_max, 1, .)
   xtbar     = J(Ti_max, p, .)
   for(t = 1; t <= Ti_max; t++){
     // ytbar[t,.] = mean( y[ ids[,1] :+ t:-1]  )
     // xtbar[t,.] = mean( x[ ids[,1] :+ t:-1, ]  )      
     ytbar[t,.] = mean( select(y, Ti :== t)  )
     xtbar[t,.] = mean( select(x, Ti :== t)  )
      
   }
   
//    209
   
   
   if (tmean ==1){
      for (i = 1; i <= n; i++) {
         y[| ids[i,1] ,. \ ids[i,2],. |]  = y[| ids[i,1] ,. \ ids[i,2],. |] - 
          ytbar[|1,. \ ids[i,3],. |]
         x[| ids[i,1] ,. \ ids[i,2],. |]  = x[| ids[i,1] ,. \ ids[i,2],. |] - 
          xtbar[|1,. \ ids[i,3],. |]
      }
   }
   
   //210
   
   y_w       = y
   x_w       = x
   
   p_Z_W     = ppZhWhGrid(Ti_max, gr0, gri, gr1)
   
   /*
   The outcome is the matrix `p_Z_W` : 3 x size of grid, where the first column is
   the pointer to the smoothing parameter used to calculate matrices
   Zh  and  Wh , pointers to which are in columns 2 and 3.
   
   Zh  and  Wh are Ti_max x Ti_max  matrices
   */
//    p_Z_W
//    2101
   
   //211
   
   L_p_i     = dim(y_w, x_w, ids, n, nt, p, Ti, Ti_max, LI, LS, 
                   p_Z_W[,1], p_Z_W[,2], p_Z_W[,3])
   /* 
   The output is a 1 x 3 vector:
      * optimal L
      * smoothing parameter when optimal L is obtained
      * order of the smoothing parameter in the grid (this is required
        for getting pointers from 'ppZhWhGrid')
   */
   //L_p_i
   //2111
   
   
   //212
   
   beta  = bcov = se_b = KSS = zz = .

   method( y_w, x_w, ids, ytbar, xtbar, n, nt, p, L_p_i[1], Ti, Ti_max, 
           p_Z_W[L_p_i[3],2], p_Z_W[L_p_i[3],3], beta, bcov, se_b, KSS, zz)

   //213
   
  mysample  = 1::nt
  EKSS      = J(nt, 1, .)
  
  if (mycost == 1){
    for (t = 1; t <= Ti_max; t++){
      mysamplet = mysample :* (Ti :== t)
      EKSS[selectindex(mysamplet)] = 
       exp(min( KSS[selectindex(mysamplet)] )  :- KSS[selectindex(mysamplet)] )
    }
  }
  else {
    for (t = 1; t <= Ti_max; t++){
      mysamplet = mysample :* (Ti :== t)
      EKSS[selectindex(mysamplet)] = 
       exp( KSS[selectindex(mysamplet)] :- max( KSS[selectindex(mysamplet)] ) )
    }
  }
  
  
  xb        = x0 * beta
  res0      = y0 - xb - KSS
  res       = y - x * beta
  ek        = res :- mean(res)

   // Diagnostics 
  ey        = y :- mean(y)
  R2        = 1 - cross(ek,ek) / cross(ey,ey)
  R2a       = 1 - (1-R2) * (nt-1) / (nt-p-n)
   
   // Score

  shat2     = variance(res0)
  shat      = sqrt(shat2)
//    998
  RSS       = cross(res0, res0)
  Cp        = RSS/shat2 - nt + 2*(p+1)
   //999
  aic       = log((nt-1)/nt*shat2)+1+2*(p+1)/nt
   //9991
	bic       = log((nt-1)/nt*shat2)+1+(p+1)*log(nt)/nt
   
   //265

	//xb        = x * bh
	bh        = beta'
	//210
	st_matrix(bsname,        bh)
	//211
   //bcov
	st_matrix(vsname,        bcov)
	//212
	st_numscalar(Lname,      L_p_i[1]);
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
	st_matrix(aname,         KSS)
	//231
	st_matrix(rezname,       res0)
	//233
	st_matrix(Eff_name,      EKSS)
	//235
	st_matrix(xbname,        xb)
	//236
	st_numscalar(shatname,   shat)
  //237
  st_numscalar(RSSname,    RSS)
  //238
  st_numscalar(Cpname,     Cp)
}


end





// begin functions for kss

/*

   This function is identical to 'vfun' in MATLAB codes except it does not 
   calculate  Zh  and  Wh , but takes pointers to where these matricies
   are stored.
   
   The output is a  nt x 1  colvector of  viz

*/
   
capture mata mata drop vfun()
mata:

real colvector vfun( real vector y_w, real matrix x_w, real matrix ids,      ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh)
{
   real colvector y_til, vi, viz

   real matrix x_til

   y_til     = J(nt, 1, .)
   x_til     = J(nt, p, .)
   for (i = 1; i <= n; i++) {
      y_til[| ids[i,1] ,. \ ids[i,2],. |] = (*Wh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       y_w[| ids[i,1] ,. \ ids[i,2],. |]
      x_til[| ids[i,1] ,. \ ids[i,2],. |] = (*Wh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       x_w[| ids[i,1] ,. \ ids[i,2],. |]
   }
   
   vi        = y_w - x_w * qrsolve(cross(x_til, x_til), cross(x_til, y_til))
   
   viz       = J(nt, 1, .)
   for (i = 1; i <= n; i++) {
      viz[| ids[i,1] ,. \ ids[i,2],. |] = (*Zh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       vi[| ids[i,1] ,. \ ids[i,2],. |]
   }
   
   return(viz)
}  
end

/*

   This function is similar to 'vfun', but  viz  for each i (i=1,...,n) is
   caculated using  beta , which is computed leaving  data for  i  out.
  
   The output is a  nt x n  matrix of  viz , in which each column has 
   missing values of viz for  i  panel.
   
   loo :=: leave-one-out
   
*/

capture mata mata drop vfunloo()
mata:

real matrix vfunloo( real vector y_w, real matrix x_w, real matrix ids,      ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh)
{
   real colvector y_til, y_wi, y_tili, vii

   real matrix x_til, x_wi, x_tili, vi, viz

   y_til     = J(nt, 1, .)
   x_til     = J(nt, p, .)
   for (i = 1; i <= n; i++) {
      y_til[| ids[i,1] ,. \ ids[i,2],. |] = (*Wh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       y_w[| ids[i,1] ,. \ ids[i,2],. |]
      x_til[| ids[i,1] ,. \ ids[i,2],. |] = (*Wh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       x_w[| ids[i,1] ,. \ ids[i,2],. |]
   }
   
   vi        = J(nt, n, .)
   // vi will have empty entries in i: leave-one-out
   mysample   = 1::nt
   for (i = 1; i <= n; i++) {
      myloosample = mysample
      myloosample[| ids[i,1],. \ ids[i,2],. |] = J(ids[i,3],1,0)
      y_wi        = select(y_w, myloosample)
      x_wi        = select(x_w, myloosample)
      y_tili      = select(y_til, myloosample)
      x_tili      = select(x_til, myloosample)
      vii         = y_wi - x_wi * 
       qrsolve(cross(x_tili, x_tili), cross(x_tili, y_tili))
      vi[selectindex(myloosample),i] = vii
   }
   
   viz       = J(nt, n, .)
   // viz will have empty entries in j: leave-one-out
   for (j = 1; j <= n; j++) {
      myloosample = mysample
      myloosample[| ids[j,1],. \ ids[j,2],. |] = J(ids[j,3],1,0)
      vizj      = J(nt, 1, .)
      for (i = 1; i <= n; i++) {
         if( i != j ){
            vizj[| ids[i,1] ,. \ ids[i,2],. |] = 
             (*Zh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
             vi[| ids[i,1] ,j \ ids[i,2],j |]       
         }
      }
      viz[selectindex(myloosample),j] = select(vizj, myloosample)
   }
    
   return(viz)
}  
end

/*

   This function is identical to 'bg2' from MATLAB codes.
   beta  and  Wh2  must be set before executing it.
   
   The outputs are modified colvector beta and  Wh2

*/

capture mata mata drop bg2()
mata:

void bg2( real vector y_w, real matrix x_w, real matrix ids,                 ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  real scalar L,                                             ///
                  real vector Ti0, real scalar Ti_max,                       ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh, ///
                  beta, Wh2)
{
   real colvector y_til, viz
   
   real rowvector vva

   real matrix x_til, temp, SIGMA, vve, gt, Zh2
   
   // Step1

   viz       = vfun(y_w, x_w, ids, n, nt, p, Zh, Wh)

   // Step2

   temp      = J(Ti_max, n, 0)
   //write into these time periods: Ti  = panelsubmatrix(Ti0, i, ids)
   for(i = 1; i <= n; i++){
      temp[panelsubmatrix(Ti0, i, ids), i] = panelsubmatrix(viz, i, ids)
   }
   
   SIGMA     = temp * temp' / n
   
   vve       = .
   vva       = .
   symeigensystem(SIGMA, vve, vva)
   // eigenvalues
   //va        = vva'
   // eigenvectors
   //vve
   // number of principals
   gt        = sqrt(Ti_max) * vve[|.,1 \ .,L|]
   
   // Step3
   
   Zh2       = gt * invsym(cross(gt,gt)) * gt'
   Wh2       = I(Ti_max) - Zh2
   
   y_til     = J(nt, 1, .)
   x_til     = J(nt, p, .)
   for (i = 1; i <= n; i++) {
      y_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
       y_w[| ids[i,1] ,. \ ids[i,2],. |]
      x_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
       x_w[| ids[i,1] ,. \ ids[i,2],. |]
   }

   beta      = qrsolve(quadcross(x_til, x_til), quadcross(x_til, y_til))
   
}  
end

/*

   This function is similar to 'bg2', but it uses 'vfunloo' and calculates 
   only a statistic  bb , required to choose optimal smoothing parameter.
   
   The output is a scalar  bb.

*/

capture mata mata drop bg2loo()
mata:

real scalar bg2loo( real vector y_w, real matrix x_w, real matrix ids,       ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  real scalar L,                                             ///
                  real vector Ti0, real scalar Ti_max,                       ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh)
{
   real colvector y_til
   
   real rowvector vva

   real matrix x_til, temp, SIGMA, vve, gt, Zh2, viz
   
   // Step1
   //0
   viz       = vfunloo(y_w, x_w, ids, n, nt, p, Zh, Wh)
   //rows(viz)
   //cols(viz)
   //1
   // now do the leave-one-out loop
   
   //ids
   //Ti_max
   
   //2
   
   mysmplCS   = 1::n
   mysample   = 1::nt
   sp         = J(nt, 1, .)
   for(j = 1; j <= n; j++){
      //11
      // Step2
      //12
      mysmplCSj = mysmplCS
      mysmplCSj[j] = 0
      //13
      temp      = J(Ti_max, n, 0)
      //write into these time periods: Ti  = panelsubmatrix(Ti0, i, ids)
      for(i = 1; i <= n; i++){
        //i
        //panelsubmatrix(viz[,j], i, ids)
        //panelsubmatrix(Ti0, i, ids)
        temp[panelsubmatrix(Ti0, i, ids), i] = panelsubmatrix(viz[,j], i, ids)
      }
      //14
      // select leave-one-out
      temp      = temp[,selectindex(mysmplCSj)]
      //15
      SIGMA     = temp * temp' / n
      //16
      vve       = .
      vva       = .
      symeigensystem(SIGMA, vve, vva)
      // eigenvalues
      //va        = vva'
      // eigenvectors
      //vve
      // number of principals
      //17
      gt        = sqrt(Ti_max) * vve[|.,1 \ .,L|]
      //18
      // Step3
   
      Zh2       = gt * invsym(cross(gt,gt)) * gt'
      //19
      Wh2       = I(Ti_max) - Zh2
      //20
      y_til     = J(nt, 1, .)
      x_til     = J(nt, p, .)
      for (i = 1; i <= n; i++) {
         y_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
          y_w[| ids[i,1] ,. \ ids[i,2],. |]
         x_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
          x_w[| ids[i,1] ,. \ ids[i,2],. |]
      }
      //21
      myloosample = mysample
      myloosample[| ids[j,1],. \ ids[j,2],. |] = J(ids[j,3],1,0)
      y_wi        = y_w[| ids[j,1],. \ ids[j,2],. |]
      x_wi        = x_w[| ids[j,1],. \ ids[j,2],. |]
      y_tili      = select(y_til, myloosample)
      x_tili      = select(x_til, myloosample)
      //22
      vic         = y_wi - x_wi * 
       qrsolve(cross(x_tili, x_tili), cross(x_tili, y_tili))
      //23
      sp[| ids[j,1],. \ ids[j,2],. |] = Wh2[| 1,1 \ ids[j,3], ids[j,3]|] * vic
      //24
   }
   //2
   bb = cross(sp, sp)
   return(bb)
}  
end

/*

   This function calculates  Ti_max x Ti_max  matrices  Zh  and  Wh  for
   different values of smoothing parameter defined by the grid values.
   
   The outcome is the matrix  3 x size of grid, where the first column is
   the pointer to the smoothing parameter used to calculate matrices
   Zh  and  Wh , pointers to which are in columns 2 and 3.

*/

capture mata mata drop ppZhWhGrid()
mata:

pointer(real matrix) matrix ppZhWhGrid( real scalar Ti_max,                 ///
               real scalar gr0,   real scalar gri, real scalar gr1)
{
   real scalar sl, ppp
   real matrix ITi_max, Zh
   real vector t1

   // Step 0
   // generate vector of pointers that contain csaps for grid points
   //1
   //601
   sl        = ceil( (gr1 - gr0) / gri + 1 )
   //sl
   ITi_max   = I(Ti_max)   
   t1        = 1 :: Ti_max
   //602
   //Zhps      = J(sl, 1, NULL)
   //Whps      = J(sl, 1, NULL)
   ZhWh      = J(sl, 3, NULL)
   //603
   for (j = 1; j <= sl; j++){
      ZhWh[j,1] = &(gr0 + (j-1)*gri)
      Zh      = J(Ti_max, Ti_max, .)
      for (t = 1; t <= Ti_max; t++){
         Zh[|1, t \ Ti_max, t|] = 
          csapseval(csaps(t1, ITi_max[,t], *ZhWh[j,1]), t1, t1)
      }
      // Wh      = ITi_max - Zh
      ZhWh[j,2]  = &(Zh * 1)
      ZhWh[j,3]  = &(ITi_max - Zh)
   }
   return(ZhWh)
}

end

/*

   The function is identical to 'dim' in MATLAB codes. The fifference:
   'dim' runs from  LI  to  LS  until either LS is reached or value  cc
   is smaller than 2.33.
   
   The output is a 1 x 3 vector:
      * optimal L
      * smoothing parameter when optimal L is obtained
      * order of the smoothing parameter in the grid (this is required
        for getting pointers from 'ppZhWhGrid')

*/

capture mata mata drop dim()
mata:

real rowvector dim( real vector y_w, real matrix x_w, real matrix ids,       ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  real vector Ti0, real scalar Ti_max,                       ///
                  real scalar LI,  real scalar LS,                           ///
                  pointer(real matrix) vector ppgrid,                        ///
                  pointer(real matrix) vector Zhps,                          ///
                  pointer(real matrix) vector Whps)
{
   //0
   //pointer(real matrix) colvector Zhps
   //pointer(real matrix) colvector Whps
   
   real scalar sl, ff, sig2, up, dw, cc, myL
   
   real matrix sgrid, ITi_max, tymch, SIGMA, gg, Pl, gt, Zh2, Wh2,           ///
    x_til, ZhPlZh
   
   real vector viz, y_til, vi, tymch2, tymch0
   
   // Step 0
   // generate vector of pointers that contain csaps for grid points
   //1
   //sl        = ceil( (gr1 - gr0) / gri + 1 )
   sl        = length(ppgrid)
   sgrid     = J(sl, 3, .)
   sgrid[,3] = 1::sl
   ITi_max   = I(Ti_max)   
   //t1        = 1 :: Ti_max
   //Zhps      = J(sl, 1, NULL)   
   //Whps      = J(sl, 1, NULL)
   for (j = 1; j <= sl; j++){
      //sgrid[j,1] = gr0 + (j-1)*gri
      sgrid[j,1] = *ppgrid[j]
      /*
      Zh      = J(Ti_max, Ti_max, .)
      for (t = 1; t <= Ti_max; t++){
         Zh[|1, t \ Ti_max, t|] = 
          csapseval(csaps(t1, ITi_max[,t], sgrid[j,1]), t1, t1)
      }
      // Wh      = ITi_max - Zh
      Zhps[j]    = &(Zh * 1)
      Whps[j]    = &(ITi_max - Zh)
      */
   }
   
   //2

   // Step1
   
   // we select L such that Cl is less than 2.33
   // "Dimensionality Test of KSS Estimator"
   //Cl        = J(LS, 2, 10)
   
  ""
  printf("{text:Select L such that Cl is less than 2.33}\n")
  printf ("{error:Please be patient!}\n")
  printf ("{input:It may take long time to compute estimates if the data size is large.}\n")
  ""
//   printf ("{text:Choosing the bandwidth which has the smallest MSE for the PSS estimator.}\n")
  printf ("{result}Going over the grid, which contains {input}%f {result}grid points\n", LS-LI+1)
//   printf ("{result}(in each grid point, {input}%f {result}bootstrap replications are used)\n", NB)
   
   stata("_dots 0")
   for (L = LI; L <= LS; L++){
     	st_numscalar("q982769726w", L - LI + 1)
			stata("_dots q982769726w 0")
	
      // 'Leave-one-out' Cross-validation
   
      for (j = 1; j <= sl; j++){
         sgrid[j,2] = 
          bg2loo(y_w, x_w, ids, n, nt, p, L, Ti0, Ti_max, Zhps[j], Whps[j])
      }
      
      sgrid  = sort(sgrid, 2)
   
      // Step2
      
      viz    = vfun(y_w, x_w, ids, n, nt, p, Zhps[sgrid[1,3]], Whps[sgrid[1,3]])
      
      // Step 3
      
      tymch   = J(Ti_max, n, 0)
      //write into these time periods: Ti  = panelsubmatrix(Ti0, i, ids)
      for (i = 1; i <= n; i++){
         tymch[panelsubmatrix(Ti0, i, ids), i] = panelsubmatrix(viz, i, ids)
      }
   
      SIGMA  = tymch * tymch' / n
   
      vve    = .
      vva    = .
      symeigensystem(SIGMA, vve, vva)
 
      //SIGMA*vve - vve*diag(vva)
      // eigenvalues
      ff     = sum(vva[| L+1 \ . |])
      gg     = vve[|.,1 \ .,L|]
      Pl     = ITi_max - gg*gg'
      gt     = sqrt(Ti_max) * gg
      
      Zh2    = gt * invsym(cross(gt,gt)) * gt'
      Wh2    = ITi_max - Zh2
   
      y_til  = J(nt, 1, .)
      x_til  = J(nt, p, .)
      for (i = 1; i <= n; i++) {
         y_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
          y_w[| ids[i,1] ,. \ ids[i,2],. |]
         x_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
          x_w[| ids[i,1] ,. \ ids[i,2],. |]
      }

      vi     = y_w - x_w * 
       qrsolve(quadcross(x_til, x_til), quadcross(x_til, y_til))
      
      viz       = J(nt, 1, .)
      for (i = 1; i <= n; i++) {
         viz[| ids[i,1] ,. \ ids[i,2],. |] = 
          (*Zhps[sgrid[1,3]])[| 1, 1 \ ids[i,3], ids[i,3] |] *
          vi[| ids[i,1] ,. \ ids[i,2],. |]
      }

      tymch2  = vi - viz
      sig2   = sum(tymch2:^2) / ((n-1) * 
       trace( (*Whps[sgrid[1,3]]) * (*Whps[sgrid[1,3]]) ))

      // Step4
      
      ZhPlZh = (*Zhps[sgrid[1,3]]) * Pl * (*Zhps[sgrid[1,3]])
      up     = n * ff - (n-1) * sig2 * trace( ZhPlZh )
      dw     = sqrt( 2 * n * sig2^2 * trace(matpowersym(ZhPlZh,2) ))
      cc     = up/dw

      // disp([up dw cc sig2])

      //Cl[L,] = cc, pp
      //L, cc
      
      myL    = L
      
      if (cc < 2.33) break
      
      if (L == Ti_max) {
         ""
         ""
         " Warning: LS must be smaller or equal to the maximum Ti. L = Ti_max"
         ""
      }
      
      if (L == LS) {
         ""
         ""
         " Warning: Cl > 2.33, consider increasing LS"
         ""
      }
   }
   // 3
   tymch0 = myL, sgrid[1,1], sgrid[1,3]
   
   return(tymch0)
   
}
end

/*

   This function is the same as 'bg' in MATLAB codes
   
   beta,gt,vi,theta,wt,y_til,x_til,va  all must be set before executing it.

*/

capture mata mata drop bg()
mata:

void bg( real vector y_w, real matrix x_w, real matrix ids,                  ///
                  real vector ytbar, real matrix xtbar,                      ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  real scalar L,                                             ///
                  real vector Ti0, real scalar Ti_max,                       ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh, ///
                  beta,gt,vi,theta,wt,y_til,x_til,va)
{
   real vector viz, vva

   real matrix tymch, tymch1, vve, Zh2, Wh2, vireshape
   
   // Step1
   //0
   viz       = vfun(y_w, x_w, ids, n, nt, p, Zh, Wh)

   // Step2
   //1
   tymch      = J(Ti_max, n, 0)
   //write into these time periods: Ti  = panelsubmatrix(Ti0, i, ids)
   for(i = 1; i <= n; i++){
     //i
      tymch[panelsubmatrix(Ti0, i, ids), i] = panelsubmatrix(viz, i, ids)
   }
   //2
   SIGMA     = tymch * tymch' / n
   //3
   vve       = .
   vva       = .
   symeigensystem(SIGMA, vve, vva)
   // eigenvalues
   va        = vva'
   // eigenvectors
   //vve
   // number of principals
   gt        = sqrt(Ti_max) * vve[|.,1 \ .,L|]
   //4
   // Step3
   
   tymch1    = invsym(cross(gt,gt)) * gt'
   Zh2       = gt * tymch1
   Wh2       = I(Ti_max) - Zh2
   //5
   y_til     = J(nt, 1, .)
   x_til     = J(nt, p, .)
   for (i = 1; i <= n; i++) {
      y_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
       y_w[| ids[i,1] ,. \ ids[i,2],. |]
      x_til[| ids[i,1] ,. \ ids[i,2],. |] = Wh2[| 1, 1 \ ids[i,3], ids[i,3] |] *
       x_w[| ids[i,1] ,. \ ids[i,2],. |]
   }
   //6
   beta      = qrsolve(quadcross(x_til, x_til), quadcross(x_til, y_til))
   //7
   vi        = y_w - x_w * beta
   //8
   wt        = (*Zh) * (ytbar - xtbar * beta)
   //9
   viz       = J(nt, 1, .)
   for (i = 1; i <= n; i++) {
      viz[| ids[i,1] ,. \ ids[i,2],. |] = (*Zh)[| 1, 1 \ ids[i,3], ids[i,3] |] *
       vi[| ids[i,1] ,. \ ids[i,2],. |]
   }
   //10
   // Step4
   
   vireshape = J(Ti_max, n, 0)
   for (i = 1; i <= n; i++) {
      vireshape[panelsubmatrix(Ti0, i, ids), i] = panelsubmatrix(vi, i, ids)
   }
   //11
   theta     = tymch1 * vireshape
   
}  
end


/*

   This function is the same as 'method' in MATLAB codes
   
   beta, bcov, se_b, KSS, zz  all must be set before executing it.

*/

capture mata mata drop method()
mata:

void method( real vector y_w, real matrix x_w, real matrix ids,                  ///
                  real vector ytbar, real matrix xtbar,                      ///
                  real scalar n,   real scalar nt, real scalar p,            ///
                  real scalar L,                                             ///
                  real vector Ti0, real scalar Ti_max,                       ///
                  pointer(real matrix) scalar Zh,                            ///
                  pointer(real matrix) scalar Wh, ///
                  beta, bcov, se_b, KSS, zz)
{
   //1110
   
   beta = gt = vi = theta = wt = y_til = x_til = va = .
   bg(y_w, x_w, ids, ytbar, xtbar, n, nt, p, L, Ti0, Ti_max, Zh, Wh, beta,gt,vi,theta,wt,y_til,x_til,va)
   
   //1111
   
   // Step5
   
   err       = y_til - x_til * beta
   bsig2     = cross(err,err) / (nt-n-p)
   bcov      = bsig2 * invsym( cross(x_til,x_til) );
   se_b      = sqrt( diagonal(bcov) )
   
   //1112
   
   KSS       = J(nt, 1, .)
   for (i = 1; i <= n; i++) {
      KSS[| ids[i,1] ,. \ ids[i,2],. |] = gt[panelsubmatrix(Ti0, i, ids),] * 
       theta[,i] + wt[panelsubmatrix(Ti0, i, ids),]
   }
   KSS       = KSS :- mean(KSS)
   
   //1113
   
   // S.e. of theta, the coefficients of factors
   
   //1114
   
   tv        = J(L, n, .)

   for (i = 1; i <= n; i++) {
      gtgtinv   = invsym( cross(gt[panelsubmatrix(Ti0, i, ids),], 
       gt[panelsubmatrix(Ti0, i, ids),]) )
      ee     = vi[| ids[i,1] ,. \ ids[i,2],. |] - 
       gt[panelsubmatrix(Ti0, i, ids),] * theta[,i]
      s2     = cross(ee,ee) / (ids[i,3]-L)
      thcov  = s2 * gtgtinv
      se_th  = sqrt( diagonal(thcov) )
      tv[,i] = theta[,i] :/ se_th
   }
   
   //1115
   
   // dimensionality test
   
   if (L == 1){
      bigI     = J(Ti_max, 1, 1)
   
      // vcov Fixed Effects Estimator
   
      y_wm     = J(n,1,.)
      x_wm     = J(n,p,.)
      y_demean = J(nt,1,.)
      x_demean = J(nt,p,.)
      for (i = 1; i <= n; i++) {
         y_wm[i, .]    = mean( panelsubmatrix(y_w,  i, ids) )
         x_wm[i, .]    = mean( panelsubmatrix(x_w,  i, ids) )
         y_demean[| ids[i,1], .\ ids[i,2], .|] = 
          panelsubmatrix(y_w,  i, ids) - J(ids[i,3], 1, y_wm[i, .])
         x_demean[| ids[i,1], .\ ids[i,2], .|] = 
          panelsubmatrix(x_w,  i, ids) - J(ids[i,3], 1, x_wm[i, .])
      }
      
      // figure collinear
      mynoncoll = 1..p
      for(j=1;j<=p;j++){
         if( variance(x_wm[,j]) < 1e-6 ) mynoncoll[,j] = 0
      }	
      //mynoncoll
      if ( sum(mynoncoll) == 0 ){
         "no variation in  X : no within estimator"
      }
      e_w1   = y_demean - x_demean * 
         qrsolve(x_demean' * x_demean, x_demean' * y_demean)
      sigf   = cross(e_w1, e_w1) / (nt - n - p + 1)
   
      ZhpItZhp = (*Zh)*(I(Ti_max) - bigI * bigI'/ Ti_max) * (*Zh)
      nume     = cross(bigI:- gt,bigI:- gt)/Ti_max :- sigf*trace(ZhpItZhp)/va[1]/n
      deno     = sigf * sqrt( 2*trace( matpowersym(ZhpItZhp,2) ) ) / va[1] / n
      zz       = nume / deno
   }
   else {
      zz     = 10
   }
   
   //1116
   
}

end


// end functions for kss


// begin functions 'csaps'

/*

   Functions setupq, chol1d, csaps, interv, and csapseval are either
   directly translated from fortran codes from 
   # http://pages.cs.wisc.edu/~deboor/pgs/ 
   or adopted for 'csaps' to work.
   
   As a result, 'csaps' is similar to MATLAB's csaps. To evaluate 
   spline on new data, use 'csapseval'.

*/

capture mata mata drop setupq()
mata:

real matrix setupq(real vector x, real vector dx,
                   real vector y, real scalar npoint){
 real colvector qty
 npm1 = npoint - 1
 v = J(npoint, 7, .)
 qty = J(npoint, 1, .)
 //1
 v[1,4] = x[2] - x[1]
 for(i = 2; i <= npm1; i++)  {
  v[i,4] = x[i+1] - x[i]
  v[i,1] = dx[i-1]/v[i-1,4]
  v[i,2] = - dx[i]/v[i,4] - dx[i]/v[i-1,4]
  v[i,3] = dx[i+1]/v[i,4]
 }
 //2
 v[npoint,1] = 0
 //3
 for(i = 2; i <= npm1; i++)  {
  v[i,5] = v[i,1]^2 + v[i,2]^2 + v[i,3]^2
 }
 //4
 if(npm1 >= 3){
  for (i = 3; i <= npm1; i++) {
   v[i-1,6] = v[i-1,2]*v[i,1] + v[i-1,3]*v[i,2]
  }
 }
 //6
 v[npm1,6] = 0
 if(npm1 >= 4){
  for (i = 4; i <= npm1; i++) {
   v[i-2,7] = v[i-2,3]*v[i,1]
  }
 }
 //8
 v[npm1-1,7] = 0
 //9
 v[npm1,7] = 0
 //10
 prev = (y[2] - y[1])/v[1,4]
 //11
 for(i = 2; i <= npm1; i++)  {
  diff = (y[i+1]-y[i])/v[i,4]
  qty[i] = diff - prev
  prev = diff
 }
 //12
 tymch = v, qty
 return(tymch)
}
end

capture mata mata drop chol1d()
mata:

real matrix chol1d( real scalar p, real matrix v, real scalar npoint ){
 //1
 qty = v[,8]
 //2
 qu = u = J(npoint, 1, .)
 npm1 = npoint - 1
 six1mp = 6*(1-p)
 twop = 2*p
 //3
 for (i = 2; i <= npm1; i++) {
  v[i,1] = six1mp*v[i,5] + twop*(v[i-1,4]+v[i,4])
  v[i,2] = six1mp*v[i,6] + p*v[i,4]
  v[i,3] = six1mp*v[i,7]
 }
 //4
 npm2 = npoint - 2
 if(npm2 < 2){
  u[1] = 0
  u[2] = qty[2]/v[2,1]
  u[3] = 0
 } else {
  // factorization
  //41
  for (i = 2; i <= npm2; i++) {
     //411
   ratio = v[i,2]/v[i,1]
   //412
   v[i+1,1] = v[i+1,1] - ratio*v[i,2]
   //413
   v[i+1,2] = v[i+1,2] - ratio*v[i,3]
   //413
   v[i,2] = ratio
   //414
   ratio = v[i,3]/v[i,1]
   //415
   v[i+2,1] = v[i+2,1] - ratio*v[i,3]
   //416
   v[i,3] = ratio
   //417
  }
  //42
  // forward substitution
  u[1] = 0
  v[1,3] = 0
  u[2] = qty[2]
  for (i = 2; i <= npm2; i++) {
   u[i+1] = qty[i+1] - v[i,2]*u[i] - v[i-1,3]*u[i-1]
  }
  //43
  // back substitution
  u[npoint] = 0
  u[npm1] = u[npm1]/v[npm1,1]
  for (i = npm2; i >= 2; i--) {
   u[i] = u[i]/v[i,1]-u[i+1]*v[i,2]-u[i+2]*v[i,3]
  }
  //44
 }
 // construct q*u
 //5
 prev = 0
 for (i = 2; i <= npoint; i++) {
  qu[i] = (u[i] - u[i-1])/v[i-1,4]
  qu[i-1] = qu[i] - prev
  prev = qu[i]
 }
 //6
 qu[npoint] = -qu[npoint]
 tymch = u, qu
 return(tymch)
}

end

// now presentable


capture mata mata drop csaps()
mata:

real matrix csaps( real vector x, real vector y, 
                   real scalar p, | real vector w){
   
   if ((n=length(x)) != length(y)) _error(3200)
   
   if ( args()==3 ) w = J(n,1,1)
   
   if (n != length(w)) _error(3200)
   
   Qmatrices = setupq(x,w,y, n)
   uQu       = chol1d(p,Qmatrices,n)

   six1mp = 6*(1-p)
   rez = J(n, 4, .)
   for (i = 1; i <= n; i++) {
      rez[i,1] = y[i] - six1mp*w[i]^2*uQu[i,2]
   }
   sixp = 6*p
   for (i = 1; i <= n; i++) {
      rez[i,3] = uQu[i,1] * sixp/2
   }
   //npm1 = n - 1
   /*
   for (i = 1; i < n; i++) {
      rez[i,4] = (rez[i+1,3]-rez[i,3])/Qmatrices[i,4]/3
      rez[i,2] = (rez[i+1,1]-rez[i,1])/Qmatrices[i,4] - 
       (rez[i,3]+rez[i,4]/3*Qmatrices[i,4])/2*Qmatrices[i,4]
   }
   // This C. de Boor eq(10), p210: does not work here
   */
   for (i = 1; i < n; i++) {
      rez[i,4] = (rez[i+1,3]-rez[i,3])/Qmatrices[i,4]/3
      rez[i,2] = (rez[i+1,1]-rez[i,1])/Qmatrices[i,4] - 
       (rez[i,3]+rez[i,4]*Qmatrices[i,4])*Qmatrices[i,4]
   }
   // # Reinsch eq(8) and eq(9) WORKS

   tymch = rez[| 1, . \ n-1,. |]
   return(tymch)
}
end

capture mata mata drop interv()
mata:

real colvector interv( real vector xdata, real vector knots ){
  n = length(knots) - 1
  nn = length(xdata)
  ind = J(nn,1,.)
  //2
  for(i = 1; i <= nn; i++){
    ind[i] = sum( xdata[i] :>= knots[|1 \ n |] )
    if(ind[i] == 0) ind[i] = 1
  }
  //3
  return(ind)
}
end

capture mata mata drop csapseval()
mata:

real colvector csapseval( real matrix  coefs, real vector oldX, real vector newX){
 // get indices (size of newX)
 //1
 indexes = interv(newX, oldX)
 //2
 // location
 xdata_loc = newX - oldX[indexes]
 //3
 // initial values
 yidata = coefs[indexes,4]
 //4
 // update
 for(i = 3; i >= 1; i--) {
   // "i begin"
    //i
    //41
    //coefs[indexes,i]
    //42
   coeffs = coefs[indexes,i]
   //43
   yidata = xdata_loc :* yidata + coeffs
    //   "i end"
    //i
 }
 //5
 return(yidata)
}
end

// end functions 'csaps'

// kss y x1 x2 x1_sd x2_sd, gr0(0.3) gr1(.7) gri(0.1) li(3) ls(7) imean tmean
// ereturn list

/*
kss y x1 x2 c.year##c.year, gr0(0.3) gr1(.7) gri(0.1) li(3) ls(7) imean tmean
ereturn list
predict xb_kss, xb
predict res_kss, residuals
predict te_kss, te
predict alpha_kss, alpha

