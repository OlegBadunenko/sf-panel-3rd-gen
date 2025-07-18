*! version 1.0.0  20Mar2020
*! version 1.1.0  24Mar2020
*! version 1.2.0  25Mar2020
*! version 1.2.1  7Apr2020
capture program drop xtsf3gpss2_p
program define xtsf3gpss2_p
	version 11
	syntax newvarname [if] [in] , [ xb te RESIDuals alpha ]
	marksample touse, novarlist
	tempname mysample
	//display "|`varlist'|"
	local case : word count `xb' `te' `residuals' `alpha'
	if `case' >1 {
		display "{err}only one statistic may be specified"
	exit 498 
	}
	if `case' == 0 {
		local n n
		display "expected te, residuals, or alpha"
	}
	generate `mysample' = e(sample)
	quietly summarize `mysample' if `mysample' == 0
	if r(N) > 0 {
		display "{text} (" as result r(N) "{text} missing values generated)"
	}
	if "`xb'" != "" {
		//quietly _predict `typlist' `varlist' if `touse' & `mysample', xb
		//label variable `varlist' "Linear prediction, PSS2 (within estimator)"
		local myvarlist `varlist'_W
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(xbW)"))
		label variable `myvarlist' "Linear prediction, PSS2 (within estimator)"
		local myvarlist `varlist'_G
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(xbG)"))
		label variable `myvarlist' "Linear prediction, PSS2 (GLS estimator)"
	}
	if "`te'" != "" {
		local myvarlist `varlist'_W
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(effW_p)"))
		label variable `myvarlist' "Efficiency, PSS2 (within estimator), time-constant"
		local myvarlist `varlist'_G
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(effG_p)"))
		label variable `myvarlist' "Efficiency, PSS2 (GLS estimator), time-constant"
	}
	if "`residuals'" != "" {
		local myvarlist `varlist'_W
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(residualsW)"))
		label variable `myvarlist' "Residuals, PSS2 (within estimator)"
		local myvarlist `varlist'_G
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(residualsG)"))
		label variable `myvarlist' "Residuals, PSS2 (GLS estimator)"
	}
	if "`alpha'" != "" {
		local myvarlist `varlist'_W
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(alphaW_p)"))
		label variable `myvarlist' "Alpha, PSS2, time-constant (within estimator)"
		local myvarlist `varlist'_G
		quietly generate `myvarlist' = .
		mata: st_store(., st_local("myvarlist"), st_local("mysample"), st_matrix("e(alphaG_p)"))
		label variable `myvarlist' "Alpha, PSS2, time-constant (GLS estimator)"
	}
end
