*! version 1.0.0  7Apr2020
*! version 1.1.0  16Apr2020
*! version 1.2.0  18Sep2020
capture program drop xtsf2gbi_p
program define xtsf2gbi_p
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
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(xb)"))
		label variable `varlist' "Linear prediction, bie"
	}
	if "`te'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(eff)"))
		label variable `varlist' "Efficiency, bie, time-varying"
	}
	if "`residuals'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(residuals)"))
		label variable `varlist' "Residuals, bie"
	}
	if "`alpha'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(alpha)"))
		label variable `varlist' "Alpha, bie, time-varying"
	}
end
