*! version 1.0.0  27Mar2020
*! version 1.1.0  7Apr2020
*! version 1.1.1  5Aug2020

capture program drop xtsf3gpss3_p
program define xtsf3gpss3_p
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
		label variable `myvarlist' "Linear prediction, PSS3"
	}
	if "`te'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(eff_p)"))
		label variable `varlist' "Efficiency, PSS3, time-constant"
	}
	if "`residuals'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(residuals)"))
		label variable `varlist' "Residuals, PSS3"
	}
	if "`alpha'" != "" {
		quietly generate `varlist' = .
		mata: st_store(., st_local("varlist"), st_local("mysample"), st_matrix("e(alpha_p)"))
		label variable `varlist' "Alpha, PSS3, time-constant"
	}
end
