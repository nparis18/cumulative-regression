**********************************************************************************
*Cumulative OLS Testing
********************************************************************************
clear all
set more off, permanently
matrix drop _all
capture log close
set scheme cleanplots

{
	if "`c(username)'" == "nicolasparis" {
		global direc "/Users/nicolasparis/Dropbox/Updated Estimation"
	}
	if "`c(username)'" == ""{
		global direc "\TradeFairs"
}

}

gl source  "${direc}/cumulative-regression"
gl data  "${direc}/data"

}


cap mkdir $data

* Load data and save as text format
sysuse auto
label drop _all
keep price mpg length weight turn
*order price length mpg 
export delimited using "$data/testData.raw", delim(",") replace nolabel

* Run command
 do "$source/cumulativelsUser.ado"

*Arguments
gl blocksize=10

*******************************************************************************
* Testing
*******************************************************************************

cumulativels length price mpg, filename("$data/testData.raw") blocksize($blocksize) cluster(turn)

reg length price mpg [aw=weight], r
