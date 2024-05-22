**********************************************************************************
*Cumulative OLS Testing
********************************************************************************
clear all
set more off, permanently
matrix drop _all
//

cap which elabel
if _rc!=0 ssc install elabel
//cap help cleanplots
//if _rc!=0 ssc install cleanplots
//set scheme cleanplots

if "`c(username)'" == "nicolasparis" {
    global direc "/Users/nicolasparis/Dropbox/Updated Estimation"
    gl source "${direc}/cumulative-regression"
    gl data  "${direc}/data"

}
if "`c(username)'" == "dcc213"| "`c(username)'" == "damian" {
    global direc "~/code"
    gl source "${direc}/cumulatitve-regression"
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

cumulativels length price mpg [aw=weight], filename("$data/testData.raw") blocksize($blocksize) r nocons
reg length price mpg [pw=weight], r nocons
