cap program drop cumulativels 

program cumulativels, eclass
version 15.0

*-------------------------------------------------------------------------------
*--- (0) Syntax diagram
*-------------------------------------------------------------------------------
#delimit;
syntax namelist(min=2),
    filename(string)
    BLOCKsize(integer)
    [
	noCONStant
	sort
	ABSorb(string)
	weight(string)
	]
;
#delimit cr

*-------------------------------------------------------------------------------
*--- (1) Unpack syntax
*-------------------------------------------------------------------------------
if "`constant'"=="noconstant" {
	local nocons = 1
}
else local nocons = 0
if "`sort'"=="sort" {
	local sort = 1
}
else local sort = 0
local K1 : word count `namelist'
local K = `K1'-1
local arglist `namelist'
local Ylist: word 1 of `namelist'
local Xlist: list namelist - Ylist

*-------------------------------------------------------------------------------
*--- (2) Run runcls and report output
*-------------------------------------------------------------------------------
#delimit ;
mata: A=runcls("`filename'","`arglist'","`weight'","`absorb'",
				`blocksize',`K',`nocons',`sort');
#delimit cr
mata: betas=*(A[1,1])'
mata: vcov=*(A[1,2])

mata: st_matrix("b", betas)
mata: st_matrix("V", vcov)

matrix colnames b  = `Xlist'
matrix colnames V  = `Xlist'
matrix rownames V  = `Xlist'

ereturn post b V, dep(`Ylist') obs(10)
ereturn display
end

*-------------------------------------------------------------------------------
*--- (3) Define auxiliary Mata details
*-------------------------------------------------------------------------------

version 15.0

mata:
struct resultsCLS {
	real scalar yy,N,W,Ybar
	real colvector Xy
	real matrix XX,Xbar
}


//------------------------------------------------------------------------------
//--- (4) Main wrapper function
//------------------------------------------------------------------------------
function runcls(string scalar filename,			// File containing data
					  string scalar arglist,	// String containing var names
					  string scalar weight,		// Weight's name variable
					  string scalar absorb,		// Fe's name variable
					  real scalar blocksize,	// Size in observations for blocks
					  real scalar K,			// Number of dep vars
					  real scalar nocons,
					  real scalar sort)
{
	cons=1-nocons
	covs=K+cons // Number of actual covs when constand is defined
	
	struct resultsCLS scalar r1
	r1.yy=r1.N=r1.W=r1.Ybar= 0
	r1.Xy		= J(covs,1,0)
	r1.XX		= J(covs,covs,0)
	r1.Xbar		= J(1,K,0)

	data=J(blocksize,K+3,.)
	datap=&data

	//ESTIMATION
	//Enter loop for clustering or repeat estimation

	groupsAccP=clusterFunOLS(datap,filename,arglist,weight,absorb,K,cons,blocksize,sort)
	Ng=length(groupsAccP) // N° groups

	// Cluster Boostrap & Multi Group OLS
	r1=*(groupsAccP[1,1])
	
	// FE ESTIMATION
	
	Xbar=&(r1.Xbar/r1.W)
	Ybar=&(r1.Ybar/r1.W)
	
	XXbar=&(quadcross(*Xbar,cons,*Xbar,cons)*r1.W)
	yybar=&(quadcross(*Ybar,*Ybar)*r1.W)
	Xybar=&(quadcross(*Xbar,cons,*Ybar,0)*r1.W)
	
	r1.XX=r1.XX-*XXbar
	r1.yy=r1.yy-*yybar
	r1.Xy=r1.Xy-*Xybar

	XX_total=&(r1.XX)
	yy_total=&(r1.yy)
	Xy_total=&(r1.Xy)
	N=r1.N
	
	for (j=2; j<=Ng; j++) {
		r1=*(groupsAccP[j,1])

		Xbar=&(r1.Xbar/r1.W)
		Ybar=&(r1.Ybar/r1.W)
		XXbar=&(quadcross(*Xbar,cons,*Xbar,cons)*r1.W)
		yybar=&(quadcross(*Ybar,*Ybar)*r1.W)
		Xybar=&(quadcross(*Xbar,cons,*Ybar,0)*r1.W)
		
		r1.XX=r1.XX-*XXbar
		r1.Xy=r1.Xy-*Xybar
		r1.yy=r1.yy-*yybar

		XX_total=&(*XX_total+r1.XX)
		yy_total=&(*yy_total+r1.yy)
		Xy_total=&(*Xy_total+r1.Xy)
			   N=N+r1.N
	}

	XXinv=&(invsym(*XX_total))
	beta=&(*XXinv**Xy_total)

	// Non-redundant degrees of freedom?
	KStd=&(colsum((*beta:!=0)))
	// 'U = (yy - 2BXy) +  B*X'X*B
	uPu= (*yy_total - 2*(*beta)'**Xy_total + (*beta)'**XX_total**beta)

	//Generate variance-covariance matrix
	// V(B)
	vcov=&(uPu/(N-*KStd-Ng)**XXinv)
// 	namesOutput=(dataNames[varPos])[2..K]

	output=(beta,vcov)
	return(output)
}

//------------------------------------------------------------------------------
//--- (X) Functions called from main
//------------------------------------------------------------------------------

struct resultsCLS scalar updateCLS(pointer X,
								   pointer y,
								   pointer w,
								   real scalar N,
								   real scalar cons,
								   struct resultsCLS scalar r)
{
	struct resultsCLS scalar a
	a.W   =r.W+quadcolsum(*w)
	a.yy  =r.yy+quadcross(*y,*y)
	a.Ybar=r.Ybar+quadcross(*w,*y)
	a.Xbar=r.Xbar+quadcross(*w,*X)
	a.Xy  =r.Xy+quadcross(*X,cons,*w,*y,0)
	a.XX  =r.XX+quadcross(*X,cons,*w,*X,cons)

	a.N =r.N+N
	return(a)
}

real colvector vec_inlistOLS(pointer B, colvector L)
{
    real colvector b, l
    real scalar minrows, answer
    
    b = J(max(*B), 1, 0)
    b[*B] = J(rows(*B), 1, 1)
    
    l = J(max(L), 1, 0)
    l[L] = J(rows(L), 1, 1)
    
    minrows = min((rows(b), rows(l)))
    answer = J(rows(b), 1, 0)
    answer[|1, 1 \ minrows, 1 |] =
        b[|1, 1 \ minrows, 1 |] :* l[|1, 1 \ minrows, 1 |]
    
    answer= answer[*B]
    return(answer)
}

struct resultsCLS matrix locationOLS(scalar Tgroup, scalar K, scalar cons) 
{
	real scalar j
	pointer(struct resultsCLS scalar) matrix Index
	covs=K+cons
	
	Index = &(resultsCLS())
	Index[1,1]->XX	= J(covs,covs,0)
	Index[1,1]->Xy	= J(covs,1,0)
	Index[1,1]->Xbar= J(1,K,0)
	Index[1,1]->yy	= 0
	Index[1,1]->Ybar= 0	
	Index[1,1]->N	= 0
	Index[1,1]->W	= 0
	if (Tgroup>1) {
		for (j=2; j<=Tgroup; j++) {
			Index = (Index \ &(resultsCLS()))

			Index[j,1]->XX  = J(covs,covs,0)
			Index[j,1]->Xy  = J(covs,1,0)
			Index[j,1]->Xbar= J(1,K,0)
			Index[j,1]->yy  = 0
			Index[j,1]->Ybar= 0	
			Index[j,1]->N   = 0
			Index[j,1]->W   = 0
		}
	}
	return(Index)
}

struct resultsCLS matrix clusterFunOLS(pointer data,
									   string scalar filename,
									   string scalar arglist,
									   string scalar weight,
									   string scalar absorb,
									   real scalar K,
									   real scalar cons,
									   real scalar blocksize,
									   real scalar sort)
{
	real scalar	i,itera
	i=itera=1
	fh = fopen(filename,"r")
	names = fget(fh)
	dataNames=tokens(subinstr(names,",", " "))	// Names from textfile
	userNames=tokens(arglist)					// Names from user
	weightName=tokens(weight)					// Name  from weight
	FeName=tokens(absorb)						// Name  from Fe

	varIdx=_aandb(dataNames,userNames)			// Match btw (index)
	wgtIdx=_aandb(dataNames,weightName)
	FeIdx=_aandb(dataNames,FeName)

	varPos=selectindex(varIdx)					// Positions
	XPos=varPos[.,2..K+1]
	Ypos=varPos[.,1]
	Wpos=selectindex(wgtIdx)
	FePos=selectindex(FeIdx)

	// Vectors that contains the list of groups and their structures. Before the loop
	groupsAccP=locationOLS(1,K,cons)
	newmatrixCols=K+1+2

	// We can save time if the data is already sorted
	if (sort==0){
		while(itera<=blocksize){
			line=fget(fh)
			vecline  = strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline
			if (i==blocksize) {
				_sort(*data,FePos)
				info=&(panelsetup(*data,FePos))
				ngroup=rows(*info) //We have max ngroup new groups						
				Gtp=&(J(ngroup,1,.))

				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*data,l,*info))
					size=rows(*portion)
					(*Gtp)[l,1]=(*portion)[1,FePos]
					groups=locationOLS(1,K,cons) // new strucutre
					groupsAccP=(groupsAccP \ groups)
					
					yp=&((*portion)[.,Ypos])
					Xp=&((*portion)[.,XPos])
					wp=&((*portion)[.,Wpos])

					pos=l+1
					groupsAccP[pos,1]=&updateCLS(Xp,yp,wp,size,cons,*(groupsAccP[pos,1]))
				}
				i=0
				(*data) = J(blocksize,newmatrixCols,.)
			}
			++i
			++itera
		}
		i=1
		groupsAccP=groupsAccP[2..ngroup+1,.]

		// Read in data 
		while((line=fget(fh)) != J(0,0,"")){
			vecline  = strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline
			if (i==blocksize) {
				_sort(*data,FePos)
				info=&(panelsetup(*data,FePos))
				ngroup=rows(*info) //We have max ngroup new groups

				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*data,l,*info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,FePos]) // group already exist (1 & 0 Vector)

					yp=&((*portion)[.,Ypos])
					Xp=&((*portion)[.,XPos])
					wp=&((*portion)[.,Wpos])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1 // position for the new added group
						Gtp=&(*Gtp \ (*portion)[1,FePos])

						groups=locationOLS(1,K,cons) // new strucutre
						groupsAccP=(groupsAccP \ groups)					
						groupsAccP[pos,1]=&updateCLS(Xp,yp,wp,size,cons,*(groupsAccP[pos,1]))
					}
					else{
						index=selectindex(add) // position of the group that existed (more info: help mata select)
						groupsAccP[index,1]=&updateCLS(Xp,yp,wp,size,cons,*(groupsAccP[index,1]))
					}
				}
				i=0
				(*data) = J(blocksize,newmatrixCols,.)
			}
			++i
		}

		if	((*data)[.,1]!=J(blocksize,1,.)){
			i=i-1
			dataP=&((*data)[1..i,.])
			_sort(*dataP,FePos)
			info=panelsetup(*dataP,FePos)
			ngroup=rows(info)

			for(l=1; l<=ngroup; l++) {
				portion=&(panelsubmatrix(*dataP,l,info))
				size=rows(*portion)
				add=vec_inlistOLS(Gtp,(*portion)[1,FePos])
				
				y=&((*portion)[.,Ypos])
				X=&((*portion)[.,XPos])
				w=&((*portion)[.,Wpos])

				if (add==J(rows(add),1,0)){
					pos=rows(*Gtp)+1
					Gtp=&(*Gtp \ (*portion)[1,FePos])

					groups=locationOLS(1,K,cons)
					groupsAccP=(groupsAccP \ groups)
					groupsAccP[pos,1]=&updateCLS(X,y,w,size,cons,*(groupsAccP[pos,1]))
				}
				else{
					index=selectindex(add)
					groupsAccP[index,1]=&updateCLS(X,y,w,size,cons,*(groupsAccP[index,1]))
				}
			}				
		};
	}
	else {
		while(itera<=blocksize){
			line=fget(fh)
			vecline  = strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline
			
			if (i==blocksize) {
				info=&(panelsetup(*data,1))
				ngroup=rows(*info) //We have max ngroup new groups						
				Gtp=&(J(ngroup,1,.))

				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*data,l,*info))
					size=rows(*portion)
					(*Gtp)[l,1]=(*portion)[1,1]
					groups=locationOLS(1,K,cons) // new strucutre
					groupsAccP=(groupsAccP \ groups)
					
					yp = &((*portion)[.,2])
					Xp = &((*portion)[.,3..XFnlPos])
					wp = &((*portion)[.,WFnlPos])
					pos=l+1
					groupsAccP[pos,1]=&updateCLS(Xp,yp,size,cons,*(groupsAccP[pos,1]),wp)
				}
				i=0
				(*data) = J(blocksize,newmatrixCols,.)
			}
			++i
			++itera
		}
		i=1
		groupsAccP=groupsAccP[2..ngroup+1,.]

		// Read in data 
		while((line=fget(fh)) != J(0,0,"")){
			vecline  = strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline
			if (i==blocksize) {
				info=&(panelsetup(*data,1))
				ngroup=rows(*info) //We have max ngroup new groups

				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*data,l,*info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,1]) // group already exist (1 & 0 Vector)
					yp = &((*portion)[.,2])
					Xp = &((*portion)[.,3..XFnlPos])
					wp = &((*portion)[.,WFnlPos])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1 // position for the new added group
						Gtp=&(*Gtp \ (*portion)[1,1])

						groups=locationOLS(1,K,cons) // new strucutre
						groupsAccP=(groupsAccP \ groups)					
						groupsAccP[pos,1]=&updateCLS(Xp,yp,size,cons,*(groupsAccP[pos,1]),wp)
					}
					else{
						index=selectindex(add) // position of the group that existed (more info: help mata select)
						groupsAccP[index,1]=&updateCLS(Xp,yp,size,cons,*(groupsAccP[index,1]),wp)
					}
				}
				i=0
				(*data) = J(blocksize,newmatrixCols,.)
			}
			++i
		}

		if	((*data)[.,1]!=J(blocksize,1,.)){
			i=i-1
			dataP=&((*data)[1..i,.])
			info=panelsetup(*dataP,1)
			ngroup=rows(info)

			for(l=1; l<=ngroup; l++) {
				portion=&(panelsubmatrix(*dataP,l,info))
				size=rows(*portion)
				add=vec_inlistOLS(Gtp,(*portion)[1,1])
				
				y = &((*portion)[.,2])
				X = &((*portion)[.,3..XFnlPos])
				w = &((*portion)[.,WFnlPos])

				if (add==J(rows(add),1,0)){
					pos=rows(*Gtp)+1
					Gtp=&(*Gtp \ (*portion)[1,1])

					groups=locationOLS(1,K,cons)
					groupsAccP=(groupsAccP \ groups)
					groupsAccP[pos,1]=&updateCLS(X,y,size,cons,*(groupsAccP[pos,1]),w)
				}
				else{
					index=selectindex(add)
					groupsAccP[index,1]=&updateCLS(X,y,size,cons,*(groupsAccP[index,1]),w)
				}
			}				
		};
	}

	// Close file, return struct j
	fclose(fh)
	return(groupsAccP)
}
end