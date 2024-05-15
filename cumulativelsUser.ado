cap program drop cumulativels 

program cumulativels, eclass
version 15.0

*-------------------------------------------------------------------------------
*--- (0) Syntax diagram
*-------------------------------------------------------------------------------
#delimit;
syntax anything,
    filename(string)
    BLOCKsize(integer)
    [
	noCONStant
	sort
	ABSorb(string)
	Robust
	CLuster(string)
	]
;
#delimit cr

*-------------------------------------------------------------------------------
*--- (1) Unpack syntax
*-------------------------------------------------------------------------------

if "`sort'"=="sort" {
	local sort = 1
}
else local sort=0

if length(`"`cluster'"')!=0 {
	local clMarker = 1
}
else local clMarker=0

if ustrregexm("`anything'", "\((.*?)\)"){
    local match "`=ustrregexs(1)'"
	local matchC : subinstr local match "=" " "
	local typeWeight = ustrregexm("fweight","`: word 1 of `matchC''")
	local weight: word 2 of `matchC'
	local wMarker = 1
	local rError = ustrregexm("pweight","`: word 1 of `matchC''")
}
else{
	local typeWeight = 0
	local wMarker = 0
	local rError = 0
}

if length(`"`robust'"')!=0 {
	local rError = 1
}

local namelist = subinstr("`anything'", "(`match')", "",.)
loc K1 : word count `namelist'
loc K = `K1'-1
loc arglist `namelist'
loc Ylist: word 1 of `namelist'
loc Xlist: list namelist - Ylist

if "`constant'"=="noconstant" {
	local nocons = 1
}
else{
    local nocons = 0
	local Xlist  `Xlist' "_cons"
} 
*-------------------------------------------------------------------------------
*--- (2) Run runcls and report output
*-------------------------------------------------------------------------------
// NP: This definitely needs some tidying up; there are too many parentheses and extra steps
#delimit ;
mata: A=runcls("`filename'","`arglist'","`weight'","`absorb'","`cluster'",`blocksize',`K',`nocons',`sort',`typeWeight',`wMarker',`clMarker',`rError');
#delimit cr

mata: betas=(*((*(A[1]))[1]))'
mata: vcov=*((*(A[1]))[2])
mata: N=*(A[2])

mata: st_matrix("b", betas)
mata: st_matrix("V", vcov)
mata: st_local("N", strofreal(N))

matrix colnames b = `Xlist'
matrix colnames V = `Xlist'
matrix rownames V = `Xlist'

ereturn post b V, dep(`Ylist') obs(`N')
ereturn display
end
*-------------------------------------------------------------------------------
*--- (3) Define auxiliary Mata details
*-------------------------------------------------------------------------------
version 15.0

mata:
struct resultsCLS 
{
	//---------------------- derived
	real scalar yy,size,Ybar
	real colvector Xy,Xgug
	real matrix XX,Xbar
}

//------------------------------------------------------------------------------
//--- (4) Main wrapper function
//------------------------------------------------------------------------------
function runcls(string scalar filename,	// File containing data
				string scalar arglist,	// String containing var names
				string scalar weight,	// Weight's name variable; Optional
				string scalar absorb,	// Fe's name variable;	 Optional
				string scalar cluster,	// Fe's name variable;	 Optional
				real scalar blocksize,	// Size in observations for blocks
				real scalar K,			// Number of dep vars
				real scalar nocons,
				real scalar sort,
				real scalar typeWeight,
				real scalar wMarker,				
				real scalar clMarker,
				real scalar rError)
{
	cons=1-nocons
	covs=K+cons // Number of actual covs when constant is defined

	// General stucture
	struct resultsCLS scalar r
	r.yy=r.size	= 0
	r.Xy		= J(covs,1,0)
	r.XX		= J(covs,covs,0)

	txtPos=inputs(filename,arglist,weight,absorb,cluster,K)

	if (absorb==""){
		// matrix calculations
		r=ols(filename,txtPos,K,cons,blocksize,wMarker,typeWeight,r)
	
		// compute results
		/* I was thinking we can create a new subrutine for computing beta and sd*/
		
		XXinv=&(invsym(r.XX))
		beta=&(*XXinv*r.Xy)
		
		// non-redundant degrees of freedom?
		KStd=&(colsum((*beta:!=0)))

		// standard error type

		if (rError==1|clMarker==1){
			struct resultsCLS scalar r2
			if (rError==1){
				r2.XX= J(covs,covs,0)
				r2=olsR(filename,txtPos,beta,K,cons,blocksize,wMarker,typeWeight,r2)
				uPuR=&(*XXinv*r2.XX**XXinv)
				vcov=&(((r.size)/(r.size-*KStd))**uPuR)
			}
			else {
				sdCluster=olsCl(filename,txtPos,beta,K,cons,blocksize,wMarker,typeWeight,sort)
				sdCluster
				Ng=length(sdCluster)-1 	 // N° groups

				r2=*(sdCluster[1,1])
				XgugT=&(transposeonly(r2.Xgug))
					
				omegaG= &(quadcross(*XgugT,*XgugT))
				omegaT=&((*omegaG)[1..covs,1..covs])

				for (j=2; j<=Ng; j++) {
					r2=*(sdCluster[j,1])
					i=j-1
					XgugT=&(transposeonly(r2.Xgug))
					omegaG= &(quadcross(*XgugT,*XgugT))
					omegaT=&(*OmegaT + (*Omegas)[1..covs,1..covs])
				}
				vcov = (N-1)/(N-*KStd)*(Ng/(Ng-1))**XXinv**omegaT**XXinv
			}
		}
		else{
			// 'U = (yy - 2BXy) +  B*X'X*B
			uPu= &(r.yy - 2*(*beta)'*r.Xy + (*beta)'*r.XX**beta)
			vcov=&(((*uPu)/(r.size-*KStd))**XXinv)
		}
	}
	outputMain=&(beta,vcov)
	outputArg=&(r.size)
	output=(outputMain,outputArg)
	return(output)
}

//------------------------------------------------------------------------------
//--- (X) Functions called from main
//------------------------------------------------------------------------------

struct resultsCLS matrix inputs(string scalar filename,
								string scalar arglist,
								string scalar weight,
								string scalar absorb,
								string scalar cluster,
								real scalar K)
{
	fh = fopen(filename,"r")
	names = fget(fh)
	filePos=&(fh)

	// Match text and user variables positions
	dataNames=tokens(subinstr(names,",", " "))	// Names from textfile
	userNames=tokens(arglist)			// Names from user
	varsPos=inlist_aposb(dataNames,userNames) // Positions

	/* 
	The folowing lines of code causes a slowdown in performance, since we got
	the positions of covs. For instance XPos might be (1,7,6,9,6). That vector enters
	as subscript in the data. Because is not a sequence there will be a slowdownx	
	*/

	Ypos=&(varsPos[.,1])
	XPos=&(varsPos[.,2..K+1])
	positions=(filePos\ Ypos \ XPos)

	if (absorb!=""){
		FeName=tokens(absorb)			// Name from Fe
		FeIdx=_aandb(dataNames,FeName)
		FePos=&selectindex(FeIdx)
		positions=(positions \ FePos)
	}
	if (weight!=""){
		weightName=tokens(weight)		// Name from weight
		wgtIdx=_aandb(dataNames,weightName)
		Wpos=&selectindex(wgtIdx)
		if (cols(*Wpos)==1){
			positions=(positions \ Wpos)
		}
		else{
            displayas("error")
            printf("%s not found\n", weightName)
            exit(601)			
		}
	}
	if (cluster!=""){
		clusterName=tokens(cluster)		// Name from cluster
		clusIdx=_aandb(dataNames,clusterName)
		Cltpos=&selectindex(clusIdx)
		if (cols(*Cltpos)==1){
			positions=(positions \ Cltpos)
		}
		else{
            displayas("error")
            printf("%s not found\n", clusterName)
            exit(601)			
		}
	}

	return(positions)
}
struct resultsCLS matrix ols(string scalar filename,
							 pointer position,
							 real scalar K,
							 real scalar cons,
							 real scalar blocksize,
							 real scalar wMarker,
							 real scalar typeWeight,
							 struct resultsCLS scalar j)
{
	// General information from the model
	real scalar	i
	i=1
	fh=*(position[1])

	//Define data
	Ypos=*(position[2])
	Xpos=*(position[3])
	Xlast=K+1
	newMtxCols=Xlast+wMarker
	data=&(J(blocksize,newMtxCols,.))

	// Read in data
	if(wMarker==0){
	    arPos=(Ypos,Xpos)
		while((line=fget(fh)) != J(0,0,"")) {
			vecline	= strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline[arPos]
			if (i==blocksize) {
				y = &((*data)[.,1])
				X = &((*data)[.,2..Xlast])
				
				j=CLS(X,y,cons,i,j)
				i=0
				(*data)=J(blocksize,newMtxCols,.)
			}
			++i
		}
		//This captures any overflow after last block (matrices aren't reset)
		if	((*data)[.,1]!=J(blocksize,1,.)){
			i=i-1
			datap=&((*data)[1..i,.])

			yp = &((*data)[.,1])
			Xp = &((*data)[.,2..Xlast])
			j=CLS(Xp,yp,cons,i,j)
		}
	}
	else{
		Wpos=*(position[4])
		arPos=(Ypos,Xpos,Wpos)
		Wlast=Xlast+1
		if(typeWeight==0){
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[arPos]

				if (i==blocksize) {
					y = &((*data)[.,1])
					X = &((*data)[.,2..Xlast])
					w = &((*data)[.,Wlast])

					j=CLS_W(X,y,w,cons,i,j)
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
			}
			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				datap=&((*data)[1..i,.])

				yp = &((*data)[.,1])
				Xp = & ((*data)[.,2..Xlast])
				wp = &((*data)[.,Wlast])

				j=CLS_W(Xp,yp,wp,cons,i,j)
			}
		}
		else{
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[arPos]

				if (i==blocksize) {
					y = &((*data)[.,1])
					X = &((*data)[.,2..Xlast])
					w = &((*data)[.,Wlast])
					iw=quadcolsum(*w)

					j=CLS_W(X,y,w,cons,iw,j)
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
			}
			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				datap=&((*data)[1..i,.])

				yp = &((*data)[.,1])
				Xp = & ((*data)[.,2..Xlast])
				wp = &((*data)[.,Wlast])

				iw=quadcolsum(*wp)
				j=CLS_W(Xp,yp,wp,cons,iw,j)
			}
		}
	}

// Close file, return struct j
fclose(fh)
return(j)
}
struct resultsCLS matrix olsR(string scalar filename,
				pointer position,
				pointer beta,
				real scalar K,
				real scalar cons,
				real scalar blocksize,
				real scalar wMarker,
				real scalar typeWeight,
				struct resultsCLS scalar j)
{
	// General information from the model
	real scalar	i
	i=1
	fh = fopen(filename,"r")
	names = fget(fh)

	//Define data
	Ypos=*(position[2])
	Xpos=*(position[3])
	Xlast=K+1
	newMtxCols=Xlast+wMarker
	data=&(J(blocksize,newMtxCols,.))

	// Read in data
	if(wMarker==0){
	    argPos=(Ypos,Xpos)
		while((line=fget(fh)) != J(0,0,"")) {
			vecline	= strtoreal(tokens(subinstr(line,",", " ")))
			(*data)[i,.]=vecline[argPos]
			if (i==blocksize) {
				y = &((*data)[.,1])
				X = &((*data)[.,2..Xlast])
				
				j=CLS_R(X,y,beta,blocksize,cons,j)
				i=0
				(*data)=J(blocksize,newMtxCols,.)
			}
			++i
		}
		//This captures any overflow after last block (matrices aren't reset)
		if	((*data)[.,1]!=J(blocksize,1,.)){
			i=i-1
			datap=&((*data)[1..i,.])

			yp = &((*data)[.,1])
			Xp = &((*data)[.,2..Xlast])
			j=CLS_R(Xp,yp,beta,blocksize,cons,j)
		}
	}
	else{
		Wpos=*(position[4])
		argPos=(Ypos,Xpos,Wpos)
		Wlast=Xlast+1
		if(typeWeight==0){
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]
				
				if (i==blocksize) {
					y = &((*data)[.,1])
					X = &((*data)[.,2..Xlast])
					w = &((*data)[.,Wlast])

					j=CLS_RAw(X,y,w,beta,blocksize,cons,j)
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
			}

			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				datap=&((*data)[1..i,.])

				yp = &((*data)[.,1])
				Xp = &((*data)[.,2..Xlast])
				wp = &((*data)[.,Wlast])

				j=CLS_RAw(Xp,yp,wp,beta,blocksize,cons,j)
			}
		}
		else{
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]
				
				if (i==blocksize) {
					y = &((*data)[.,1])
					X = &((*data)[.,2..Xlast])
					w = &((*data)[.,Wlast])

					j=CLS_RFw(X,y,w,beta,blocksize,cons,j)
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
			}

			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				datap=&((*data)[1..i,.])

				yp = &((*data)[.,1])
				Xp = &((*data)[.,2..Xlast])
				wp = &((*data)[.,Wlast])

				j=CLS_RFw(Xp,yp,wp,beta,blocksize,cons,j)
			}			
		}
	}

// Close file, return struct j
fclose(fh)
return(j)
}
struct resultsCLS matrix olsCl(string scalar filename,
							   pointer position,
							   pointer beta,
							   real scalar K,
							   real scalar cons,
							   real scalar blocksize,
							   real scalar wMarker,
							   real scalar typeWeight,
							   real scalar sort)
{
	// General information from the model
	real scalar	i,itera
	i=itera=1
	fh = fopen(filename,"r")
	names = fget(fh)

	//Define data
	Ypos=*(position[2])
	Xpos=*(position[3])
	Xlast=K+2
	newMtxCols=Xlast+wMarker
	data=&(J(blocksize,newMtxCols,.))
	clustersAccP=matStructCl(1,K) // Vectors containing the list of groups and their structures, prior to the loop

	// Read in data
	if (sort==0){
		if(wMarker==0){
			Clpos=*(position[4])
			argPos=(Clpos,Ypos,Xpos)
			while(itera<=blocksize){
				line=fget(fh)

				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]
				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) // we have max ngroup new groups						
					Gtp=&(J(ngroup,1,.))

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						(*Gtp)[l,1]=(*portion)[1,1]
						groups=matStructCl(1,K) // new strucutre
						clustersAccP=(clustersAccP \ groups)

						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						pos=l+1
						clustersAccP[pos,1]=&CLS_Cl(X,y,beta,size,cons,*(clustersAccP[pos,1]))
					}
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
				++itera
			}
			i=1
			clustersAccP=clustersAccP[2..ngroup+1,.]

			// Read in data 
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) //We have max ngroup new groups

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						add=vec_inlistOLS(Gtp,(*portion)[1,1]) // group already exist (1 & 0 Vector)

						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])

						if (add==J(rows(add),1,0)){
							pos=rows(*Gtp)+1 // position for the new added group
							pos
							Gtp=&(*Gtp \ (*portion)[1,1])

							groups=matStructCl(1,K) // new strucutre
							clustersAccP=(clustersAccP \ groups)
							clustersAccP[pos,1]=&CLS_Cl(X,y,beta,cons,size,*(clustersAccP[pos,1]))
						}
						else{
							index=selectindex(add) // position of the group that existed (more info: help mata select)
							groupsAccP[index,1]=&CLS_Cl(X,y,beta,cons,size,*(clustersAccP[index,1]))
						}
					}
					i=0
					i
					(*data) = J(blocksize,newmatrixCols,.)
				}
				++i
			}

			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				dataP=&((*data)[1..i,.])
				_sort(*dataP,1)
				info=panelsetup(*dataP,1)
				ngroup=rows(info)
				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*dataP,l,info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,1])

					yp = &((*portion)[.,1])
					Xp = &((*portion)[.,2..Xlast])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1
						Gtp=&(*Gtp \ (*portion)[1,1])

						groups=matStructCl(1,K)
						clustersAccP=(groupsAccP \ groups)
						clustersAccP[pos,1]=&CLS_Cl(Xp,yp,beta,size,cons,*(clustersAccP[pos,1]))
					}
					else{
						index=selectindex(add)
						clustersAccP[index,1]=&CLS_Cl(Xp,yp,beta,size,cons,*(clustersAccP[index,1]))
					}
				}
			}
		}
		else{
			Wpos=*(position[4])
			argPos=(Clpos,Ypos,Xpos,Wpos)
			Wlast=Xlast+1
			while(itera<=blocksize){
				line=fget(fh)

				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) // we have max ngroup new groups						
					Gtp=&(J(ngroup,1,.))

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						(*Gtp)[l,1]=(*portion)[1,1]
						groups=matStructCl(1,K) // new strucutre
						clustersAccP=(clustersAccP \ groups)

						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						w = &((*portion)[.,Wlast])

						pos=l+1
						clustersAccP[pos,1]=&CLS_ClW(X,y,w,beta,size,cons,*(clustersAccP[pos,1]))
					}
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
				++itera
			}
			i=1
			clustersAccP=clustersAccP[2..ngroup+1,.]

			// Read in data 
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) //We have max ngroup new groups

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						add=vec_inlistOLS(Gtp,(*portion)[1,1]) // group already exist (1 & 0 Vector)
						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						w = &((*portion)[.,Wlast])

						if (add==J(rows(add),1,0)){
							pos=rows(*Gtp)+1 // position for the new added group
							Gtp=&(*Gtp \ (*portion)[1,1])

							groups=matStructCl(1,K) // new strucutre
							clustersAccP=(clustersAccP \ groups)					
							clustersAccP[pos,1]=&CLS_ClW(Xp,yp,beta,cons,size,*(clustersAccP[pos,1]),w)
						}
						else{
							index=selectindex(add) // position of the group that existed (more info: help mata select)
							clustersAccP[index,1]=&CLS_ClW(Xp,yp,beta,cons,size,*(clustersAccP[index,1]),w)
						}
					}
					i=0
					(*data) = J(blocksize,newmatrixCols,.)
				}
				++i
			}

			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				dataP=&((*data)[1..i,.])
				_sort(*dataP,1)
				info=panelsetup(*dataP,1)
				ngroup=rows(info)
				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*dataP,l,info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,1])

					yp = &((*portion)[.,2])
					Xp = &((*portion)[.,3..Xlast])
					wp = &((*portion)[.,Wlast])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1
						Gtp=&(*Gtp \ (*portion)[1,1])

						groups=matStructCl(1,K)
						clustersAccP=(clustersAccP \ groups)
						clustersAccP[pos,1]=&CLS_ClW(Xp,yp,beta,size,cons,*(clustersAccP[pos,1]),wp)
					}
					else{
						index=selectindex(add)
						clustersAccP[index,1]=&CLS_ClW(Xp,yp,beta,size,cons,*(clustersAccP[index,1]),wp)
					}
				}
			}
		}
	}
	else{
		if(typeWeight==0){
			argPos=(Clpos,Ypos,Xpos)
			while(itera<=blocksize){
				line=fget(fh)

				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) // we have max ngroup new groups						
					Gtp=&(J(ngroup,1,.))

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						(*Gtp)[l,1]=(*portion)[1,1]
						groups=matStructCl(1,K) // new strucutre
						clustersAccP=(clustersAccP \ groups)

						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						pos=l+1
						clustersAccP[pos,1]=&CLS_Cl(X,y,beta,size,cons,*(clustersAccP[pos,1]))
					}
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
				++itera
			}
			i=1
			clustersAccP=clustersAccP[2..ngroup+1,.]

			// Read in data 
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) //We have max ngroup new groups

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						add=vec_inlistOLS(Gtp,(*portion)[1,1]) // group already exist (1 & 0 Vector)
						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])

						if (add==J(rows(add),1,0)){
							pos=rows(*Gtp)+1 // position for the new added group
							Gtp=&(*Gtp \ (*portion)[1,1])

							groups=matStructCl(1,K) // new strucutre
							clustersAccP=(clustersAccP \ groups)					
							clustersAccP[pos,1]=&CLS_Cl(Xp,yp,beta,cons,size,*(groupsAccP[pos,1]))
						}
						else{
							index=selectindex(add) // position of the group that existed (more info: help mata select)
							groupsAccP[index,1]=&CLS_Cl(Xp,yp,beta,cons,size,*(groupsAccP[index,1]))
						}
					}
					i=0
					(*data) = J(blocksize,newmatrixCols,.)
				}
				++i
			}

			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				dataP=&((*data)[1..i,.])
				_sort(*dataP,1)
				info=panelsetup(*dataP,1)
				ngroup=rows(info)
				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*dataP,l,info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,1])

					yp = &((*portion)[.,1])
					Xp = &((*portion)[.,2..Xlast])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1
						Gtp=&(*Gtp \ (*portion)[1,1])

						groups=matStructCl(1,K)
						clustersAccP=(groupsAccP \ groups)
						clustersAccP[pos,1]=&CLS_Cl(Xp,yp,beta,size,cons,*(clustersAccP[pos,1]))
					}
					else{
						index=selectindex(add)
						clustersAccP[index,1]=&CLS_Cl(Xp,yp,beta,size,cons,*(clustersAccP[index,1]))
					}
				}
			}
		}
		else{
			Wpos=*(position[4])
			argPos=(Clpos,Ypos,Xpos,Wpos)
			Wlast=Xlast+1
			while(itera<=blocksize){
				line=fget(fh)

				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) // we have max ngroup new groups						
					Gtp=&(J(ngroup,1,.))

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						(*Gtp)[l,1]=(*portion)[1,1]
						groups=matStructCl(1,K) // new strucutre
						clustersAccP=(clustersAccP \ groups)

						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						w = &((*portion)[.,Wlast])

						pos=l+1
						clustersAccP[pos,1]=&CLS_ClW(X,y,w,beta,size,cons,*(clustersAccP[pos,1]))
					}
					i=0
					(*data)=J(blocksize,newMtxCols,.)
				}
				++i
				++itera
			}
			i=1
			clustersAccP=clustersAccP[2..ngroup+1,.]

			// Read in data 
			while((line=fget(fh)) != J(0,0,"")) {
				vecline	= strtoreal(tokens(subinstr(line,",", " ")))
				(*data)[i,.]=vecline[argPos]

				if (i==blocksize) {
					_sort(*data,1)
					info=&(panelsetup(*data,1))
					ngroup=rows(*info) //We have max ngroup new groups

					for(l=1; l<=ngroup; l++) {
						portion=&(panelsubmatrix(*data,l,*info))
						size=rows(*portion)
						add=vec_inlistOLS(Gtp,(*portion)[1,1]) // group already exist (1 & 0 Vector)
						y = &((*portion)[.,2])
						X = &((*portion)[.,3..Xlast])
						w = &((*portion)[.,Wlast])

						if (add==J(rows(add),1,0)){
							pos=rows(*Gtp)+1 // position for the new added group
							Gtp=&(*Gtp \ (*portion)[1,1])

							groups=matStructCl(1,K) // new strucutre
							clustersAccP=(clustersAccP \ groups)					
							clustersAccP[pos,1]=&CLS_ClW(Xp,yp,beta,cons,size,*(clustersAccP[pos,1]),w)
						}
						else{
							index=selectindex(add) // position of the group that existed (more info: help mata select)
							clustersAccP[index,1]=&CLS_ClW(Xp,yp,beta,cons,size,*(clustersAccP[index,1]),w)
						}
					}
					i=0
					(*data) = J(blocksize,newmatrixCols,.)
				}
				++i
			}
			//This captures any overflow after last block (matrices aren't reset)
			if	((*data)[.,1]!=J(blocksize,1,.)){
				i=i-1
				dataP=&((*data)[1..i,.])
				_sort(*dataP,1)
				info=panelsetup(*dataP,1)
				ngroup=rows(info)
				for(l=1; l<=ngroup; l++) {
					portion=&(panelsubmatrix(*dataP,l,info))
					size=rows(*portion)
					add=vec_inlistOLS(Gtp,(*portion)[1,1])

					yp = &((*portion)[.,2])
					Xp = &((*portion)[.,3..Xlast])
					wp = &((*portion)[.,Wlast])

					if (add==J(rows(add),1,0)){
						pos=rows(*Gtp)+1
						Gtp=&(*Gtp \ (*portion)[1,1])

						groups=matStructCl(1,K)
						clustersAccP=(clustersAccP \ groups)
						clustersAccP[pos,1]=&CLS_ClW(Xp,yp,beta,size,cons,*(clustersAccP[pos,1]),wp)
					}
					else{
						index=selectindex(add)
						clustersAccP[index,1]=&CLS_ClW(Xp,yp,beta,size,cons,*(clustersAccP[index,1]),wp)
					}
				}
			}
		}
	}
// Close file, return struct j
fclose(fh)
return(j)
}
struct resultsCLS scalar CLS(pointer X,
							 pointer y,
							 real scalar cons,
							 real scalar N,
							 struct resultsCLS scalar r)
{
	r.Xy  =r.Xy+quadcross(*X,cons,*y,0)
	r.XX  =r.XX+quadcross(*X,cons,*X,cons)
	r.yy  =r.yy+quadcross(*y,*y)
	r.size=r.size+N
	return(r)
}
struct resultsCLS scalar CLS_W(pointer X,
								pointer y,
								pointer w,
								real scalar cons,
								real scalar N,
								struct resultsCLS scalar r)
{
	r.Xy  =r.Xy+quadcross(*X,cons,*w,*y,0)
	r.XX  =r.XX+quadcross(*X,cons,*w,*X,cons)
	r.yy  =r.yy+quadcross(*y,*w,*y)
	r.size=r.size+N
	return(r)
}
struct resultsCLS scalar CLS_R(pointer X,
							   pointer y,
							   pointer beta,
							   real scalar blocksize,
							   real scalar cons,
							   struct resultsCLS scalar r)
{
	Xtemp =	&(*X,J(blocksize,cons,1))
	uR = &(*y - *Xtemp*(*beta))
	uR2= &(*uR:^2)
	r.XX = r.XX + quadcross(*Xtemp,*uR2,*Xtemp)
	return(r)
}
struct resultsCLS scalar CLS_RFw(pointer X,
							    pointer y,
							    pointer w,
								pointer beta,
								real scalar blocksize,
							    real scalar cons,
							    struct resultsCLS scalar r)
{
	Xtemp =	&(*X,J(blocksize,cons,1))
	XW = &(sqrt(*w):**Xtemp)
	uR = &(sqrt(*w):**y - *XW*(*beta))
	uR2= &(*uR:^2)
	r.XX = r.XX + quadcross(*Xtemp,*uR2,*Xtemp)
	return(r)
}
struct resultsCLS scalar CLS_RAw(pointer X,
							    pointer y,
							    pointer w,
								pointer beta,
								real scalar blocksize,
							    real scalar cons,
							    struct resultsCLS scalar r)
{
	Xtemp =	&(*X,J(blocksize,cons,1))
	XW = &(sqrt(*w):**Xtemp)
	uR = &(sqrt(*w):**y - *XW*(*beta))
	uR2= &(*uR:^2)
	r.XX = r.XX + quadcross(*XW,*uR2,*XW)
	return(r)
}
struct resultsCLS scalar CLS_Cl(pointer X,
								pointer y,
								pointer beta,
								real scalar blocksize,
								real scalar cons,
								struct resultsCLS scalar r)
{
	Xtemp =	&(*X,J(blocksize,cons,1))
	uR = &(*y - *Xtemp*(*beta))
	uR2= &(*uR:^2)
	r.XX = r.XX + quadcross(*Xtemp,*uR2,*Xtemp)
	r.size=r.size+blocksize
	return(r)
}
struct resultsCLS scalar CLS_ClW(pointer X,
								 pointer y,
								 pointer w,
								 pointer beta,
								 real scalar blocksize,
								 real scalar cons,
								 struct resultsCLS scalar r)
{
	Xtemp =	&(*X,J(blocksize,cons,1))
	XW = &(sqrt(*w):**Xtemp)
	uR = &(sqrt(*w):**y - *XW*(*beta))
	uR2= &(*uR:^2)
	r.XX = r.XX + quadcross(*XW,*uR2,*XW)
	r.size=r.size+blocksize
	return(r)
}
real rowvector inlist_aposb(rowvector B, rowvector A)
{
    real scalar                                     i
    real rowvector                                  R

    R = J(1,cols(A),0)

    for (i=1; i<=cols(A); i++) {
        pos= selectindex(_aposb(A[1,i],B))
        if (cols(pos)==1){
            R[1,i] = pos         
        }
        else{
            displayas("error")
            printf("variable %s not found\n", A[1,i])
            exit(601)
        }
    }
    return(R)
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
struct resultsCLS matrix matStructCl(scalar Tgroup, scalar K) 
{
	real scalar j
	pointer(struct resultsCLS scalar) matrix Index
	covs=K+1
	
	Index = &(resultsCLS())
	Index[1,1]->XX	 = J(covs,covs,0)
	Index[1,1]->size = 0
	if (Tgroup>1) {
		for (j=2; j<=Tgroup; j++) {
			Index = (Index \ &(resultsCLS()))
			Index[j,1]->XX   = J(covs,covs,0)
			Index[j,1]->size = 0
		}
	}
	return(Index)
}
end
