#include "foo_mrmr.h"


SEXP mrnet_adapted2(SEXP Rdata, SEXP Rnamat, SEXP Rprior, SEXP Rprior_weight, SEXP Rmaxparents, SEXP Rnvar, SEXP Rnsample, SEXP Rpredn, SEXP Rnpredn, SEXP Rthreshold){     
	double *data, *prior, *prior_weight;
	double *rel, *red, *res, *mim, score=1,*threshold, *tmp, *res_final;
	
	const int *maxparents, *nvar, *nsample, *npredn;
	int *predn, *ind,*namat;
	
	
	unsigned int n, jmax=0;
	
	PROTECT(Rdata = AS_NUMERIC(Rdata));
	PROTECT(Rnamat = AS_INTEGER(Rnamat));
	PROTECT(Rprior = AS_NUMERIC(Rprior));
	PROTECT(Rprior_weight = AS_NUMERIC(Rprior_weight));
	PROTECT(Rmaxparents= AS_INTEGER(Rmaxparents));
	PROTECT(Rnvar= AS_INTEGER(Rnvar));
	PROTECT(Rnsample= AS_INTEGER(Rnsample));
	PROTECT(Rpredn= AS_INTEGER(Rpredn));
	PROTECT(Rnpredn= AS_INTEGER(Rnpredn));
	PROTECT(Rthreshold= AS_NUMERIC(Rthreshold));
	
	
	data = NUMERIC_POINTER(Rdata);
	namat=INTEGER_POINTER(Rnamat);
	prior = NUMERIC_POINTER(Rprior);
	prior_weight = NUMERIC_POINTER(Rprior_weight);
	nvar = INTEGER_POINTER(Rnvar);
	nsample = INTEGER_POINTER(Rnsample);
	predn = INTEGER_POINTER(Rpredn);
	npredn = INTEGER_POINTER(Rnpredn);
	threshold = NUMERIC_POINTER(Rthreshold);
	maxparents = INTEGER_POINTER(Rmaxparents);
	n = *nvar;

	//new variables
	SEXP Rmim, Rres,Rred,Rrel,Rind,Rtmp,Rres_final;


	PROTECT(Rmim = NEW_NUMERIC(n*n));
	PROTECT(Rres = NEW_NUMERIC(n*n));
	PROTECT(Rres_final = NEW_NUMERIC(n*n));
	PROTECT(Rrel = NEW_NUMERIC(n));
	PROTECT(Rred = NEW_NUMERIC(n)); 
	PROTECT(Rind = NEW_INTEGER(*nsample));
	PROTECT(Rtmp = NEW_NUMERIC(n)); 
	
	tmp= NUMERIC_POINTER(Rtmp);	
	ind = INTEGER_POINTER(Rind);	
	mim = NUMERIC_POINTER(Rmim);
	res = NUMERIC_POINTER(Rres);
	rel = NUMERIC_POINTER(Rrel);
	red = NUMERIC_POINTER(Rred);
	res_final = NUMERIC_POINTER(Rres_final);
	
/*	if(*maxparents==*nvar){
		std::cout<<"warning: the maximum number of parents can be maximally "<<*nvar-1<<std::endl;
	}
 */
	for(unsigned int i=0;i < *nsample; ++i){
		ind[i]=i;
	}
		
	build_mim_subset(mim, data, namat, n, *nsample, ind, *nsample);

	
/*	for( unsigned int i=0; i< n; ++i ){
		for( unsigned int j=0; j<n; ++j ){ 
			mim[i*n+j] = (1- (*prior_weight))*mim[i*n+j] + (*prior_weight) * prior[i*n+j];
		}
	}
*/	
	
	// initialize res and res_final to threshold (usually something like -1000)
	for( unsigned int i=0; i< n; ++i ){
		for( unsigned int j=0; j<n; ++j ){ 
			res[i*n+j]=*threshold;
			res_final[i*n+j]=*threshold;
		}
	}
	
	// outer loop: every variable is target exactly once
	for(unsigned int i=0; i<n; ++i ) {
		for( unsigned int j=0; j<n; ++j ) {
			//first variable to be selected will be that with highest mutual information
			// init rel and red and select first
			rel[j]= mim[i*n+j];
			red[j]=0;
			if( rel[j] > rel[jmax])
				jmax = j;
		}
				
		score = rel[jmax];
		if( res[i*n+jmax] < score ) {
			res[jmax*n+i] = score;
			res[i*n+jmax] = score;
		}
		
		rel[jmax]=-1000; 
		
		// redundancy for each variable with the previously selected variable is increased by their MI
		for(unsigned int l=0; l<n; ++l )  
			red[l] += mim[l*n+jmax];
				
				//select others
		for(unsigned int k=1; k < n-1; k++ ) { 
			jmax = 0;
			for(unsigned int j=1; j < n; ++j ) {   
				if( (rel[j] - red[j]/k) > (rel[jmax] - red[jmax]/k) ) 
					jmax = j;
			}      
			score = (rel[jmax] - (red[jmax]/k));
			if( res[i*n+jmax] < score ) {
				res[i*n+jmax] = score;
			}
								
					//update rel and red
			rel[jmax]=-1000; 
			for( int l=0; l<n; ++l )  
				red[l] += mim[l*n+jmax];
							
					// stop criterion
		//          if( score < 0 ) k=n;
				if( score < *threshold ) k=n;
		}
	}

		// force diagonal to zero
	for( unsigned int i=0; i< n; ++i ) {
		for( unsigned int j=0; j<n; ++j ){
			if(j==i){
				res[i*n+j]=-1000;//diagonal
			}
		}
	}



	double maxmrmr=*threshold;
	// force symmetry
	for( unsigned int i=0; i< n; ++i ) {
		for( unsigned int j=i+1; j<n; ++j ){
			if(res[i*n+j]>res[j*n+i]){
				res[j*n+i]=res[i*n+j];
			}else{
				res[i*n+j]=res[j*n+i];
			}
			if(abs(res[j*n+i])>maxmrmr){
				maxmrmr=abs(res[j*n+i]);
			}
		}
	}

	//normalize with maximum absolute value
	for( unsigned int i=0; i< n; ++i ) {
		for( unsigned int j=i+1; j<n; ++j ){
			res[i*n+j]=res[i*n+j]/maxmrmr;
			res[j*n+i]=res[j*n+i]/maxmrmr;
		}
	}
	cout << "normalized mrmr matrix"<<endl;
	for( unsigned int i=0; i< n; ++i ) {
		cout << res[i]<< " ";
	}
	cout<<endl;
	for( unsigned int i=0; i< n; ++i ){
		for( unsigned int j=0; j<n; ++j ){ 
			res[i*n+j] = (1- (*prior_weight))*res[i*n+j] + (*prior_weight) * prior[i*n+j];
		}
	}
	// force symmetry after prior integration
/*	for( unsigned int i=0; i< n; ++i ) {
		for( unsigned int j=i+1; j<n; ++j ){
			if(res[i*n+j]>res[j*n+i]){
				res[j*n+i]=res[i*n+j];
			}else{
				res[i*n+j]=res[j*n+i];
			}
		}
	}*/
	
	//i corresponds to column; j to row in mrmr matrix; for columns which actually should be predicted
	for(unsigned int j=0; j< *npredn; ++j ) {
		for(unsigned int i=0; i< n; ++i ) {
			tmp[i]= res[i+(predn[j]-1)*n];
		}
	
		sort(tmp,tmp+n);
		for(unsigned int i=0; i< n; ++i ) {
			if(res[i+(predn[j]-1)*n]>tmp[n-(*maxparents) - 1]){
				res_final[i+(predn[j]-1)*n]=res[i+(predn[j]-1)*n];
			}
		}
	}

	UNPROTECT(17);

	return Rres_final;
}

