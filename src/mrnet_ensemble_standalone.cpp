#include "foo_mrmr.h"

double mrnet_onegene(double mim [], int size, int nbvar,int *var_ind,int var_target, int var_interest){
	// mim:			mutual information matrix
	// size:		total number of variables
	// nbvar:		number of previously selected variables (not the target)
	// var_ind:		the indices of the previously selected variables as vector
	// var_target:	the index of the target gene
	// var_interest: the variable for which the mrmr score with the target has to be computed; will be used for bootstrapping it

	unsigned int jmax;
	double rel, red,res;
	double max_val=-1000;

	jmax=var_target-1;
	//initialize the remaining entries to zero
	red=0;
	// the relevance for variable of interest with the target is simply its mutual information with it
	rel=mim[(var_target-1)*size+var_interest-1];
	if(nbvar > 0){
		// in case other variables have been previously selected; compute their redundancy with the variable of interest
		for(unsigned int j=0;j< nbvar; j++){
			red+=mim[(var_ind[j]-1)*size+var_interest-1];
		}		
		res=rel-red/nbvar ;
	}else{
		res=rel;
	}
	return res;
}


void bootstrap_mrmr(double &mean, double &sd, double data[],int size, int rep_boot, int size_boot,int nsamples, int var_target, int var_interest, int nprev_sel,int* var_ind)
{
	//mean
	//sd
	//data
	//size
	//rep_boot
	//size_boot
	//nsamples
	//var_target
	//var_interest
	//nprev_sel
	//var_ind

	int *ind;
	double *mim, *boot_val;
	
	ind=new int[size_boot];	
	mim = new double [size*size];
	boot_val=new double[rep_boot];
	
	
	for(unsigned int k=0; k< rep_boot; ++k){
	//in total there will be rep_boot times the mrmr sampled
		//determine the subset of samples that should be used (in total size_boot samples will be selected)
		for(unsigned int i=1;i<= size_boot;++i){
			ind[i-1]=rand () %nsamples;
		}
		// compute mi matrix for the subset
		build_mim_subset(mim, data, size, nsamples, ind, size_boot);
		boot_val[k]=mrnet_onegene( mim, size, nprev_sel, var_ind, var_target, var_interest);	
	}

	// determine mean and variance of bootstrapped values
	for(unsigned int i=0;i< rep_boot;++i){
		if(boot_val[i]==boot_val[i]){
			mean+=boot_val[i];
		}
	}
	mean=mean/rep_boot;

	for(unsigned int i=0;i< rep_boot;++i){
		if(boot_val[i]==boot_val[i]){
			sd+=(boot_val[i]-mean) * (boot_val[i]-mean);
		}
	}
	sd=sqrt(sd/rep_boot);
	
	delete [] ind;
	delete [] mim;
	delete [] boot_val;
}	

void mrmr_ensemble_one_gene (tree<int>& res, tree<int>::pre_order_iterator one, double data[], int nsamples,int n , int max_elements, int predn , int rep_boot){
	
	//predn:			index of target node

	int *nsub, *prev_sel; // nsub: the variables which have been previously selected + target; prev_sel=nsub-target
	
	int nsamples_boot=nsamples, nprev_sel=0;
	double *res_vec, *vec_mean, *vec_sd,  *vec_local_max_mean, *vec_local_max_sd;;
	
	double local_max_mean=0.0, local_max_sd=0.0;
	double max_mrmr;

	int cnt=0;
	int max_mrmr_ind;

	res_vec=new double [n];
	vec_mean=new double [n];
	vec_sd=new double [n];
	vec_local_max_mean=new double [max_elements];
	vec_local_max_sd=new double [max_elements];

	prev_sel=new int[max_elements];
	nsub = new int [max_elements];
	
	for(unsigned int j=0;j< max_elements;++j){
		vec_local_max_mean[j]=-1000;
	}
	
	
/*	tree<int> res;
	tree<int>::iterator top,one;*/
	
	tree<double> res_mean;
	tree<double>::iterator top_mean,one_mean;
	
	//initialize trees: res_mean for mean of the bootstrap
//	top=res.begin();
	top_mean=res_mean.begin();
	
	prev_sel[0]=0;
	
//has to be removed later
		//srand (54321);
////	
		
	//one=res.insert(top, predn);
	one_mean=res_mean.insert(top_mean, predn);
	prev_sel[0]=0;
	
	for(unsigned int k=0;k< n;++k){
		vec_mean[k]=0;
		vec_sd[k]=0;
	}
		
	//mrmr score should not be predicted to the target node
	vec_mean[predn-1]=-1000;
	vec_sd[predn-1]=-1000;
	
	for(unsigned int k=1;k<= n;++k){
		if(vec_mean[k-1]!= (-1000)){
			bootstrap_mrmr(vec_mean[k-1],vec_sd[k-1], data,n, rep_boot,nsamples_boot,nsamples,predn,k, nprev_sel,prev_sel);
		}
	}

	max_mrmr=-1000.0;
	max_mrmr_ind=-1;
	
	for(unsigned int j=0;j<n;++j){
		if(j!=(predn-1) && vec_mean[j]>max_mrmr){
			max_mrmr=vec_mean[j];
			max_mrmr_ind=j;
		}
	}
	if(max_mrmr>0){
		res.append_child(one,max_mrmr_ind+1);
		res_mean.append_child(one_mean,max_mrmr);
		for(unsigned int k=0;k< n;++k){
			//add equivalent nodes
			 if(k!=max_mrmr_ind && vec_mean[k]< vec_mean[max_mrmr_ind]+vec_sd[max_mrmr_ind] && vec_mean[k]> vec_mean[max_mrmr_ind]-vec_sd[max_mrmr_ind]){
				res.append_child(one,k+1);
				res_mean.append_child(one_mean,vec_mean[k]);
			 }
		}
	}
	
	//max mean and sd for this tree level
	vec_local_max_mean[cnt]=max_mrmr;
	vec_local_max_sd[cnt]=vec_sd[max_mrmr_ind];
	
	
	// in the first level (only relevance is computed), the first one is the maximum and only the equivalent nodes are added
	// starting from the second level, the different sub-trees might not be equivalent anymore -> has to be checked
	// -> has to be checked and lower one have to be removed
	max_mrmr=-1000;
	max_mrmr_ind=-1;
	int tmp_leaf;
	int max_elements_tmp=2;
	bool notincremented=true;
	
	while (max_elements_tmp<=max_elements) {
//has to be removed later
	//	srand (54321);
	////	
		//initialize different iterators for the two trees
		tree<int>::pre_order_iterator it=res.begin(), it_tmp=res.begin(), it_tmp2=res.begin();
		tree<double>::pre_order_iterator it_mean=res_mean.begin(), it_mean_tmp=res_mean.begin(), it_mean_tmp2=res_mean.begin();

		tree<int>::leaf_iterator li;
		tree<double>::leaf_iterator li_mean;
			
		if(!res.is_valid(it)) return;
	
		int rootdepth=res.depth(it);
		int cnt_lowlev=0;
		while(it!=res.end()) {
			if(!notincremented or cnt==0){	
				li=res.begin_leaf(it);
				li_mean=res_mean.begin_leaf(it_mean);
				while(li_mean!=res_mean.end_leaf(it_mean_tmp)) {
					//delete leafs for which the local mrmr score is not in the confidence interval
					if((*li_mean) < (vec_local_max_mean [res_mean.depth(li_mean)-1]- vec_local_max_sd [res_mean.depth(li_mean)-1]) ){
						res.erase(li);
						res_mean.erase(li_mean);
						li=res.begin_leaf(it);
						li_mean=res_mean.begin_leaf(it_mean);
					}else{
						++li;
						++li_mean;
					}
				}
				
				notincremented=true;
				cnt++;
					
				if(cnt!=0){ // delete elements which have no children in the next level
					it_tmp2=res.begin(); it_mean_tmp2=res_mean.begin();
						
					int rootdepth=res.depth(it_tmp2);
					while(it_tmp2!=res.end()) {
						if(res.depth(it_tmp2)==(max_elements_tmp-2) && res.number_of_children(it_tmp2)==0){
							it_tmp=res.parent(it_tmp2);it_mean_tmp=res_mean.parent(it_mean_tmp2);
							res.erase(it_tmp2);
							res_mean.erase(it_mean_tmp2);
							it_tmp2=it_tmp; it_mean_tmp2=it_mean_tmp;
						}else{
							++it_tmp2;++it_mean_tmp2;
						}							
					}
				}
			}
			
			//if actually in the next level
			if(res.depth(it)-rootdepth==max_elements_tmp-1){
				nsub[res.depth(it)-1]=*res.parent(it);   //previous node did not have any children, should be replaced by actual parent; same parent is the same for all siblings
				nsub[res.depth(it)]=(*it) ;

				//add all previously selected variables to prev_sel
				for (unsigned int i=0;i<max_elements_tmp-1;++i){	
					prev_sel[i]=nsub[i+1];					
				}
				
				max_mrmr=-1000.0;
				for(unsigned int k=0;k< n;++k){
					vec_mean[k]=0;
					vec_sd[k]=0;
				}
				//mrmr score previously selected variables should not be computed
				for(unsigned int k=0;k<(max_elements_tmp-1);++k){
					vec_mean[prev_sel[k]-1]=-1000;
					vec_sd[prev_sel[k]-1]=-1000;
				}
				
				vec_mean[nsub[0]-1]=-1000;
				vec_sd[nsub[0]-1]=-1000;
				
				for(unsigned int k=1;k<= n;++k){
					max_mrmr=-1000.0;
					if(vec_mean[k-1]!= (-1000)){
						bootstrap_mrmr(vec_mean[k-1],vec_sd[k-1], data,n, rep_boot,nsamples_boot,nsamples,nsub[0],k, max_elements_tmp-1,prev_sel);				
					}	
				}
				
				max_mrmr_ind=-1;
				for(unsigned int j=0;j<n;++j){
					if(j!=(nsub[0]-1) && vec_mean[j]>max_mrmr){
						max_mrmr=vec_mean[j];
						max_mrmr_ind=j;
					}
				}

				if(max_mrmr>0){
					if(local_max_mean==0){
						local_max_mean=max_mrmr;
						local_max_sd=vec_sd[max_mrmr_ind];
					}
					if( vec_local_max_mean[res.depth(it) ]== (-1000) ){
						vec_local_max_mean[res.depth(it)]=max_mrmr;
						vec_local_max_sd[res.depth(it)]=vec_sd[max_mrmr_ind];
					}else if( max_mrmr > vec_local_max_mean[res.depth(it)]){
						vec_local_max_mean[res.depth(it)]=max_mrmr;
						vec_local_max_sd[res.depth(it)]=vec_sd[max_mrmr_ind];
					}
					
					res.append_child(it,max_mrmr_ind+1);
					res_mean.append_child(it_mean,max_mrmr);
					
					//append equivalently good nodes to the two trees
					for(unsigned int k=0;k< n;++k){
						if(k!=max_mrmr_ind && vec_mean[k]< vec_mean[max_mrmr_ind]+vec_sd[max_mrmr_ind] && vec_mean[k]> vec_mean[max_mrmr_ind]-vec_sd[max_mrmr_ind]){
							res.append_child(it,k+1);
							res_mean.append_child(it_mean,vec_mean[k]);
						}
					}
				}
			}else{
				nsub[res.depth(it)]=(*it) ;
			}
			local_max_mean=0.0;
			++it;
			++it_mean;
			li=res.begin_leaf(it);
			li_mean=res_mean.begin_leaf(it_mean);
		}
		cnt++;
		max_elements_tmp++;
		notincremented=false;
		
	}
	///////////////////////
	//last level
	///////////////////////
	//remove horizontally the non-equivalent nodes
	tree<int>::pre_order_iterator it=res.begin();
	tree<double>::pre_order_iterator it_mean=res_mean.begin();
	tree<int>::leaf_iterator li;
	tree<double>::leaf_iterator li_mean;
	
	li=res.begin_leaf(it);
	li_mean=res_mean.begin_leaf(it_mean);
	while(li_mean!=res_mean.end_leaf(it_mean)) {
		if((*li_mean) < (vec_local_max_mean [res_mean.depth(li_mean)-1]- vec_local_max_sd [res_mean.depth(li_mean)-1]) ){
			res.erase(li);
			res_mean.erase(li_mean);
			li=res.begin_leaf(it);
			li_mean=res_mean.begin_leaf(it_mean);
		}else{
		++li;
		++li_mean;
		}
	}
	//remove the childless node in the next-to-last level
	tree<int>::pre_order_iterator it_tmp,it_back;
	tree<double>::pre_order_iterator it_mean_tmp,it_mean_back;
	if(cnt!=0){ // delete elements which have no children in the next level
		it=res.begin();
		it_mean=res_mean.begin();
		
		int rootdepth=res.depth(it);

		while(it!=res.end() ) {
			if(res.depth(it)<=(max_elements_tmp-2) && res.number_of_children(it)==0){
				it_tmp=res.parent(it);it_mean_tmp=res_mean.parent(it_mean);
				bool found_child=false;
				bool multiple=false;
				while(!found_child && (it_tmp!=res.begin() || res.number_of_children(it_tmp)>1)){ //end loop if there is a node with more than one child or if back at top and number of children==1 for top node
					if(res.number_of_children(it_tmp)==1){
						it_back=it_tmp;it_mean_back=it_mean_tmp; //if this was the last level for which there is only one child
						multiple=true;
						if(it_tmp!=res.begin()){
							it_tmp=res.parent(it_tmp);
							it_mean_tmp=res_mean.parent(it_mean_tmp);
						}else{//in case of having arrived at the top node
							res.erase(it_back);
							res_mean.erase(it_mean_back);
							found_child=true;
						}

					}else{
						if(multiple){
							res.erase(it_back);
							res_mean.erase(it_mean_back);
						}else{
							res.erase(it);
							res_mean.erase(it_mean);
						}
						found_child=true;
					}
				}
				it=it_tmp;
				it_mean=it_mean_tmp;
			}else{
				++it;
				++it_mean;
			}							
		}
	}	
}

SEXP mrmr_ensemble( SEXP Rdata, SEXP Rmaxparents, SEXP Rnvar, SEXP Rnsample, SEXP Rpredn, SEXP Rnpredn, SEXP Rrep_boot, SEXP Rmaxnsol){
	// Rdata:		data should be passed as vector, variable-wise appended
	// Rmaxparents:	number of maximum number of parents
	// Rnvar:		number of variables in the dataset
	// Rnsample:	number of samples in the dataset
	// Rpredn:		vector of target genes to consider
	// Rnpredn:		number of target genes (number of elements in Rpredn)
	// Rrep_boot:	how many bootstrap iterations
	// Rmaxnsol:	maximum number of children for each node at each step
	
	double *data;
	const int* maxparents, * nvar, *nsample, *maxnsol;

	int *predn, *rep_boot,*res,*res_all,*res_all2;
	int vec_tmp;
	const int *npredn;
	
	SEXP Rres;

	srand (54321);
	
	PROTECT(Rdata = AS_NUMERIC(Rdata));
	PROTECT(Rmaxparents= AS_INTEGER(Rmaxparents));
	PROTECT(Rnvar= AS_INTEGER(Rnvar));	
	PROTECT(Rnsample= AS_INTEGER(Rnsample));
	PROTECT(Rpredn = AS_INTEGER(Rpredn));
	PROTECT(Rnpredn = AS_INTEGER(Rnpredn));
	PROTECT(Rrep_boot = AS_INTEGER(Rrep_boot));
	PROTECT(Rmaxnsol= AS_INTEGER(Rmaxnsol));
	
	data=NUMERIC_POINTER(Rdata);
	maxparents = INTEGER_POINTER(Rmaxparents);
	nvar= INTEGER_POINTER(Rnvar);
	nsample= INTEGER_POINTER(Rnsample);
	predn= INTEGER_POINTER(Rpredn);
	npredn= INTEGER_POINTER(Rnpredn);
	rep_boot= INTEGER_POINTER(Rrep_boot);	
	maxnsol= INTEGER_POINTER(Rmaxnsol);

	tree<int> res_tree;
	tree<int>::iterator top,one;
	tree<int>::breadth_first_queued_iterator  it_final;
	top=res_tree.begin();
	int length_res=0;
	int length_res_old;
	for(unsigned int i=0;i< *npredn;++i){
		one=res_tree.insert(top, predn[i]);
		mrmr_ensemble_one_gene(res_tree, one, data,*nsample,*nvar,*maxparents,predn[i],*rep_boot);
		
		//std::cout<< std::endl ;
		//print_tree_int(res_tree,res_tree.begin(),res_tree.end());
		//std::cout<< std::endl ;
		////////////////////////
		//convert tree to vector
		////////////////////////
		int *tmp_nchildren,*res_tmp;
		res_tmp=new int [2*(res_tree.size())+1];
		tmp_nchildren= new int [(res_tree.size())];

		it_final=res_tree.begin_breadth_first();
		int cnt=1,cnt2=0;
		
		res_tmp[0]=res_tree.size();
		int rootdepth=res_tree.depth(it_final);
		while(it_final!=res_tree.end_breadth_first()) {
			res_tmp[cnt]=*it_final;
			tmp_nchildren[cnt-1]=(res_tree.number_of_children(it_final));
			cnt++;
			++it_final;		
		}
		
		////////////////
		//save in final result vector
		////////////////
		length_res_old=length_res;
		length_res+=2*(res_tree.size())+1;
		int *res_all, *res_old;
		int ind=0;
		res_all=new int[length_res];
		if(length_res_old>0){
			for(unsigned int k=0;k<length_res_old;k++){
				res_all[k]=res_old[k];
			}
		}
			
		for(unsigned int k=0;k<=res_tree.size();k++){
				res_all[length_res_old+k]=res_tmp[k];
		}
		for(unsigned int k=0;k<res_tree.size();k++){
			res_all[length_res_old+k+res_tree.size()+1]=tmp_nchildren[k];
		}
		
		delete [] res_old;
		res_old=new int[length_res];
		for(unsigned int k=0;k<length_res;k++){
			res_old[k]=res_all[k];
		}

		delete [] res_all;
		if(i==(*npredn-1)){
			PROTECT(Rres = NEW_INTEGER(length_res));
			res = INTEGER_POINTER(Rres);
			for(unsigned int k=0;k<length_res;k++){
				res[k]=res_old[k];
			}
			delete [] res_old;
		}
		////////////////
		//erase old tree
		////////////////

		delete [] tmp_nchildren;
		delete [] res_tmp;
		res_tree.erase(one);
	}

	UNPROTECT(9);

	return Rres;
}
