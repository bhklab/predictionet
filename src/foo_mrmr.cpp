#include "foo_mrmr.h"

double get_correlation(double data [],int ind_x, int ind_y, int size){
	//compute correlation of two variables;
	//data: contains all data in a vector; variable-wise appended
	//ind_x: starting index of first variable in data
	//ind_y: starting index of second variable in data
	//size: number of samples for both variables
	double mean_data_x=0.0,mean_data_y=0.0;
	double correlation_nom=0.0,correlation_den_x=0.0,correlation_den_y=0.0;
	
	for( unsigned int i=0; i< size; ++i ) {
		mean_data_x+=data[ind_x+i];
		mean_data_y+=data[ind_y+i];
	}
	
	mean_data_x=mean_data_x/size;
	mean_data_y=mean_data_y/size;
	
	for( unsigned int i=0; i< size; ++i ) {		
		correlation_nom+=(data[ind_x+i]-mean_data_x)*(data[ind_y+i]-mean_data_y);
		correlation_den_x+=(data[ind_x+i]-mean_data_x)*(data[ind_x+i]-mean_data_x);
		correlation_den_y+=(data[ind_y+i]-mean_data_y)*(data[ind_y+i]-mean_data_y);
	}
	return correlation_nom/(sqrt(correlation_den_x*correlation_den_y));
}


void build_mim_subset(double mim[],double data[],int nvar,int nsample, int subset [],int size_subset){
	//compute mutual information matrix
	//mim:			matrix (stored as vector) in which the mi values will be stored
	//data:			contains all data in a vector; variable-wise appended
	//nvar:			number of variables
	//nsample:		number of samples in dataset
	//subset:		indices of samples to be included in the bootstrapping data
	//size_subset:	number of variables in the bootstrapped dataset
	
	double tmp;
	double *data_x;
	
	data_x = new double [nvar*size_subset];
	
	for(unsigned int i=0; i< size_subset; ++i){
		for(unsigned int j=0; j< nvar; ++j){
			data_x[size_subset*j+i]=data[(subset[i])+nsample*j];
		}
	}
	
	for(unsigned int i=0; i< nvar; ++i){
		for(unsigned int j=i+1; j< nvar; ++j){
			tmp=get_correlation(data_x,i*size_subset,j*size_subset,size_subset);
			mim[j*nvar+i]= -0.5* log (1-tmp*tmp);
			mim[i*nvar+j]=mim[j*nvar+i];
		}
	}
	delete [] data_x;
}