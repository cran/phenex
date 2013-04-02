#include<stdio.h>
#include<math.h>
#include<R.h>
#include "sigmoid.h"
#include "quicksort.h"
#include "condition.h"

//sigmoid function creation
//F .. base
//G+F .. maximum
double *sigmoid(int *rdays, double *ndvi, double *rF, double *rG,double *model){
	int i=0, j=0, count=0;
	
	//initialize parameters
	int days=rdays[0];
	double F=rF[0];
	double G=rG[0];
	double GF=0;

	//initialize parameter vectors
	double kvec[100];
	double dvec[8];

	//initialize other variables
	double k;
	double d;

	double convec[100][8];
	double minconvec[100*8];
	int minconvecpos[100*8];
	for (i=0; i<(100*8); i++){
		minconvecpos[i]=i;
	}
	
	//initialize sigmoid function vector
	double *sfct;
	sfct = (double *) calloc(days,sizeof(double));

	//fill kvec
	for (i=0; i<100; i++){
		kvec[i]=0.01+0.01*((double) i);
	}

	//fill dvec
	for (i=0; i<8;i++){
		GF = (G/F-1.0 < 0) ? (1.0-G/F) : (G/F-1.0);
		dvec[i]= pow(10,(double) i)* GF;
	}
	
	//conditioning
	for (i=0; i<100; i++){
		for (j=0; j<8; j++){
			convec[i][j] = 0;
			convec[i][j] = condition(ndvi, sigmoidfct(F,G, kvec[i], dvec[j], days, sfct), days);
		}
	}
	
	//search minima
	for (j=0; j<8; j++){
		for (i=0; i<100; i++){
			minconvec[j*100+i]=convec[i][j];
		}
	}
	sort(minconvec, minconvecpos, 100*8);

	count = 0;
	do{
		//get optimized parameters
		k = kvec[minconvecpos[count] % 100];
		d = dvec[minconvecpos[count]/100];
	
		sfct = sigmoidfct(F,G,k,d,days,sfct);

		count++;
		if (count >= 100){
			k = kvec[minconvecpos[0] % 100];
			d = dvec[minconvecpos[0]/100];
			sfct = sigmoidfct(F,G,k,d,days,sfct);
			break;
		}
	} while (maximum(sfct, days) < ( (G+F)-0.01*(G+F) ) );

	for (i=0; i<days; i++){
		model[i]=sfct[i];
	}

	free(sfct);

	return(model);
}

//sigmoid function
double *sigmoidfct(double F, double G, double k, double d, int length, double *sfct){
	int i=0;
	for (i=0; i<length; i++){
		sfct[i]=0;
		sfct[i]= G/(1 + d*exp(-k*G*((double) i))) + F;
	}
	return(sfct);
}

//returns the maximum of a vector (type double)
double maximum(double *values, int length){
	int i=0;
	double max=values[0];
	for (i=0; i<length; i++){
		if (values[i] > max){
			max=values[i];
		}
	}
	return(max);
}
