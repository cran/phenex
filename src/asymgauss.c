#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include "asymgauss.h"
#include "quicksort.h"
#include "condition.h"

//fits a gaussian function at ndvi-values
double *asymgauss(int *rdays, double *ndvi, int *mustart, double *sigstart, double *rbase,double *model){
	int i=0, j=0, k=0, m=0;
	
	//initialize parameters
	int days = rdays[0];
	int mu = mustart[0];
	double sig = sigstart[0];
	double base = rbase[0];
	double scal = 1;
	double sig1, sig2;
	
	//initialize parameter vectors
	int muvec[8];
	double sigvec[20], scalvec[10];

	//initialize other variables
	int *sorted;
	sorted = (int *) calloc(days, sizeof(int));
	
	double convec;
	double minconvec;
	int minconvecpos[4];
	
	//initialize gaussian function vector
	double *gfct;
	gfct = (double *) calloc(days,sizeof(double));

	//fill muvec
	sort(ndvi, sorted, days);
	for (i=0; i<5; i++){
		muvec[i]=sorted[days-1-i];
	}

	//fill sigvec
	for (i=0; i < 10; i++){
		sigvec[i]= 0.001 + (sig/5)*((double) (i+1));
	}
	
	//fill scalvec
	for (i=0; i<10; i++){
		scalvec[i]= 1 + i*0.5;
	}

	//conditioning
	minconvec=0;
	minconvecpos[0]=0;
	minconvecpos[1]=0;
	minconvecpos[2]=0;
	minconvecpos[3]=0;

	for (i=0; i < 10; i++){
		for (j=0; j < 5; j++){
			for (k=0; k < 10; k++){
				for (m=0; m<10;m++){
					convec=condition(ndvi, asymgaussianfct(muvec[j], sigvec[i], sigvec[m], scalvec[k], base, days, gfct), days);
					if ((i==0)&&(j==0)&&(k==0)){
						minconvec = convec;
						minconvecpos[0]=i;
						minconvecpos[1]=m;
						minconvecpos[2]=j;
						minconvecpos[3]=k;
					} else {
						if (convec < minconvec){
							minconvec = convec;
							minconvecpos[0]=i;
							minconvecpos[1]=m;
							minconvecpos[2]=j;
							minconvecpos[3]=k;
						}
					}
				}
			}
		}
	}

	//get optimized parameters
	sig1 = sigvec[minconvecpos[0]];
	sig2 = sigvec[minconvecpos[1]];	
	mu = muvec[minconvecpos[2]];
	scal = scalvec[minconvecpos[3]];
	
	gfct = asymgaussianfct(mu, sig1, sig2, scal, base, days, gfct);
	
	for (i=0; i<days; i++){
		model[i]=gfct[i];
	}

	free(sorted);
	free(gfct);
	
	return(model);
	
}

//creates the gaussian function
double *asymgaussianfct(int mu, double sig1, double sig2, double scal, double base, int length, double *gfct){
	int i=0;
	double maxfirst=0, max=0;
	double sig=sig1;
	
	for (i=0; i<length; i++){
		if (i > mu){
			//second half with other sig for asymmetric Gaussian function
			sig = sig2;
		}
		gfct[i] = 0;
		gfct[i] = pow((1.0/(sig*sqrt(2*M_PI))),(-0.5*pow((((((double) i)/((double) length))-((double) mu)/ ((double) length))/sig),2))) / scal;
		maxfirst = (gfct[i] > maxfirst) ? (gfct[i]) : (maxfirst);
		gfct[i] += base;
		max = (gfct[i] > max) ? (gfct[i]) : (max);
	}
	for (i=0; i < length; i++){
		gfct[i] = (gfct[i]/max)*maxfirst;
	}	

	return(gfct);
}
