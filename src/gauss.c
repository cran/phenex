#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include "gauss.h"
#include "quicksort.h"
#include "condition.h"

//fits a gaussian function at ndvi-values
double *gauss(int *rdays, double *ndvi, int *mustart, double *sigstart, double *rbase,double *model){
	int i=0, j=0, k=0;
	
	//initialize parameters
	int days = rdays[0];
	int mu = mustart[0];
	double sig = sigstart[0];
	double base = rbase[0];
	double scal = 1;
	
	//initialize parameter vectors
	int muvec[8];
	double sigvec[20], scalvec[10];

	//initialize other variables
	int *sorted;
	sorted = (int *) calloc(days, sizeof(int));
	
	double convec;
	double minconvec;
	int minconvecpos[3];
	
	//initialize gaussian function vector
	double *gfct;
	gfct = (double *) calloc(days,sizeof(double));

	//fill muvec
	sort(ndvi, sorted, days);
	for (i=0; i<8; i++){
		muvec[i]=sorted[days-1-i];
	}

	//fill sigvec
	for (i=0; i < 20; i++){
		sigvec[i]= 0.001 + (sig/10)*((double) (i+1));
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

	for (i=0; i < 20; i++){
		for (j=0; j < 8; j++){
			for (k=0; k < 10; k++){
				convec=condition(ndvi, gaussianfct(muvec[j], sigvec[i], scalvec[k], base, days, gfct), days);
				if ((i==0)&&(j==0)&&(k==0)){
					minconvec = convec;
					minconvecpos[0]=i;
					minconvecpos[1]=j;
					minconvecpos[2]=k;
				} else {
					if (convec < minconvec){
						minconvec = convec;
						minconvecpos[0]=i;
						minconvecpos[1]=j;
						minconvecpos[2]=k;
					}
				}
			}
		}
	}

	//get optimized parameters
	sig = sigvec[minconvecpos[0]];	
	mu = muvec[minconvecpos[1]];
	scal = scalvec[minconvecpos[2]];
	
	gfct = gaussianfct(mu, sig, scal, base, days, gfct);
	
	for (i=0; i<days; i++){
		model[i]=gfct[i];
	}

	free(sorted);
	free(gfct);
	
	return(model);
	
}

//creates the gaussian function
double *gaussianfct(int mu, double sig, double scal, double base, int length, double *gfct){
	int i=0;
	double maxfirst=0, max=0;
	
	for (i=0; i<length; i++){
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
