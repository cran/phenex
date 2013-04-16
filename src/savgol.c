#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include "savgol.h"

double *savGol(int *rdays, double *ndvi, double *coef, int *numberofCoef, double *corrected_ndvi){
	int days, h=21, i, j, count, step, position, errorflag=0;
	int nCoef, middle;
	double *factor;
	double *temp_ndvi, *nothing_found;	
	double tmp;

	//copy values from pointer variables to integer variables
	nCoef = numberofCoef[0];
	middle = (int) nCoef/2;
	days = rdays[0];

	//allocate memory
	temp_ndvi = (double *) calloc(3*days,sizeof(double));
	factor = (double *) calloc(nCoef,sizeof(double));
	nothing_found = (double *) calloc(days, sizeof(double));

	//initialize factor
	for (i=0; i < nCoef; i++){
		factor[i] = 0;
	}

	//initialize nothing_found
	for (i=0; i<days; i++){
		nothing_found[i] = -1;
	}

	//Triple data to avoid inconsistency of start and end of year
	for (i=0; i < days;i++ ){
		temp_ndvi[i] = ndvi[i];
		temp_ndvi[i+days] = ndvi[i];
		temp_ndvi[i+2*days] = ndvi[i];
	}
	
	//run Savitzky-Golay smoothing
	for (i=days; i < (2*days-1); i++){
		count=0;
		step=1;
		//search next valid value (through alternation)
		while (temp_ndvi[(int) round(pow(-1, count)*step+i)] <= 0){
			count++;
			if ((count % 2)==0){
				step++;
			}
		}
		position = (int) round(pow(-1, count)*step+i);
		if ((position < 0)||(position > (3*days-1))){
			for (j=0; j < days; j++)
				corrected_ndvi[j] = nothing_found[j];
			break;
		}
		factor[middle] = temp_ndvi[position];
		
		//search valid factors
		//search left
		step=1;
		for (j=0; j < middle; j++){
			if ((position-step) < 0) {
				errorflag=1;
				break;
			}
			while(temp_ndvi[position-step] <=0){
				step++;
				if (step >= days){
					factor[j]=0;
				}
			}
			factor[j]=temp_ndvi[position-step];
			step++;
		}
		//search right
		step=1;
		for (j=0; j < middle; j++){
			if ((position+step) > (3*days-1)){
				errorflag=1;
				break;

			}
			while(temp_ndvi[position+step] <=0){
				step++;
				if (step >= days){
					factor[middle + j+1]=0;
					break;
				}
			}
			factor[middle + j+1]=temp_ndvi[position+step];
			step++;
		}
		if (errorflag == 1){
			for (j=0; j < days; j++)
				corrected_ndvi[j] = nothing_found[j];
			break;
		}
		//calculate new value at position i
		tmp=0;
		for (j=0; j < nCoef; j++){
			tmp = tmp+coef[j]*factor[j];
		}
		corrected_ndvi[i-days] = tmp;
	}
	
	//free the allocated memory
	free(temp_ndvi);
	free(factor);
	free(nothing_found);

	return corrected_ndvi;
}
