#include<stdio.h>
#include<R.h>
#include<math.h>
#include "condition.h"

//conditioning function, mean of the differences between the input and the model
double condition(double *ndvi, double *model, int length){
	int i=0,count=0;
	double condition=0;
	for (i=0; i < length; i++){
		if (ndvi[i]>0){
			count++;
			condition += (ndvi[i] > model[i]) ? pow(ndvi[i]-model[i],2) : pow(model[i]-ndvi[i],2);
		}
	}
	condition /= (double) count;
	return(condition);
}



