#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include "ravg.h"

//running average
double *runAVG(int *rdays, double *ndvi, int *rwindow, double *corrected_ndvi){
	int days = rdays[0];
	int window = rwindow[0];
	int i=0, j=0, count=0, position=0;
	int noNAlength=0;
	double avg;
	
	double *temp_corrected_ndvi;
	double *temp_ndvi;
	int *noNA;

	temp_corrected_ndvi = (double *) calloc(3*days,sizeof(double));
	temp_ndvi = (double *) calloc(3*days,sizeof(double));

	//Triple data to avoid inconsistency of start and end of year
	for (i=0; i < days;i++ ){
		temp_ndvi[i] = ndvi[i];
		temp_ndvi[i+days] = ndvi[i];
		temp_ndvi[i+2*days] = ndvi[i];
		if (ndvi[i] >= 0){
			noNAlength++;
		}
	}

	//find no-NA positions	
	noNAlength = noNAlength * 3;
	noNA = (int *) calloc(noNAlength, sizeof(int));

	count=0;
	for (i=0; i < (3*days); i++){
		if (temp_ndvi[i] >= 0){
			noNA[count] = i;
			count++;
		}
		if (count > (noNAlength-1)){
			break;
		}
	}

	//make window odd
	if ((window % 2) == 0){
		window++;
	}

	//smooth data
	for(i=days; i < (2*days-1); i++){
		if (temp_ndvi[i] < 0){
			temp_corrected_ndvi[i] = -1;
		} else {
			avg=0;
			//search position
			count=0;
			position=0;
			while (i != noNA[count]){
				count++;
				if (count > noNAlength){
					break;
				}
			}
			if (i == noNA[count]){
				position = count;
			} else {
				//invalid position
				temp_corrected_ndvi[i] = -1;
				continue;
			}

			//search valid values
			if (((position-((window-1)/2)) >= 0)&&((position+((window-1)/2)) < noNAlength)){
				for (j=(position-((window-1)/2)); j <=(position+((window-1)/2)); j++){
					avg = avg + temp_ndvi[noNA[j]];
				}
				avg = avg/ ((double) window);
			} else {
				//not enough values
				temp_corrected_ndvi[i] = -1;
				continue;
			}
			
			temp_corrected_ndvi[i]=avg;
		}
	}
	
	//cut results
	for ( i=days; i<(2*days-1); i++){
		corrected_ndvi[i-days] = temp_corrected_ndvi[i];
	}

	free(temp_corrected_ndvi);
	free(temp_ndvi);
	free(noNA);

	return corrected_ndvi;
}
