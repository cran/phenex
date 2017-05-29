#include<stdio.h>
#include<R.h>
#include "quicksort.h"

void sort(double *values, int *sorted, int length){
	int i=0;
	for (i=0; i<length; i++){
		sorted[i]=i;
	}
	quicksort(0, length-1, values, sorted);
}

void quicksort(int left, int right, double *values, int *sorted){
	int divide=0;
	if (left < right){
		divide = partition(left, right, values, sorted);
		quicksort(left, divide-1, values, sorted);
		quicksort(divide+1, right, values, sorted);
	}
}

int partition(int left, int right, double *values, int *sorted){
	int i=left;
	int j=right-1;
	double pivot = values[sorted[right]];

	do {
		while ((values[sorted[i]] <= pivot)&&(i < right)){
			i++;
		}
		while ((values[sorted[j]] >= pivot)&&(j > left)){
			j--;
		}
		if (i < j){
			swap(i,j,sorted);
		}
	} while (i < j);
	
	if (values[sorted[i]] > pivot) {
		swap(i, right, sorted);
	}
	return(i);
}

void swap(int first, int second, int *sorted){
	int tmp=0;
	tmp = sorted[second];
	sorted[second]=sorted[first];
	sorted[first]=tmp; 
}

