#include<stdio.h>
#include<R.h>
#include "analyzeBits.h"

char * readBits(short number){
	static char out[16];
	int i=0;

	//initialize with zeros
	memset( out, '\0',sizeof(out));
	for (i=15; i>=0; i--){
		if ((number & 1) == 0)
			out[15-i] = '0';
		else
			out[15-i] = '1';

		number = number >> 1;
	}

	return out;
}

int * getBits(int * number, int * pos, int * m){
	char * bits;
	char c;
	int erg;
	bits = readBits(number[0]);
	c = bits[pos[0]];	
	*m = atoi(&c);
	return m;
}

int * getLandWater(int * number, int * lw){
	char * bits;
	char c;
	int erg;
	bits = readBits(number[0]);
	c = bits[15];	
	*lw = atoi(&c);
	return lw;
}

int * getCloudmask(int * number, int *cm){
	char * bits;
	char c;
	bits = readBits(number[0]);
	c = bits[14];
	*cm = atoi(&c);
	return cm;
}

int * getDayNo(int *number, int *day){
	char * bits;
	char cday[4+1];
	int dual=0, dec=0;
	bits = readBits(number[0]);
	cday[0]=bits[0];
	cday[1]=bits[1];
	cday[2]=bits[2];
	cday[3]=bits[3];
	cday[4]='\0';

	dual = atoi(cday);

	dec = (dual % 10)*1;
	dual = dual/10;
	dec = dec + (dual % 10)*2;
	dual = dual/10;
	dec = dec + (dual % 10)*4;
	dual = dual/10;
	dec = dec + (dual % 10)*8;
	*day = dec;

	return day;
}

