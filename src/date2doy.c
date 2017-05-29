#include<stdio.h>
#include<R.h>
#include "date2doy.h"

//rdate as YYMMDD
//returns day of year
int *dateToDoy(int *rdate, int *doy)
{
	int date = rdate[0];
	int day = 0, month = 0, year = 0;
	int daysofmonth[12];
	int i=0;
	
	*doy=0;

	day = date % 100;
	date = (int) ((double) date / 100.0);
	month = date % 100;
	date = (int) ((double) date / 100.0);
	year = date;

	year = getFullYear(year);

	//printf("Year: %d, Month: %d, Day: %d, LeapYear: %d\n", year, month, day, isLeap(year));

	for (i=0; i<12; i++){
		if ((i==3) || (i==5) || (i==8) || (i==10)){
			daysofmonth[i] = 30;
		} else {
			if (i==1){
				daysofmonth[i] = 28+isLeap(year);
			} else {
				daysofmonth[i] = 31;
			}
		}
	}
	for (i=0; i < month-1; i++)
		*doy = *doy+daysofmonth[i];
	*doy = *doy + day;

	return doy;
}

int *getDaysCount(int *rdate, int *count){
	int date = rdate[0];
	date = (int) ((double) date / 10000.0);
	if (isLeap(getFullYear(date))==0)
		*count=365;
	else
		*count=366;
	return count;
} 

//is year YYYY a leap year ?
int isLeap(int year) {
	if (year % 4 == 0)
		if (!(year % 100 == 0))
			return 1;
		else
			if (year % 400 == 0)
				return 1;
	return 0;
}

//Parameter YY, returns YYYY
//if YY is greater then 80 - 20th Century
//else 21th Century
int getFullYear(int year){
	if (year > 80)
		return (year + 1900);
	else
		return (year + 2000);
}

