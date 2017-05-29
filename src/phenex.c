#include <R_ext/Rdynload.h>
#include "analyzeBits.h"
#include "date2doy.h"
#include "ravg.h"
#include "savgol.h"
#include "sigmoid.h"
#include "quicksort.h"
#include "condition.h"
#define NULL 0

static const R_CMethodDef cMethods[] = {
        {"sigmoid",  (DL_FUNC) &sigmoid, 5},
        //{"sigmoidfct",  (DL_FUNC) &sigmoidfct, 6},
        //{"maximum",  (DL_FUNC) &maximum, 2},
	{"savGol",  (DL_FUNC) &savGol, 5},
	{"getBits",  (DL_FUNC) &getBits, 3}, 
	{"getLandWater", (DL_FUNC) &getLandWater, 2},
	//{"readBits", (DL_FUNC) &readBits, 1},
	{"getCloudmask", (DL_FUNC) &getCloudmask, 2},	
	{"getDayNo", (DL_FUNC) &getDayNo, 2},
	{"dateToDoy",  (DL_FUNC) &dateToDoy, 2},
        {"getDaysCount", (DL_FUNC) &getDaysCount, 2},
        {"isLeap", (DL_FUNC) &isLeap, 1},
        {"getFullYear", (DL_FUNC) &getFullYear, 1},
	//{"sort",  (DL_FUNC) &sort, 3},
        {"quicksort",  (DL_FUNC) &quicksort, 4},
        //{"partition",  (DL_FUNC) &partition, 4},
        //{"swap",  (DL_FUNC) &swap, 3},
	//{"condition",  (DL_FUNC) &condition, 3},
	{"runAVG",  (DL_FUNC) &runAVG, 4},
        {NULL, NULL, 0}
};

void R_init_phenex(DllInfo *info)                                                                                                                                                                             
{
        R_registerRoutines(info,
                cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);                                                                                                                                                                       

}
