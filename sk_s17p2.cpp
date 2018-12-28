
#define _CRT_SECURE_NO_DEPRECATE
#define SEARCH17SOL
/* program organisation
	main is the standard frame including the basic brute force 
*/

//#define MODE66_ON
#define GTEST17_ON 1
#define UALIMSIZE 20
#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
#define TUA64_12SIZE 2000
//============================================== 

#define MAXN5 51520
#define MAXN6 237770 

#define G17MORESIZE 64

#define G17TESTUASGUASLIMITS 1

#include <sys/timeb.h>
#include "main.h"  // main and main tables and basic brute force
#include "go_17sol_tables.h"     
#include "Zh1b2b.h"  // brute force 2 bands  
//_________________ brute force handling 1 2 3 bands 
extern ZHOU    zhou[50],zhou_i;// , zhou_i, zhou_solve;
extern ZH_GLOBAL zh_g;
extern SGO sgo;

extern ZH2B_GLOBAL   zh2b_g;
extern ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands

extern ZH2B zh2b[40], zh2b_i, zh2b_i1;
extern ZH2B5 zh2b5[10]; // solved digit per digit
extern ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row

extern ZHONE_GLOBAL   zh1b_g;
extern ZHONE zhone[20];
extern ZHONE zhone_i;



//extern GENUAS_1B genuas1b;

ofstream  fout1, fout2;
#include "go_17sol.h"    
//XY_EXPAND xye6[2][MAXN6], xye5[2][MAXN5];
G17XY g17xy;
G17CHUNK g17chunk;
G17INDEXSTEP indexstep;
G17B g17b;
uint64_t p_cptg[40], p_cpt1g[20], p_cpt2g[20];
uint64_t p_cpt[40], p_cpt1[20], p_cpt2[20];
int a_17_found_here;

FINPUT finput;

#include "go_17_bands.h"  
#include "go_17_genb12.h"     
#include "go_17sol_tables.h"

//STD_B3 mytband3[256];
STD_B1_2 myband1, myband2;
GENUAS_B12 genuasb12;
GEN_BANDES_12 genb12;
#include "go_17_bands_cpp.h"  

#include "go_17_genb12_cpp.h"     
#include "go_17sol_bs_cpp.h"     
#include "go_17sol_zx_cpp.h"  
#include "go_17sol_commands_cpp.h"



void Go_0() {
	// open  outputs files 1.txt
	if (sgo.foutput_name) {
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_file1.txt");
		fout1.open(zn);
	}
	if (sgo.command >= 100) {// input file expected
		finput.open(sgo.finput_name);
		if (!finput.is_open()) {
			cerr << "error open file " << sgo.finput_name << endl;
			return;
		}
	}
	cerr << "running command " << sgo.command << endl;
	switch (sgo.command) {
	case 0: Go_c17_00(); break; // search one band1

	case 80: Go_c17_80(); break; // enumeration test 
	case 90: Go_c17_90(); break; // regressive test uas one band	
	case 91: Go_c17_91(); break; // test uas collector 2 bands
	case 92: Go_c17_92(); break;// regressive test uas 5 6  expand
	}
	cerr << "go_0 return" << endl;
}

