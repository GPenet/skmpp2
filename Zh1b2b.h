
/*
ZhouSolver.h
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/

/* the BIT_SET_30 pattern is defined in a band for a digit
   it is a map of candidates plus 3 bits for unsolved rows
   bits 0-8 9-17 18-26 for the rows 27-29 for unsolved rows
*/
//#include "t_128GP.h"
// tables specific to the brute force located in zh4_tables
const extern int TblRowMask[8];// rows where single  found  000 to 111
const extern int  Tblstartblock[27]; // (i%3) * 27 in zhou brute force
const extern int TblShrinkMask[512];// existing minirows 000 to 111
const extern int TblComplexMask[512]; // keep mini rows still valid optimised process
const extern int TblMaskSingle[512]; // kill in other blocks locked column /box
const extern int TblMaskDouble[512];// kill for locked in box / column
const extern int TblColumnSingle[512]; // single in column applied to shrinked bloc
const extern int TblShrinkSingle[512]; // keep only rows with single
const extern int TblRowUniq[512]; // 1 is row not defined in block  mode  to 111
const extern T128 AssignMask_Digit[81];
//const extern T128 AssignMask_OtherDigits[81];
// tables for UA collectors
extern int floors_2d[36];
extern int floors_3d[84]; // 84 3d  
extern int floors_4d[126];
extern int subsets_2d_3d[84][3];
extern int subsets_3d_4d[126][4];
extern int subsets_2d_4d[126][6];
extern int perm_0_8_size3[84][3];
extern int perm_0_8_size4[126][4];
//int * floors_3d_4d = floors_2d; // 84 3d   126 4d/5d


struct ZHOU;




/* class encapsulating the brute force with one known band
remaining clues are given in a 0-53 "cells" space
the 2 bands are located in a 64 bits field
unknown rows per digit are using 2x32 bits
first 5x6 bits for digits 0-4
second 4x6bits for digits 5-8

*/
struct ZH2B_GLOBAL { // global variables for the game table
	int nsol, lim, icount, ntsol, single_applied, new_single_in_Update,
		rdigit, nctlg, go_back,
		test;
	uint64_t * digsols; // pointer to solution grid per digit
	uint64_t ua_ret;
	BF64 val_init1_81, pairs, triplets;
	BF64 Digit_cell_Assigned_init[9];
	char * zsol, out54[55];// bands 12 in output mode
	char zdebug[82];

	// band UA collection active band pointers and UA table to build
	int modeguess;
	uint64_t *tuaold,tua;// old floors uas /  fresh uas
	uint32_t nuaold, nua,ndigits;;
	int  puz0[54], gangster[9];
	BF64 fd_sols[2][9];//start puzzle/ solution
	BF64 fdsw[3][9];//morphed digits puzzle/ solution rev
	BF64 previous_ua_status[6];// index 0 is for digit 3
	// bands sols handling
	BF64 sols_buffer[3000], ua_buffer[3000],*tsolw,*tuaw,
		mysol,mystart,andsol,myandsol;
	int nsolw;

	ZH2B_GLOBAL();
	inline void InitsGetSols(int i,int ibuf){
		mystart = fdsw[1][i];
		mysol = fdsw[0][i];
		tsolw =&sols_buffer[ibuf];
		tuaw = &ua_buffer[ibuf];
	}
	void GetBands(int * g1, int * g2);
	void Genuas2();

};
/* 2 BF 128 per digit
	equivalent to F and Fcomp in Zhou
	Last 32 bits in FD[0] is a  bits field for unknown rows
	Last 32 bits in FD[1]] contain  digit mapped
*/
// class encapsulating the brute force 
struct ZH2B {// size 32 bytes 
	BF64 FD[9], CompFD[9], cells_unsolved, rows_unsolved;

	void Init_std_bands(); // after getbands in zh2b_g
	void DebugSol();


	inline void Copy(ZH2B & o) { *this = o; }
	inline void Assign(int digit, int cell, int xcell) {
		FD[digit] &= AssignMask_Digit[cell].u64[0];
		cells_unsolved.Clear(xcell);
		int ddig = 6 * digit;
		if (digit > 4) ddig += 2;// second bloc of 32 bits
		rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	}
	inline int Unsolved_Count() { return rows_unsolved.Count(); }
	inline void ComputeNext() { if (FullUpdate())Guess(); }
	inline void ComputeNextFalse() { if (FullUpdate())GuessFalse(); }
	int Isvalid();
	int IsvalidNoUpdate(int debug = 0); // usually after init 2 steps
	uint64_t GetUa();
	void Init_x_(GINT64 t, int n);
	uint64_t Init_y_(GINT64 t, int n);
	void Init_xy(ZH2B & o, GINT64 tx, int nx, GINT64 ty, int ny);

	int Update();
	int FullUpdate();
	int FullUpdateNoGuess();
	void Guess();
	void GuessFalse();
	int ApplySingleOrEmptyCells();
	uint64_t CheckUa(uint64_t ua);
	char * SetKnown(char * zs);
	int Seta(int digit, int xcell);
	inline int SetaC(int digit, int cell) { return Seta(digit, C_To128[cell]); }
	inline void ClearCandidate(int dig, int xcell) { FD[dig].Clear(xcell); }
	int GetFreeDigits(int cell);
	int GetAllDigits(int cell);
	void Debug(int all=0);
	void DebugDigit(int digit);
	void ImageCandidats();
	int DebugCheckUa(uint64_t ua);



// located in go_17sol_zx_UaCollector_cpp  generation of UAs on a multi floor
	void Guess2(); void GenUas3();
//void InitGenUas(int * zpuz);
	//void InitUas2();
	//void GenUas(int floors);
	//void GenUasBands12LessMinirow(int floors, int c1, int c2, int c3);
	//void GenUAs2_minirowpair(int c1, int c2, int dig1, int dig2, int i_81);
	//int GenMoreGUAs2(int c1, int c2, int dig1, int dig2);
	//void GenUAs2_minirowtriplet(int c1, int c2, int c3, int dig1, int dig2, int dig3, int i_81);
	//void GenUas2();		void GenUas4();	void GenUas5();
	//int Update2();	int Update3();	int Update4();	int Update5();
	//	void Guess3();	void Guess4();	void Guess5();
	//void CollectFinal(TUA64 & mt);
	//int CollectFinalua2s(uint64_t *td, int maxt, int n0);
	//void DebugFloors(int nfloors, int all = 0);



	
	/*
 
 // located in zh_uacollector.cpp  generation of UAs on a multi floor
	void InitGenUas(char * zpuz);
	void GenUas( int floors);
	void GenUasBands12(int floors);
	void GenUas6();
	void GenMinirows2();
	void GenMinirows3();
	int Update6();
	void Guess6();
	void Guess5_4();
	void Guess5_3();
	int CollectFinal(BF128 *td,int &lim10);
   */
 };
 /* class encapsulating the brute force for one band
initial is a solved band giving the clues in each column
each bit map per digit is a 32 bit field
as in the original code, unknown rows per digit are located in bits 28-29-30 of the bit map
remaining clues are given in a 0-26 "cells" space

*/
struct ZHONE_GLOBAL { // global variables for the game table
	int diag;
	int nsol, lim, icount, ntsol, single_applied,
		new_single_in_Update,
		type,//gocatmode,    
		rdigit, nctlg, go_back;
	int pairs, triplets;
	char * zsol, out27[28];// band1 in output mode 

	// band UA collection active band pointers and UA table to build
	uint32_t *tua, nua;//   maximum 81  
	int * band0,*gangster;
	int floors_mini_row, digmap[9];// to adjust pm
	uint32_t fd_sols[2][9];//start puzzle/ solution
	uint32_t fdsw[3][9];//morphed digits puzzle/ solution rev
	uint32_t ndigits;
	uint32_t previous_ua_status[6];// index 0 is for digit 3

	ZHONE_GLOBAL();
	void GetBand(int * b, uint32_t * t);
	inline void SetPat(char * pat, char * zsol, char * puzfinal) {
		pat = pat; zsol = zsol; puzfinal = puzfinal;
	};
	inline void InitIsvalid() { // usually after init 2 steps
		nsol = go_back = type = 0; lim = 1;
	}
	void ValidPuzzle(uint32_t * sol);
	void AddUA(uint32_t ua);
	void PrintTua();
	void FindMissingUAs();
};
struct ZH2B_1D {// size 32 bytes
	BF64 FD, CompFD; 
	int GetSols( int ru);
	inline void Assign( int cell) {
		FD &= AssignMask_Digit[cell].u64[0];
	}
	void ComputeNext(int ru);
	int Update(int &ru);
	void Guess(int ru);
	int IsValid(uint64_t v);
};

 struct ZHONE {
	 uint32_t FD[9], CompFD[9], cells_unsolved;

	 // standard calls and basic process
	 inline void Assign(int digit, int cell) {
		 FD[digit] &= AssignMask_Digit[cell].u32[0];
		 cells_unsolved &= ~(1 << cell);
		 FD[digit] &= ~(1 << (27 + C_row[cell]));//9*digit + row
	 }
	 inline int Unsolved_Count() { return __popcnt(cells_unsolved); }
	 int ApplySingleOrEmptyCells();
	 void InitOne_std_band(); // after getband in zh1b_g
	 int InitSudokux(GINT * t, int n);
	 void AddMissingUAs(int * tcells,int ncells);
	 int Update();
	 int FullUpdate();
	 void Guess();
	 inline void ComputeNext() { if (FullUpdate())	Guess(); }
	 char * SetKnown(char * zs);
	 void Seta(int digit, int xcell);
	 int Isvalid();
	 int GetAllDigits(int cell);
	 void Debug(int all = 0);
	 void DebugDigit(int digit);
	 void ImageCandidats();
	 //==================================================
	 // uacollector 
	 void Checkstart();
	 int  Start_nFloors(int floors);
	 void Start_Uas_Mini(int floors, int floors_mini_row);
	 void ApplyGangsterChanges(int * g0, int * g1);
	 void InitGuess();
	 int UpdateDigit(int digit);
	 void Guess2();	void Guess3();	void Guess4();
	 void Guess5(); void Guess6();  void Guess7();
	 void Set2(int cell);
 };

