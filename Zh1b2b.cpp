/*
Based on code posted to <http://forum.enjoysudoku.com/3-77us-solver-2-8g-cpu-testcase-17sodoku-t30470.html>
by user zhouyundong_2012.
The copyright is not specified.
*/
#include "main.h"   
#include "Zh1b2b.h"

ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row

//========================= global variable and working areas 
ZH2B zh2b[50]; // must host main brute force plus minimality analysis and recursive generation
ZH2B zh2b_i, zh2b_i1;
ZH2B_GLOBAL   zh2b_g;   // global variables for the class GAME
ZHONE zhone[50]; // must host main brute force plus minimality analysis and recursive generation
ZHONE zhone_i;//================== initial for a GAME
ZHONE_GLOBAL zh1b_g;
extern ZH_GLOBAL zh_g;
#include "go_17sol_zx_UaCollector_cpp.h"
//GENUAS_1Bx genuas1b;
//============================= ZH_GLOBAL code
ZH2B_GLOBAL::ZH2B_GLOBAL(){
	zsol = 0; // no solution unless required buy the user
	nctlg = test = 0;
	val_init1_81.bf.u64 = BIT_SET_2X;
	for (int i = 0; i < 9; i++) {
		zh2b_i.FD[i] = val_init1_81;
		zh2b_i.CompFD[i] = val_init1_81;
	}
	zh2b_i.cells_unsolved = val_init1_81;
	zh2b_i.rows_unsolved.bf.u32[0] = BIT_SET_30; //6x5
	zh2b_i.rows_unsolved.bf.u32[1] = 077777777;  //6x4
	for (int i = 0; i < 50; i++) {// for X+Y+27 format 50 blocs
		zh2b[i] = zh2b_i;
	}

}
void ZH2B_GLOBAL::GetBands(int * g1, int * g2) {
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = puz0[cell];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[0] |= 1 << cell;
			dig = puz0[cell+27];
			fd_sols[1][dig].bf.u32[0] |= Zhoucol << i;
			fd_sols[1][dig].bf.u32[1] |= Zhoucol << i;
			fd_sols[0][dig].bf.u32[1] |= 1 << cell;
		}
	}
}

void ZH2B_GLOBAL::Genuas2() {// 2 digits

}

//============================================ ZH2B code
void ZH2B::Init_std_bands() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	*this = zh2b_i;
	memcpy(FD, zh2b_g.fd_sols[1], sizeof FD);
	memset(CompFD, 0, sizeof CompFD);
}
void ZH2B::DebugSol() {//init after zh1b_g getband
	zh2b_g.ndigits = 9;
	memset(this, 0, sizeof this);
	memcpy(FD, zh2b_g.fd_sols[0], sizeof FD);
	ImageCandidats();
}

#define UPDN(I,J)A=FD[I].bf.u32[J];\
Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
FD[I].bf.u32[1 - J] &= TblMaskSingle[S]; \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
CompFD[I].bf.u32[J] = FD[I].bf.u32[J] = A;

#define UPWCL(I,P,Q,R,T,U,V,W,X)cl = ~(A & TblRowMask[S]); \
cells_unsolved.bf.u32[I] &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;

#define UPDN1(I,P,Q,R,T,U,V,W,X)Shrink = (TblShrinkMask[A & 0x1FF] | \
TblShrinkMask[ (A>>9) & 0x1FF]<<3 | \
TblShrinkMask[ (A>>18) & 0x1FF]<<6);\
if ((A &=TblComplexMask[Shrink]) ==0)  return 0; \
S = ((A | (A >> 9) | (A >> 18)) & 0x1FF); \
S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]]; \
if ((A>>27) != S){\
cl = ~(A & TblRowMask[S]); \
A=(A&BIT_SET_27) | (S<<27);\
cells_unsolved &= cl; \
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;\
wcl[U]&= cl;wcl[V]&= cl;wcl[W]&= cl;wcl[X]&= cl;}\
CompFD[I] = FD[I] = A



//================================ZH2B code
int ZH2B::Isvalid() { // usually after init 2 steps
	zh2b_g.nsol = zh2b_g.go_back = zh2b_g.modeguess = 0;
	zh2b_g.lim = 1; ComputeNext();  return zh2b_g.nsol;
}
int ZH2B::IsvalidNoUpdate(int debug ) { // usually after init 2 steps
	zh2b_g.nsol = zh2b_g.go_back = 0; zh2b_g.lim = 1; Guess();  return zh2b_g.nsol;
}
uint64_t ZH2B::GetUa() {
	zh2b_g.ua_ret = 0;
	zh2b_g.go_back = 0;
	ComputeNextFalse();
	return zh2b_g.ua_ret;
}

void ZH2B::Init_x_(GINT64 t, int n) {// cells must be 0-53
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	*this= zh2b_i;
	for (int icell = 0; icell < n; icell++) {
		int cell = t.u8[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | zh2b_g.Digit_cell_Assigned_init[i];
}
uint64_t ZH2B::Init_y_(GINT64 t, int n) {// cells must be 0-53
	BF64 Digit_cell_Assigned[9];
	memcpy(Digit_cell_Assigned, zh2b_g.Digit_cell_Assigned_init, sizeof Digit_cell_Assigned);
	Copy(zh2b_i1);
	for (int icell = 0; icell < n; icell++) {
		int cell = t.u8[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
//		if (FD[digit].Off(xcell))  return 1;// no check for not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | Digit_cell_Assigned[i];
	//	Debug(1);
	zh2b_g.ua_ret = 0;
	zh2b_g.go_back = 0;
	//ImageCandidats();
	// find now first UA false or return 0
	ComputeNextFalse();
	return zh2b_g.ua_ret;
}
void ZH2B::Init_xy(ZH2B & o, GINT64 tx, int nx, GINT64 ty, int ny) {// cells must be 0-53
	memset(zh2b_g.Digit_cell_Assigned_init, 0,
		sizeof zh2b_g.Digit_cell_Assigned_init);
	*this = o;
	for (int icell = 0; icell < nx; icell++) {
		int cell = tx.u8[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int icell = 0; icell < ny; icell++) {
		int cell = ty.u8[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | zh2b_g.Digit_cell_Assigned_init[i];
}

uint64_t ZH2B::CheckUa(uint64_t ua) {// purely debugging sequence
	BF64 Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0,
		sizeof Digit_cell_Assigned);
	Copy(zh2b_i);
	{
		register uint64_t R = 1;
		for (int cell = 0; cell < 54; cell++) {
			if (!(R&ua)) {// assign all cells not UA
				int digit = zh2b_g.puz0[cell], xcell = C_To128[cell];
				Assign(digit, cell, xcell);
				Digit_cell_Assigned[digit].Set(xcell);
			}
		}
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | Digit_cell_Assigned[i];
	zh2b_g.ua_ret = 0;
	zh2b_g.go_back = 0;
	ComputeNextFalse();
	return zh2b_g.ua_ret;
}

#define NAKED(X) 	Map=map[X];R3|=R2&Map;R2|=R1&Map;R1|=Map;
#define NAKED2(X) 	Map=map[X];R4|=R3&Map;R3|=R2&Map;R2|=R1&Map;R1|=Map;

int ZH2B::ApplySingleOrEmptyCells() {
	zh2b_g.single_applied = 0;
	unsigned _int64 * map = &FD[0].bf.u64;
	unsigned _int64 unsolved = cells_unsolved.bf.u64;
	register unsigned _int64 R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map, R3, R4;// digits 12
	Map = map[2]; R3 = R2 & Map; R2 |= R1 & Map; R1 |= Map;
	Map = map[3]; R4 = R3 & Map; R3 |= R2 & Map; R2 |= R1 & Map; R1 |= Map;
	//NAKED(3)	looking for triplets in the same search
	NAKED2(4)	NAKED2(5)	NAKED2(6)	NAKED2(7)	NAKED2(8) // digits 3-9
		if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (R1) zh2b_g.single_applied = 1;
	else {
		zh2b_g.pairs.bf.u64 = R2 & (~R3);
		zh2b_g.triplets.bf.u64 = R3 & (~R4);
		return 0;
	}
	while (R1) {// usually a very small number of cells to assign
		unsigned long res;
		if (!_BitScanForward64(&res, R1)) break;
		unsigned _int64 bit = 1; bit <<= res; // switch to the bit value
		//		cout << "naked cell" << From_128_To_81[res]<< endl;
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit
		for (int idig = 0; idig < 9; idig++) {
			if (map[idig] & bit) {// this is the digit
				//				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				//				Seta(idig, res);
				int cell = From_128_To_81[res];
				Assign(idig, cell, res);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assugn
	nextr1: {}
	}
	return 0;// not locked
}

int ZH2B::Seta(int digit, int xcell) { // single in cell
	int cell = From_128_To_81[xcell],
		block = TblBoard_Block[cell];
	if (FD[digit].Off(xcell)) return 1; // not valid
	Assign(digit, cell, xcell);
	BF64 *Fd = &FD[digit];
	*Fd &= AssignMask_Digit[cell].u64[0];
	int ddig = 6 * digit;
	if (digit > 4) ddig += 2;// second bloc of 32 bits
	rows_unsolved.Clear(ddig + C_row[cell]);//6*digit + row
	cells_unsolved.Clear(xcell);
	BF64 * RF = &FD[8];
	for (; RF >= FD; RF--)RF->Clear(xcell);
	Fd->Set(xcell); // restore bit for digit assigned
	return 0;
}



int ZH2B::Update(){
	int Shrink = 1;
	register int S, A;
	register unsigned int cl, *wcl = FD[0].bf.u32;
	while (Shrink ){
	  Shrink = 0;
	  if (!rows_unsolved.bf.u32[0])goto digit5;

	  {register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2,3,4
	  if (!(AR & 077))goto digit1;

//=digit 0
	  if (FD[0].bf.u32[0] ==  CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0,0)	if ((AR & 7) != S){
				AR &= 07777777770 | S;	UPWCL(0, 2,4,6,8,10,12,14,16)	}

digit0b:if (FD[0].bf.u32[1] ==CompFD[0].bf.u32[1])goto digit1;
		UPDN(0,1)	if (((AR >> 3) & 7) != S){
				AR &= 07777777707 | (S << 3);	UPWCL(1, 3,5,7,9,11,13,15,17)	}

digit1:	if (!(AR  & 07700))goto digit2;

		if (FD[1].bf.u32[0] ==	CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1,0)	if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);UPWCL(0, 0, 4, 6, 8, 10, 12, 14, 16)	}

digit1b:if (FD[1].bf.u32[1] ==	CompFD[1].bf.u32[1])goto digit2;
		UPDN(1,1)		if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9);UPWCL(1, 1, 5, 7, 9, 11, 13, 15, 17)	}

digit2:	if (!(AR  & 0770000))goto digit3;

		if (FD[2].bf.u32[0] ==	CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2,0)	if (((AR >> 12) & 7) != S){
			AR &= 07777707777 | (S << 12);	UPWCL(0, 0, 2, 6, 8, 10, 12, 14, 16)}

digit2b:if (FD[2].bf.u32[1] ==	CompFD[2].bf.u32[1])goto digit3;
		UPDN(2,1)	if (((AR >> 15) & 7) != S){
			AR &= 07777077777 | (S << 15);	UPWCL(1, 1, 3, 7, 9, 11, 13, 15, 17)	}

digit3: if (!(AR & 077000000))goto digit4;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		  UPDN(3,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18); UPWCL(0, 0, 2, 4, 8, 10, 12, 14, 16)	 }

digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		  UPDN(3, 1)if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21); UPWCL(1, 1, 3, 5, 9, 11, 13, 15, 17)	}

digit4:if (!(AR & 07700000000))goto end01234;

		if (FD[4].bf.u32[0] ==	CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S){
			AR &= 07077777777 | (S << 24);  UPWCL(0, 0, 2, 4, 6, 10, 12, 14, 16)	}

digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto end01234;
		UPDN(4, 1)if (((AR >> 27) & 7) != S){
			AR &= 0777777777 | (S << 27);  UPWCL(1, 1, 3, 5, 7, 11, 13, 15, 17)	}

end01234: rows_unsolved.bf.u32[0] = AR;
		}// end of validity for AR 01234


digit5:
	  if (!rows_unsolved.bf.u32[1])continue; // second lot  4 digits

	  {register unsigned int  AR = rows_unsolved.bf.u32[1];// valid for digits 5,6,7,8
	  if (!(AR & 077))goto digit6;

	  if (FD[5].bf.u32[0] ==	CompFD[5].bf.u32[0])goto digit5b;
		UPDN(5, 0) if ((AR & 7) != S){
			AR &= 07777777770 | S;	UPWCL(0, 0, 2, 4, 6, 8, 12, 14, 16)		}

digit5b:if (FD[5].bf.u32[1] == CompFD[5].bf.u32[1])goto digit6;
		UPDN(5, 1) 	if (((AR >> 3) & 7) != S){
			AR &= 07777777707 | (S << 3); UPWCL(1, 1, 3, 5, 7, 9, 13, 15, 17)	}

digit6:  if (!(AR & 07700))goto digit7;

		if (FD[6].bf.u32[0] ==  CompFD[6].bf.u32[0])goto digit6b;
		UPDN(6, 0) if (((AR >> 6) & 7) != S){
			AR &= 07777777077 | (S << 6);	UPWCL(0, 0, 2, 4, 6, 8, 10, 14, 16)	  }

digit6b: if (FD[6].bf.u32[1] == CompFD[6].bf.u32[1])goto digit7;
		UPDN(6, 1) if (((AR >> 9) & 7) != S){
			AR &= 07777770777 | (S << 9); UPWCL(1, 1, 3, 5, 7, 9, 11, 15, 17)  }

digit7:   if (!(AR  & 0770000))goto digit8;

		if (FD[7].bf.u32[0] ==			  CompFD[7].bf.u32[0])goto digit7b;
		  UPDN(7, 0)if (((AR >> 12) & 7) != S){
			  AR &= 07777707777 | (S << 12);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 16)  }

digit7b:  if (FD[7].bf.u32[1] ==  CompFD[7].bf.u32[1])goto digit8;
		  UPDN(7, 1) if (((AR >> 15) & 7) != S){
			  AR &= 07777077777 | (S << 15);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 17)  }

digit8:  if (!(AR  & 077000000))goto end5678;

		  if (FD[8].bf.u32[0] ==  CompFD[8].bf.u32[0])goto digit8b;
		  UPDN(8,0)	  if (((AR >> 18) & 7) != S){
			  AR &= 07770777777 | (S << 18);  UPWCL(0, 0, 2, 4, 6, 8, 10, 12, 14)	  }

digit8b: if (FD[8].bf.u32[1] == CompFD[8].bf.u32[1])goto end5678;
		  UPDN(8,1)	  if (((AR >> 21) & 7) != S){
			  AR &= 07707777777 | (S << 21);  UPWCL(1, 1, 3, 5, 7, 9, 11, 13, 15)  }

end5678:rows_unsolved.bf.u32[1] = AR;
	  }// end of validity for AR
  }// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
  return 1;
}
char * ZH2B::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	int tdig[9];// build the table of digits
	register int A = rows_unsolved.bf.u32[0];
	tdig[0] = A & 077; A >>= 6; tdig[1] = A & 077; A >>= 6; tdig[2] = A & 077; A >>= 6;
	tdig[3] = A & 077; A >>= 6; tdig[4] = A & 077;
	A = rows_unsolved.bf.u32[1];
	tdig[5] = A & 077; A >>= 6; tdig[6] = A & 077; A >>= 6;
	tdig[7] = A & 077; A >>= 6; tdig[8] = A & 077;

	for (int idig = 0; idig < 9; idig++) {// one digit
		int arowsj = tdig[idig];// 6 rows
		if (arowsj == 077) continue;
		for (int ib = 0; ib < 2; ib++) {// 3 blocs per digit
			int arows = (arowsj >> (3 * ib)) & 7;
			if (arows == 7) continue; // not assigned
			unsigned int band = FD[idig].bf.u32[ib];
			for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
				int row = (band >> TblMult9[j]) & 0x1ff;
				unsigned long  irow;
				_BitScanForward(&irow, row);
				int	cell = Tblstartblock[ib] + TblMult9[j] + irow;
				zs[cell] = idig + '1';
			}
		}
	}
	return zs;
}
int ZH2B::FullUpdate() {
	if (zh2b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh2b_g.single_applied) {
			continue;
		}
		break;
	}
	return 1;
}
int ZH2B::FullUpdateNoGuess() {// called if partial puzzle, just solve
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh2b_g.single_applied == 0)	break;
	}
	return 1;
}
void ZH2B::Guess() {
	if (rows_unsolved.isEmpty()) {
		if (zh2b_g.zsol && (!zh2b_g.nsol)) SetKnown(zh2b_g.zsol);// store the first solution
		zh2b_g.nsol++;
		if (zh2b_g.nsol > zh2b_g.lim) zh2b_g.go_back = 1;
		return;
	}

	unsigned long res;
	int ndig = 2;
	register unsigned _int64 R3 = zh2b_g.pairs.bf.u64;
	if (_BitScanForward64(&res, R3))goto gogo; // first pair is ok to go
	ndig = 3;
	R3 = zh2b_g.triplets.bf.u64;
	if (_BitScanForward64(&res, R3))goto gogo; // first triplet is ok to go
	R3 = cells_unsolved.bf.u64;
	if (!_BitScanForward64(&res, R3))return; // should not be
	ndig = 6;//unknown, do it  always in a new area


gogo:
	int cell = From_128_To_81[res];
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig].On(res)) {
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->Seta(idig, res);
			mynext->ComputeNext();
		}
}
/*

the guess sequence should select first cells in the "X" band
or may be the "x" band shloud be always band 1 to have cell order in the same
sequence as the guess priority
Target, be forced to refresh the table as much as possble for a new "x"
*/
void ZH2B::GuessFalse() {// try to find a false
	if (rows_unsolved.isEmpty()) {// stop at first false
		register uint64_t Rs = 0, Rw;
		for (int i = 0; i < 9; i++) {
			Rw = FD[i].bf.u64;
			Rw &= (~zh2b_g.digsols[i]);
			Rs |= Rw;
		}
		if (Rs) {
			zh2b_g.ua_ret = Rs;
			zh2b_g.go_back = 1;
		}
		return;
	}
	unsigned long res;
	register unsigned _int64 R3 = zh2b_g.pairs.bf.u64;
	if (_BitScanForward64(&res, R3))goto gogo; // first pair is ok to go
	R3 = zh2b_g.triplets.bf.u64;
	if (_BitScanForward64(&res, R3))goto gogo; // first triplet is ok to go
	R3 = cells_unsolved.bf.u64;
	if (!_BitScanForward64(&res, R3))return; // should not be

gogo:
	int cell = From_128_To_81[res], digit = zh2b_g.puz0[cell];
	if (FD[digit].On(res)) {// start with the good digit if possible
		ZH2B * mynext = this + 1; // start next guess
		mynext->Copy(*this);
		mynext->Seta(digit, res);
		mynext->ComputeNextFalse();
		if (zh2b_g.go_back) return;
	}

	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;// finish with the ok
		if (FD[idig].On(res)) {// one valid digit		
			if (0)cout << "brute force guess  cell " << cell
				<< " dig" << idig + 1 << endl;
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->Seta(idig, res);
			mynext->ComputeNextFalse();
			if (zh2b_g.go_back) return;
		}
	}
}
/*
//==================== ZH2B ua collector code
void ZH2B::InitGenUas(int * zpuz) {// so far nothing to do except checking the index (done)
	genuasb12.Initgen(zpuz);// maximum is the size -1 for ztemp ; ntemp
	zh2b_g.test = 0;
}
void ZH2B::InitUas2() { genuasb12.InitUas2(); }
void ZH2B::CollectFinal(TUA64 & mt) {
	register uint64_t b1 = 0xffffffff, b2 = b1 << 32;
	for (int i = 4; i <= genuasb12.limsize; i++) {
		if (genuasb12.nftemp[i])
			for (int j = 0; j < genuasb12.nftemp[i]; j++) {
				register uint64_t w = genuasb12.zftemp[i][j];
				if ((!(w&b1)) || (!(w&b2)))continue;	// forget if purely band 1 or band2
				mt.wt = w;
				mt.AddUA();
			}
	}
}
int ZH2B::CollectFinalua2s(uint64_t *td, int maxt, int n0) {// must clear subsets
	// done for a specific
	int n = n0, istart;
	for (int i = 0; i <= genuasb12.limsize; i++) {
		istart = n;
		if (genuasb12.nftemp[i])
			for (int j = 0; j < genuasb12.nftemp[i]; j++) {
				if (n >= maxt)return maxt;
				register uint64_t w = genuasb12.zftemp[i][j], wn = ~w;
				for (int k = n0; k < istart; k++) {// clear subsets
					if (!(wn&td[k])) goto next;// subset found
				}
				td[n++] = w;
			next: {}
			}
	}
	memset(genuasb12.nftemp, 0, sizeof genuasb12.nftemp);
	return n;
}
void ZH2B::GenUas(int floors) {
	genuasb12.OrderFloors(floors);// first relabel sol to have floors as lower digits
	// active floors are now located in pseudo digits 0-...
	// init brute force for active floors to empty
	__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
	cells_unsolved = 0;
	for (int i = 0; i < genuasb12.nfloors; i++) {
		cells_unsolved |= genuasb12.ordered_digsols[i];
		FD[i] = genuasb12.ordered_digpm[i];
	}
	rows_unsolved.bf.u64 = maskLSB[6 * genuasb12.nfloors].u64[0];
	for (int i = 0; i < genuasb12.nfloors; i++) {
		FD[i] &= cells_unsolved;
	}
}
void ZH2B::GenUAs2_minirowpair(int col1, int col2, int dig1, int dig2, int i_81) {
	genuasb12.OrderFloorsPair(dig1, dig2);
	__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
	cells_unsolved = 0;
	for (int i = 0; i < 2; i++) {
		cells_unsolved |= genuasb12.ordered_digsols[i];
		FD[i] = genuasb12.ordered_digpm[i];
	}
	rows_unsolved.bf.u64 = maskLSB[6 * 2].u64[0];
	for (int i = 0; i < 2; i++) {
		FD[i] &= cells_unsolved;
	}
	// add and assign now col1 dig1 and col2 dig2
	// where was col1 dig2 and col2 dig1
	unsigned long cx1, cx2;
	{
		register uint64_t RC = Zhoucol, RU = cells_unsolved.bf.u64;
		RC |= RC << 32; // now col1 six rows
		register uint64_t Rw = (RC << col1)&RU;// one cell
		_BitScanForward64(&cx1, Rw); // now cx1 is the cell mode 2x32
		Rw = (RC << col2)&RU;
		Rw &= FD[0].bf.u64; // must be one cell and only one
		_BitScanForward64(&cx2, Rw); // now cx1 is the cell mode 2x32
	}

	FD[0] |= (uint64_t)1 << cx1;
	FD[1] |= (uint64_t)1 << cx2;
	Seta(0, cx1); Seta(1, cx2);
	if (Update2()) {
		if (!(rows_unsolved.bf.u32[0] & 07777)) {
			register uint64_t w0 = ~genuasb12.ordered_digsols[0];
			w0 &= FD[0].bf.u64;
			register uint64_t w1 = ~genuasb12.ordered_digsols[1];
			w1 &= FD[1].bf.u64;
			w0 |= w1;
			int n = (int)_popcnt64(w0);
			genuasb12.AddUA(w0, n);
		}
		else Guess2();
		genuasb12.Collectua2s(i_81);
	}
	//=====================================================
	//cout << "try now 3 digits" << endl;
	int flbase = (1 << dig1) | (1 << dig2);
	genuasb12.nfloors = 3;
	genuasb12.uax[3] = 0;
	{// next step is to try all 2+1 digits patterns
		for (int i = 0; i < 84; i++) {// try now all 4 digits
			int fl = floors_3d[i], last_floor = fl ^ flbase;
			unsigned long  diglast;
			if (_popcnt32(last_floor) != 1)continue;// must include the 2 digits in minirow
			genuasb12.BuildTsubsets3(i);
			//cout << "try digits 0" << oct << fl << dec << endl;
			_BitScanForward(&diglast, last_floor);// this is the last dig
			genuasb12.ordered_digpm[2] = genuasb12.digpm[diglast];
			genuasb12.ordered_digsols[2] = genuasb12.digsols[diglast];


			__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
			cells_unsolved = 0;
			for (int i = 0; i < 3; i++) {
				cells_unsolved |= genuasb12.ordered_digsols[i];
				FD[i] = genuasb12.ordered_digpm[i];
			}
			rows_unsolved.bf.u64 = maskLSB[6 * 3].u64[0];
			for (int i = 0; i < 3; i++) 	FD[i] &= cells_unsolved;
			//ImageCandidats();
			// add and assign now col1 dig1 and col2 dig2
			// where was col1 dig2/diglast and col2 dig1/diglast
			unsigned long cx1, cx2;
			{
				register uint64_t RC = Zhoucol, RU = cells_unsolved.bf.u64;
				RC |= RC << 32; // now col1 six rows
				register uint64_t Rw1 = (RC << col1)&RU;
				ZH2B rzh = *this;
				while (_BitScanForward64(&cx1, Rw1)) { // can be one or 2 cells
					Rw1 ^= (uint64_t)1 << cx1;// clear bit
					register uint64_t Rw2 = (RC << col2)&RU;
					while (_BitScanForward64(&cx2, Rw2)) {// can be one or 2 cells
						Rw2 ^= (uint64_t)1 << cx2; // clear bit
						Copy(rzh);// restaure the start point
						FD[0] |= (uint64_t)1 << cx1;
						FD[1] |= (uint64_t)1 << cx2;
						Seta(0, cx1); Seta(1, cx2);
						if (Update3()) Guess3();
						genuasb12.Collectua2s(i_81);
					}
				}
			}
		}
		//==================================
		genuasb12.nfloors = 4;
		genuasb12.uax[4] = 0;
		{// next step is to try all 2+1 digits patterns
			for (int i = 0; i < 126; i++) {// try now all 4 digits
				int fl = floors_4d[i], last_floors = fl ^ flbase;
				unsigned long dig3, dig4;
				if (_popcnt32(last_floors) != 2)continue;// must include the 2 digits in minirow
				genuasb12.BuildTsubsets4(i);
				//cout << "try digits 0" << oct << fl << dec << endl;
				_BitScanForward(&dig3, last_floors);
				genuasb12.ordered_digpm[2] = genuasb12.digpm[dig3];
				genuasb12.ordered_digsols[2] = genuasb12.digsols[dig3];
				last_floors ^= 1 << dig3;
				_BitScanForward(&dig4, last_floors);
				genuasb12.ordered_digpm[3] = genuasb12.digpm[dig4];
				genuasb12.ordered_digsols[3] = genuasb12.digsols[dig4];

				__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
				cells_unsolved = 0;
				for (int i = 0; i < 4; i++) {
					cells_unsolved |= genuasb12.ordered_digsols[i];
					FD[i] = genuasb12.ordered_digpm[i];
				}
				rows_unsolved.bf.u64 = maskLSB[6 * 4].u64[0];
				for (int i = 0; i < 4; i++) 	FD[i] &= cells_unsolved;
				//ImageCandidats();
				// add and assign now col1 dig1 and col2 dig2
				// where was col1 dig2/diglast and col2 dig1/diglast
				unsigned long cx1, cx2;
				{
					register uint64_t RC = Zhoucol, RU = cells_unsolved.bf.u64;
					RC |= RC << 32; // now col1 six rows
					register uint64_t Rw1 = (RC << col1)&RU;
					ZH2B rzh = *this;
					while (_BitScanForward64(&cx1, Rw1)) { // can be one or 2 cells
						Rw1 ^= (uint64_t)1 << cx1;// clear bit
						register uint64_t Rw2 = (RC << col2)&RU;
						while (_BitScanForward64(&cx2, Rw2)) {// can be one or 2 cells
							Rw2 ^= (uint64_t)1 << cx2;// clear bit
							Copy(rzh);// restaure the start point
							FD[0] |= (uint64_t)1 << cx1;
							FD[1] |= (uint64_t)1 << cx2;
							Seta(0, cx1); Seta(1, cx2);
							if (Update4()) Guess4();
							genuasb12.Collectua2s(i_81);
						}
					}
				}
			}
		}
	}

}
void ZH2B::GenUAs2_minirowtriplet(int col1, int col2, int col3, int dig1, int dig2, int dig3, int i_81) {
	int flbase = (1 << dig1) | (1 << dig2) | (1 << dig3);
	for (int i = 0; i < 84; i++) if (floors_3d[i] == flbase) {
		genuasb12.BuildTsubsets3(i); break;
	}
	genuasb12.OrderFloorsTriplet(dig1, dig2, dig3);
	__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
	cells_unsolved = 0;
	for (int i = 0; i < 3; i++) {
		cells_unsolved |= genuasb12.ordered_digsols[i];
		FD[i] = genuasb12.ordered_digpm[i];
	}
	rows_unsolved.bf.u64 = maskLSB[6 * 3].u64[0];
	for (int i = 0; i < 3; i++) {
		FD[i] &= cells_unsolved;
	}
	unsigned long cx1, cx2, cx3;
	{
		register uint64_t RC = Zhoucol, RU = cells_unsolved.bf.u64;
		RC |= RC << 32; // now col1 six rows
		register uint64_t Rw1 = (RC << col1)&RU;
		ZH2B rzh = *this;
		while (_BitScanForward64(&cx1, Rw1)) { // can be one or 2 cells
			Rw1 ^= (uint64_t)1 << cx1;// clear bit
			register uint64_t Rw2 = (RC << col2)&RU;
			while (_BitScanForward64(&cx2, Rw2)) {// can be one or 2 cells
				Rw2 ^= (uint64_t)1 << cx2;// clear bit
				register uint64_t Rw3 = (RC << col3)&RU;
				while (_BitScanForward64(&cx3, Rw3)) {
					Rw3 ^= (uint64_t)1 << cx3;// clear bit
					Copy(rzh);// restaure the start point
					FD[0] |= (uint64_t)1 << cx1;
					FD[1] |= (uint64_t)1 << cx2;
					FD[2] |= (uint64_t)1 << cx3;
					Seta(0, cx1); Seta(1, cx2); Seta(2, cx3);
					if (Update3()) {
						if (!(rows_unsolved.bf.u32[0] & 0777777)) {
							register uint64_t w0 = ~genuasb12.ordered_digsols[0];
							w0 &= FD[0].bf.u64;
							register uint64_t w1 = ~genuasb12.ordered_digsols[1];
							w1 &= FD[1].bf.u64;
							w0 |= w1;
							w1 = ~genuasb12.ordered_digsols[2];
							w1 &= FD[2].bf.u64;
							w0 |= w1;
							int n = (int)_popcnt64(w0);
							genuasb12.AddUA(w0, n);
						}
						else Guess3();
					}
				}
			}
		}
	}
	genuasb12.Collectua2s(i_81);
	//================================ try one more floor
	//if (1) return;
	genuasb12.nfloors = 4;
	genuasb12.uax[4] = 0;
	{// next step is to try all 2+1 digits patterns
		for (int i = 0; i < 126; i++) {// try now all 4 digits
			int fl = floors_4d[i], last_floors = fl ^ flbase;
			unsigned long dig4;
			if (_popcnt32(last_floors) != 1)continue;// must include the 2 digits in minirow
			genuasb12.BuildTsubsets4(i);
			//cout << "try digits 0" << oct << fl << dec << endl;
			_BitScanForward(&dig4, last_floors);
			genuasb12.ordered_digpm[3] = genuasb12.digpm[dig4];
			genuasb12.ordered_digsols[3] = genuasb12.digsols[dig4];

			__stosq(&FD[0].bf.u64, 0, 9); __stosq(&CompFD[0].bf.u64, 0, 9);
			cells_unsolved = 0;
			for (int i = 0; i < 4; i++) {
				cells_unsolved |= genuasb12.ordered_digsols[i];
				FD[i] = genuasb12.ordered_digpm[i];
			}
			rows_unsolved.bf.u64 = maskLSB[6 * 4].u64[0];
			for (int i = 0; i < 4; i++) 	FD[i] &= cells_unsolved;
			//ImageCandidats();
			unsigned long cx1, cx2, cx3;
			{
				register uint64_t RC = Zhoucol, RU = cells_unsolved.bf.u64;
				RC |= RC << 32; // now col1 six rows
				register uint64_t Rw1 = (RC << col1) &RU;
				ZH2B rzh = *this;
				while (_BitScanForward64(&cx1, Rw1)) {
					Rw1 ^= (uint64_t)1 << cx1;// clear bit
					register uint64_t Rw2 = (RC << col2) & RU;
					while (_BitScanForward64(&cx2, Rw2)) {
						Rw2 ^= (uint64_t)1 << cx2; // clear bit
						register uint64_t Rw3 = (RC << col3) & RU;
						while (_BitScanForward64(&cx3, Rw3)) {
							Rw3 ^= (uint64_t)1 << cx3;// clear bit
							Copy(rzh);// restaure the start point
							FD[0] |= (uint64_t)1 << cx1;
							FD[1] |= (uint64_t)1 << cx2;
							FD[2] |= (uint64_t)1 << cx3;
							Seta(0, cx1); Seta(1, cx2); Seta(2, cx3);
							if (Update4()) {
								if (!(rows_unsolved.bf.u32[0] & 077777777)) {
									register uint64_t w1 = ~genuasb12.ordered_digsols[3];
									w1 &= FD[3].bf.u64;
									if (w1) {
										register uint64_t w0 = ~genuasb12.ordered_digsols[0];
										w0 &= FD[0].bf.u64;
										w0 |= w1;
										w1 = ~genuasb12.ordered_digsols[1];
										w1 &= FD[1].bf.u64;
										w0 |= w1;
										w1 = ~genuasb12.ordered_digsols[2];
										w1 &= FD[2].bf.u64;
										w0 |= w1;
										int n = (int)_popcnt64(w0);
										genuasb12.AddUA(w0, n);
									}
								}
								else Guess4();
							}
						}

					}
				}
			}
			genuasb12.Collectua2s(i_81);
		}
	}

}

void ZH2B::GenUas2() {
	genuasb12.ntsubsets = 0;
	genuasb12.uax[2] = 0;
	for (int i = 0; i < 36; i++) {// try now all 2 digits
		GenUas(floors_2d[i]);
		if (Update2()) {
			if (!(rows_unsolved.bf.u32[0] & 077)) continue;
			Guess2();
			genuasb12.Collect(genuasb12.u2d[i], genuasb12.nnu2[i]);
		}
	}

}
void ZH2B::GenUas3() {
	genuasb12.uax[3] = 0;
	for (int i = 0; i < 84; i++) {// try now all 2 digits
		genuasb12.BuildTsubsets3(i);
		GenUas(floors_3d[i]);
		if (Update3()) {
			if (!(rows_unsolved.bf.u32[0] & 0777777)) continue;
			Guess3();
			genuasb12.Collect(genuasb12.u3d[i], genuasb12.nnu3[i]);

		}
	}
}
void ZH2B::GenUas4() {
	genuasb12.uax[4] = 0;
	for (int i = 0; i < 126; i++) {// try now all 4 digits
		genuasb12.BuildTsubsets4(i);
		GenUas(floors_4d[i]);
		if (Update4()) {
			if (!(rows_unsolved.bf.u32[0] & 077777777)) continue;
			Guess4();
			genuasb12.Collect(genuasb12.u4d[i], genuasb12.nnu4[i]);
		}
	}
}
void ZH2B::GenUas5() {
	zh2b_g.test = 0;
	genuasb12.limstep = 16;
	for (int i = 0; i < 126; i++) {// try now all54 digits
		int floors = 0x1ff ^ floors_4d[i];
		genuasb12.BuildTsubsets5(floors);
		GenUas(floors);
		Update5();			Guess5();
		genuasb12.Collect5();
	}
}

#define UPWCL2(I,P)cl = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[I] &= cl;\
wcl[P]&= cl;

#define UPWCL3(I,P,Q)cl = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[I] &= cl;\
wcl[P]&= cl;wcl[Q]&= cl;

#define UPWCL4(I,P,Q,R)cl = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[I] &= cl;\
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;

#define UPWCL5(I,P,Q,R,T)cl = ~(A & TblRowMask[S]);\
cells_unsolved.bf.u32[I] &= cl;\
wcl[P]&= cl;wcl[Q]&= cl;wcl[R]&= cl;wcl[T]&= cl;


int ZH2B::Update2() {
	register unsigned int Shrink = 1, S, A, cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2
		if (!(AR & 07777)) return 1;
		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL2(0, 2)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL2(1, 3)
		}
	digit1:	if (!(AR & 07700))goto loop;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL2(0, 0)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto loop;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL2(1, 1)
		}


	loop: rows_unsolved.bf.u32[0] = AR;
	}// end while
	return 1;
}
void ZH2B::Guess2() {
	register uint32_t Rr = rows_unsolved.bf.u32[0] & 077;
	// note if all 1 are settled, the puzzle is solved
	if (Rr) {// select next unsolved row for digit 1 and try each possible cell
		unsigned long row;
		_BitScanForward(&row, Rr);
		int  band = C_box[row], r_row = C_box_rel[row];
		int x = (FD[0].bf.u32[band] >> TblMult9[r_row]) & 0x1ff;
		while (x) {
			unsigned long rcell;
			_BitScanForward(&rcell, x);
			int  drow = 27 * band + TblMult9[r_row];
			x ^= (1 << rcell);
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaC(0, rcell + drow);
			if (mynext->Update2())	mynext->Guess2();
		}
		return;
	}
	// all digit 0 solved this is a solution
	register uint64_t w0 = ~genuasb12.ordered_digsols[0];
	w0 &= FD[0].bf.u64;
	if (!w0) return;
	register uint64_t w1 = ~genuasb12.ordered_digsols[1];
	w1 &= FD[1].bf.u64;
	if (!w1) return;
	w0 |= w1;
	if (genuasb12.nfloors == 2) {
		int n = (int)_popcnt64(w0);
		genuasb12.AddUA(w0, n);
		return;
	}
	w0 |= genuasb12.uax[2];
	int n = (int)_popcnt64(w0);
	if (n > genuasb12.limstep) return;
	genuasb12.AddUA(w0, n);
}

int ZH2B::Update3() {
	register unsigned int Shrink = 1, S, A, cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2
		if (!(AR & 0777777)) return 1;
		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL3(0, 2, 4)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL3(1, 3, 5)
		}
	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL3(0, 0, 4)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL3(1, 1, 5)
		}
	digit2:	if (!(AR & 0770000))goto loop;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL3(0, 0, 2)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto loop;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL3(1, 1, 3)
		}


	loop: rows_unsolved.bf.u32[0] = AR;
	}// end while
	return 1;
}
void ZH2B::Guess3() {// three floors empty at start solve floor 3
	register uint32_t Rr = (rows_unsolved.bf.u32[0] >> 12) & 077;
	// note if all 1 are settled, the puzzle is solved
	if (Rr) {// select next unsolved row for digit 2 and try each possible cell
		unsigned long row;
		_BitScanForward(&row, Rr);
		int  band = C_box[row], r_row = C_box_rel[row];
		int x = (FD[2].bf.u32[band] >> TblMult9[r_row]) & 0x1ff;
		while (x) {
			unsigned long rcell;
			_BitScanForward(&rcell, x);
			int  drow = 27 * band + TblMult9[r_row];
			x ^= (1 << rcell);
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaC(2, rcell + drow);// digit 2
			if (mynext->Update3())	mynext->Guess3();
		}
		return;
	}
	{
		// all digit 3 solved go to guess2 if not a zero deviation
		register uint64_t w0 = ~genuasb12.ordered_digsols[2];
		w0 &= FD[2].bf.u64;
		if (!w0) return;
		w0 |= genuasb12.uax[3];
		if (!(rows_unsolved.bf.u32[0] & 0777777)) {
			if (FD[0].bf.u64 == genuasb12.ordered_digsols[0] ||
				FD[1].bf.u64 == genuasb12.ordered_digsols[1]) return;
		}
		genuasb12.uax[2] = w0;
	}
	Guess2();
}

int ZH2B::Update4() {
	register unsigned int Shrink = 1, S, A, cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2
		if (!(AR & 077777777)) return 1;
		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL4(0, 2, 4, 6)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL4(1, 3, 5, 7)
		}
	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL4(0, 0, 4, 6)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL4(1, 1, 5, 7)
		}
	digit2:	if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL4(0, 0, 2, 6)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL4(1, 1, 3, 7)
		}

	digit3: if (!(AR & 077000000))goto loop;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL4(0, 0, 2, 4)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto loop;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL4(1, 1, 3, 5)
		}

	loop: rows_unsolved.bf.u32[0] = AR;
	}// end while
	return 1;
}
void ZH2B::Guess4() {// three floors empty at start solve floor 3
	register uint32_t Rr = (rows_unsolved.bf.u32[0] >> 18) & 077;
	// note if all 1 are settled, the puzzle is solved
	if (Rr) {// select next unsolved row for digit 3 and try each possible cell
		unsigned long row;
		_BitScanForward(&row, Rr);
		int  band = C_box[row], r_row = C_box_rel[row];
		int x = (FD[3].bf.u32[band] >> TblMult9[r_row]) & 0x1ff;
		while (x) {
			unsigned long rcell;
			_BitScanForward(&rcell, x);
			int  drow = 27 * band + TblMult9[r_row];
			x ^= (1 << rcell);
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaC(3, rcell + drow);// digit 3
			if (mynext->Update4())	mynext->Guess4();
			//  if not  game locked ?? continue;
		}
		return;
	}
	{
		// all digit 3 solved go to guess3 if not a zero deviation
		register uint64_t w0 = ~genuasb12.ordered_digsols[3];
		w0 &= FD[3].bf.u64;
		if (!w0) return;
		w0 |= genuasb12.uax[4];
		if (!(rows_unsolved.bf.u32[0] & 077777777)) {
			if (FD[0].bf.u64 == genuasb12.ordered_digsols[0] ||
				FD[1].bf.u64 == genuasb12.ordered_digsols[1] ||
				FD[2].bf.u64 == genuasb12.ordered_digsols[2]) return;
		}
		genuasb12.uax[3] = w0;
	}
	Guess3();
}

int ZH2B::Update5() {
	register unsigned int Shrink = 1, S, A, cl, *wcl = FD[0].bf.u32;
	while (Shrink) {
		Shrink = 0;
		register unsigned int  AR = rows_unsolved.bf.u32[0];// valid for digits 0,1,2
		if (!(AR & 077777777)) return 1;
		//=digit 0
		if (FD[0].bf.u32[0] == CompFD[0].bf.u32[0])goto digit0b;
		UPDN(0, 0)	if ((AR & 7) != S) {
			AR &= 07777777770 | S;	UPWCL5(0, 2, 4, 6, 8)
		}

	digit0b:if (FD[0].bf.u32[1] == CompFD[0].bf.u32[1])goto digit1;
		UPDN(0, 1)	if (((AR >> 3) & 7) != S) {
			AR &= 07777777707 | (S << 3);	UPWCL5(1, 3, 5, 7, 9)
		}
	digit1:	if (!(AR & 07700))goto digit2;

		if (FD[1].bf.u32[0] == CompFD[1].bf.u32[0])goto digit1b;
		UPDN(1, 0)	if (((AR >> 6) & 7) != S) {
			AR &= 07777777077 | (S << 6); UPWCL5(0, 0, 4, 6, 8)
		}

	digit1b:if (FD[1].bf.u32[1] == CompFD[1].bf.u32[1])goto digit2;
		UPDN(1, 1)		if (((AR >> 9) & 7) != S) {
			AR &= 07777770777 | (S << 9); UPWCL5(1, 1, 5, 7, 9)
		}
	digit2:	if (!(AR & 0770000))goto digit3;

		if (FD[2].bf.u32[0] == CompFD[2].bf.u32[0])goto digit2b;
		UPDN(2, 0)	if (((AR >> 12) & 7) != S) {
			AR &= 07777707777 | (S << 12);	UPWCL5(0, 0, 2, 6, 8)
		}

	digit2b:if (FD[2].bf.u32[1] == CompFD[2].bf.u32[1])goto digit3;
		UPDN(2, 1)	if (((AR >> 15) & 7) != S) {
			AR &= 07777077777 | (S << 15);	UPWCL5(1, 1, 3, 7, 9)
		}

	digit3: if (!(AR & 077000000))goto loop;

		if (FD[3].bf.u32[0] == CompFD[3].bf.u32[0])goto digit3b;
		UPDN(3, 0)	  if (((AR >> 18) & 7) != S) {
			AR &= 07770777777 | (S << 18); UPWCL5(0, 0, 2, 4, 8)
		}

	digit3b:  if (FD[3].bf.u32[1] == CompFD[3].bf.u32[1])goto digit4;
		UPDN(3, 1)if (((AR >> 21) & 7) != S) {
			AR &= 07707777777 | (S << 21); UPWCL5(1, 1, 3, 5, 9)
		}
	digit4:if (!(AR & 07700000000))goto loop;

		if (FD[4].bf.u32[0] == CompFD[4].bf.u32[0])goto digit4b;
		UPDN(4, 0)if (((AR >> 24) & 7) != S) {
			AR &= 07077777777 | (S << 24);  UPWCL5(0, 0, 2, 4, 6)
		}

	digit4b:if (FD[4].bf.u32[1] == CompFD[4].bf.u32[1])goto loop;
		UPDN(4, 1)if (((AR >> 27) & 7) != S) {
			AR &= 0777777777 | (S << 27);  UPWCL5(1, 1, 3, 5, 7)
		}

	loop: rows_unsolved.bf.u32[0] = AR;
	}// end while
	return 1;
}
void ZH2B::Guess5() {// three floors empty at start solve floor 3.0
	register uint32_t Rr = (rows_unsolved.bf.u32[0] >> 24) & 077;
	// note if all 1 are settled, the puzzle is solved
	if (Rr) {// select next unsolved row for digit 3 and try each possible cell
		unsigned long row;
		_BitScanForward(&row, Rr);
		int  band = C_box[row], r_row = C_box_rel[row];
		int x = (FD[4].bf.u32[band] >> TblMult9[r_row]) & 0x1ff;
		while (x) {
			unsigned long rcell;
			_BitScanForward(&rcell, x);
			int  drow = 27 * band + TblMult9[r_row];
			x ^= (1 << rcell);
			ZH2B * mynext = this + 1; // start next guess
			mynext->Copy(*this);
			mynext->SetaC(4, rcell + drow);// digit 4
			if (mynext->Update5())	mynext->Guess5();
		}
		return;
	}
	{
		//  digit 4 solved go to guess4 if not a zero deviation
		register uint64_t w0 = ~genuasb12.ordered_digsols[4];
		w0 &= FD[4].bf.u64;
		if (!w0) return;
		int n = (int)_popcnt64(w0);
		if (n > 4) return;
		if (!(rows_unsolved.bf.u32[0] & 07777777777)) {
			if (FD[0].bf.u64 == genuasb12.ordered_digsols[0] ||
				FD[1].bf.u64 == genuasb12.ordered_digsols[1] ||
				FD[2].bf.u64 == genuasb12.ordered_digsols[2] ||
				FD[3].bf.u64 == genuasb12.ordered_digsols[3]) return;
		}
		genuasb12.uax[4] = w0;
	}
	Guess4();
}*/

void ZH2B::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh2b_g.nsol << " unsolved=" << Unsolved_Count() << endl;
	//	cout << zh1b_g.out27 << " band1 "<<endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_6 digits " << endl;
	if (!all) return;

	cout << "map per digit bands 2 3" << endl;
	for (int ib = 0; ib < 2; ib++) {
		for (int ir = 0; ir < 3; ir++) {
			for (int idig = 0; idig < 9; idig++) {
				unsigned vf = FD[idig].bf.u32[ib];
				unsigned wr = (vf >> (9 * ir)) & 0x1ff;
				for (int k = 0; k < 9; k++) {
					if (wr & (1 << k))		cout << idig + 1;
					else 		cout << ".";
					if (k == 2 || k == 5) 	cout << " ";
				}
				cout << "  ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
void ZH2B::DebugDigit(int digit) {
	cout << "DEBUG  digit=" << digit + 1 << endl;
	for (int ib = 0; ib < 2; ib++) {
		for (int ir = 0; ir < 3; ir++) {
			unsigned vf = FD[digit].bf.u32[ib];
			unsigned  wr = (vf >> (9 * ir)) & 0x1ff;
			for (int k = 0; k < 9; k++) {
				if (wr & (1 << k))		cout << digit + 1;
				else 		cout << ".";
				if (k == 2 || k == 5) 	cout << " ";
			}
			cout << endl; //end of row
		}
		cout << endl; // end of block
	}
	cout << endl; // end of map per digit

}
int ZH2B::GetAllDigits(int cell) {
	int ir = 0, xcell = C_To128[cell];;
	for (int i = 0; i < 9; i++) if (FD[i].On(xcell))ir |= (1 << i);
	return ir;
}
void ZH2B::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 54; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // attention ici i indice colonne
		lcol[i] = 2;    // 2  mini tous chiffres imposés
		for (j = 0; j < 6; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 6; i++) { // maintenant indice ligne
		if ((i == 3) || (i == 6)) {
			for (int ix = 0; ix < (tcol + 10); ix++)       cout << (char)'-';
			cout << endl;
		}
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < 9; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}

int  ZH2B::DebugCheckUa(uint64_t ua) {
	Init_std_bands();
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	// assign all cells not ua
	{
		register uint64_t Rua = ua;
		for (int i = 0; i < 54; i++) {
			int xi = C_To128[i];
			uint64_t bit = (uint64_t)1 << xi;
			if (Rua&bit) continue;
			int digit = zh2b_g.puz0[i];
			Assign(digit, i, xi);
			zh2b_g.Digit_cell_Assigned_init[digit].Set(xi);
		}
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | zh2b_g.Digit_cell_Assigned_init[i];
	//ImageCandidats();
	return Isvalid();
}
/*
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	*this= zh2b_i;
	for (int icell = 0; icell < n; icell++) {
		int cell = t.u8[icell], digit = zh2b_g.puz0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved | zh2b_g.Digit_cell_Assigned_init[i];
void ZH2B::Guess2() {	// note if all 1 are settled, the puzzle is solved
	int v = FD[0], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown rows
	//cout << Char27out(v) << " v guess2" << endl;
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZH2B * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Set2(cell);
		if (mynext->UpdateDigit(0)) {// update digit 1/0
			if (mynext->FD[0] >> 27) {
				mynext->Guess2();
			}
			else {// all digit 0 solved this is a solution
				//cout << Char27out(mynext->FD[0]) << " solution digit 0 " << endl;
				register uint32_t w0 = mynext->FD[0] & zh1b_g.fdsw[2][0],
					w1 = mynext->FD[1] & zh1b_g.fdsw[2][1] & mynext->cells_unsolved;
				if ((!w0) || (!w1)) continue;
				w0 |= (w1 | zh1b_g.previous_ua_status[0]);
				w0 |= _popcnt32(w0) << 27;
				if (zh1b_g.diag)
					cout << Char27out(w0) << " ua added guess2" << endl;
				zh1b_g.AddUA(w0);
			}
		}
	}
}


*/


//======================== ZH2B1B one digit 2 bands table of solutions
int ZH2B_1D::GetSols( int ru) {
	CompFD.bf.u64 = 0;
	FD = zh2b_g.mystart &zh2b_g.myandsol;
	//cout << Char2Xout(FD.bf.u64) << " start get sols" << endl;
	zh2b_g.nsolw = 0;
	//zh2b_g.andsol.bf.u64 = BIT_SET_2X;
	ComputeNext(ru);
	return zh2b_g.nsolw;
}
void ZH2B_1D::ComputeNext(int ru) {
	if (Update(ru)) {
		if (!ru) {// new sol put it in table
			if (FD != zh2b_g.mysol) {
				//cout << Char2Xout(FD.bf.u64) << "store it" << endl;
				//cout << Char2Xout((FD - zh2b_g.mysol).bf.u64) << "store ua" << endl;

				zh2b_g.tuaw[zh2b_g.nsolw] = FD- zh2b_g.mysol;// partial ua
				zh2b_g.tsolw[zh2b_g.nsolw++] = FD;// partial solution
				zh2b_g.myandsol &= FD;
			}
		}
		else Guess(ru);
	}
}
int ZH2B_1D::Update(int &ru) {// solve as much as possible 
	//cout << "update ru=0" <<oct<< ru <<dec<< endl;
	int Shrink = 1;
	register int S, A;
	while (Shrink) {
		Shrink = 0;
		if (!(ru & 7)) goto band2;
		A = FD.bf.u32[0];
		if (A == CompFD.bf.u32[0]) goto band2;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[1] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[0] = 
			FD.bf.u32[0] = A;
		ru = (ru & 070) | S;
	band2:;
		if (!(ru & 070)) goto exit;
		A = FD.bf.u32[1];
		if (A == CompFD.bf.u32[1]) goto exit;
		Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		FD.bf.u32[0] &= TblMaskSingle[S];
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		CompFD.bf.u32[1] = 
			FD.bf.u32[1] = A;
		ru = (ru & 7) | (S<<3);
	exit:;
	}
	//cout << Char2Xout(FD.bf.u64) << " sortie update ru=" 
	//	<< oct << ru << dec << endl;

	return 1;
}
void ZH2B_1D::Guess(int ru) {//
	//cout << Char2Xout(FD.bf.u64) << "guess ru =0" << oct << ru << dec << endl;
	int dcell = 0,v,ruw=ru;
	{
		register uint32_t a = FD.bf.u32[0], b = FD.bf.u32[1];
		if (!(ru & 7)) goto guess_b2;
		if ((ru & 070) && (_popcnt32(b) < _popcnt32(a))) goto guess_b2;
		// guess in band1
		ruw = (-ru) &ru;
		v = FD.bf.u32[0] & TblRowUnsolved[ruw];// unknown in last unknown row
		goto gogo;
	guess_b2:	// guess in band2
		dcell = 27;
		ruw=ru >> 3;
		ruw = (-ruw) &ruw;
		v = FD.bf.u32[0] & TblRowUnsolved[ruw];// unknown in last unknown row
		ruw <<= 3; //relocate bit to the right place
	}
gogo:;
	uint32_t cell;
	ru ^= ruw;// clear the row for next steps
	while (bitscanforward(cell, v)) {
		v ^= 1 << cell;
		cell += dcell;
		ZH2B_1D * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(cell);
		mynext->ComputeNext(ru);
	}
}

int ZH2B_1D::IsValid(uint64_t v) {
	CompFD.bf.u64 = 0;
	FD.bf.u64 = v;
	int ru = 077;
	return Update(ru);
}



//================================================
//===================== now ZHone code 
ZHONE_GLOBAL::ZHONE_GLOBAL() {
	zsol  = 0; // no solution unless required buy the user
}
void ZHONE_GLOBAL::GetBand(int * b, uint32_t * t) {
	band0 = b; tua = t;  nua = 0;
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = b[cell];
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}
}

void ZHONE_GLOBAL::AddUA(uint32_t ua) {// add if and clear supersets
	// ua 27 bit + 5 bit length
	if (nua >= 100) return; // size given in STD_B416 go_17_bands
	//if(diag)cout << "add in tua " << tua << " nua="<< nua<< endl;
	register uint32_t ua27 = ua & BIT_SET_27;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint32_t R = tua[iua];
		if (R < ua) {// is it subset
			R&= BIT_SET_27;
			if ((R&ua27) == R) {
				//if (diag)cout << Char27out(R) << " subset of (iua=)"<<iua << endl
				//	<< Char27out(ua27) << endl;
				return;// we have a subset
			}
		}
		else {
			for (uint32_t jua = nua; jua > iua; jua--)tua[jua] = tua[jua - 1];// to insert the new
			tua[iua] = ua;// inserted
			nua++;
			//if (diag)cout << Char27out(ua27) << " inserted nua=" << nua << endl;
			// is it a subset of a previous entry
			for (iua++; iua < nua; iua++) {
				if ((tua[iua] &ua27) == ua27) {// we have a subset
					//if (diag)cout << Char27out(tua[iua]) << " superset cleared nua=" << nua << endl;
					for (uint32_t k = iua + 1; k < nua; k++)tua[k - 1] = tua[k];
					nua--;
					iua--; //continue same position
				}
			}
			return;
		}
	}
	tua[nua++] = ua;// added
	//if (diag) {
		//cout << Char27out(ua27) << " added nua=" << nua << endl;
		//cout << Char27out(tua[nua-1]) << " added stored"  << endl;
	//}
}
void ZHONE_GLOBAL::FindMissingUAs() {
	struct SPB3 {// spots to find band 3 missing uas
		int  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[12], *s3, *sn3;
	s3 = spb3;
	memset(s3, 0, sizeof spb3[0]);
	s3->active_cells = BIT_SET_27;// all cells active
	s3->possible_cells = tua[0]& BIT_SET_27;
	int tcells[10], ispot;
	cout << "search nua=" << nua << endl;
	//____________________  here start the search
next:
	uint32_t cell;
	ispot = (int)(s3 - spb3); // be sure to have always the right one
	if (ispot > 10) return;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = (uint64_t)1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		s3->active_cells ^= bit;
		sn3 = s3 + 1;
		*sn3 = *s3;
		sn3->all_previous_cells |= bit;
		register int filter = sn3->all_previous_cells,
			ac = s3->active_cells;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nua; i++) {
			if (tua[i] & filter)continue;
			//cout << Char27out(tua[i]) << "next ua i=" <<i << endl;
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}
		// no more ua is it valid or take more uas
		zhone[0].AddMissingUAs(tcells, ispot + 1);
		goto next;
	}
	// going back, for a non empty index, count it back
back:
	//if (1) return;
	if (--s3 >= spb3)goto next;

}
void ZHONE_GLOBAL::ValidPuzzle(uint32_t * sol) {//depending on type
	//cout << "entry valid puzzle" << endl;
	switch (type) {
	case 1: {// mode add uas
		uint32_t or_sol = 0;
		for (int i = 0; i < 9; i++) {// 9 digits
			or_sol |= sol[i] & ~(fd_sols[0][i]);// invalid cells
		}
		if (or_sol) {
			//cout << Char27out(or_sol) << " ua found "<< _popcnt32(or_sol) << endl;
			AddUA(or_sol | _popcnt32(or_sol) << 27);
		}
	}
		break;
	}
	return;
}
/* old code for a valid puzzle
		if (zh1b_g.type == 2) {//compare to start puzzle
			memcpy(zh1b_g.fd_sols[1], FD, sizeof FD); //store the solution
			for (int i = 0; i < 9; i++) {
				if (FD[i] - zh1b_g.fd_sols[0][i]) {
					zh1b_g.nsol = 1;// force exit after that one
					break;
				}
			}
		}
		else if (zh1b_g.type == 1) {
			if (zh1b_g.nsol < 2) memcpy(zh1b_g.fd_sols[zh1b_g.nsol], FD, sizeof FD);
		}
*/

void ZHONE_GLOBAL::PrintTua() {
	cout << "band status" << endl;
	for (uint32_t i = 0; i < nua; i++) {
		uint32_t ua = tua[i], ua27 = ua & BIT_SET_27, nd = ua >> 27;
		cout << i << "\t" << Char27out(ua27) << "\t" << nd << endl;
	}
}
//============================================ ZHONE
void ZHONE::InitOne_std_band() {//init after zh1b_g getband
	cells_unsolved = BIT_SET_27;
	memcpy(FD, zh1b_g.fd_sols[1], sizeof FD);
	for (int i = 0; i < 9; i++)	FD[i] |= 7 << 27;//set unknown rows
	zh1b_g.ndigits = 9;
}
int ZHONE::InitSudokux(GINT * t, int n) {// cells must be 0-26
	int Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
	*this=zhone_i;
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell].u8[0], digit = t[icell].u8[1], bit = 1 << cell;
		if ((FD[digit] & bit) == 0)  return 1;// check not valid entry
		Assign(digit, cell);
		Digit_cell_Assigned[digit] |= bit;
	}
	for (int i = 0; i < 9; i++)  FD[i] &=
		cells_unsolved |
		Digit_cell_Assigned[i] |
		07000000000;// don't touch rows unsolved

//	Debug(1);
//	ImageCandidats();
	return 0;
}
int ZHONE::Isvalid() { // usually after init 2 steps
	zh1b_g.InitIsvalid();
	ComputeNext();
	return zh1b_g.nsol;
}

#define NAKED(X) 	Map=map[X];R3|=R2&Map;R2|=R1&Map;R1|=Map;
int ZHONE::ApplySingleOrEmptyCells() {
	zh1b_g.single_applied = 0;
	uint32_t * map = FD, unsolved = cells_unsolved;
	register int R2 = map[0] & map[1],
		R1 = (map[0] | map[1]), Map, R3;// digits 12
	Map = map[2]; R3 = R2 & Map; R2 |= R1 & Map; R1 |= Map;
	NAKED(3) NAKED(4) NAKED(5) NAKED(6)	NAKED(7) NAKED(8) // digits 3-9
		if (unsolved & (~R1)) return 1; // locked
	R1 &= ~R2;
	R1 &= unsolved; // forget solved seen as singles
	if (R1) zh1b_g.single_applied = 1;
	else {
		zh1b_g.pairs = R2 & (~R3);
		zh1b_g.triplets = R3;
		return 0;
	}
	while (R1) {// usually a very small number of cells to assign
		unsigned long res;
		if (!_BitScanForward(&res, R1)) break;
		int cell = res, bit = 1 << cell; // switch to the bit value
		//		cout << "naked cell" << From_128_To_81[res]<< endl;
		R1 &= ~bit;  // clear the bit
		// call Seta(int digit, int xcell) so find digit 
		for (int idig = 0; idig < 9; idig++) {
			if (map[idig] & bit) {// this is the digit
				//				if (FD[idig].Off(res))  return 1; // invalid, gane locked
				//				Seta(idig, res);
				int cell = res;
				Assign(idig, cell);
				goto nextr1;// finished for that cell
			}
		}
		return 1; //conflict with a previous cell assign
	nextr1: {}
	}
	return 0;// not locked 
}
void ZHONE::Seta(int digit, int cell) { // single in cell
	Assign(digit, cell);
	register  int  bit= 1 << cell,cl = ~bit;
	FD[0]&=cl; FD[1] &= cl; FD[2] &= cl; FD[3] &= cl; FD[4] &= cl;
	FD[5] &= cl; FD[6] &= cl; FD[7] &= cl; FD[8] &= cl;
	FD[digit] |= bit;// restore digit
}

int ZHONE::Update() {
	int Shrink = 1;
	register int S, A, cl;
	register uint32_t *wcl = FD;
	while (Shrink) {
		Shrink = 0;
		if ((A = FD[0]) - CompFD[0]) { UPDN1(0, 1, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[1]) - CompFD[1]) { UPDN1(1, 0, 2, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[2]) - CompFD[2]) { UPDN1(2, 0, 1, 3, 4, 5, 6, 7, 8); }
		if ((A = FD[3]) - CompFD[3]) { UPDN1(3, 0, 1, 2, 4, 5, 6, 7, 8); }
		if ((A = FD[4]) - CompFD[4]) { UPDN1(4, 0, 1, 2, 3, 5, 6, 7, 8); }
		if ((A = FD[5]) - CompFD[5]) { UPDN1(5, 0, 1, 2, 3, 4, 6, 7, 8); }
		if ((A = FD[6]) - CompFD[6]) { UPDN1(6, 0, 1, 2, 3, 4, 5, 7, 8); }
		if ((A = FD[7]) - CompFD[7]) { UPDN1(7, 0, 1, 2, 3, 4, 5, 6, 8); }
		if ((A = FD[8]) - CompFD[8]) { UPDN1(8, 0, 1, 2, 3, 4, 5, 6, 7); }
		//	  Debug(1);
	}// end while
#ifdef DIAG
	cout << "end update cycle" << endl;
	Debug(1);
	ImageCandidats();
#endif
	return 1;
}
void ZHONE::Guess() {
	if (!cells_unsolved) {
		if (zh1b_g.type) {
			//ImageCandidats();
			zh1b_g.ValidPuzzle(FD); return;	}
		else if (zh1b_g.zsol && (!zh1b_g.nsol)) SetKnown(zh1b_g.zsol);// store the first solution
		zh1b_g.nsol++;
		if (zh1b_g.nsol > zh1b_g.lim) zh1b_g.go_back = 1;
		return;
	}
	uint32_t res;
	int ndig = 2;
	register int R3 = zh1b_g.pairs;
	if (bitscanforward(res, R3))goto gogo; // first pair is ok to go
	ndig = 3;
	R3 = zh1b_g.triplets;
	bitscanforward(res, R3); // first triplet is ok to go
gogo:
	int cell = res, bit = 1 << cell;
	//cout << "brute force guess  cell " << cellsFixedData[cell].pt << endl;
	for (int idig = 0; idig < 9; idig++) {
		if (FD[idig] & bit) {// one valid digit		
			if (--ndig) {
				ZHONE * mynext = this + 1; // start next guess
				*mynext=*this;
				mynext->Seta(idig, res);
				mynext->ComputeNext();
				if (zh1b_g.go_back) return;
			}
			else {// this is the last do it in the same place
				Seta(idig, res);
				ComputeNext();
				return;
			}
		}
	}
}
int ZHONE::FullUpdate() {
	if (zh1b_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (ApplySingleOrEmptyCells())			return 0; // locked empty cell or conflict singles in cells
		if (zh1b_g.single_applied) {
			//			cout << "after singles" << endl; Debug();
			continue;
		}
		break;
	}
	return 1;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
char * ZHONE::SetKnown(char * zs) {
	strcpy(zs, &empty_puzzle[27]);
	zs[27] = 0;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {// one digit
		int arows = FD[idig] >> 27;// 3 rows
		if (arows == 7) continue;
		int band = FD[idig];
		for (int j = 0; j < 3; j++) if (!(arows & (1 << j))) {
			int row = (band >> TblMult9[j]) & 0x1ff;
			uint32_t  irow;
			bitscanforward(irow, row);
			int	cell = TblMult9[j] + irow;
			zs[cell] = idig + '1';
		}
	}
	return zs;
}
void ZHONE::Debug(int all) {
	cout << "DEBUG  nbsol=" << zh1b_g.nsol << " unsolved=" << Unsolved_Count() 
		<< " index to zhone="<<this-zhone<< endl;
	//cout << zh1b_g.out27 << " band1 " << endl;
	char zi[82];  SetKnown(zi);
	cout << zi << " known rows 1_3 digits " << endl;
	if (!all) return;

	cout << "map per digit one band" << endl;
	for (int ir = 0; ir < 3; ir++) {
		for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++) {
			int vf = FD[idig], wr = (vf >> (9 * ir)) & 0x1ff;
			for (int k = 0; k < 9; k++) {
				if (wr & (1 << k))		cout << idig + 1;
				else 		cout << ".";
				if (k == 2 || k == 5) 	cout << " ";
			}
			cout << "  ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of map per digit

}
void ZHONE::DebugDigit(int digit) {
	cout << "DEBUG  digit=" << digit + 1 << endl;
	for (int ir = 0; ir < 3; ir++) {
		int vf = FD[digit], wr = (vf >> (9 * ir)) & 0x1ff;
		for (int k = 0; k < 9; k++) {
			if (wr & (1 << k))		cout << digit + 1;
			else 		cout << ".";
			if (k == 2 || k == 5) 	cout << " ";
		}
		cout << endl; //end of row
	}
	cout << endl; // end of block

}
int ZHONE::GetAllDigits(int cell) {
	int ir = 0, bit = 1 << cell;
	for (int i = 0; i < (int)zh1b_g.ndigits; i++) if (FD[i] & bit)ir |= (1 << i);
	return ir;
}
void ZHONE::ImageCandidats() {
	int dig_cells[81]; for (int i = 0; i < 27; i++) dig_cells[i] = GetAllDigits(i);
	int i, j, l, lcol[9], tcol = 0, ncand = 0;
	cout << "PM map " << endl << endl;
	for (i = 0; i < 9; i++) {  // column
		lcol[i] = 2;    // 2  mini 
		for (j = 0; j < 3; j++) {
			l = _popcnt32(dig_cells[9 * j + i]);
			if (l > lcol[i])       lcol[i] = l;
		}
		tcol += lcol[i];
	}
	for (i = 0; i < 9; i++) {
		if ((i == 3) || (i == 6))cout << "|";
		cout << (char)('A' + i) << Blancs(lcol[i], 1);
	}
	cout << endl;
	for (i = 0; i < 3; i++) { // now row 
		for (j = 0; j < 9; j++) {
			if ((j == 3) || (j == 6))cout << "|";
			int cell = 9 * i + j, digs = dig_cells[cell], 
				ndigs = _popcnt32(digs);
			ncand += ndigs;
			for (int id = 0; id < (int)zh1b_g.ndigits; id++)if (digs & (1 << id))
				cout << id + 1;
			cout << Blancs(lcol[j] + 1 - ndigs, 1);
		} // end for j
		cout << endl;
	} // end for i
	cout << endl;

}


//============== collector mode
void ZHONE::Checkstart() {
	zh1b_g.ndigits = 9;
	memcpy(FD, zh1b_g.fd_sols[1], sizeof FD);
	cells_unsolved = BIT_SET_27 ;
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	ImageCandidats();
}
int ZHONE::Start_nFloors(int floors) {
	uint32_t solved_cells = 0,nd=0;
	for (int idig = 0,bit=1; idig < 9; idig++,bit<<=1) {
		if (floors & bit) {
			zh1b_g.fdsw[0][nd]= zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig];
		}
		else solved_cells|= zh1b_g.fd_sols[0][idig];
	}
	zh1b_g.ndigits = nd;
	cells_unsolved=BIT_SET_27^ solved_cells;
	for (uint32_t i = 0; i < nd; i++) {
		FD[i] &= cells_unsolved;
		if (_popcnt32(FD[i]) == 3)return 1;
	}
	return 0;
}
void ZHONE::Start_Uas_Mini(int floors, int floors_mini_row) {
	//cout << "start uas mini lfoors0" << oct << floors << " digits0" << floors_mini_row << dec << endl;
	zh1b_g.floors_mini_row = floors_mini_row; //used to check gangster changes 
	uint32_t solved_cells = 0, nd = 0;
	// first assigned digits
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (!(floors & bit))  solved_cells |= zh1b_g.fd_sols[0][idig];
	}
	cells_unsolved = BIT_SET_27 ^ solved_cells;
	// then digits not exchanged
	floors &= ~floors_mini_row;
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (floors & bit) {
			zh1b_g.fdsw[0][nd] = zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig]& cells_unsolved;
		}
	}
	floors =floors_mini_row;// then digits  exchanged
	for (int idig = 0, bit = 1; idig < 9; idig++, bit <<= 1) {
		if (floors & bit) {
			zh1b_g.digmap[idig] = nd; // used in gangster changes
			zh1b_g.fdsw[0][nd] = zh1b_g.fd_sols[0][idig];
			zh1b_g.fdsw[2][nd] = (~zh1b_g.fdsw[0][nd]) & BIT_SET_27;
			FD[nd++] = zh1b_g.fd_sols[1][idig] & cells_unsolved;
		}
	}
	zh1b_g.ndigits = nd;
	//cout << "end start uas nd=" << nd << endl;
	//ImageCandidats();
}
void ZHONE::ApplyGangsterChanges(int * g0, int * g1) {
	for (int i = 0,col=Zhoucol; i < 9; i++,col<<=1) {
		if (g0[i] == g1[i])continue;
		int changes = g0[i] ^ g1[i]; // one added one cleared
		// safety temp control, must be 2 digits exchanged
		if (changes & (~zh1b_g.floors_mini_row))return;
		//if (zh1b_g.diag > 1)cout << "apply change gang go i=" << i 
		//	<<" change 0"<<oct<<changes<<dec<< endl;
		for (int d = 0,bit=1; d < 9; d++,bit<<=1) {// check digits
			if (!(changes & bit)) continue;
			int digit = zh1b_g.digmap[d];// digit rank in the brute force
			if (g0[i] & bit)FD[digit] &= ~col; 
			else FD[digit] |= col & cells_unsolved;
		}
	}
}
void ZHONE::InitGuess() {// after start floors before guess
	for (uint32_t idig = 0; idig < zh1b_g.ndigits; idig++)
		FD[idig] |= 7 << 27;
	memset(zh1b_g.previous_ua_status, 0,
		sizeof zh1b_g.previous_ua_status); //no previous status
}
int ZHONE::UpdateDigit(int digit) {
	register int S, A = FD[digit], cl;
	if (A  != CompFD[digit]) {
		int Shrink = (TblShrinkMask[A & 0x1FF] |
			TblShrinkMask[(A >> 9) & 0x1FF] << 3 |
			TblShrinkMask[(A >> 18) & 0x1FF] << 6);
		if ((A &= TblComplexMask[Shrink]) == 0)  return 0;
		S = ((A | (A >> 9) | (A >> 18)) & 0x1FF);
		S = TblRowUniq[TblShrinkSingle[Shrink] & TblColumnSingle[S]];
		if ((A >> 27) != S) {
			cl = ~(A & TblRowMask[S]);
			A = (A&BIT_SET_27) | (S << 27);
			cells_unsolved &= cl;
		}
		FD[digit] = CompFD[digit] = A;
	}// end if
	return 1;
}
void ZHONE::Set2(int cell) { // single in cell
	Assign(0, cell);
	register  int  cl = ~(1 << cell);
	FD[1] &= cl;
}
void ZHONE::Guess2() {	// note if all 1 are settled, the puzzle is solved
	int v = FD[0], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown rows
	//cout << Char27out(v) << " v guess2" << endl;
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Set2(cell);
		if (mynext->UpdateDigit(0)) {// update digit 1/0
			if (mynext->FD[0] >> 27) {
				mynext->Guess2();
			}
			else {// all digit 0 solved this is a solution
				//cout << Char27out(mynext->FD[0]) << " solution digit 0 " << endl;
				register uint32_t w0 = mynext->FD[0] & zh1b_g.fdsw[2][0],
					w1 = mynext->FD[1] & zh1b_g.fdsw[2][1]& mynext->cells_unsolved;
				if ((!w0) || (!w1)) continue;
				w0 |= (w1  | 	zh1b_g.previous_ua_status[0]);
				w0 |= _popcnt32(w0) << 27;
				if (zh1b_g.diag)
					cout << Char27out(w0) << " ua added guess2 + count "<<(w0>>27)  << endl;
				zh1b_g.AddUA(w0);
			}
		}
	}
}
void ZHONE::Guess3() {	// solve digit 3/2
	int v = FD[2], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(2, cell);
		if (mynext->UpdateDigit(2)) {// update digit 3/2
			if (mynext->FD[2] >> 27) {
				mynext->Guess3();
			}
			else {// all digit 3 solved go to guess2 if not a zero deviation
				if (zh1b_g.diag > 1)
					cout << Char27out(mynext->FD[2]) << " guess2 " << endl;
				uint32_t w3 = mynext->FD[2] & zh1b_g.fdsw[2][2];
				if (!w3) continue;
				zh1b_g.previous_ua_status[0] = w3 | zh1b_g.previous_ua_status[1];
				mynext->FD[1] &= (mynext->cells_unsolved|7<<27);
				mynext->FD[0] &= (mynext->cells_unsolved | 7 << 27);
				mynext->Guess2();
			}
		}
	}
}
void ZHONE::Guess4() {	//  solve digit 4/3
	int v = FD[3], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(3, cell);
		if (mynext->UpdateDigit(3)) {// update digit 4/3
			if (mynext->FD[3] >> 27) {
				mynext->Guess4();
			}
			else {// all digit 4 solved go to guess3 if not a zero deviation
				if (zh1b_g.diag > 1)
					cout << Char27out(mynext->FD[3]) << " guess3 " << endl;
				uint32_t w4 = mynext->FD[3] & zh1b_g.fdsw[2][3];
				if (!w4) continue;
				zh1b_g.previous_ua_status[1] = w4 | zh1b_g.previous_ua_status[2];
				mynext->FD[2] &= (mynext->cells_unsolved | 7 << 27);
				mynext->Guess3();
			}
		}
	}
}
void ZHONE::Guess5() {	//  solve digit 5/4
	int v = FD[4], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(4, cell);
		if (mynext->UpdateDigit(4)) {// update digit 4/3
			if (mynext->FD[4] >> 27) {
				mynext->Guess5();
			}
			else {// all digit 5 solved go to guess4 if not a zero deviation
				if (zh1b_g.diag > 1)
					cout << Char27out(mynext->FD[4]) << " guess5 " << endl;
				uint32_t w5 = mynext->FD[4] & zh1b_g.fdsw[2][4];
				if (!w5) continue;
				zh1b_g.previous_ua_status[2] = w5 | zh1b_g.previous_ua_status[3];
				mynext->FD[3] &= (mynext->cells_unsolved | 7 << 27);
				mynext->Guess4();
			}
		}
	}
}
void ZHONE::Guess6() {	//  solve digit 6/5
	int v = FD[5], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(5, cell);
		if (mynext->UpdateDigit(5)) {// update digit 4/3
			if (mynext->FD[5] >> 27) {
				mynext->Guess6();
			}
			else {// all digit 6 solved go to guess5 if not a zero deviation
				if (zh1b_g.diag>1)
					cout << Char27out(mynext->FD[5]) << " guess6 " << endl;
				uint32_t w6 = mynext->FD[5] & zh1b_g.fdsw[2][5];
				if (!w6) continue;
				zh1b_g.previous_ua_status[3] = w6 | zh1b_g.previous_ua_status[4];
				mynext->FD[4] &= (mynext->cells_unsolved | 7 << 27);
				mynext->Guess5();
			}
		}
	}
}
void ZHONE::Guess7() {	//  solve digit 7/6
	int v = FD[6], vr = v >> 27;
	vr = (-vr)&vr; // catch last bit
	v &= TblRowUnsolved[vr];// unknown in last unknown row
	uint32_t cell;
	while (bitscanforward(cell, v)) {// first cell
		v ^= 1 << cell; // clear bit
		ZHONE * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->Assign(6, cell);
		if (mynext->UpdateDigit(6)) {// update digit 4/3
			if (mynext->FD[6] >> 27) {
				mynext->Guess7();
			}
			else {// all digit 7 solved go to guess6 if not a zero deviation
				uint32_t w7 = mynext->FD[6] & zh1b_g.fdsw[2][6];
				if (!w7) continue;
				zh1b_g.previous_ua_status[4] = w7 | zh1b_g.previous_ua_status[5];
				mynext->FD[5] &= (mynext->cells_unsolved | 7 << 27);
				mynext->Guess6();
			}
		}
	}
}
void ZHONE::AddMissingUAs(int * tcells, int ncells) {
	zh1b_g.type = 1;
	zh1b_g.ndigits = 9;
	zh1b_g.go_back = 0;
	int Digit_cell_Assigned[9];
	memset(Digit_cell_Assigned, 0, sizeof Digit_cell_Assigned);
	memset(this, 0, sizeof this);
	InitOne_std_band();
	for (int icell = 0; icell < ncells; icell++) {
		int cell = tcells[icell],
			digit = zh1b_g.band0[cell],
			bit = 1 << cell;
		if ((FD[digit] & bit) == 0)  return;// check not valid entry
		Assign(digit, cell);
		Digit_cell_Assigned[digit] |= bit;
	}
	for (int i = 0; i < 9; i++)  FD[i] &=
		cells_unsolved |
		Digit_cell_Assigned[i] |
		07000000000;// don't touch rows unsolved
	//Debug(1);
	//ImageCandidats();
	if (FullUpdate()) Guess();
}

