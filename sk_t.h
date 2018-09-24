// all these tables are used in many functions as
#pragma once
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

//platform-specific settings
#ifdef _MSC_VER
inline char * stpcpy(char * d, const char * o)
	{strcpy(d,o); return (d+strlen(d));}
#else
//stpcpy is part of POSIX: http://pubs.opengroup.org/onlinepubs/9699919799/functions/stpcpy.html
#endif

#ifndef _MSC_VER
//assume every compiler but MS is C99 compliant and has inttypes
#include <inttypes.h>

#else
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;
#endif


#ifdef   _MSC_VER
#include <intrin.h>
#else
#include <immintrin.h>
#endif
#include <limits.h>

#ifdef   _MSC_VER
#define _popcnt64(a) __popcnt64(a)
#define _popcnt32(a) __popcnt(a)
//__movsb is builtin
//__movsd is builtin
//__movsq is builtin
//__stosd is builtin
//__stosq is builtin
//_bittestandset64 is builtin
//_bittestandreset64 is builtin
//_bittest64 is builtin
//_BitScanForward64 is builtin
//_BitScanForward is builtin
//_BitScanReverse64 is builtin
//_BitScanReverse is builtin
//strcpy_s is implemented
//strncpy_s is implemented
//_timeb has the same fields as timeb
//_ftime64_s does the same as ftime64
#else
#define _popcnt64(a) __builtin_popcountll(a)
#define _popcnt32(a) __builtin_popcount(a)
#define __movsb(dst,src,count) memcpy(dst,src,count)
#define __movsd(dst,src,count) memcpy(dst,src,count*4)
#define __movsq(dst,src,count) memcpy(dst,src,count*8)
#define __stosd(dst, c, N) \
   __asm__ __volatile__( \
       "rep stosl %%eax, (%%rdi)\n\t" \
       : : "D"((dst)), "a"((c)), "c"((N)) : "memory");
//#define __stosq(a,b,c) memset(a,b,c*8)
#define __stosq(dst, c, N) \
   __asm__ __volatile__( \
       "rep stosq %%rax, (%%rdi)\n\t" \
       : : "D"((dst)), "a"((c)), "c"((N)) : "memory");
#define _bittestandset64(dest,offset) (offset < 64 ? ((uint64_t*)dest)[0] |= ((uint64_t)1 << offset) : ((uint64_t*)dest)[1] |= ((uint64_t)1 << (offset - 64)))
#define _bittestandreset64(dest,offset) (offset < 64 ? ((uint64_t*)dest)[0] &= !((uint64_t)1 << offset) : ((uint64_t*)dest)[1] &= !((uint64_t)1 << (offset - 64)))
#define _bittest64(a, b) (((*((uint64_t*)a)) >> (b)) & 1)
#define _BitScanForward64(res, src) (*res = __builtin_ctzll(src))
#define _BitScanForward(res, src) (*res = __builtin_ctz(src))
#define _BitScanReverse64(res, src) (*res = __builtin_ffsll(src))
#define _BitScanReverse(res, src) (*res = __builtin_ffs(src))
//strcpy_s isn't implemented
#define strcpy_s(dest, size, src) (strncpy(dest, src, size))
//strncpy_s isn't implemented
#define strncpy_s(dest, size, src, n) {strncpy(dest, src, size < n ? size : n); dest[n] = 0;}
#define _timeb timeb
#define _ftime64_s ftime
#endif

#ifndef __AVX2__
#define _blsi_u32(v) ((v | (v - 1)) ^ (v - 1))
#endif


//=======
#define BIT_SET_27         0777777777
#define BIT_SET_30         07777777777
#define BIT_SET_64         0xffffffffffffffff
#define Zhoucol 01001001
#define Zhoubox 07007007

typedef unsigned char byte;
typedef unsigned short word;
typedef unsigned short int  USHORT;
typedef unsigned int UINT;
typedef unsigned long ULONG;
typedef unsigned char UCHAR;


// candidate 7bits cells + digit * 128 
typedef unsigned short int  UCAND;
typedef unsigned short int  SCAND;// UCAND + status off bit 12 set to 1
typedef unsigned short int  RCAND;// UCAND + (<<11)  unit or 27 if cell


typedef union GINT16_t {
	uint16_t   u16;
	uint8_t    u8[2];
} GINT16;

typedef union GINT_t {
	uint32_t   u32;
	uint16_t   u16[2];
	uint8_t    u8[4];
} GINT;

// in GINT mode a candidate is {cell + digit<<8}   or {cell=u8[0],digit=u8[1}}

typedef union GINT64_t {
	uint64_t   u64;
	uint32_t   u32[2];
	uint16_t   u16[4];
	uint8_t    u8[8];
} GINT64;



typedef union T128 {// new definition closer to GINTx
	uint64_t    u64[2];
	uint8_t     u8[16];
	uint16_t    u16[8];
	uint32_t    u32[4];
	__m128i		u128;
} T128;

typedef union p9x9 {
	char puz[82];
	char gr[9][9];
} p9_9;
//==============================tables and functions  in tab0 and tab0b

extern char *  empty_puzzle;
extern char *  puzstart;

// printing and debugging 
extern char * Blancs(int n, int no);
extern void Coutg9(char * zl, int endl = 1);
extern void CoutGrid(char * zc);
extern char * Char9out(int w);
extern char * Char27out(int w);
extern char * Char54out(uint64_t v);
extern char * Char64out(uint64_t v);
extern char * Char2Xout(uint64_t v);
extern char * CoutGintPuzzle(GINT * t, int n);
extern char * CoutGint64Puzzle(GINT * t, int n);
// general correspondance row column box band pattern 

extern int TblMult3[9];  // 3*i
extern int TblMult9[9];  // 9*i
extern int C_minirow[81];
extern int * Tbldiv3; //i/3
extern int C_row[81];
extern int * Tblofftorow;  // i/9 cell to row
extern int C_col[81];
extern int * TblMod9; //  i%9 cell to column
extern int C_box[81];
extern int * Tblofftobox; // cell in band to relative box in band
extern int C_stack[81];
extern int C_box_rel[81]; // relative position in the box
extern int * TblCelltoBox;

extern int  C_transpose_c[81];
extern int  C_transpose_d[81];
extern int  C_transpose_d2[81];
extern int  C_transpose_r[81];// rotate col 9 to row 1
extern int  C_transpose_rm[81]; // rotate row1 to col9

extern int  C_mod27[81];
extern int * Tbloffset; // i%27  giving cell to offset in a band
extern int C_div27[81];
extern int * TblBoard_Block; //  i/27  giving cell -> band 
extern int C_To128[81];
extern int From_128_To_81[128];
extern int box_col_to_row[512];// Flash transpose a box
extern int tband_box[6], tband_row[3], tband_box_skrink[3];

extern byte cellsInGroup[27][9];
extern byte miniline[54][3]; // row order
extern byte box_miniline[54][3]; // box order
extern uint64_t tbl_2x_cols[9]; // column bit field in 2x mode(2x27 in 2x32)
extern uint64_t tbl_2x_stacks[3];
extern uint64_t tbl_2x_stacksn[3];
extern int tperm6[6][3];
extern int tperm3[3][3]; // here all boxes are once in pos 1

class BUILDSTRING;
#include "sk_bitfields.h"
//  here included a merged version of fsss2 and Gp 128 bits base


//extern const t_128 bitSet[128];
extern const T128 maskLSB[129]; // not in fss2
extern const T128 maskffff; // not in fss2
extern T128 cellsInBandBM[6];
extern T128 cellsInHouseBM[27];
extern T128 band3xBM[6];
extern T128 units3xBM[27];

//!Small class containing the permanent data for a cell
class CELL_FIX {
public:
	USHORT
		i8,		///< cell index (0-80)
		el,     ///< row index (0-8)
		pl,     ///< column index (0-8)
		plu,     /// same 9 17
		eb,     ///< box index(0-8) 
		ebu,     ///< box index(18_��) 
		pb;     ///< relative position in the box (0-8)
	char pt[5];	///< printing string like R4C9 with \0
	T128 z;     ///< list of the  20 cells controled by a cell in a bit field

	int ObjCommun(const CELL_FIX *p8) const {
		return((el == p8->el) || (pl == p8->pl) || eb == p8->eb);
	}
	inline void GetRegions(int * tt){ tt[0] = el; tt[1] = plu; tt[2] = ebu; }
	int GetTableRegions(USHORT * tt, CELL_FIX & cell2);
};
extern  CELL_FIX  cellsFixedData[81];
extern T128 cell_z3x[81];
extern int C_rbc27[81];
/* COMBINE is a small module giving combination
for a given size   (maths C n;p )

the goal is to have a simpler code, not to improve the performance

*/
class COMBINE
{
public:
	UINT inds[15],  // internal table size p
		p,      // maxi 15
		lim,    // set to p-1 final index at the end
		n;      // numer of objects in the table
	USHORT * entry, *output;
	void Sort(); // internal function formating the output

public:
	void First(int ne, int pe, USHORT * d, USHORT * f); // initial
	int Next(); // get next group of index
};

//========================== various old struct or UA collection
struct PUZ {
	p9x9 gg;
	void OutGrid();
};

struct TUA64 {// table outside the struct UA handling 
	uint64_t * t, wt;
	int nua, tsize;
	inline void Init(uint64_t * e_t, int n) { t = e_t; tsize = n; nua = 0; }
	void AddUA(int debug = 0);
	void Print(int start = 0, int end = 0);
	void Print(char * lib);
};


class GG {
public:
	char g[9][9],	///<The grid
		filler,		///<a null to terminate 81 char string
		*pg;		///<a pointer to the beginning of the string

	GG(){ pg = g[0]; }		// constructor

	inline void Copie(const GG & ge) {
		strcpy_s(pg, 82, ge.pg);
	}

	inline void Copie(const char * pge) {
		strcpy_s(pg, 82, pge);
	}

	inline void operator =(const GG &ge) {
		strcpy_s(pg, 82, ge.pg);
	}
	int NBlancs() const;
	int Nfix() const;
	//void Image(FSR * EE, char * lib) const;
};

/* class VV9 is a small class
extracting a row or a column   from a GG described puzzle

to see later another use of VV9
*/
class VV9{
public: char v[10];
		VV9(){ v[9] = 0; };
		//   void base(){for(USHORT i=0;i<9;i++)v[i]=(char)i;}
		void row(char * puz, USHORT r){
			for (USHORT i = 0, ii = 9 * r; i<9; i++)
				v[i] = puz[ii++];
		}
		void col(char * puz, USHORT c) {
			for (USHORT i = 0, ii = c; i<9; i++, ii += 9)
				v[i] = puz[ii];
		}

};

//================= symmetry of given 

extern USHORT sh_36[36][2];
extern USHORT sv_36[36][2];
extern USHORT sd1_36[36][2];
extern USHORT sd2_36[36][2];
extern USHORT sst_36[36][2];
extern USHORT sc_40[40][2];
extern USHORT sr90_20[20][4];
extern USHORT sym_81[5][81];
extern USHORT sym_f9[3][9];
extern USHORT  sym_tcor[3][9];
