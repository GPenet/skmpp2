

#define GTEST17_ON 1
#define UALIMSIZE 20
#define GUALIMSIZE 18
#define UA32_10 0xffc00000
#define UA64_54 0x3fffffffffffff
#define TUA64_12SIZE 2000
//============================================== 

#define MAXN5 51520
#define MAXN6 237770 

#define G17MORESIZE 32
#include "go_17sol_tables.h"
#include "Zh1b2b.h"  // brute force 2 bands  
extern ZH2B_GLOBAL   zh2b_g;
extern ZH2B5_GLOBAL   zh2b5_g;   // 2_5 digits 2 bands
extern ZH2B_1D_GLOBAL zh2b1d_g;  // one digit 2 bands

extern ZH2B zh2b[40], zh2b_i, zh2b_i1;
extern ZH2B5 zh2b5[10]; // solved digit per digit
extern ZH2B_1D zh2b1d[6]; // maxi 6 guesses, one per row
extern ZHONE_GLOBAL   zh1b_g;
extern ZHONE zhone[20];
extern ZHONE zhone_i;
const char * libs_c17_00_cpt2g[40] = {
	"0 bands1+2 processed entries M10",//0
	"1 total bands 3",//1
	"2 no more ua b2 ",//2
	"3 valid brute force b2 ",//3
	"4 max band3 ",//4
	"5 7 clues b2",//5
	"6 8 clues b2 ",//6
	"7 expand 7",//7
	"8 final check",//8
	"9 valid after final check",//9
	"10 critical band3",//10
	"11 not critical 1",//11
	"12 not critical 234",//12
	"13 not critical >4",//13
	"14 using socket 3",//14
	"15 entry band3 handler excluding critical+ua outfield",//15
	"16 entry critical ",//16
	"17 add 1 from active",//17
	"18 n uas at start",//18
	"19 n gua2s at start  ",//19
	"20 n gua3s at start  ",//20
	"21 n sockets2",//21
	"22 n sockets1",//22
	"23 max bands 3",//23
	"24/4 min", "25/5 min", "26/6 min","27/7 min","28/8 min","29/9 min","30/10 min",

};

// standard first band (or unique band)
struct STD_B416 {
	char band[28];
	int band0[27], i416, gangster[9], map[27], dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i);
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas();
	void InitC10(int i);
	void InitG12(int i);
	void InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
		, int iband = 1);
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	// 12 115 maximum see band 28
	int index1[30][3], index2[135][3], index3[2000][2],
		n5, n6, nind[3];// bitfiedl,current index 5 current index 6
	//XY_EXPAND xye6[MAXN6], xye5[MAXN5];
	// row solution pattern in digit
	int mini_digs[9], mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27], nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb);
	void DoExpandBand(int dband);// dband 0/27
	void DebugIndex2();
	void Debug_2_3();
	void PrintShortStatus();
} myband1, myband2;
struct STD_B3 :STD_B416 {// data specific to bands 3
	struct GUAs {
		BF128 isguasocket2, isguasocket3, isguasocket4;// active i81
		int triplet[9];//same gua3s
		int triplet_imini[81];
		int ua_pair[81], ua_triplet[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81];
		int ua2_i27[81];
	}guas;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	//BF128 tbands_UA4_6s, tbands_pairs, tbands_triplets;
	//int tuas46[81];
	//_______________________
	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);
	int IsGua(int i81);
	int IsGua3(int i81);
	int GetI81_2(int bf){
		for (int i = 0; i < 81; i++) if (guas.ua_pair[i] == bf) return i;
		return -1;
	}
	int GetI81_3(int bf) {
		for (int i = 0; i < 81; i++) if (guas.ua_triplet[i] == bf) return i;
		return -1;
	}
	void PrintB3Status();
};
//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
		revised_gangbf[9],// same revised UA2s UA3s ***
		mini_digs[9], mini_pairs[27], // UA2s UA3  ***
		//valid_pairs, //  27 bits valid sockets UA2s ***
		nfloors, limstep, map[9], cptdebug, modemore;
	BF128 valid_sockets;

	//=============== uas collector 
	int limsize, floors;
	uint64_t  tuaold[1000],// previous non hit uas infinal table of uas for bands 1+2
		tua[TUA64_12SIZE],// 
		tuab1b2[200];// collecting bands uas in 2x mode
	uint32_t nuaold, nua, nuab1b2,
		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	int ib, digp;
	uint64_t w0, ua;
	//_____________________ functions collect UAs bands 1+2
	int Initgen();
	void BuildFloorsAndCollectOlds(int fl);
	//int AddUA64(uint64_t * t, uint32_t & nt);
	inline void AddUA(uint64_t v) {
		ua = v; AddUA64(tua, nua, ua);
	}
	inline void AddUACheck(uint64_t v) {
		if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
		ua = v; AddUA64(tua, nua, ua);
	}
	void BuilOldUAs(uint32_t r0);
	int CheckOld();
	int CheckMain(uint64_t wua);
	void CollectMore();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	//_____________________ functions collect UA2s UAs3 socket 

	void ProcessSocket2(int i81);
	int DebugUas();
}genuasb12;

#define SIZETGUA 150
#define GUAREDSIZE 100
struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[290];
	int modeb12, go_back, diagmore,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	int   *gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// actice cols for a given digit
	//____________structs hosting the 81 GUA entries
	struct SGUA2 {// 81 possible UA2 sockets
		// permanent data
		uint64_t tua[SIZETGUA];
		int col1, col2;// columns of the socket
		int i_81; // index 0_80 for this 
		int i9;// complementary column in minirow
		int id1, id2; // index of digits in gang 27 
		// Current band1+2 data
		int digs, dig1, dig2;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		int gangcols[9];// revised gangster
		uint32_t nua;// nua_start, nua_end;
		void Debug(const char * lib);

	}tsgua2[81];
	struct SGUA3 {// 81 possible UA3 sockets
		// permanent data
		uint64_t tua[SIZETGUA];
		int col1;// first columns 0-9 
		int i_81, imini; // index 0_80 for this 
		int id1, id2, id3; // index of digits in gang 27 
		// Current band1+2 data
		int  dig1, dig2, dig3;// depending on gang27 status
		int valid, // valid if guas 
			validuas,// gua2s found
			used;// if needed in bands3
		uint32_t nua;// nua_start, nua_end, nua;
		void Debug(const char * lib);
	}tsgua3[81];
	// __________________________  primary UAs tables and creation of such tables
	uint64_t  // tua3x[3000],// dynamic sub tables
		*ptua2;// pointer to current table cycle search 2/3
	uint32_t tuasb2[500], nuasb2;
	uint32_t  ntua2, ntua3, nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int tactive2[81], nactive2, tactive3[81], nactive3;
	int   tcolok[2], ncolok;
	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//__________________________ secondary guas table 
	uint32_t tuar2[81][GUAREDSIZE], tuar3[81][GUAREDSIZE],
		ntuar2[81], ntuar3[81];
	uint32_t guar2i81[80],guar3i81[80], nguared_2, nguared_3;
	BF128 forced81_2, forced81_3,final81_2,final81_3;
	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	//=== expand band 2 and process 
	int cptclues[16];
	uint64_t ispot ;
	//_____________ one b12 set to check
	uint32_t tcluesb12[12], ncluesb3;
	int ncluesb12,ntb3,nmiss;
	uint32_t mini_bf1 , mini_bf2 , mini_bf3 ,pairsbf , pairs27 , mini_triplet ;
	uint32_t all_used_minis, mincount;
	uint32_t uasb3_1[2000], uasb3_2[500], nuasb3_1, nuasb3_2;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		gang27 = gang[0];
		InitialSockets2Setup();
		InitialSockets3Setup();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void Build_CheckUAs_Subsets_Table();
	void Build_tuacheck(int fl);
	int Have_tuacheck_subset();
	void SecondSockets2Setup();// band1+2 level
	void SecondSockets2MoreUAs();// band1+2 level
	void GuaCollectMore();
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl, int diag = 0);
	//================================= 
	void GetStartB2(int i); // one of the 20 starts 
	void NewBand1();
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	void Go_Sol_Band1();
	void ExpandBand2();
	void ExpandBand3();
	void FinalBand2(uint32_t filter);// final for band2
	void ValidInitGang();
	void Find_band3B();
	//====================== debugging
	int DebugFreshUA(uint64_t ua);
	void GUAs2_Status(int all=0) {// print reduction status
		char ws[82];
		cout << forced81_2.String3X(ws) << " forced 2; n81=" << nguared_2 << endl;
		/*
		cout << "reduction GUAs status nguas reduit=" << ntguab22 << endl;
		if(all)
		for (uint32_t i = 0; i < nguared_2; i++) {
			cout << w.i81 << "\tnua=\t" << w.nua << endl;
			//for (uint32_t j = 0; j < w.nua; j++)
			//	cout << Char27out(w.tua[j]) << endl;
		}
		if(all)
		for (uint32_t i = 0; i < nguared_3; i++) {
			GUA_RED & w = guared_3[i];
			cout << w.i81 << "\tnua=\t" << w.nua << endl;
			//for (uint32_t j = 0; j < w.nua; j++)
			//	cout << Char27out(w.tua[j]) << endl;
		}

*/
	}
	void Go511();

}genb12;
void GEN_BANDES_12::InitialSockets2Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 & w = tsgua2[i];
		w.i_81 = i;
		register int i9 = i / 9;
		int tpcol[3][2] = { {1,2},{0,2},{0,1} };
		int rdcol = i9 % 3, dcol = 3 * (i9 / 3), *p = tpcol[rdcol];
		w.i9 = i9;
		w.col1 = dcol + p[0];
		w.col2 = dcol + p[1];
		int rd1 = i - 9 * w.i9;//relative d1;d2		
		w.id1 = rd1 / 3 + 3 * w.col1;
		w.id2 = rd1 % 3 + 3 * w.col2;
		if (0) {
			cout << "sua2=" << w.i_81 << " i9=" << w.i9 + 1
				<< " cols12 " << w.col1 + 1 << w.col2 + 1
				<< " id1;id2 " << w.id1 << ";" << w.id2 << endl;
		}
	}
}
void GEN_BANDES_12::InitialSockets3Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 3
		SGUA3 & w = tsgua3[i];
		int minir(i / 27);
		w.i_81 = i;
		w.imini = minir;
		w.col1 = 3 * minir;// minirow first column in gangster
		int dp = i - 27 * minir,// the 9 perms of the gangster minirow
			dg = 9 * minir;// gangster 27 start
		 // pointers to gangster digits
		w.id1 = dp / 9;
		dp -= 9 * w.id1;
		w.id1 += dg;
		dg += 3;// next gangster column inthe minirow
		w.id2 = dp / 3;
		dp -= 3 * w.id2;
		w.id2 += dg;
		dg += 3;// last gangster column in the minirow
		w.id3 = dp + dg;
		if (0) {
			cout << "sua2=" << w.i_81 << " col1=" << w.col1 + 1
				<< " id1;id2,id3 " << w.id1
				<< ";" << w.id2 << ";" << w.id3 << endl;
		}
	}
}
void GEN_BANDES_12::SecondSockets2Setup() {
	ntua2 = 0; nactive2 = 0;
	Build_CheckUAs_Subsets_Table();
	for (i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA2 & w = tsgua2[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.digs = (1 << w.dig1) | (1 << w.dig2);
		w.valid = 0;
		//skip if no band3 uses it gua2 gua4 gua6 / build guas 
		for (int iband3 = 0; iband3 < nband3; iband3++) {
			if(bands3[iband3].IsGua(i81))// setup guas in band
				w.valid=1;
		}
		if (!w.valid)continue;
		// build revised gangster
		memcpy(w.gangcols, gangb12, sizeof gangb12);
		w.gangcols[w.col1] ^= w.digs;
		w.gangcols[w.col2] ^= w.digs;
		// find guas of the socket
		zh2b_g.InitGangster(gangb12, w.gangcols);
		zh2b5_g.sizef5 = GUALIMSIZE;
		zh2b5_g.modevalid = 0;
		w.nua = 0;
		//w.nua_start = ntua2;
		ptua2 = w.tua;
		nua2 = 0;

		//================== GUA collector 2 bands 
		GuaCollect(w.digs );
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if ((fl & w.digs) == w.digs)	GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			int fl = floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		for (int i = 0; i < 126; i++) {// find UAs 5digits
			int fl = 0x1ff ^ floors_4d[i];
			if ((fl & w.digs) == w.digs) GuaCollect(fl);
		}
		SecondSockets2MoreUAs();
		if (nua2) {
			tactive2[nactive2++]=i81;
			ntua2 += nua2;
		}
		w.nua = nua2;
	}
}
void GEN_BANDES_12::SecondSockets3Setup() {
	ntua3 = 0; nactive3 = 0;
	for (int i81 = 0; i81 < 81; i81++) {// initial socket 2
		SGUA3 &w = tsgua3[i81];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.dig3 = gang27[w.id3];
		w.valid = 0;

		//skip if no band3 uses it   
		for (int iband3 = 0; iband3 < nband3; iband3++) {
			if (bands3[iband3].IsGua3(i81))// setup guas in band
				w.valid = 1;
		}
		if (!w.valid)continue;
		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = bita | bitb | bitc;
		int triplet_perms[2][3];

		int * p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		ptua2 = w.tua;
		nua2 = 0;

		// build revised gangster
		//w.gangcols[w.col2] ^= w.digs;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, gangb12, sizeof gangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];

			// find guas of the socket
			zh2b_g.InitGangster(gangb12, rgangcols);
			zh2b5_g.sizef5 = UALIMSIZE - 1;//; - 3;
			zh2b5_g.modevalid = 0;
			//w.nua_start = ntua3;

			//================== GUA collector 2 bands 
			GuaCollect(digs);// find UAs 3 digits
			for (int i = 0; i < 126; i++) {// find UAs 4 digits
				int fl = floors_4d[i];
				if ((fl& digs) == digs)	GuaCollect(fl);
			}
			for (int i = 0; i < 126; i++) {// find UAs 5digits
				int fl = 0x1ff ^ floors_4d[i];
				if ((fl & digs) == digs) GuaCollect(fl);
			}
		}
		if (nua2) {
			tactive3[nactive3++] = i81;
			ntua3 += nua2;
		}

		//w.nua_end = ntua3;// store the final count
		w.nua = nua2;// store the final count
	}


}
void GEN_BANDES_12::Build_CheckUAs_Subsets_Table() {
	nua_for_check = 0;
	for (uint32_t i = 0; i < genuasb12.nuab1b2; i++) {// uAs bands
		register uint64_t R = genuasb12.tuab1b2[i]& BIT_SET_2X;
		uint64_t cc = _popcnt64(R);
		if (cc > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R| (cc<<59);
	}
	for (uint32_t i = 0; i < genuasb12.nua; i++){
		register uint64_t R = genuasb12.tua[i];
		if ((R >> 59) > 10) continue;// nearly no chance to be subset
		tua_for_check[nua_for_check++] = R;
	}
	// and build the digit pattern of the uas in tua_for_check
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = 0,res;
		register uint64_t R = tua_for_check[i]&BIT_SET_2X;
		while (bitscanforward64(res, R)) {
			R ^= (uint64_t)1 << res;
			int cell = From_128_To_81[res],digit=grid0[cell];
			digs |= 1 << digit;
		}
		uadigs_for_check[i] = digs;
		if (0) {
			cout << Char2Xout(tua_for_check[i]) << " ua for check digs 0"
				<< oct << digs << dec <<" ndigs "<<_popcnt32(digs)<< endl;
		}
	}

}
void GEN_BANDES_12::Build_tuacheck(int fl) {
	nua_check=0;
	int ndigs = 0x1ff ^ fl;
	for (uint32_t i = 0; i < nua_for_check; i++) {
		uint32_t digs = uadigs_for_check[i];
		if (digs & ndigs) continue;
		tuacheck[nua_check++]= tua_for_check[i];
	}
}
int GEN_BANDES_12::Have_tuacheck_subset(  ) {// ua 2x27  + 5 bit length
	uint64_t * t = tuacheck;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nua_check; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}
void GEN_BANDES_12::SecondSockets2MoreUAs() {
	diagmore = 0;
	SGUA2 s = tsgua2[i81];
	int rx0[2], bx0[2], dx3[2];
	rx0[0] = zh2b_g.fd_sols[0][s.dig1].SolRow(s.col2);
	dx3[0] = genb12.grid0[9 * rx0[0] + s.col1];
	rx0[1] = zh2b_g.fd_sols[0][s.dig2].SolRow(s.col1);
	dx3[1] = genb12.grid0[9 * rx0[1] + s.col2];
	bx0[0] = rx0[0] / 3;
	bx0[1] = rx0[1] / 3;

	if (bx0[0] == bx0[1]) return; // want digitsin 2 bands
	if (dx3[0] == dx3[1]) return; // ua 6 cells
	for (int idig = 0; idig < 2; idig++) {
		int dx[2], cx[2],row,dig3;
		if (idig) {
			dx[0] = s.dig1; dx[1] = s.dig2;
			cx[0] = s.col1; cx[1] = s.col2;
			row = rx0[0]; dig3 = dx3[0];
		}
		else {
			dx[0] = s.dig2; dx[1] = s.dig1;
			cx[0] = s.col2; cx[1] = s.col1;
			row = rx0[1]; dig3 = dx3[1];
		}
		band=row/3;
		// dig3 must be possible in 
		if (!((1 << dig3) & s.gangcols[cx[1]])) {
			//cout << " not a valid pattern dig3 col" << endl;
			continue;
		}
		// now fill band bx[0] except the 2 cells
		zh2b[0].Init_gang();
		BF64 cells_to_fill;
		cells_to_fill.bf.u64=BIT_SET_2X;// all cells
		cells_to_fill.bf.u32[1-band] = 0;// nothing in the other band
		c1 = row * 9 + s.col1;
		c2 = row * 9 + s.col2;
		if (row >= 3) { c1 += 5; c2 += 5; }// must be 2X mode
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c1;
		cells_to_fill.bf.u64 ^= (uint64_t)1 << c2;
		zh2b[0].Init_2digits_banda(cells_to_fill);
		zh2b[0].FullUpdate();
		if(diagmore)zh2b[0].ImageCandidats();
		// at this point,we have only one band with a gangster 
		//not equal to the start and a partial GUA 2cells in the other band
		//we  prepare ZHONE to find partial uas in one band
		uint32_t tua[300]; // supply a tua to ZHONE
		zh1b_g.tua = tua;
		for (int i = 0; i < 9; i++) {// load the gangster pm
			zh1b_g.fd_sols[1][i] = zh2b[0].FD[i].bf.u32[1 - band];
		}
		STD_B416 * myb = &myband1;
		if(!band) myb= &myband2; // to catch the solution
		memcpy(zh1b_g.fd_sols[0], myb->fd_sols[0], sizeof zh1b_g.fd_sols[0]);
		zh1b_g.band0 = myb->band0;
		digp = s.digs | (1 << dig3);// must be in the floor
		//==============  now find more uas 6/7 digits
		ngua6_7 = 6;
		wua0 = ((uint64_t)1 << c1) | ((uint64_t)1 << c2);
		for (int i6 = 0; i6 < 84; i6++) {
			int fl3 = floors_3d[i6];// , fl6 = 0x1ff ^ fl3;
			if (fl3&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl3;
			if (diagmore)cout << "try 6 floors 0" << oct << floors << dec << endl;
			GuaCollectMore();
		}
		ngua6_7 = 7;
		for (int i7 = 0; i7 < 36; i7++) {
			int fl2 = floors_2d[i7];// , fl7 = 0x1ff ^ fl2;
			if (fl2&digp) continue;// digits must be in the multi floors
			floors = 0x1ff ^ fl2;
			if (diagmore)cout << "try 7 floors 0" << oct << floors << dec << endl;
			GuaCollectMore();
		}
	}
}
void GEN_BANDES_12::GuaCollectMore() {
	zh1b_g.nua = 0;
	zhone[0].Start_nFloors(floors);
	zhone[0].InitGuess();
	if (diagmore > 1) zhone[0].ImageCandidats();
	else zh1b_g.diag = 0;
	if (ngua6_7 == 6)zhone[0].Guess6();
	else zhone[0].Guess7();
	if (zh1b_g.nua) {
		Build_tuacheck(floors);
		SGUA2 s = tsgua2[i81];
		for (uint32_t i = 0; i < zh1b_g.nua; i++) {
			uint32_t w = zh1b_g.tua[i]&BIT_SET_27,	cc=_popcnt32(w),cell;
			//cout << Char27out(w) << endl;;
			if (cc < 8) continue;
			if (cc < 12 && ngua6_7 == 7) continue;
			if (cc > GUALIMSIZE - 2) continue;
			int digpw = digp;// assigned outside the band
			register uint32_t R = w;
			while (bitscanforward(cell, R)) {
				R ^= 1 << cell;
				int dig = zh1b_g.band0[cell]; // digit changed
				digpw |= 1 << dig;
			}
			if (_popcnt32(digpw) != ngua6_7) continue;
			//cout << "valid ua to add to the table" << endl;
			ua = wua0;
			ua |= (uint64_t)w << 32 * (1 - band);
			if (Have_tuacheck_subset()) {
				//cout << Char2Xout(ua) << " ua has subset" << endl;
				continue;
			}
			ua |= (uint64_t)(cc + 2) << 59;
			if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
			int ir = 		AddUA64(ptua2, nua2, ua);
			if(ir && diagmore)		cout << Char2Xout(ua) << " wua added cc=" << cc + 2 <<" nua2="<<nua2<< endl;

		}
	}
	// try to add the ua to the file
}
void GEN_BANDES_12::GuaCollect(int fl,int diag) {//use revised gangster
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 1);
	if (!solved_cells) return;// one digit solved true
	if (diag) cout << "return from zh2b5_g.FindUAsInit(fl, 1);" << endl;
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	if (diag) cout << "zh2b5_g.CollectUas5();" << endl;
	if (!zh2b5_g.nuaf5) return;
	Build_tuacheck( fl);

	// check subsets and add to main table
	for (uint32_t i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64;
		//genuasb12.ua = ua;
		if (diag)cout << Char2Xout(ua) << " to add if new2 nua=" << nua2 << endl;
		if (Have_tuacheck_subset()) continue;// superset of a previous ua
		//		if (genuasb12.CheckOld()) continue;// superset of a previous ua
		uint64_t cc = _popcnt64(ua&BIT_SET_2X);
		ua |= cc << 59;
		if (nua2 >= SIZETGUA)nua2 = SIZETGUA - 1; // guess it will be a smaller
		int ir = AddUA64(ptua2, nua2,ua);
		if (ir) {
			if (diag)cout << Char2Xout(genuasb12.ua) <<" ir="<<ir << " added nua2="<<nua2 << endl;
		}
	}
}





void GEN_BANDES_12::GetStartB2(int ip) {//set  rows 3_9 column 1
	char const *tp[20] = {// 3 out of 6 ordered 
		"012345", "345012", "013245", "245013", "014235", "235014", "015234", "234015",
		"023145", "145023", "024135", "135024", "025134", "134025",
		"034125", "125034", "035124", "124035", "045123", "123045"
	};
	char tpw[7];	strcpy( tpw , tp[ip]);
	for (int i = 0; i < 6; i++)boxd[i] = 0x1ff;
	for (int j = 0, jc = 27; j < 6; j++, jc += 9) {
		int ic = tpw[j] - '0', c0 = tc[ic], bit = 1 << c0;
		grid0[jc] = c0;
		zsol[jc] = c0+'1';
		rowd[j] = 0x1ff ^ bit;
		if (j < 3)boxd[0] ^= bit; else  boxd[3] ^= bit;
	}
}
uint32_t t2clues_band1[9] = {
0400000001,0200000004,02000020,
04000040,010000100,01000010,
020000200,040000400,0100000002,
};
void GEN_BANDES_12::NewBand1() {
	go_back = 0;
	modeb12 = 0;
	myband1.Initstd();
	zsol[81] = 0;
	nb12 = 0;
	i1t16 = 415;	it16 = 28;
	myband1.InitG12(28);
	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	memcpy(gcheck, myband1.band0, sizeof myband1.band0);
	strcpy(zsol, myband1.band);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	for (int i = 0; i < 9; i++) // init columns status
		cold[i] = 0x1ff ^ myband1.gangster[i];
	zsol[27] = 0;
	cout << "i1t16=" << i1t16 << " it16=" << it16 
		<< " n auto morphs=" << n_auto_b1 << endl;
/*
	uint32_t *tua=myband1.tua, nua=myband1.nua;//   maximum 81
	//myband1.PrintStatus();

	int tua1[20], ntua1 = 0;
	BitsInTable32(tua1, ntua1, tua[0]);// first ua
	cout << " ua1 size " << ntua1 << endl;
	for (int ic = 0; ic < ntua1; ic++) {
		uint32_t cell = tua1[ic], bit = 1 << cell,and=BIT_SET_27;
		for (uint32_t iu = 1; iu < nua; iu++) {
			register uint32_t R = tua[iu];
			if (R&bit) continue;
			and &= R;
		}
		cout << Char27out(bit) << endl;
		cout << Char27out(and)<< "and" << endl << endl;
		if (and) cout << oct << (bit | and) << dec << endl;
	}
*/
	ntc = 0;
	BitsInTable32(tc, ntc, cold[0]);// first col 6 digits in table
	for (int ip = 0; ip < 20; ip++) {//0;ip<20 setup initial values for rows columns
		GetStartB2(ip);
		Find_band2B();
		if (go_back)return;
	}
}


void Go_c500() {// find 18 clues puzzles one band 2 clues
	cout << "entry 500 18 with band 2 clues" << endl;
	zh_g2.grid0=genb12.grid0;
	genb12.NewBand1();
	cout << "print final stats" << endl;
	for (int i = 0; i < 40; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
	cout << "compte clues bande 2" << endl;
	for (int i = 0; i < 15; i++) if (p_cpt1g[i]) {
		cout << i << "\t" << p_cpt1g[i] << endl;
	}
	cout << "exit" << endl;
	cerr << "exit" << endl;
}

void GEN_BANDES_12::Find_band2B() {
	int * zs0 = &grid0[27];
	register int  *crcb, bit;
	int *rd = rowd, *cd = cold, *bd = boxd; // to updates rows cols boxes
	char * zs = zsol;
	// now loop over the 24 cells not yet assigned in the band to fill the band (use relative cell) 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{	crcb = tgen_band_cat[ii];//cell_r_c_b  24 cells to fill
	register int fr0 = cd[crcb[2]] & bd[crcb[3]], fr = rd[crcb[1]] & fr0;
	if (crcb[4])if (_popcnt32(fr0) < 3) goto back; // 3 clues needed here
	if (!fr)goto back;
	free[ii] = fr;
	}
	goto next_first;
next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 27] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if lexically minimale 
		int ir = bandminlex.Getmin(zs0, &pband2, 0);
		if (ir < 0) {//would be bug  did not come in enumeration
			cerr << "Find B2 invalid return  Getmin" << endl;
			return;
		}
		pcheck2 = pband2;
		it16_2 = pband2.i416;
		ib2check = i2t16 = t416_to_n6[it16_2];
		{
			memcpy(&gcheck[27], zs0, 27 * sizeof gcheck[0]);
			if (Band2Check())goto next;// do nothing if p2b
		}

		nb12++;
		p_cpt2g[0] ++;// compte band2
		if (p_cpt2g[0] < sgo.vx[0]) goto next;
		myband2.InitBand2_3(it16_2, &zsol[27], pband2);
		//myband2.PrintStatus();

		//__________________________ bands 3
		Find_band3B();
		if (!nband3) goto next;
		for (int i = 0; i < nband3; i++) {
			STD_B3 wb3= bands3[i];
			int min= t16_min_clues[wb3.i416];
			//cout << wb3.band << "\t" << wb3.i416 << "\t" << min << endl;
		}
		p_cpt2g[1] += nband3;// compte ba32
		if (nband3 > p_cpt2g[4]) p_cpt2g[4] = nband3;
		//if (1) goto next;
		// setup gangster
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
			gangb12[i] = 0x1ff ^ cold[i];
		}
		memcpy(gangcols, cold, sizeof gangcols);

		//=========================== collect UAs  GUAs 
		zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more

		if (genuasb12.Initgen()) return;
		if (1) {
			BuildGang9x3();
			zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
			genb12.SecondSockets2Setup();// collect GUA2s
			genb12.SecondSockets3Setup();// collect GUA3s
			cout << "==========" << p_cpt2g[0] << " b2id="<<myband2.i416 << endl;
			if (p_cpt2g[0] == sgo.vx[3]) {
				cout << "n bands3      \t" << genb12.nband3 << endl;
				cout << "ua bands1+2   \t" << genuasb12.nua << endl;
				cout << "guas socket2  \t" << genb12.ntua2 << endl;
				cout << "guas socket3  \t" << genb12.ntua3 << endl;
				cout << "active socket2\t" << genb12.nactive2 << endl;
				cout << "active socket3\t" << genb12.nactive3 << endl;

			}
		}
		//if (1) return;
		//________________loop on the 9 band1 
		Go_Sol_Band1();
	
		if (p_cpt2g[0]>sgo.vx[2]) return; 
		goto next;
	}
	/*
//_______________________ std process
	if ((nb12 >> 6) < skip) return 0;// here restart value, kept untouched if no band 3 found
	Find_band3B();
	{// print a restart point every 64 bands 1+2 seen
		uint64_t w = genb12.nb12, w1 = w >> 6;
		w &= 63;
		if (w == 0) {
			long tfin = GetTimeMillis();
			cout << "next skip value to use=\t" << w1 <<"\t"<<(tfin-sgo.tdeb)/1000<<"\t"<< p_cpt2g[0]<< endl;
		}
	}
	if ((nb12 >> 6) >= last)return 1;
	return 0;
*/
back:
	if (--ii >= 0) goto next;
}
void GEN_BANDES_12::Find_band3B() {
	register int  *crcb, bit;
	nband3 = 0;
	int *rd = rowdb3, *cd = cold, *bd = boxdb3; // to updates rows cols boxes
	char * zs = zsol;
	int * zs0 = &grid0[54];
	memcpy(boxdb3, &boxd[3], sizeof boxdb3);
	memcpy(rowdb3, &rowd[3], sizeof rowdb3);
	// now loop over the 24 cells not yet assigned in the band to fill the band use relative cell 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{
		crcb = tgen_band_cat[ii];//cell_row_col_box one of the 24 cells to fill
		register int fr = cd[crcb[2]] & bd[crcb[3]] & rd[crcb[1]];
		if (!fr)goto back;
		free[ii] = fr;
	}
	goto next_first;

next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 54] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if canonical 
		int ir = bandminlex.Getmin(zs0, &pband3, 0);
		if (ir < 0) {//would be bug  did not come in enumeration
			cerr << "gen band 3 invalid return Getmin" << endl;
			return;
		}
		int it16_3 = pband3.i416;
		ib3check = i3t16 = t416_to_n6[it16_3];
		if (Band3Check())goto next;
		//if(ib3check==ib2check)//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		bands3[nband3++].InitBand3(it16_3, &zs[54], pband3);
		goto next;
	}
back:
	if (--ii >= 0) goto next;

	//cout << "nb3=" << nband3 << endl;
}

void GEN_BANDES_12::Go_Sol_Band1() {
	uint32_t nua = genuasb12.nua;
	uint64_t  *tua = genuasb12.tua;
	for (int i = 0; i < 9; i++) {
		register uint32_t b1 = t2clues_band1[i];
		ncluesb12 = 0;
		BitsInTable32((int *)tcluesb12,  ncluesb12, b1);
		nuasb2 = 0;
		for (uint32_t iua = 0; iua < nua; iua++) {
			register uint64_t ua = tua[iua];
			if (ua&b1) continue; // hit by the 2 clues in b1
			ua >>= 32; // keep only b2
			register uint32_t b2 = (uint32_t)ua & BIT_SET_27;
			b2 |= (_popcnt32(b2) << 27);
			if (nuasb2 == 500)nuasb2--;
			AddUA32(tuasb2, nuasb2, b2);
		}
		// add uas specific to band2
		{
			uint32_t  *tua = myband2.tua, nua = myband2.nua;
			for (uint32_t iua = 0; iua < nua; iua++) {
				register uint32_t b2 = tua[iua] & BIT_SET_27;
				b2 |= (_popcnt32(b2) << 27);
				AddUA32(tuasb2, nuasb2, b2);
			}

		}
		for (uint32_t j = 0; j < nuasb2; j++)	tuasb2[j] &= BIT_SET_27;// kill count
		// find a minimum number of clues through disjoint uas	
		// assuming most are minirow uas
		uint32_t disj_uas = tuasb2[0], ndisj = 1;
		for (uint32_t iua = 1; iua < nuasb2; iua++) {
			if (tuasb2[iua] & disj_uas) continue;
			disj_uas |= tuasb2[iua];
			ndisj++;
		}
		if (p_cpt2g[0]==sgo.vx[3]) {
			cout << Char27out(b1) << " b1 2 clues" << endl;
			cout << "minband2 = " << ndisj << endl;
			cout << "reduced table nuas= " << nuasb2 << endl;
			for (uint32_t j = 0; j < nuasb2; j++)
				cout << Char27out(tuasb2[j]) << " s=" << (tuasb2[j] >> 27) << " i=" << j << endl;
		}
		if (ndisj > 8) continue; // nothing to do for a 18 clues
		p_cpt2g[20 + ndisj]++;
		if (0) {
			ExpandBand2();
			continue;
		}
		//_______________ guas filter => already killed/forced plus sub table
		// build new subtables still active not fixed 
		nguared_2 = 0;
		memset(ntuar2, 0, sizeof ntuar2);
		forced81_2.SetAll_0(); forced81_3.SetAll_0();
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			uint32_t * tuasw = tuar2[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t Ru = w.tua[iua];
				if (Ru&b1) continue; // hit by the 2 clues in b1
				Ru >>= 32; // keep only b2
				register uint32_t b2 = (uint32_t)Ru & BIT_SET_27;
				if (!b2) {// empty in band 2 stop and force
					forced81_2.Set_c(i81);
					nt = 0;
					goto nexti81_2;
				}
				b2 |= (_popcnt32(b2) << 27);
				AddUA32(tuasw, nt, b2);
			}
			if (nt) {
				ntuar2[i81] = nt;
				guar2i81[nguared_2++] = i81;
			}
		nexti81_2:;
		}
		nguared_3 = 0;
		memset(ntuar3, 0, sizeof ntuar3);
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			uint32_t * tuasw = tuar3[i81], nt = 0;
			for (uint32_t iua = 0; iua < w.nua; iua++) {
				register uint64_t Ru = w.tua[iua];
				if (Ru&b1) continue; // hit by the 2 clues in b1
				Ru >>= 32; // keep only b2
				register uint32_t b2 = (uint32_t)Ru & BIT_SET_27;
				if (!b2) {// empty in band 2 stop and force
					forced81_3.Set_c(i81);
					nt = 0;
					goto nexti81_3;
				}
				b2 |= (_popcnt32(b2) << 27);
				AddUA32(tuasw, nt, b2);
			}
			if (nt) {
				ntuar3[i81] = nt;
				guar3i81[nguared_3++] = i81;
			}
		nexti81_3:;
		}
		//GUAs2_Status(0);
		ExpandBand2();
	}
}
/*

}
void G17INDEXSTEP::ShrinkGuas() {// shrink table for gangster uas2 ua3s
// group all gua2 gua2 soskets in one table
	ntgua = ntgua_raw = nactive2 = nactive3 = 0;
	for (int i = 0; i < genb12.nactive2; i++) {
		int i81 = genb12.tactive2[i];
		GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
		int n = ShrinkGuasLoad(w.tua,  w.nua);
		if (n) {// still valid guas for this i81 gua2
			tactive2_end[nactive2] = ntgua;
			tactive2_start[nactive2] = ntgua - n;
			tactive2[nactive2++] = i81;
		}
	}
	for (int i = 0; i < genb12.nactive3; i++) {
		int i81 = genb12.tactive3[i];
		GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
		int n = ShrinkGuasLoad(w.tua,  w.nua);
		if (n) {// still valid guas for this i81 gua3
			tactive3_end[nactive3] = ntgua;
			tactive3_start[nactive3] = ntgua - n;
			tactive3[nactive3++] = i81;

		}
	}
	if (g17b.debug17>1)cout << "ntgua = " << ntgua << " limit " <<128 * G17BLOCGSUA << endl;
	if (ntgua > 128 * G17BLOCGSUA) ntgua = 128 * G17BLOCGSUA;
	//working with a limited number of GUAs for an index chunk
	n128vgua = (ntgua + 127) >> 7;
	{ // initial vectors to appropriate value and actives vector
		memset(v256guas, 255,sizeof v256guas);
		memset(&v256ga, 0, sizeof v256ga);
		register BF128 * Ra = v256ga.v;
		register int Rn = ntgua;
		while (Rn > 0) {
			if (Rn >= 128){
				Ra->SetAll_1();
				Rn -= 128;
				Ra++;
			}
			else {
				*Ra = maskLSB[Rn];
				break;
			}
		}
	}
	uint32_t cc;
	for (int i = 0; i < ntgua; i++) {// set uas
		register int bloc = i >> 7, ir = i & 127;
		register uint64_t Rw = tgua[i], biti = (uint64_t)1 << ir;
		Rw &= BIT_SET_2X;
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			Rw ^= bit2;// clear bit
			v256guas[From_128_To_81[cc] ].v[bloc].clearBit(ir);
		}
	}

	// ____________________________________create id81 vectors of corresponding guas
//	uint64_t * vx = vid81;
	for (int i = 0; i < nactive2; i++)
		SetV(vid81s2[tactive2[i]], tactive2_start[i], tactive2_end[i]);
	for (int i = 0; i < nactive3; i++)
		SetV(vid81s3[tactive3[i]], tactive3_start[i], tactive3_end[i]);

}
void G17INDEXSTEP::SetV(V256_GUAS & vid, int i1, int i2) {
	memset(&vid, 0, sizeof vid);
	BF128 * v = vid.v;
	while (1) {// *v initial value set to 0 outside the call
		if (i1 >= 128) {
			v++;
			i1 -= 128;
			i2 -= 128;
			continue;
		}
		if (i1) {// last bloc i1
			if (i2 > 128) {
				v->SetAll_1(); // active bloc initial value to all '1'
				*v-= maskLSB[i1];
				i1 = 0;
				i2 -= 128;
				v++;
			}
			else {// last bloc
				*v= maskLSB[i2];
				*v-= maskLSB[i1];
				return;
			}
		}
		// it is always the last bloc i1=0 i2 residual count
		*v -= maskLSB[i2];
		return;
	}
}
*/
void GEN_BANDES_12::ExpandBand2() {
	memset(cptclues, 0, sizeof cptclues);
	// expand first all uas size <=3 (disjoints except 2/3 pairs same mini row )
	// usually 6/7 clues
	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[15], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tuasb2[0];
	//cout << Char27out(s3->possible_cells)<<" depart" << endl;
	int tcells[15];
	//____________________  here start the search
next:
	ispot = s3 - spb3;
	//cout <<ispot<<" ispot"<<endl;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		if (ispot ==6)p_cpt2g[7] ++;

		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nuasb2; i++) {
			if (tuasb2[i] & filter)continue;
			if (ispot >=7)goto next;
			sn3->iuab3 = i;
			sn3->possible_cells = tuasb2[i] & ac;
			s3 = sn3; // switch to next spot
			if (0)			
				cout << Char27out(s3->all_previous_cells) << " next ispot="<<ispot << endl;
			goto next;
		}

	}
	// no more ua
	p_cpt2g[2] ++;
	cptclues[ispot+1]++;
	if (ispot >7)goto next;
	// check if this is a valid band 1+2
	ncluesb12 = 2;// keep band 1
	BitsInTable32((int *)tcluesb12, ncluesb12, sn3->all_previous_cells,27);
	register uint64_t myua = zh2b[0].ValidXY(tcluesb12, ncluesb12);
	if (myua) {
		//if (p_cpt2g[3] < 100)
			//cout << Char2Xout(myua) << "new ua" << endl;
		if (genuasb12.nua < TUA64_12SIZE) genuasb12.AddUA(myua);
		if (nuasb2 < 500) {// 500 is the limit for tuasb2
			tuasb2[nuasb2++] = myua >> 32;
		}
		if (ispot < 7) {// must use the UA as new spot
			sn3->iuab3 = nuasb2-1;
			sn3->possible_cells = tuasb2[sn3->iuab3] & sn3->active_cells;
			s3 = sn3; // switch to next spot
		}
		goto next;
	}
	p_cpt2g[3] ++;// valid bands 1+2
	p_cpt2g[ispot - 1] ++;
	FinalBand2(sn3->all_previous_cells);// do final	
	if (ispot < 7) {// if below 8 one more step all active cells
		sn3->possible_cells =  sn3->active_cells;
		s3 = sn3; // switch to next spot
	}
	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
	for (int i = 0; i < 15; i++) if (cptclues[i]) {
		//cout << i << "\t" << cptclues[i] << endl;
		p_cpt1g[i] += cptclues[i];
	}

}
void GEN_BANDES_12::ExpandBand3() {
	struct SPB3 {// spots to find band 3 minimum valid solutions
		uint32_t  possible_cells, all_previous_cells, active_cells, iuab3;
	}spb3[10], *s3, *sn3;
	uint32_t * tua = uasb3_2, nua = nuasb3_2;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tua[0] & BIT_SET_27;
	int tcells[15];
	//____________________  here start the search
next:
	ispot = s3 - spb3;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = 1 << cell;
		tcells[ispot] = cell;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop
		for (uint32_t i = s3->iuab3 + 1; i < nua; i++) {
			if (tua[i] & filter)continue;
			if (ispot >= ncluesb3 - 1)goto next;
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}
	}	// no more ua
	// check if this is a valid band 1+2+3
	int ir = zhou[1].CallMultipleB3(zhou[0], sn3->all_previous_cells, 0);
	if (ir) {//consider store the fresh ua b3
		if (nua < 500) {// 500 is the limit for tuasb2
			tua[nua++] = zh_g2.cells_assigned.bf.u32[2];
		}
	}
	goto next;
back:
	if (--s3 >= spb3)goto next;
}
struct G17B3HANDLER {
	int known_b3, rknown_b3, active_b3,ib3,nb3,
		active_sub, ndead, wactive0, nmiss, ncritical,
		irloop, *uasb3, nuasb3, wua;
	uint32_t mini_bf1, mini_bf2, mini_bf3, pairsbf, pairs27, mini_triplet;
	int diagh;
	// ================== entry in the proces
	void Init(int i);
	int IsMultiple(int bf);
	int ShrinkUas1(int * to, int no);
	void Go();
	//=============== process critical
	void CriticalAssignCell(int Ru);
	void Critical2pairs();
	void Go_Critical();
	void CriticalEntryFromSubcritical(int *to, int no);
	void CriticalLoop(int *to, int no);
	void CriticalExitLoop(int *uasb3, int nuasb3);
	void Critical_0_UA();
	void CriticalFinalCheck();
	//===================== process not critical
	void Go_Not_Critical_missn();

	//==================== process subcritical no cell added outside the GUAs field
	void SubMini( int M, int mask);
	void Go_Subcritical(int docheck = 1);
	void Go_SubcriticalMiniRow();
	void Go_SubcriticalMiniRow_End(int stack);
	//===============  debugging 
	void PrintStatus(int mode=0);
	void PrintStatusAfterFirstLoop(int ir, int * uasb3, int nuasb3);
	void Print_Ok3();
	void Print_Ok3GUAs();
	void Print_Ok3_Not_Critical(int ractive_outfield);
	void PrintB3UasOut(int wac);
	void DNC(int wua, int rawua, char * lib);

};

void GEN_BANDES_12::FinalBand2(uint32_t filter) {// reaching a valid bands 1+2
	//cout << "finalb2" << endl;
	final81_2 = forced81_2;
	for (uint32_t i = 0; i < nguared_2; i++) {
		int i81 = guar2i81[i];
		uint32_t * tua = tuar2[i81];
		for (uint32_t iua = 0; iua < ntuar2[i81]; iua++) {
			register uint32_t Ru = tua[iua];
			if (Ru&filter)continue;
			final81_2.Set_c(i81);
			break;// one ua not hit is enough here
		}
	}
	final81_3 = forced81_3;
	for (uint32_t i = 0; i < nguared_3; i++) {
		int i81 = guar3i81[i];
		uint32_t * tua = tuar3[i81];
		for (uint32_t iua = 0; iua < ntuar3[i81]; iua++) {
			register uint32_t Ru = tua[iua];
			if (Ru&filter)continue;
			final81_3.Set_c(i81);
			break;// one ua not hit is enough here
		}
	}
	ntb3 = 0;
	ncluesb3 = 16 - _popcnt32(filter);
	int zinitdone = 0;
	if (p_cpt2g[0] == sgo.vx[3])
		cout<< Char27out(filter) << " new valid bands 1+2 check bands 3<<<<<<<<" << endl;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3::GUAs & myb = genb12.bands3[ib3].guas;
		BF128 ws2 = final81_2 & myb.isguasocket2;
		BF128 ws3 = final81_3 & myb.isguasocket3;
		// switch to mini rows patterns
		int tix[81], ntix = ws2.Table3X27(tix);
		mini_bf1 = mini_bf2 = mini_bf3 = pairsbf = pairs27 = mini_triplet = 0;
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int i81 = tix[i], imini = myb.ua2_imini[i81], bit = 1 << imini;
			if (mini_bf2&bit) mini_bf3 |= bit;
			if (mini_bf1&bit) mini_bf2 |= bit;
			mini_bf1 |= bit;
			pairsbf |= myb.ua_pair[i81];
			pairs27 |= 1 << myb.ua2_i27[i81];
		}
		ntix = ws3.Table3X27(tix);// now triplets to mini rows
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int imini = myb.ua3_imini[tix[i]], bit = 1 << imini;
			mini_triplet |= bit;
		}
		//___________________________ prepare a new band to process 
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(mini_bf1) + _popcnt32(mini_bf3) + _popcnt32(mini_triplet);
		if (mincount > ncluesb3) continue; // too many clues 
		nmiss = ncluesb3 - mincount;
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs
		// set up pair + triplet bitfield
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					pairsbf |= field;
		}

		//============= collect Gua46 and uas b3 for the band split them "in-field" "out-field"
		nuasb3_1 = nuasb3_2 = 0;
		register int  Rfilt = pairsbf;
		// first GUA46 usually shorter than UAs band3
		BF128  socket4 = genb12.bands3[ib3].guas.isguasocket4;// i81 3X
		socket4 &= final81_2;
		int * ua_46 = genb12.bands3[ib3].guas.ua_pair; // ua pattern
		int i81;
		while ((i81 = socket4.getFirstCell()) >= 0) {
			socket4.Clear_c(i81);// clear bit
			register uint32_t Ru = ua_46[i81] & BIT_SET_27;
			Ru |= _popcnt32(Ru) << 27;
			if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
			else if (!nmiss) goto nextib3;// critical + outfield uas
			else AddUA32(uasb3_2, nuasb3_2, Ru);
		}
		uint32_t * to = genb12.bands3[ib3].tua;
		for (uint32_t i = 0; i < genb12.bands3[ib3].nua; i++) {
			register uint32_t Ru = to[i] & BIT_SET_27;
			Ru |= _popcnt32(Ru) << 27;
			if (Ru & Rfilt)	AddUA32(uasb3_1, nuasb3_1, Ru);
			else if (!nmiss) goto nextib3;// critical + outfield uas
			else AddUA32(uasb3_2, nuasb3_2, Ru);
		}
		//cout << " nmiss=" << nmiss <<" nuas 1="<< nuasb3_1 << " nuas 2=" << nuasb3_2 << endl;
		if (nmiss < 2)p_cpt2g[10 + nmiss] ++;
		else if (nmiss < 5)p_cpt2g[12] ++;
		else p_cpt2g[13] ++;
		p_cpt2g[15] ++;
		memcpy(&genb12.grid0[54], genb12.bands3[ib3].band0, 4 * 27);
		if (!zinitdone) {
			zinitdone = 1;
			if (zhou[0].PartialInitSearch17(tcluesb12, ncluesb12))
				return;// would be  bug
		}
		G17B3HANDLER hh0; hh0.Init(ib3);
		//if (nmiss) goto nextib3;//<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (!mincount) ExpandBand3();
		else if (!nmiss)hh0.Go_Critical();
		else hh0.Go_Not_Critical_missn();
	nextib3:;
	}
}

void G17B3HANDLER::Init(int i) {
	ib3 = i;
	known_b3 = rknown_b3 = 0;
	ndead = BIT_SET_27;
	active_b3 = genb12.pairsbf;// active cells in field
	wactive0 = BIT_SET_27 ^ active_b3;//  active cells out field
	nmiss = genb12.nmiss;
	ncritical = genb12.mincount;
	nb3 = nmiss + ncritical;
	mini_bf1 = genb12.mini_bf1;
	mini_bf2 = genb12.mini_bf2;
	mini_bf3 = genb12.mini_bf3;
	mini_triplet = genb12.mini_triplet;
	diagh = 0;
}

void G17B3HANDLER::CriticalAssignCell(int Ru){// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & mini_bf3){// cell is in minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		mini_bf3 ^= bit; // now only a pair to hit
		mini_bf1 |= bit;
	}
	else{// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_b3 &= (~Mask); // kill the minirow as active
		mini_bf1 &= ~bit;
		mini_triplet &= ~bit;
	}
}
int G17B3HANDLER::IsMultiple(int bf) {
	if (_popcnt32(bf) > 25) return 0;
	// check first if all tuab3 is hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf, diagh);
	if (ir) {//consider store the fresh ua b3
		STD_B3 & myb3 = genb12.bands3[ib3];

		if (p_cpt2g[0] == sgo.vx[3]) {
			for (int i = 0; i < 54; i++)cout << genb12.grid0[i] + 1;
			cout << endl;
			for (int i = 54; i < 81; i++)cout << genb12.grid0[i] + 1;
			cout << endl;
			//zhou[1].ImageCandidats();
			cout << Char27out(zh_g2.cells_assigned.bf.u32[2]) << " ua" << endl;
		}
		BF128 wua = zh_g2.cells_assigned;
		uint32_t ua = wua.bf.u32[2];
		int cc = _popcnt32(ua);
		if (cc < 4) {
			if (cc == 2) {// fresh gua2
				int i81 = myb3.GetI81_2(ua);
				if (i81 >= 0) {
					genb12.final81_2.Set_c(i81);// new valid gua2 for other bands
					if (!genb12.ntuar2[i81]) genb12.guar2i81[genb12.nguared_2++] = i81;;
					if(genb12.ntuar2[i81]< GUAREDSIZE) 
						genb12.tuar2[i81][genb12.ntuar2[i81]++]= wua.bf.u32[1];
					GEN_BANDES_12::SGUA2 & sg = genb12.tsgua2[i81];
					if (sg.nua >= SIZETGUA)sg.nua = SIZETGUA - 1;
					AddUA64(sg.tua,sg.nua,wua.bf.u64[0]);
				}
				if (p_cpt2g[0] == sgo.vx[3])cout << Char27out(ua) << "gua2 added" << endl;
			}
			if (cc == 3) {// fresh gua2
				int i81 = myb3.GetI81_3(ua);
				if (i81 >= 0) {
					genb12.final81_3.Set_c(i81);// new valid gua2 for other bands
					if (!genb12.ntuar3[i81]) genb12.guar3i81[genb12.nguared_3++] = i81;;
					if (genb12.ntuar3[i81] < GUAREDSIZE)
						genb12.tuar3[i81][genb12.ntuar3[i81]++] = wua.bf.u32[1];
					GEN_BANDES_12::SGUA3 & sg = genb12.tsgua3[i81];
					if (sg.nua >= SIZETGUA)sg.nua = SIZETGUA - 1;
					AddUA64(sg.tua, sg.nua, wua.bf.u64[0]);
				}
				if (p_cpt2g[0] == sgo.vx[3])cout << Char27out(ua) << "gua3 added" << endl;
			}
		}
	}
	return ir;
}

void G17B3HANDLER::CriticalFinalCheck(){// no more ua is it a valid solution
	p_cpt2g[8]++;
	if (p_cpt2g[8] < 50 ||p_cpt2g[0] == sgo.vx[3]) {
		for (int i = 54; i < 81; i++)cout << genb12.grid0[i] + 1;
		cout << endl;
		cout << Char27out(known_b3) << " known b3 final" << endl;
	}
	int ncl = _popcnt32(known_b3);
	if (ncl != genb12.ncluesb3) return; // should be seen earlier if possible
	register int ir = IsMultiple(known_b3); 
	if (!ir){// one sol to print
		p_cpt2g[9]++;
		cout << "one sol to print" << endl;
		char ws[82];
		strcpy(ws, empty_puzzle);
		for (int i = 0; i < genb12.ncluesb12; i++) {
			int cell = genb12.tcluesb12[i];
			ws[cell] = genb12.grid0[cell] + '1';
		}
		for(int i=0,bit=1;i<27;i++,bit<<=1)if(known_b3 & bit)
			ws[54+i] = genb12.grid0[54+i] + '1';
		fout1 << ws << endl;
	}
	else if(p_cpt2g[8]<50){
			if (p_cpt2g[0] == sgo.vx[3])zhou[1].ImageCandidats();
			//char ws[82];
			cout <<Char27out(zh_g2.cells_assigned.bf.u32[2]) << " ua" << endl;
			//cout << zh_g2.cells_assigned.String3X(ws)<<" ua" << endl;
	}
}
void G17B3HANDLER::Go_Not_Critical_missn() {
	if (1) {
		cout << "entry not_critical miss " << nmiss << endl;
	}
	int wua = wactive0, ncells = 27, rawua = 0;
	{  // select ua to use
		register int Ra = wactive0, Rfilt = known_b3;
		for (uint32_t iua = 0; iua < genb12.nuasb3_2; iua++) {
			register int Ru = genb12.uasb3_2[iua];
			if (Ru & Rfilt) continue;
			Ru &= Ra;
			if (!Ru)return; // empty need at least one outfield
			if (nmiss == 1) Ra = Ru;
			else {
				register int cc = _popcnt32(Ru);
				if (cc < ncells) { ncells = cc; wua = Ru; rawua = genb12.uasb3_2[iua]; }
				if (cc > 6 && ncells < 5) break;
			}
		}
	}
	if (1) {
		cout <<Char27out(wua)<< "wua to use "  << endl;
	}
	{ // apply first UA to use or all out field cells 
		uint32_t res;
		int x = wua;
		while (bitscanforward(res, x)) {
			int bit = 1 << res; x ^= bit; wactive0 ^= bit;
			G17B3HANDLER hn = *this; hn.nmiss--; hn.known_b3 |= bit;
			if (hn.nmiss) hn.Go_Not_Critical_missn();
			else 	hn.Go_Critical();
		}
	}
	if (ncells == 27)	Go_Subcritical();// finish in Subcritical if no ua
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (mini_bf2) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[genb12.mini_bf2 & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~genb12.pairs27);// and set the common cell as assigned
			}
		}
		mini_bf2 = 0;
	}
}
int G17B3HANDLER::ShrinkUas1(int * to, int no) {
	irloop = 0;
	uasb3 = &to[no];
	nuasb3 = 0;
	for (int iua = 0; iua < no; iua++) {
		register int Ru = to[iua];
		if (Ru & known_b3) continue;// already hit, forget it
		Ru &= active_b3;
		if (!Ru) return 1;// dead branch
		if (_popcnt32(Ru) == 1) {// assign it and reduce the active cells
			CriticalAssignCell(Ru);
			irloop = 1;// should loop for new singles
		}
		else uasb3[nuasb3++] = Ru;
	}
	if (!nuasb3) irloop = 0;// no need to loop again
	return 0;

}
void G17B3HANDLER::Go_Critical() {// critical situation all clues in pairs tripl:ets
	if (0)cout << Char27out(active_b3) << " active b3 entry critical" << endl;
	Critical2pairs();// assign 2 pairs in minirow to common cell
	rknown_b3 = known_b3 | active_b3;
	if (!active_b3) {
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (IsMultiple(rknown_b3)) {

		return;
	}
	if (p_cpt2g[0] == sgo.vx[3]) {
		for (int i = 54; i < 81; i++)cout << genb12.grid0[i] + 1;
		cout << endl;
		cout << Char27out(known_b3) << " known b3 after crit 2 pairs" << endl;
		cout << Char27out(active_b3) << " active_b3 after crit 2 pairs" << endl;
	}
	if (ShrinkUas1((int *)genb12.uasb3_1, genb12.nuasb3_1)) return;// dead branch
	int wknown = known_b3 | active_b3;
	if (!active_b3) {
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (rknown_b3 != wknown)
		if (IsMultiple(wknown))	return;// not valid using all cells
	rknown_b3 = wknown;
	if (p_cpt2g[0] == sgo.vx[3]) {
		zhou[0].ImageCandidats();
	}
	if (irloop)		CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalLoop(int *to, int no){
	if (ShrinkUas1(to, no)) return;
	if (irloop)CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalExitLoop(int *uasb3, int nuasb3){
	if (p_cpt2g[0] == sgo.vx[3]) {
		cout << "exit critical loop" << endl;
		cout << Char27out(known_b3) << " known b3 exit loop" << endl;
		cout << Char27out(active_b3) << " active_b3 exit loop" << endl;
	}
	int nmissb = genb12.ncluesb3 - _popcnt32(known_b3);// missing clues
	if (nmissb < 0)return;
	if (!active_b3){// nothing more to assign
		if (nuasb3)return; // still not hit uas
		if (nmissb)return;// dead cell in a mini row 3 pairs
		CriticalFinalCheck();
		return;
	}
	// check known + active with brute force
	int wknown = known_b3 | active_b3;
	if (rknown_b3 != wknown) {
		if (IsMultiple(wknown))	return;// not valid using all cells
		rknown_b3 = wknown;
	}
	if (p_cpt2g[0] == sgo.vx[3]) {
		cout << "exit loop nuasb3=" << nuasb3 << " nmissb=" << nmissb << endl;
		PrintStatus(1);
	}
	if (nuasb3){		// find the smallest ua and apply it
		int wua = 0, sizeua = 27;
		uint32_t cell;
		if (nmissb == 1) {//most frequent case
			register int and_uas = active_b3;
			for (int i = 0; i < nuasb3; i++) {
				and_uas &= uasb3[i];
			}
			if (!and_uas) return; // no possibility
			wua = and_uas;
		}
		else if (mini_bf1){	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, mini_bf1);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3&mask;// catch the minirow
		}
		else{
			for (int i = 0; i < nuasb3; i++){
				register int ua = uasb3[i], cc = _popcnt32(ua);
				if (cc < sizeua){ wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum
			}
			if (sizeua >= 2 && mini_triplet){// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, mini_triplet);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 &mask;// catch the minirow

			}
		}
		if (p_cpt2g[0] == sgo.vx[3])		cout << Char27out(wua) << " wua to use" << endl;

		while (bitscanforward(cell, wua)){
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bit);
			if (p_cpt2g[0] == sgo.vx[3])hn.PrintStatus(1);
			hn.CriticalLoop(uasb3, nuasb3);
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void G17B3HANDLER::Critical_0_UA(){
	int nmissb = genb12.ncluesb3 - _popcnt32(known_b3);// missing clues
	if (nmissb < 0)return;
	if (!nmissb){// nothing more to assign (granted at first call in a branch)
		CriticalFinalCheck();
		return;
	}
	if (mini_bf3)	{// in active minirows with 3 pairs, assign 2
		while (mini_bf3){
			uint32_t mini;
			bitscanforward(mini, mini_bf3);
			int shift = 3 * mini, bit = 1 << shift;
			mini_bf3 ^= 1 << mini; //clear bit the mini row is always killed
			active_b3 &= ~(7 << shift); // clear also the bitfield of active cells
			int tp[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
			for (int i = 0; i < 3; i++){
				int * tpi = tp[i];
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bit << tpi[0]);
				hn.CriticalAssignCell(bit << tpi[1]);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	if (mini_bf1){// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, mini_bf1);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		int x = active_b3&mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		mini_bf1 ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++){
			int bb = bit << i;
			if (x&bb){
				G17B3HANDLER hn = *this;
				hn.CriticalAssignCell(bb);
				hn.Critical_0_UA();
			}
		}
		return;
	}
	// now must be active triplet in minirow
	if (mini_triplet){// safety control should always be
		uint32_t mini;
		bitscanforward(mini, mini_triplet);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		active_b3 &= ~mask;// clear the minirow
		mini_triplet ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++){
			int bb = bit << i;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}


void  G17B3HANDLER::Go_Subcritical(int docheck) {
	active_b3 = active_sub = genb12.pairsbf;
	// check first if a global solution  is still possible (likely of no use here)
	if (docheck)if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	int cct = _popcnt32(genb12.pairs27) - ncritical;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	ndead = 0;
	uasb3 =(int*) genb12.uasb3_1;
	nuasb3 = genb12.nuasb3_1;
	Go_SubcriticalMiniRow();// find the first miss
}
void G17B3HANDLER::PrintStatus(int mode) {
	cout << "remaining uas table" << endl;
	if(mode&1)for (int i = 0; i < nuasb3; i++) 		cout << Char27out(uasb3[i]) << endl;
	cout << Char27out(active_b3) << " active_b3" << endl;
	//cout << Char27out(genb12.pairs27) << "genb12.pairs27" << endl;
	cout << Char9out(mini_bf1) << "mini bf1" << endl;
	cout << Char9out(mini_bf2) << "mini bf2" << endl;
	cout << Char9out(mini_bf3) << "mini bf3" << endl;
	cout << Char9out(mini_triplet) << "mini triplet" << endl;

}

void G17B3HANDLER::Go_SubcriticalMiniRow(){
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	int bit = 1<<ndead, mask = 7<<(3*ndead);
	for (int i = ndead; i < 9; i++, bit <<= 1, mask <<= 3){
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & mini_bf1) {// gua2 pair assign both
			G17B3HANDLER hn = *this;
			hn.mini_bf1 ^= bit;
			hn.SubMini( M, mask);
		}
		else if (bit & mini_bf2) {// 2 gua2 pairs assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				G17B3HANDLER hn = *this;
				hn.mini_bf2 ^= bit;
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini( M, mask);
			}
		}
		else if (bit & mini_bf3) {// 3 gua2 pairs assign all
			G17B3HANDLER hn = *this;
			hn.mini_bf3 ^= bit;
			hn.SubMini( M, mask);
		}
		else if (bit & mini_triplet) {// gua3 assign 2 out of 3
			for (int j = 0; j < 3; j++) {
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				G17B3HANDLER hn = *this;
				hn.mini_triplet ^= bit;
				hn.SubMini( M, mask);
			}
		}
		else { // second add in the mini row one residual cell take it
			G17B3HANDLER hn = *this;
			hn.SubMini( M, mask);
		}
	}
}
void G17B3HANDLER::SubMini( int M, int mask){
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added
	active_b3 &= ~mask;
	active_sub ^= M;
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else{	// leave sub critical mode and enter the critical mode
		Critical2pairs();// assign 2 pairs in minirow to common cell
		if (ShrinkUas1((int *)genb12.uasb3_1, genb12.nuasb3_1)) return;// dead branch
		rknown_b3 = known_b3 | active_b3;
		if (IsMultiple(rknown_b3))return;// not valid using all cells
		if (irloop)		CriticalLoop(uasb3, nuasb3);
		else CriticalExitLoop(uasb3, nuasb3);
	}
}
int GEN_BANDES_12::Band2Check() {
	int * zs0 = &gcheck[27];
	n_auto_b1b2 = 0;
	if (n_auto_b1) {
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			int band[27];// morph the band
			for (int i = 0; i < 9; i++) {
				band[i] = p.map[zs0[p.cols[i]]];
				band[i + 9] = p.map[zs0[p.cols[i] + 9]];
				band[i + 18] = p.map[zs0[p.cols[i] + 18]];
			}
			int ir = G17ComparedOrderedBand(zs0, band);
			if (ir == 1)				return 1;
			else if (!ir) {// auto morph b1 b2 store it for later
				t_auto_b1b2[n_auto_b1b2++] = p;
			}
		}
	}

	n_auto_b2b1 = 0;// possible automorph after perm b1b2
	if (i1t16 == ib2check) {// must try perm bands 12 auto morphs
		int b23[3][9];
		for (int i = 0; i < 3; i++) {// morph band1 to band2 minlex
			register int * rrd = b23[i], *rro = &gcheck[9 * i];
			for (int j = 0; j < 9; j++)
				rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
		}
		int ir = G17ComparedOrderedBand(zs0, b23[0]);// is it same as base
		if (ir == 1) 			return 1;
		else if (!ir)// auto morph b1 b2 store it for later
			t_auto_b2b1[n_auto_b2b1++].InitBase(ib2check);
		// must also test all auto morphs b2b1
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {// same automorphs b1 b2
			BANDMINLEX::PERM &pp = t_auto_b1[imorph];
			int b23_a[3][9];
			for (int i = 0; i < 3; i++) {
				register int * rrd = b23_a[i], *rro = b23[i];
				for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
			}
			int ir = G17ComparedOrderedBand(zs0, b23_a[0]);
			if (ir == 1)return 1;
			else if (!ir)// auto morph b1 b2 store it for later
				t_auto_b2b1[n_auto_b2b1++] = pp;
		}
	}
	return 0;
}
int GEN_BANDES_12::Band3Check() {
	{
		//========================== morphs on b1b2 base test
		if (n_auto_b1b2) {// still direct automorphism b1b2
			for (int imorph = 0; imorph < n_auto_b1b2; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1b2[imorph];
				int b23[3][9];
				// direct
				for (int i = 0; i < 3; i++) {// band 3 only
					register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)
					return 1;
			}
		}
		//=========================== perm b1b2 and base test (b1=b2)
		if (n_auto_b2b1) {// possible lower band3 with a perm band1 band2
			int b23[3][9];//first morph to band 2 min lexical
			for (int i = 0; i < 3; i++) {// rows 4 to 9 as of band 2 perm
				register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
				for (int j = 0; j < 9; j++)
					rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
			}
			for (int imorph = 0; imorph < n_auto_b2b1; imorph++) {// then apply auto morphs 
				BANDMINLEX::PERM &pp = t_auto_b2b1[imorph];
				int b23_a[3][9];
				for (int i = 0; i < 3; i++) {
					register int * rrd = b23_a[i], *rro = b23[i];
					for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23_a[0]) == 1)
					return 1;
			}
		}
		if (1) return 0;// forget redundancy b3=b2
		//========================= (b2=b3)#b1  perm b2b3 to consider (direct done)
		if (ib3check == ib2check) {// check b3b2 on  auto morphs b1
			if (gcheck[27] - 1) return 1; // must be '2' in r4c1
			for (int imorph = 0; imorph < n_auto_b1; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1[imorph];
				int b23[6][9];
				for (int i = 0; i < 6; i++) {// rows 4 to 9 from band 2
					register int * rrd = b23[i], *rro = &gcheck[27 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				int ir = G17ComparedOrderedBand(&gcheck[27], b23[3]);
				if (ir == 1)return 1;
				if (ir < 1 && G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)return 1;
			}
		}
	}
	return 0;
}

//==================================

void STD_B416::Initstd() {
	strcpy(band, "12345678945");
	for (int i = 0; i < 11; i++) band0[i] = band[i] - '1';
	for (int i = 0; i < 27; i++)map[i] = i;// no mapping
	dband = 0;
}
void STD_B416::GetBandTable(int i) {
	i416 = i;
	strncpy(&band[11], t416[i], 16);
	band[27] = 0;
	for (int i = 11; i < 27; i++) band0[i] = band[i] - '1';
}
void STD_B416::SetGangster() {
	memset(gangster, 0, sizeof gangster);
	for (int ir = 0, i = 0; ir < 3; ir++)
		for (int ic = 0; ic < 9; ic++, i++)
			gangster[ic] |= 1 << band0[i];
	// build sol per digit and pm per digit at start
	memset(fd_sols, 0, sizeof fd_sols);
	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 3; j++) {
			int cell = 9 * j + i, dig = band0[cell];
			fd_sols[1][dig] |= Zhoucol << i; // add candidates in the column
			fd_sols[0][dig] |= 1 << cell;
		}
	}


}
void STD_B416::InitC10(int i) {
	GetBandTable(i); SetGangster();
	zh1b_g.GetBand(band0, tua);// set zhone_i
}
void STD_B416::InitG12(int i) {
	GetBandTable(i); SetGangster(); GetUAs();
}
void STD_B416::MorphUas() {
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i], ua = 0;
		register uint32_t cc;
		while (bitscanforward(cc, uao)) {
			uao ^= 1 << cc;
			ua |= 1 << map[cc];
		}
		tua[i] = ua;
	}
}

void STD_B416::InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p
	, int iband) {
	i416 = i16;
	dband = 27 * iband;
	GetUAs();
	strncpy(band, ze, 27);
	for (int i = 0; i < 27; i++) band0[i] = band[i] - '1';
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
			map[vr0 + j] = vr + p.cols[j];
	}
	MorphUas();// morph all uas

	SetGangster();
}
void STD_B416::PrintStatus() {
	cout << "band status i=" << i416 << "\tstart=" << dband << endl << "map ";
	for (int i = 0; i < 27; i++)cout << map[i] << " ";
	cout << endl;
	cout << band << endl << "gangster status" << endl;;
	zh1b_g.GetBand(band0, tua);// set zhone_i
	zhone[0].InitOne_std_band();
	zh1b_g.ndigits = 9;
	zhone[0].ImageCandidats(); // gangster status
	cout << "UAs table" << endl;
	for (uint32_t i = 0; i < nua; i++)
		cout << Char27out(tua[i]) << endl;
}
void STD_B1_2::FillMiniDigsMiniPairs(STD_B1_2 & bb) {
	if (0) {
		cout << "gangster other band" << oct << endl;
		for (int i = 0; i < 9; i++)cout << bb.gangster[i] << " ";
		cout << endl << "my gangster  " << endl;
		for (int i = 0; i < 9; i++)cout << gangster[i] << " ";
		cout << dec << endl;
		cout << band << " my band" << endl;
	}

	nvpairs = 0;
	for (int i = 0, j = 0; i < 9; i++, j += 3) {
		int a = (1 << band0[j]), b = (1 << band0[j + 1]), c = (1 << band0[j + 2]);
		mini_digs[i] = a | b | c;
		mini_pairs[j] = b | c;// missing a  relative columns 2,3
		mini_pairs[j + 1] = a | c;// missing b
		mini_pairs[j + 2] = a | b;// missing c 
		int jcol = j % 9;// start col for the mini row
		int * gg = bb.gangster;
		if ((gg[jcol + 1] & c) && (gg[jcol + 2] & b))
			tv_pairs[nvpairs++] = j;
		if ((gg[jcol] & c) && (gg[jcol + 2] & a))
			tv_pairs[nvpairs++] = j + 1;
		if ((gg[jcol] & b) && (gg[jcol + 1] & a))
			tv_pairs[nvpairs++] = j + 2;
	}
	if (0) {
		cout << " valid pairs table ";
		for (int i = 0; i < nvpairs; i++)
			cout << tv_pairs[i] << " ";
		cout << endl;
		cout << " valid pairs pairs " << oct;
		for (int i = 0; i < nvpairs; i++)
			cout << mini_pairs[tv_pairs[i]] << " ";
		cout << dec << endl;
	}
}
int STD_B1_2::ReviseG_triplet(int imini, int ip, STD_B1_2 * bb) {
	int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
	int dcell = 3 * imini, dcol = dcell % 9,
		*cols = &bb->gangster[dcol],
		*myp = tp3f[ip];
	int digit[3], digit_bit[3], pdigit[3], pdigit_bit[3], digp = 0;
	for (int i = 0; i < 3; i++) {// collect digits
		digit[i] = band0[dcell + i];
		int bit = 1 << digit[i];
		digit_bit[i] = bit;
		digp |= bit;
	}
	for (int i = 0; i < 3; i++) {// collect digits perm
		pdigit[i] = digit[myp[i]];
		pdigit_bit[i] = 1 << pdigit[i];
	}
	int *g12 = &zh2b_g.gangster[dcol];
	if (!(pdigit_bit[0] & g12[0]) ||
		!(pdigit_bit[1] & g12[1]) ||
		!(pdigit_bit[2] & g12[2])) return 0;
	for (int ic = 0; ic < 3; ic++)
		bb->revised_g[dcol + ic] ^= (pdigit_bit[ic] | digit_bit[ic]);
	return digp;
}

uint32_t  STD_B1_2::GetMiniData(int index, uint32_t & bcells, STD_B1_2 *bb) {
	//index is cell 0-26 in the band assumed free in the mini-row
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	uint32_t tcells[3] = { 6,5,3 };// corresponding pairs
	int imini = index / 3, dmini = 3 * imini, dcol = dmini % 9,
		perm = index % 3, *pcol = tpcol[perm];
	bcells = tcells[perm] << (3 * imini);
	uint32_t digs = mini_pairs[index];
	//bb->InitRevisedg();// must be done by the caller
	bb->revised_g[dcol + pcol[0]] ^= digs;
	bb->revised_g[dcol + pcol[1]] ^= digs;
	return digs;
}

void STD_B3::InitBand3(int i16, char * ze, BANDMINLEX::PERM & p) {
	InitBand2_3(i16, ze, p);
	memset(&guas, 0, sizeof guas);
	// setup minirows bit fields
	for (int i = 0; i < 9; i++) {
		minirows_bf[i] = 0;
		int * p = &band0[3 * i];
		for (int j = 0; j < 3; j++)
			minirows_bf[i] |= 1 << p[j];
	}

}
/* GUA4 GUA6
.x. .y.
..y .x.

.x. ... .y.   c2a=c2b c3a=c3b
..y .x. ...   3 rows
... .y. .x.

.x. .y. .z.   2 rowssame third
..y .z. .x.
*/
int STD_B3::IsGua(int i81) {
	GEN_BANDES_12::SGUA2 w81 = genb12.tsgua2[i81];
	int * d10 = &band0[w81.col1], *d20 = &band0[w81.col2];
	int d1 = w81.dig1, d2 = w81.dig2, r1, r2, r3, ua, cell1, cell2;
	for (int irow = 0; irow < 3; irow++) {
		int c1 = d10[9 * irow], c2 = d20[9 * irow];
		if (c1 == d1 && c2 == d2) {// gua2
			r1 = r2 = irow;		cell1 = 9 * irow + w81.col1;
			cell2 = 9 * irow + w81.col2;			break;
		}
		if (c1 == d1) { r1 = irow; cell1 = 9 * irow + w81.col1; }
		else if (c2 == d2) { r2 = irow; cell2 = 9 * irow + w81.col2; }
		else r3 = irow;
	}
	ua = (1 << cell1) | (1 << cell2);
	if (r1 == r2) {// gua2
		guas.ua_pair[i81] = ua;
		int i27 = 9 * r1 + w81.i9;// index 0-26 of the pair
		guas.ua2_i27[i81] = i27;
		guas.isguasocket2.Set_c(i81);
		guas.ua2_imini[i81] = 3 * r1 + w81.i9 / 3;
		//cout << "set gua2 i81=" << i81 << " imini=" << 3 * r1 + w81.i9 / 3 
		//	<<	" guas.ua2_i27="<< guas.ua2_i27[i81] << endl;;
		return 1;
	}
	// is it a gua4 gua6 catch data to use
	int tc[2], ntc = 0, digs = w81.digs;//colums with the 2 digits
	int col1, col2, *p1 = &band0[9 * r1], *p2 = &band0[9 * r2];
	for (int i = 0; i < 9; i++) if ((gangster[i] & digs) == digs)
		tc[ntc++] = i;
	for (int icol = 0; icol < 9; icol++, p1++, p2++) {
		if (*p1 == d2)col1 = icol;
		if (*p2 == d1)col2 = icol;
	}


	if (ntc == 2) {// gua6 first type find and store ua
		int cella = 9 * r1 + col1, cellb = 9 * r2 + col2,
			cellc = 9 * r3 + col1, celld = 9 * r3 + col2;
		ua |= (1 << cella) | (1 << cellb) | (1 << cellc) | (1 << celld);
		guas.ua_pair[i81] = ua;
		guas.isguasocket4.Set_c(i81);
		return 4;
	}
	if (ntc) {
		int c = band0[9 * r3 + tc[0]];
		if (c != d1 && c != d2) {// gua4 
			int cella = 9 * r1 + col1, cellb = 9 * r2 + col2;
			ua |= (1 << cella) | (1 << cellb);
			guas.ua_pair[i81] = ua;
			guas.isguasocket4.Set_c(i81);
			return 2;
		}
	}
	// last  is gua6 with d1,r2 ; d2,r1
	if (band0[9 * r1 + col2] == band0[9 * r2 + col1]) {// gua6 second type
		int cella = 9 * r1 + col1, cellb = 9 * r2 + col2,
			cellc = 9 * r2 + col1, celld = 9 * r1 + col2;
		ua |= (1 << cella) | (1 << cellb) | (1 << cellc) | (1 << celld);
		guas.ua_pair[i81] = ua;
		guas.isguasocket4.Set_c(i81);
		return 8;
	}
	return 0;
}
int STD_B3::IsGua3(int i81) {
	GEN_BANDES_12::SGUA3 w81 = genb12.tsgua3[i81];
	int *g = genb12.gang27,//the gangster 
		d1 = g[w81.id1], d2 = g[w81.id2], d3 = g[w81.id3];
	//catch the minirow pattern needed
	int bita = 1 << d1, bitb = 1 << d2, bitc = 1 << d3;
	int mrpat = bita | bitb | bitc;
	int stack = w81.col1 / 3;
	// minirow of the stack must fit
	for (int irow = 0; irow < 3; irow++) {
		int imini = stack + 3 * irow,
			*pmini = &band0[9 * irow + 3 * stack];
		if (mrpat != minirows_bf[imini])continue;
		// possible triplet, must be right digit in right place
		if (d1 != pmini[0] || d2 != pmini[1])continue;
		guas.triplet[imini] = i81;// valid triplet
		guas.triplet_imini[i81] = imini;// valid triplet
		guas.isguasocket3.Set_c(i81);
		guas.ua_triplet[i81] = 7 << (3 * imini);
		guas.ua3_imini[i81] = imini;
		return 1;
	}
	return 0;
}

//==================== collect UAs 2 bands

int GENUAS_B12::Initgen() {// buil start myband1 myband2
	limsize = UALIMSIZE;
	zh2b5_g.sizef5 = UALIMSIZE;
	zh2b5_g.modevalid = 0;
	// prepare zh2b_g___________________________________________
	memcpy(zh2b_g.puz0, myband1.band0, sizeof myband1.band0);
	memcpy(&zh2b_g.puz0[27], myband2.band0, sizeof myband2.band0);
	for (int i = 0; i < 9; i++)
		zh2b_g.gangster[i] = myband1.gangster[i] | myband2.gangster[i];
	zh2b_g.GetBands(myband1.gangster, myband2.gangster);// set sol/pm
	//_______________________________________________________
	nua = 0;// final table of uas bands 12 empty at start

	nuab1b2 = 0;// switch uas band1 and uas band2 to 2X mode
	for (uint32_t i = 0; i < myband1.nua; i++) // collect band 1
		tuab1b2[nuab1b2++] = myband1.tua[i] & BIT_SET_27;
	for (uint32_t i = 0; i < myband2.nua; i++) {// collect band 2
		register uint64_t  R = myband2.tua[i] & BIT_SET_27;
		R <<= 32;
		tuab1b2[nuab1b2++] = R;
	}
	//___________________________ Start collection of uas
	zh2b_g.nua = 0;// new uas 
	for (int i = 0; i < 36; i++) BuildFloorsAndCollectOlds(floors_2d[i]);
	for (int i = 0; i < 84; i++) BuildFloorsAndCollectOlds(floors_3d[i]);
	for (int i = 0; i < 126; i++)BuildFloorsAndCollectOlds(floors_4d[i]);
	for (int i = 0; i < 126; i++) BuildFloorsAndCollectOlds(0x1ff ^ floors_4d[i]);
	//if(g17b.debug17&&g17b.GodebugCheckUas("standars uas")) return 1;
	//==================== collect more uas 6/7 digit 
	//cout << "initial status for UAS bands 1+2 nua="<<nua << endl;
	CollectMore();
	//cout << "after collect more nua=" << nua << endl;
	CollectTriplets();
	//cout << "after collect triplets nua=" << nua << endl;
	CollectMore2minirows();
	//cout << "after collect 2 mini rows nua=" << nua << endl;
	//if (g17b.debug17&&g17b.GodebugCheckUas(" all us")) return 1;
	///DebugUas();
	return 0; // ok
}
int GENUAS_B12::DebugUas() {
	cout << "  debug uas" << endl;

	for (uint32_t i = 0; i < nua; i++) {
		uint64_t w = tua[i];
		int cc = (w >> 59);
		//if(cc==13)
		cout << Char2Xout(w) << " " << cc << " i=" << i << endl;
	}
	cout << " end debug uas" << endl;
	return 0;
}

void GENUAS_B12::BuildFloorsAndCollectOlds(int fl) {
	int diag = 0;
	if (0 && fl == 055) {
		diag = 1;
		zh2b5_g.diag = 1;
	}
	else zh2b5_g.diag = 0;
	if (diag == 1) {
		cout << "entry GENUAS_B12::BuildFloorsAndCollectOlds(int fl)" << endl;
		cout << "floors 0" << oct << fl << dec << endl;
	}
	floors = fl;// for debugging only
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 0);
	//cout << Char2Xout(solved_cells) << " solved cells" << endl;
	if (!solved_cells) return;// one digit solved true
	// now collect UAs not hit by solved cells  
	nuaold = 0;
	{	register uint64_t R = solved_cells;
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	for (uint32_t i = 0; i < nua; i++)
		if (!(R & tua[i])) tuaold[nuaold++] = tua[i];
	}
	if (diag)cout << " nuaold=" << nuaold << endl;
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	//if (1) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< first test
	// check subsets and add to main table
	for (uint32_t i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64&BIT_SET_2X;// be sure to use only relevant bits
		uint64_t cc = _popcnt64(ua);
		if (cc > limsize) continue;
		if (diag) cout << Char2Xout(ua) << "verif cc=" << cc << endl;
		if (CheckOld()) continue;// superset of a previous ua
		if (diag)cout << "goadd nua=" << nua << endl;
		ua |= cc << 59;
		AddUA64(tua, nua, ua);
		if (diag)cout << "retour add nua=" << nua << endl;
	}
	if (diag)cout << " nua =" << nua << endl;

}
void GENUAS_B12::BuilOldUAs(uint32_t r0) {
	//====extract uas not hit in the band where is the mini row
	nuaold = 0;
	register uint64_t R = BIT_SET_27 ^ r0;
	if (ib)R <<= 32;
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	for (uint32_t i = 0; i < nua; i++)
		if (!(R & tua[i])) tuaold[nuaold++] = tua[i];
	//cout<<Char2Xout(R) << "BuilOldUAs nuas=" << nuaold << endl;
}
int GENUAS_B12::CheckOld() {// ua 2x27  + 5 bit length
	uint64_t * t = tuaold;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nuaold; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}
int GENUAS_B12::CheckMain(uint64_t wua) {//subset in the main table
	register uint64_t ua2x = wua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nua; iua++) {
		register uint64_t R = tua[iua];
		if (R > wua) return 0; // no subset
		if (R < wua) {// is it subset
			R &= BIT_SET_2X;
			if ((R&ua2x) == R) return 1;// subset discard
		}
		else if (R == ua) return 1;// same discard
	}
	return 0;
}
void GENUAS_B12::CollectMore() {// special 6 7 digits minirow
	//zh1b_g.tua = tuamore;
	modemore = 2;
	myband1.FillMiniDigsMiniPairs(myband2);
	myband2.FillMiniDigsMiniPairs(myband1);
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (ib = 0; ib < 2; ib++) {
		ba = mybx[ib];
		bb = mybx[1 - ib];
		zh1b_g.GetBand(bb->band0, tuamore);
		int npairs = ba->nvpairs;
		//cout << " CollectMore()vpairs status nvpairs=" << npairs << endl;
		uint32_t bcells1;
		for (int i = 0; i < npairs; i++) {
			int cell1 = ba->tv_pairs[i];
			bb->InitRevisedg();
			digp = ba->GetMiniData(cell1, bcells1, bb);
			uint32_t  R0 = bcells1;
			w0 = R0;// w0 is tha ua part located in band ba
			if (ib) w0 <<= 32; // 2 cells each ua
			BuilOldUAs(R0);
			CollectMoreTry6_7();// then go 6/7
		}
	}
}
void GENUAS_B12::CollectMoreTry6_7() {
	nfloors = 6;// for debugging
	//____________ try 6 digits unsolved in band a
	for (int i6 = 0; i6 < 84; i6++) {
		int fl3 = floors_3d[i6];// , fl6 = 0x1ff ^ fl3;
		if (fl3&digp) continue;// digits must be in the multi floors
		floors = 0x1ff ^ fl3;
		if (zh1b_g.diag)
			cout << "start floors 0" << oct << floors << dec << endl;
		//zhone[0].InitOne_std_band();
		zhone[0].Start_Uas_Mini(floors, digp);
		if (0 && modemore == 4) {
			cout << "gangsters at call apply gangster changes" << oct << endl;
			for (int i = 0; i < 9; i++) cout << bb->gangster[i] << "\t";
			cout << endl;
			for (int i = 0; i < 9; i++) cout << bb->revised_g[i] << "\t";
			cout << dec << endl;
			if (floors == 0374)zh1b_g.diag = 1;
			else 		zh1b_g.diag = 0;
		}
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		//if (modemore == 4) zhone[0].ImageCandidats();
		if (zhone[0].Update6())
			zhone[0].Guess6();
		if (zh1b_g.nua) EndCollectMoreStep();
	}

	nfloors = 7;// for debugging
	//____________ try now 7 digits unsolved in band a
	for (int i7 = 0; i7 < 36; i7++) {
		int fl2 = floors_2d[i7];// , fl7 = 0x1ff ^ fl2;
		if (fl2&digp) continue;// digits must be in the multi floors
		floors = 0x1ff ^ fl2;
		if (zh1b_g.diag)
			cout << "start floors 0" << oct << floors << dec << endl;
		//zhone[0].InitOne_std_band();
		zhone[0].Start_Uas_Mini(floors, digp);
		if (0 && modemore == 4) {
			cout << "gangsters at call apply gangster changes" << oct << endl;
			for (int i = 0; i < 9; i++) cout << bb->gangster[i] << "\t";
			cout << endl;
			for (int i = 0; i < 9; i++) cout << bb->revised_g[i] << "\t";
			cout << dec << endl;
		}
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		//if (modemore == 4) zhone[0].ImageCandidats();
		if (zhone[0].Update7())		zhone[0].Guess7();
		if (zh1b_g.nua) EndCollectMoreStep();
	}
}
void GENUAS_B12::EndCollectMoreStep() {
	int diag = 0;
	//if (zh1b_g.diag) 
	if (diag)	cout << "end collect step nua=" << zh1b_g.nua << endl;
	for (uint32_t i = 0; i < zh1b_g.nua; i++) {
		ua = zh1b_g.tua[i] &= BIT_SET_27;
		if (!ib) ua <<= 32;
		ua |= w0;
		ua &= BIT_SET_2X;// be sure to use only relevant bits
		uint64_t cc = _popcnt64(ua);
		if (diag)cout << Char2Xout(ua) << "\t " << cc << "  ua to check" << endl;
		if (cc > limsize) continue;
		if (CheckOld()) continue;
		cc <<= 59;
		ua |= cc;
		if (diag)cout << "try add" << endl;
		AddUA64(tua, nua, ua);
	}
}
void GENUAS_B12::CollectTriplets() {// special 6 7 digits full minirow
	modemore = 3;
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (ib = 0; ib < 2; ib++) {
		//check possible  triplets  must be three possible digits
		ba = mybx[ib];
		bb = mybx[1 - ib];
		// init the brute force 
		zh1b_g.GetBand(bb->band0, tuamore);
		for (int imini = 0, cell = 0; imini < 9; imini++, cell += 3) {
			for (int ip = 0; ip < 2; ip++) {// 2 false triplets in mini row
				bb->InitRevisedg();
				digp = ba->ReviseG_triplet(imini, ip, bb);
				if (!digp)continue;// not a possible perm
				// ______________________need the cells to assign in band bb
				uint32_t R0 = 7 << (3 * imini);
				w0 = R0;
				if (ib) w0 <<= 32; // 2 cells each ua
				BuilOldUAs(R0);
				CollectMoreTry6_7();// then go 6/7
			}
		}
	}
}
void GENUAS_B12::CollectMore2minirows() {
	modemore = 4;
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (ib = 0; ib < 2; ib++) {
		//if(ib)	zh1b_g.diag = 1;
		ba = mybx[ib];
		bb = mybx[1 - ib];
		zh1b_g.GetBand(bb->band0, tuamore);
		if (zh1b_g.diag) {
			cout << "new 2 minis ib=" << ib << endl;
			zhone[0].CheckSolPerDigit();
		}
		int npairs = ba->nvpairs;
		//cout << " CollectMore2minirows vpairs status nvpairs=" << npairs <<" ib="<<ib<< endl;
		uint32_t bcells1, bcells2;
		for (int i1 = 0; i1 < npairs - 1; i1++) {
			int cell1 = ba->tv_pairs[i1],
				box1 = cellsFixedData[cell1].eb;
			for (int i2 = i1 + 1; i2 < npairs; i2++) {
				int cell2 = ba->tv_pairs[i2],
					box2 = cellsFixedData[cell2].eb;
				if (box1 == box2) continue;
				// 2 mini rrows
				bb->InitRevisedg();
				digp = ba->GetMiniData(cell1, bcells1, bb);
				digp |= ba->GetMiniData(cell2, bcells2, bb);
				uint32_t  R0 = bcells1 | bcells2;
				w0 = R0;// w0 is tha ua part located in band ba
				if (ib) w0 <<= 32; // 2 cells each ua
				BuilOldUAs(R0);
				if (0) {
					cout << Char2Xout(w0) << " try this as 2 pairs digp=0"
						<< oct << digp << dec << " i1=" << i1 << " i2=" << i2 << " cell1=" << cell1
						<< " cell2=" << cell2 << endl;
					cout << "gangsters" << oct << endl;
					for (int i = 0; i < 9; i++) cout << bb->gangster[i] << "\t";
					cout << endl;
					for (int i = 0; i < 9; i++) cout << bb->revised_g[i] << "\t";
					cout << dec << endl;

				}
				//zh1b_g.diag = 1;
				CollectMoreTry6_7();// then go 6/7
				//zh1b_g.diag = 0;
				//_____ get 1  pair + 1 triplet
			}
			//_______________________ get 1  pair + 1 triplet
			for (int ibox2 = 0; ibox2 < 3; ibox2++) {
				if (box1 == ibox2) continue;
				for (int iminirow = 0; iminirow < 3; iminirow++) {
					for (int ip = 0; ip < 2; ip++) {// 2 false triplets in mini row
						int imini = 3 * iminirow + ibox2;
						bb->InitRevisedg();
						digp = ba->ReviseG_triplet(imini, ip, bb);
						if (!digp)continue;// not a possible perm
						digp |= ba->GetMiniData(cell1, bcells1, bb);
						uint32_t R0 = 7 << (3 * imini) | bcells1;
						w0 = R0;// w0 is the ua part located in band ba
						if (ib) w0 <<= 32; // 2 cells each ua
						BuilOldUAs(R0);
						CollectMoreTry6_7();
					}
				}
			}
		}
	}

}

void GEN_BANDES_12::BuildGang9x3() {
	gang27 = gang[0];
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int * d = gang[i], c = gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

}


void ZH2B::Init_2digits_banda(BF64  cellsbf) {
	int txcells[50], nxcells = 0;
	BitsInTable64(txcells, nxcells, cellsbf.bf.u64);
	memset(zh2b_g.Digit_cell_Assigned_init, 0, sizeof zh2b_g.Digit_cell_Assigned_init);
	for (int i = 0; i < nxcells; i++) {
		int xcell = txcells[i], cell = From_128_To_81[xcell],
			digit = zh2b_g.puz0[cell];
		Assign(digit, cell, xcell);
		zh2b_g.Digit_cell_Assigned_init[digit].Set(xcell);
	}
	for (int i = 0; i < 9; i++)  FD[i] &= cells_unsolved |
		zh2b_g.Digit_cell_Assigned_init[i];
}

int ZHOU::CallMultipleB3(ZHOU & o, uint32_t bf, int diag) {
	*this = o;
	BF128 dca[9];
	int digitsbf = zh_g2.digitsbf;
	memcpy(dca, zh_g2.Digit_cell_Assigned, sizeof dca);
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = zh_g2.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 0;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	//cout << "digitsbf 0"<<oct << digitsbf<<dec << endl;
	if (_popcnt32(digitsbf < 8)) return 1;// can not be one solution
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	//ImageCandidats();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0, diag);

	return zh_g.nsol;
}

int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4 = R3 & FD[3][0];
	R3 |= R2 & FD[3][0]; R2 |= R1 & FD[3][0]; R1 |= FD[3][0];
	BF128 R5 = R4 & FD[4][0]; R4 |= R3 & FD[4][0];
	R3 |= R2 & FD[4][0]; R2 |= R1 & FD[4][0]; R1 |= FD[4][0];
	R5 |= R4 & FD[5][0]; R4 |= R3 & FD[5][0];
	R3 |= R2 & FD[5][0]; R2 |= R1 & FD[5][0]; R1 |= FD[5][0];
	R5 |= R4 & FD[5][6]; R4 |= R3 & FD[6][0];
	R3 |= R2 & FD[6][0]; R2 |= R1 & FD[6][0]; R1 |= FD[6][0];
	R5 |= R4 & FD[7][0]; R4 |= R3 & FD[7][0];
	R3 |= R2 & FD[7][0]; R2 |= R1 & FD[7][0]; R1 |= FD[7][0];
	R5 |= R4 & FD[8][0]; R4 |= R3 & FD[8][0];
	R3 |= R2 & FD[8][0]; R2 |= R1 & FD[8][0]; R1 |= FD[8][0];
	if ((cells_unsolved - R1).isNotEmpty()) 	return 1; // empty cells
	R1 -= R2; // now true singles
	R1 &= cells_unsolved; // these are new singles
	if (R1.isEmpty()) {// no single store pairs and more
		zh_g.pairs = R2 - R3;
		zh_g2.triplets = R3 - R4;
		zh_g2.quads = R4 - R5;
		return 0;
	}
	int tcells[80], ntcells = R1.Table3X27(tcells);
	for (int i = 0; i < ntcells; i++) {
		int cell = tcells[i];
		for (int idig = 0; idig < 9; idig++) {
			if (FD[idig][0].On_c(cell)) {
				Assign(idig, cell, C_To128[cell]);
				goto nextr1;
			}
		}
		return 1; // conflict with previous assign within this lot
	nextr1:;
	}
	zh_g.single_applied = 1;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		if (!Unsolved_Count()) return 2;
		if (Apply17SingleOrEmptyCells())	return 0; // locked empty cell or conflict singles in cells
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Guess17(int index, int diag) {
	if (zh_g2.s17_b3_mini) {// look once for mini rows
		zh_g2.s17_b3_mini = 0;
		uint32_t b3=cells_unsolved.bf.u32[2],aig=0;
		if ( p_cpt2g[0] == sgo.vx[3]) {
			cout << "new call" << endl;
			cout << Char27out(b3) << "start mini study" << endl;

		}
		if (!b3)return; // if b3 solved, all is solved
		for (uint32_t i = 0,mask=7; i < 9; i++,mask<<=3) {
			uint32_t mini = b3 & mask;
			if (_popcnt32(mini) < 2) continue;
			aig = 1;
			//cout << Char27out(mini) << " mini " << endl;
			// try this mini row as unsolved
			ZHOU * mynext = (this + 1);
			*mynext = *this;
			for (uint32_t j = 0,bit=1; j < 27; j++,bit<<=1) {
				if (bit&mini)continue;
				if (b3&bit) {
					int cell = j + 54,	digit = zh_g2.grid0[cell];
					mynext->Seta_c(digit, cell);
				}
			}
			//cout << Char27out(mynext->cells_unsolved.bf.u32[2]) << " unsolved b3 " << endl;
			//cout << "appel suite zh_g.go_back="<< zh_g.go_back << endl;
			int ir = mynext->Full17Update();// solve as much as possible
			if (ir==2) continue;
			uint32_t b3_n = mynext->cells_unsolved.bf.u32[2];
			if (!b3_n) continue; // now solved
			//ImageCandidats();
			if (p_cpt2g[0] == sgo.vx[3])mynext->ImageCandidats(); // pour voir
			mynext->Guess17(0,0);
			if (zh_g.nsol) return;
			zh_g.go_back = 0;// see why it is 1
		}
		if (p_cpt2g[0] == sgo.vx[3])cout << "exit mini aig=" << aig << endl;
		if (aig) {
			int ir=Apply17SingleOrEmptyCells();// restore zh_g
			if (p_cpt2g[0] == sgo.vx[3]) {
				ImageCandidats();
				cout << "retour apply17single " << ir <<" zh_g.go_back= "<< zh_g.go_back<< endl;
				char ws[82];
				cout << zh_g.pairs.String3X(ws) << " paires" << endl;
				cout << zh_g2.triplets.String3X(ws) << " triplets" << endl;
			}
			zh_g.go_back = 0;// see why it is 1
		}
	}

	BF128 w = zh_g.pairs;
	if (w.isEmpty())w = zh_g2.triplets;
	if (w.isEmpty())w = zh_g2.quads;
	if (w.isEmpty())w = cells_unsolved;
	{ // select band with more unsolved cells 
		uint32_t nfreecells = 0, nw;
		if (w.bf.u32[0]) {
			nfreecells = _popcnt32(cells_unsolved.bf.u32[0]);
		}
		if (w.bf.u32[1]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[1]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
				}
			}
			else	nfreecells = _popcnt32(cells_unsolved.bf.u32[1]);
		}
		if (w.bf.u32[2]) {
			if (nfreecells) {
				nw = _popcnt32(cells_unsolved.bf.u32[2]);
				if (nw > nfreecells) {
					nfreecells = nw;
					w.bf.u32[0] = 0;
					w.bf.u32[1] = 0;
				}
			}
		}
	}
	int xcell = w.getFirst128(),
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell],
		tdig[10], ndig = 0;
	if (p_cpt2g[0] == sgo.vx[3]) {
		char ws[82];
		cout << w.String3X(ws) << " w for index "<<index << endl;
	}
	// if first step try first false
	if (!index)	ClearCandidate_c(digit, cell);// force false
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig][0].On(xcell))tdig[ndig++] = idig;
	for (int idig = 0; idig < ndig; idig++) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(tdig[idig], cell, xcell);
		mynext->Compute17Next(index + 1, diag);
		if (zh_g.go_back) return;

	}
	if (!index) {
		FD[digit]->Set_c(cell);// restore the candidate
		SetaCom(digit, cell, xcell);
		Compute17Next(index, diag);

	}
}

void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128 & wua = zh_g2.cells_assigned;
			int * sol = genb12.grid0;
			wua.SetAll_0();;
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			zh_g.nsol++;
		}
		zh_g.go_back = 1;// closed anyway
		return;
	}
	Guess17(index, diag);// continue the process
}
void Msp_ReorderBand(char * ze, char * zep )
{
	char temp[9];
	if (ze[0] > ze[9]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[9], sizeof temp);
		memmove(&ze[9], temp, sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[9], sizeof temp);
			memmove(&zep[9], temp, sizeof temp);

		}
	}
	if (ze[0] > ze[18]) {
		memmove(temp, ze, sizeof temp);
		memmove(ze, &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, zep, sizeof temp);
			memmove(zep, &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
	if (ze[9] > ze[18]) {
		memmove(temp, &ze[9], sizeof temp);
		memmove(&ze[9], &ze[18], sizeof temp);
		memmove(&ze[18], temp, sizeof temp);
		if (zep) {
			memmove(temp, &zep[9], sizeof temp);
			memmove(&zep[9], &zep[18], sizeof temp);
			memmove(&zep[18], temp, sizeof temp);
		}
	}
}
void Go_c510() {
	if (!sgo.finput_name) {
		cerr << "missing input file name" << sgo.finput_name << endl; return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze, zout[82], zdiag[82], zediag[82], zes[82], zzs[82], zesf[82];
	zout[81] = zdiag[81] = zediag[81] = zes[81] = zzs[81] = 0;
	zh_g2.zsol = zout;
	while (finput.GetLigne()) {
		if (zhou[0].CheckValidityQuick(ze) != 1) continue;
		for (int i = 0; i < 81; i++) {
			zdiag[C_transpose_d[i]] = zout[i];
			zediag[C_transpose_d[i]] = ze[i];
		}
		int count[6], i2;
		memset(count, 0, sizeof count);
		for (int i = 0; i < 81; i++) {
			if (ze[i] - '.') {
				int band = i / 27, stack = C_stack[i] + 3;
				++count[band];
				++count[stack];
			}
		}
		for (i2 = 0; i2 < 6; i2++) if (count[i2] == 2) goto morph1;
		continue;// not a puzzle with a band 2 clues
	morph1:
		if (i2 > 2) {
			strcpy(ze, zediag);
			strcpy(zout, zdiag);
			i2 -= 3;
			memcpy(count, &count[3], 12);
		}
		switch (i2) {
		case 0:
			memcpy(zes, &ze[0], 27); memcpy(zzs, &zout[0], 27);
			if (count[1] <= count[2]) {
				memcpy(&zes[27], &ze[27], 27); memcpy(&zzs[27], &zout[27], 27);
				memcpy(&zes[54], &ze[54], 27); memcpy(&zzs[54], &zout[54], 27);
			}
			else {
				memcpy(&zes[27], &ze[54], 27); memcpy(&zzs[27], &zout[54], 27);
				memcpy(&zes[54], &ze[27], 27); memcpy(&zzs[54], &zout[27], 27);
			}
			break;
		case 1:
			memcpy(zes, &ze[27], 27); memcpy(zzs, &zout[27], 27);
			if (count[0] <= count[2]) {
				memcpy(&zes[27], &ze[0], 27); memcpy(&zzs[27], &zout[0], 27);
				memcpy(&zes[54], &ze[54], 27); memcpy(&zzs[54], &zout[54], 27);
			}
			else {
				memcpy(&zes[27], &ze[54], 27); memcpy(&zzs[27], &zout[54], 27);
				memcpy(&zes[54], &ze[0], 27); memcpy(&zzs[54], &zout[0], 27);
			}
			break;
		case 2:;
			memcpy(zes, &ze[54], 27); memcpy(zzs, &zout[54], 27);
			if (count[0] <= count[1]) {
				memcpy(&zes[27], &ze[0], 27); memcpy(&zzs[27], &zout[0], 27);
				memcpy(&zes[54], &ze[27], 27); memcpy(&zzs[54], &zout[27], 27);
			}
			else {
				memcpy(&zes[27], &ze[27], 27); memcpy(&zzs[27], &zout[27], 27);
				memcpy(&zes[54], &ze[0], 27); memcpy(&zzs[54], &zout[0], 27);
			}
			break;
		}

		cout << zzs << endl;
		cout << zes << endl;

		// morph it to canonical 
		BANDMINLEX::PERM pr;
		int g0[81];
		for (int i = 0; i < 81; i++)g0[i] = zzs[i] - '1';
		bandminlex.Getmin(g0, &pr);
		char zs[82]; zs[81] = 0;
		cout << pr.i416 << " must be 28" << endl;
		strcpy(zs, "12345678945");
		strncpy(&zs[11], t416[pr.i416], 16);
		strcpy(zesf, empty_puzzle);
		// relabel all
		for (int i = 3, ij = 27; i < 9; i++)
			for (int j = 0; j < 9; j++, ij++) {
				int is = 9 * i + pr.cols[j];
				int c = zzs[is] - '1';
				zs[ij] = pr.map[c] + '1';
				if (zes[is] != '.')zesf[ij] = zs[ij];
			}
		// reorder band1 p18
		for (int i = 0, ij = 0; i < 3; i++)
			for (int j = 0; j < 9; j++, ij++) {
				int is = 9 * pr.rows[i] + pr.cols[j],
					c = zes[is];
				if (c != '.')zesf[ij] = zs[ij];
			}
		Msp_ReorderBand(&zs[27], &zesf[27]);
		Msp_ReorderBand(&zs[54], &zesf[54]);
		fout1 << zs << ";" << zesf << endl;
		fout2 << zesf << endl;
		//return;
	}
}
void Go_c511() {// process an existing file
	if (!sgo.finput_name) {
		cerr << "missing input file name" << sgo.finput_name << endl; return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze;
	int * zs0 = genb12.grid0, npuz = 0;
	while (finput.GetLigne()) {
		zh_g.modevalid = 1;
		zh_g2.grid0 = genb12.grid0;
		zh_g2.zsol = zh_g2.stdfirstsol;
		strncpy(zh_g2.stdfirstsol, ze, 81);
		// search 17 using a file having known  as entry and one 17 given 6 6 5
		npuz++;
		cout << ze << " to process  n=" << npuz << "\t";
		long tdeb = GetTimeMillis();
		for (int i = 0; i < 81; i++)zs0[i] = ze[i] - '1';
		BANDMINLEX::PERM perm_ret;
		bandminlex.Getmin(zs0, &perm_ret);
		if (perm_ret.i416 != 28) continue;
		myband1.InitBand2_3(perm_ret.i416, ze, perm_ret, 0);
		bandminlex.Getmin(&zs0[27], &perm_ret);
		myband2.InitBand2_3(perm_ret.i416, &ze[27], perm_ret, 1);
		bandminlex.Getmin(&zs0[54], &perm_ret);
		genb12.bands3[0].InitBand3(perm_ret.i416, &ze[54], perm_ret);
		genb12.nband3 = 1;
		cout << "band2 N°" << myband2.i416 << "\tband3 N° " << perm_ret.i416 << endl;
		//myband1.DoExpandBand(0);// expand band1
		ze[81] = 0;
		genb12.Go511();
		//return;
	}
	cout << "print final stats" << endl;
	for (int i = 0; i < 40; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}

}

void GEN_BANDES_12::Go511() {
	for (int i = 0; i < 9; i++) {// init columns status
		cold[i] = 0x1ff;
		for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		gangb12[i] = 0x1ff ^ cold[i];
	}
	memcpy(gangcols, cold, sizeof gangcols);

	//=========================== collect UAs  GUAs
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more

	if (genuasb12.Initgen()) return;
	if (0) {
		for (uint32_t i = 0; i < genuasb12.nua; i++)
			cout << Char2Xout(genuasb12.tua[i]) << endl;
			return;
	}
	BuildGang9x3();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s
	genb12.SecondSockets3Setup();// collect GUA3s
	cout << "=======================" << p_cpt2g[0] << endl;
	if (1) {
		cout << "ua bands1+2   \t" << genuasb12.nua << endl;
		cout << "guas socket2  \t" << genb12.ntua2 << endl;
		cout << "guas socket3  \t" << genb12.ntua3 << endl;
		cout << "active socket2\t" << genb12.nactive2 << endl;
		cout << "active socket3\t" << genb12.nactive3 << endl;
	}
	//________________loop on the 9 band1
	p_cpt2g[0]++;
	Go_Sol_Band1();

}