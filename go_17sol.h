
struct XY_EXPAND {// storing valid puzzles in a band and main data  size 32 bytes
	GINT64 cells;// up to 6 cells sorted in the entry order
	uint64_t cellsbf; //cells bit field in 54 bits mode
	GINT64 stacks;// stack count 
	inline void Create5(int * tcells) {
		stacks.u64 = cells.u64 = cellsbf = 0;
		for (int i = 0; i < 5; i++) {
			register int c = tcells[i];
			cellsbf |= (uint64_t)1 << c;
			cells.u8[i] = c;
			stacks.u16[C_stack[c]]++;
		}

	}
	inline void Create6(int * tcells) {
		stacks.u64 = cells.u64 = cellsbf = 0;
		for (int i = 0; i < 6; i++) {
			register int c = tcells[i];
			cellsbf |= (uint64_t)1 << c;
			cells.u8[i] = c;
			stacks.u16[C_stack[c]]++;
		}
	}

} ;

struct SPOT17_NOUAS2_3{//hosting a classical search with stack limit control
	int known, active;
	GINT64  stacks;
	int * tua, nua;
	void newspot(int * oua, int onua, SPOT17_NOUAS2_3 * sp, int cell, int bit = 0);
};

struct G17TMORE{// FIFO table of more for bands 1+2
	uint64_t  t[G17MORESIZE];
	int nt, maxt, curt;
	inline void Init(){ maxt = G17MORESIZE; nt = 0; }
	inline void Add(uint64_t v){//add a new more in FIFO 
		if (nt < maxt){// use next location
			curt = nt;
			t[nt++] = v;
		}
		else{// replace the oldest
			curt++;
			if (curt == maxt)curt = 0;
			t[curt] = v;
		}

	}
	inline void Add_If_New(uint64_t v){// check oldest first
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// from curt to 0
		if ((*Rt) == V)return;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) goto add;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)if ((*Rt) == V)return;
	add:
		Add(v);
	}
	inline int Check(uint64_t v){// check oldest first
		if (!nt) return 0;
		register uint64_t V = v, *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		if (!((*Rt) & V))return 1;
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return 0;
		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)if (!((*Rt) & V))return 1;
		return 0;
	}
	void Print(int modegua){
		register uint64_t *Rt = &t[curt], *Rtl;
	loop1:// form curt to 0
		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}
		if (--Rt >= t)goto loop1;
		if (nt < maxt) return;

		Rtl = &t[curt];
		Rt = &t[maxt];
		while (--Rt>Rtl)		{
			register uint64_t w = *Rt;
			if (modegua)			cout << Char54out(w) << " " << (w >> 56) << endl;
			else cout << Char2Xout(w) << endl;
		}

	}
	void PrintGuasDirect(){
		for (int i = 0; i < nt; i++){
			register uint64_t w = t[i];
			cout << Char54out(w) << " " << (w >> 56) << endl;
		}
	}

	void PrintUasDirect(){
		for (int i = 0; i < nt; i++){
			register uint64_t w = t[i];
			cout << Char2Xout(w) << endl;
		}
	}

}g17more, g17morebig, g17moreg[3], g17moreg0[3];

struct G17TB3GO{
	int ib3, minirows_triplets;
	GINT64 pairs, // 27 bits  pairs holes 27 bits pairs cells(later pairs+triplets)
		triplets,
		count, // 4 items 9 bits/9 minirows 1;2;3 pairs any number of pairs
		countsum,// 4 values (4x16 bits) min clues total and per stack
		countstack;// count of clues per stack {bands 12 + countsum)
	void Debug();
}g17tb3go[256];
struct G17B3HANDLER{
	G17TB3GO wg3;
	int known_b3,rknown_b3, active_b3, active_sub, ndead, wactive0, nmiss, ncritical,
		irloop, *uasb3, nuasb3, wua;
	int diagh;
	// ================== entry in the proces
	void Init(G17TB3GO & wgo);
	void Not_Critical_wactive();
	int IsMultiple(int bf);
	inline int ShrinkUas1(int * to, int no){
		irloop = 0;
		uasb3 = &to[no];
		nuasb3 = 0;
		for (int iua = 0; iua < no; iua++){
			register int Ru = to[iua];
			if (Ru & known_b3) continue;// already hit, forget it
			Ru &= active_b3;
			if (!Ru) return 1;// dead branch
			if (__popcnt(Ru) == 1){// assign it and reduce the active cells
				CriticalAssignCell(Ru);
				irloop = 1;// should loop for new singles
			}
			else uasb3[nuasb3++] = Ru;
		}
		if (!nuasb3) irloop = 0;// no need to loop again
		return 0;
	}
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
	void Go_Not_Critical_miss5();
	void Go_Not_Critical_miss4();
	void Go_Not_Critical_miss3();
	void Go_Not_Critical_miss2();
	void Go_Not_Critical_miss1();
	//==================== process subcritical no cell added outside the GUAs field
	void SubMini(G17B3HANDLER * o, int M, int mask, int stack, int imode, int bit);
	void Go_Subcritical(int docheck=1);
	void Go_SubcriticalMiniRow();
	void Go_SubcriticalMiniRow_End(int stack);
	//===============  debugging 
	void PrintStatus();
	void PrintStatusAfterFirstLoop(int ir, int * uasb3, int nuasb3);
	void Print_Ok3();
	void Print_Ok3GUAs();
	void Print_Ok3_Not_Critical(int ractive_outfield);
	void PrintB3UasOut(int wac);
	void DNC(int wua, int rawua, char * lib);
	// documentation sequence for the preprocessor code
	void FindWua();
	void ApplyWua();
	int SubLoop(int * to, int no);

};


struct G17XY{
	XY_EXPAND xxe, yye;
	uint64_t vxg[G17BLOCGSUA], vyg[G17BLOCGSUA], vw[G17BLOCGSUA+1];// one dummy in vw
	uint64_t vm1xg[3], vm1yg[3], vm2xg[3], vm2yg[3];
	BF128 bands_active_pairs, bands_active_triplets;
	BF128 more_active_pairs, more_active_triplets;
	GINT64 stacksxy;
	int ntb3, b3lim,// number of active bands after Guas and stack filters
		nbx, nby,more3;// band x and band y count
	int tuab3[40], ntuab3;// UAs band3 found inside a ba,d3 cycle
	G17TB3GO wg3;// current band 
	G17B3HANDLER g17hh0;
	int uasb3_1[2000], nuasb3_1, uasb3_2[2000], nuasb3_2;
	//============= process
	void Go_0();// check stacks, collect Guas status
	void Go_Guas_collect();
	void Go_Guas_collectMore();
	uint64_t GetPairs(int ib3);
	int GetTriplets(int ib3);
	void BuildActiveBands();
	int CheckValidBand12();
	int Rebuild();
	void Go_ValideB12();
	void FoutValid17(int bf3, int ib3);
};
struct G17CHUNK{
	int nxc, nyc, n64vua, n64vgua, ix, iy;
	uint64_t vxc[G17CHKX * G17BLOCSUA], vxgc[G17CHKX * G17BLOCGSUA]; // X chunk vectors UAs and GUAs
	uint64_t vyc[G17CHKY * G17BLOCSUA], vygc[G17CHKY * G17BLOCGSUA]; // Y chunk vectors UAs and GUAs
	uint16_t tstore[256 * 256]; // storing XY passing filter 1 will never reach the limit
	XY_EXPAND tyec[G17CHKY], txec[G17CHKX];// current tables
	// room for more vectors guas
	uint64_t vm1xgc[G17CHKX * 3], vm1ygc[G17CHKY * 3]; //X  Y chunk vectors more 1

	void GoChunk();
	void GoChunk_2();
};

struct G17INDEXSTEP{ // one pair 2 clues band1 2 clues band2
	XY_EXPAND xyew1[MAXNIND6], xyew2[MAXNIND6]; //12072 as seen maximum
	uint64_t tvuas[MAXNIND6 * G17BLOCSUA], tvguas[MAXNIND6 * G17BLOCGSUA],
		vcellsuas[54][G17BLOCSUA], vcellsguas[54][G17BLOCGSUA],
		vauas[G17BLOCSUA], vaguas[G17BLOCGSUA], vid81[162][G17BLOCGSUA],
		tua[TUA64_12SIZE], tgua[5000];
	// moreguas partiel vectors
	uint64_t tmgua[2][192], vmcellsguas[54][3], vmid81[162][3],vmuas[3];
	uint64_t tvm1guas[MAXNIND6 * 3],	tvm2guas[MAXNIND6 * 3];


	uint64_t  b1b2_2Xbf, b1b2_54bf, nb1b2_2Xbf, nb1b2_54bf, searched17_54;
	int searched17_band3, searched_nclues1, searched_nclues2;
	BF128 pairsok, tripletsok;// found one empty gua 
	int ntua, ntgua, n64vua, n64vgua,added_more;
	int * t1, *t2; // index used in band1 and band2
	int xfirstdead, x2dead, yfirstdead, y2dead, oldx1, oldy1, oldx2;
	uint64_t dead, dead54;
	int n51, n52, n61, n62, nw1, nw2, ncluesw2;
	//======================================= process
	void Debug();
	void CheckVectors();// debugging sequence
	void CheckUas();// debugging sequence
	// initial for a new step
	inline void StartDead(){ xfirstdead = 0; oldx1 = 0; }
	void GoAddNewUas_sub_step();
	void GoAddNewGUas_sub_step();
	void Init(int i1, int i2);
	int ShrinkUas();// shrink table of UAs for bands1+2
	void ShrinkGuas();// shrink table for gangster uas2 ua3s
	// process summary
	void Do66(); // do 6 clues b1 6 clues b2  
	void Do65(); // do 6 clues b1 5 clues b2  
	void Do56(); // do 5 clues b1 6 clues b2
	void Do_Common();// after xye w1('Y')  and w2('X') have been loaded
	void Do_Common_2();// build  Y tables of vectors
	void Do_Common_3_BuildXvectors();// in fact, in the chunk for 256 X maximum
	void Do_Common_3();// launch chunks 256 x256 in G17CHUNK
};


struct G17B{// hosting the search in 6 6 5 mode combining bands solutions
	int xyexpand_size,debug17,p17;	
	G17B();

	//=====================process
	void GoM10();// end preparation for a given band1+band2+ table band3 pass 656 566
	void Go();
	void GoAddNewUas();
	void GoAddNewGUas();
	void InsertGua(uint64_t * tu, int & ntu, uint64_t gu);
	//================ debugging code
	void PrintEndPuz();
	int Godebug17_1(int i1);
	int Godebug17_2(int i2);
	void Godebug17_3();
	void Godebug17_4(int i1,int i2);

};

//#include "g17_debuggingcode.cpp"
