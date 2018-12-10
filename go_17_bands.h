

// standard first band (or unique band)
struct STD_B416 {
	char band[28];
	int band0[27],i416,gangster[9],map[27],dband;
	uint32_t tua[100], nua;//   maximum 81  
	void Initstd();
	void GetBandTable(int i) ;
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void InitC10(int i);
	void InitG12(int i) ;
	void InitBand2_3(int i16,char * ze, BANDMINLEX::PERM & p) ;
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	// 12 115 maximum see band 28
	int index1[30][3], index2[135][3],
		n5,n6,nind[2];// bitfiedl,current index 5 current index 6
	XY_EXPAND xye6[MAXN6], xye5[MAXN5];
	// row solution pattern in digit
	int mini_digs[9],mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  valid_pairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs();
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetCellsBf(int box, int imini, int icase);
	uint32_t GetMiniData(int index,  uint32_t & bcells, STD_B1_2 *bb);
	void DoExpandBand(int dband);// dband 0/27
	void DebugIndex(int ind6 = 0);
};

//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
		nfloors, limstep,map[9], cptdebug;
	//=============== uas collector cluding bands
	int limsize,floors;
	uint64_t  tuaold[1000],// previous non hit uas infinal table of uas for bands 1+2
		tua[TUA64_12SIZE],// 
		tuab1b2[200];// collecting bands uas in 2x mode
	uint32_t nuaold,nua,nuab1b2,
		tuamore[500];
	//_______________uamore control
	STD_B1_2 *ba, *bb;
	int ib,digp;
	uint64_t w0, ua;
	//_____________________ functions collect UAs bands 1+2
	void Initgen();
	int BuildFloorsAndCollectOlds(int fl);
	int CollectNews();
	void EndCollectNewUAs();
	int AddUA64(uint64_t * t, uint32_t & nt);
	inline void AddUA(uint64_t v) {
		ua = v; AddUA64(tua, nua);
	}
	void BuilOldUAs( uint32_t r0);
	int CheckOld();
	void CollectMore();
	void CollectMoreTry6_7();
	void EndCollectMoreStep();
	void CollectTriplets();
	void CollectMore2minirows();
	//_____________________ functions collect UA2s UAs3 socket 

	//	genb12.CollectUA2s();// collect GUA2s
	//	genb12.CollectUA3s();//collect GUA3s

};


