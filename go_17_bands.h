

// standard first band (or unique band)
struct STD_B416 {
	char band[28];
	int band0[27],i416,gangster[9],map[27],dband;
	uint32_t tua[100], nua;//   maximum 81  
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();
	void GetBandTable(int i) ;
	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas();
	void InitC10(int i);
	void InitG12(int i) ;
	void InitBand2_3(int i16,char * ze, BANDMINLEX::PERM & p
		,int iband=1) ;
	void PrintStatus();
};
struct STD_B1_2 :STD_B416 {
	// 12 115 maximum see band 28
	int index1[30][3], index2[135][3],index3 [2000][2],
		n5,n6,nind[3];// bitfiedl,current index 5 current index 6
	XY_EXPAND xye6[MAXN6], xye5[MAXN5];
	// row solution pattern in digit
	int mini_digs[9],mini_pairs[27],
		revised_g[9];// after false forced in 1/2 minirows
	int  tv_pairs[27],nvpairs; //9 or 27 bits 
	void FillMiniDigsMiniPairs(STD_B1_2 & bb);
	inline void InitRevisedg() {
		memcpy(revised_g, gangster, sizeof gangster);
	}
	int ReviseG_triplet(int imini, int ip, STD_B1_2 * bb);
	uint32_t GetMiniData(int index,  uint32_t & bcells, STD_B1_2 *bb);
	void DoExpandBand(int dband);// dband 0/27
	void DebugIndex2();
	void Debug_2_3();
	void PrintShortStatus();
};

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
	void PrintB3Status();
};

//================== UA collector 2 bands 

struct GENUAS_B12 {// for uas collection in bands 1+2 using brute force 
	int dig_cells[9][9],
		gangbf[9],// columns for band 3 in bit field
		revised_gangbf[9],// same revised UA2s UA3s ***
		mini_digs[9], mini_pairs[27], // UA2s UA3  ***
		//valid_pairs, //  27 bits valid sockets UA2s ***
		nfloors, limstep,map[9], cptdebug,modemore;
	BF128 valid_sockets;

	//=============== uas collector 
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
	int Initgen();
	void BuildFloorsAndCollectOlds(int fl);
	//int AddUA64(uint64_t * t, uint32_t & nt);
	inline void AddUA(uint64_t v) {
		ua = v; AddUA64(tua, nua,ua);
	}
	inline void AddUACheck(uint64_t v) {
		if (nua >= TUA64_12SIZE) nua = TUA64_12SIZE - 1;
		ua = v; AddUA64(tua, nua,ua);
	}
	void BuilOldUAs( uint32_t r0);
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
};


