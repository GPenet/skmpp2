#define SIZETGUA 150
struct GEN_BANDES_12{// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12,go_back,diagmore,
		it16, it16_2, imin16_1, imin16_2, imin16_3; 
	int i1t16, i2t16, i3t16,maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82],tc[6],ntc;
	int gcheck [82], ib2check, ib3check;
	int skip,last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2,pband3,pcheck2,pcheck3; 
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], rowd[6], boxd[6],rowdb3[3],boxdb3[3]; //free digits 
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
		int i_81,imini; // index 0_80 for this 
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
	uint32_t  ntua2, ntua3, nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int tactive2[81], nactive2, tactive3[81], nactive3;
	int   tcolok[2], ncolok;
	BF128 bands_pairs, tbands_pairs[256];// 81 bits in 3x27 mode
	BF128 tbands_UA4_6s[256];// 81 bits in 3x27 mode
	int tbands_UA4_6s_pat[256][81];
	BF128 bands_triplets, tbands_triplets[256];// 81 bits in 3x27 mode
	GINT64 tipairs[256][96];
	int tindexUA4s[256][96];// pair id_81 (3x27) to bit 0_8 in the band
	int tindextriplets[256][96];// triplet id_81 (3x27) to bit 0_8 in the band
	int pairs_cols_digits[81][4];
	int triplets_mini_digits[81][4];
	int ngua6_7, c1, c2, band,floors,digp,i81;
	uint64_t wua0,ua;// partial gua ua to check
	uint64_t tuacheck[100],tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12(){ 
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
	void GuaCollect(int fl,int diag=0);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode=0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10=1);
	int DebugFreshUA(uint64_t ua);
	int Debug17(SGUA2 & w);
	//int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3


	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);

};
