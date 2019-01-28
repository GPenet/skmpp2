#define SIZETGUA 150
struct GEN_BANDES_12{// encapsulating global data 
	STD_B3 bands3[256];
	//STD_B1_2 band1s, band2s;
	int modeb12,aigstop, ndiag, 
		it16, it16_2, imin16_1, imin16_2, imin16_3; 
	int i1t16, i2t16, i3t16,maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[82],tc[6],ntc;
	int skip,last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2,pband3; 
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
	G17TMORE tmore_sockets2[81];
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
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12(){ 
		gang27 = gang[0];
		InitialSockets2Setup(); 
		InitialSockets3Setup();	
	}
	void Init_tmore_sockets2() {
		for (int i = 0; i < 81; i++)
			tmore_sockets2[i].Init();
	}
	void InitialSockets2Setup();// batch level
	void InitialSockets3Setup();// batch level
	void BuildGang9x3();
	void SecondSockets2Setup();// band1+2 level
	void Sockets2SetupForB12(uint64_t cluesbf);
	void SecondSockets3Setup();// band1+2 level
	void GuaCollect(int fl,int diag=0);
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode=0);
	void NewBand1(int iw);
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10=1);
	int DebugFreshUA(uint64_t ua);
	int Debug17(SGUA2 & w);
	//int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3

	//================= UA and GUAs collection
	/*
int Get_I_81x(int c1, int d1, int c2, int d2){// c1,c2 same box
		int stack = C_box[c1], col1 = c1 - 3 * stack, col2 = c2 - 3 * stack;
		int p1 = 0, p2 = 0;
		for (int i = 0; i < 3; i++) if (gang[c1][i] == d1){ p1 = i; break; }
		for (int i = 0; i < 3; i++) if (gang[c2][i] == d2){ p2 = i; break; }
		int i_81 = 27 * stack + 3 * p1 + p2;
		if (col1)i_81 += 18;
		else i_81 += 9 * (col2 - 1);
		return i_81;
	}

	inline int Where3x(int * t, int v){// extract digit gangster index
		if (*t == v) return 0;
		t++;
		if (*t == v) return 1;
		return 2;
	}

	int GetSocketx(int bf, int i3){// UA; band 3 index
		register int cc = _popcnt32(bf);
		if (cc > 3) return -1; // can not be a mini row
		if (cc <2) return -1; // can not be a mini row
		register int mask = 7;
		for (int i = 0; i < 9; i++, mask <<= 3){
			if ((bf&mask) == bf){
				int ibox = i % 3, irow = i / 3;
				int * mybox = &gang27[9 * ibox],
					*tcol1 = mybox, *tcol2 = mybox + 3, *tcol3 = mybox + 6;
				int *bb = bands3[i3].band0, v1, v2, v3, cell1;
				uint32_t cell;
				register int x = bf;
				bitscanforward(cell, x);		x ^= 1 << cell;		v1 = bb[cell];
				cell1 = cell % 3;
				bitscanforward(cell, x);		x ^= 1 << cell;		v2 = bb[cell];
				if (cc == 3){// triplet
					bitscanforward(cell, x);		v3 = bb[cell];
					return 27 * ibox + 9 * Where3(tcol1, v1) +
						3 * Where3(tcol2, v2) + Where3(tcol3, v3);
				}
				// now a pair
				if (cell1)// it is relative cols 1 2
					return 27 * ibox + 18
					+ 3 * Where3(tcol2, v1) + Where3(tcol3, v2);
				if ((cell % 3)<2)// it is relative cols 0 1
					return 27 * ibox + 3 * Where3(tcol1, v1) + Where3(tcol2, v2);
				return 27 * ibox + 9  // it is relative cols 0 2
					+ 3 * Where3(tcol1, v1) + Where3(tcol3, v2);
			}
		}
		return-1;// not a guA2 GUA3 socket
	}
	*/
	

	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);

};
