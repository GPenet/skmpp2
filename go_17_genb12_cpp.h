
#include "go_17sol_tables.h"

struct GEN_BANDES_12{// encapsulating global data 
	STD_B416 bands3[256];
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
	//================== bands 3 and gangster band 3 analysis
	int tband3[256][27], nband3;
	int gangcols[9], gang[9][3], *gang27, gang_digits_cols[9][3],
		tcolok[2], ncolok;
	BF128 bands_pairs, tbands_pairs[256];// 81 bits in 3x27 mode
	BF128 tbands_UA4_6s[256];// 81 bits in 3x27 mode
	int tbands_UA4_6s_pat[256][81];
	BF128 bands_triplets, tbands_triplets[256];// 81 bits in 3x27 mode
	GINT64 tipairs[256][96];
	int tindexUA4s[256][96];// pair id_81 (3x27) to bit 0_8 in the band
	int tindextriplets[256][96];// triplet id_81 (3x27) to bit 0_8 in the band
	int pairs_cols_digits[81][4];
	int triplets_mini_digits[81][4];
	// __________________________  primary UAs tables and creation of such tables
	//TUA64 mtua;// using tua as table
	uint64_t  tua2[3000], tua3[2000];
	int  ntua2, ntua3;// , nyuas_2y;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode=0);
	void NewBand1(int iw);
	void Find_band2B();
	int ValidBand2();
	void M10Find_band3B(int m10=1);
	int FindBand3Unique();//test or  debugging code see the corresponding file
	//================ B creating a catalogue for the 17 search 
	//same as A exchanging bands 2/3

	//================= UA and GUAs collection
	int Get_I_81(int c1, int d1, int c2, int d2){// c1,c2 same box
		int stack = C_box[c1], col1 = c1 - 3 * stack, col2 = c2 - 3 * stack;
		int p1 = 0, p2 = 0;
		for (int i = 0; i < 3; i++) if (gang[c1][i] == d1){ p1 = i; break; }
		for (int i = 0; i < 3; i++) if (gang[c2][i] == d2){ p2 = i; break; }
		int i_81 = 27 * stack + 3 * p1 + p2;
		if (col1)i_81 += 18;
		else i_81 += 9 * (col2 - 1);
		return i_81;
	}

	inline void BuildGang9x3(){
		gang27 = gang[0];
		for (int i = 0; i < 9; i++){// 9 cols to build out of gangcols
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
	inline int Where3(int * t, int v){// extract digit gangster index
		if (*t == v) return 0;
		t++;
		if (*t == v) return 1;
		return 2;
	}

	int GetSocket(int bf, int i3){// UA; band 3 index
		register int cc = __popcnt(bf);
		if (cc > 3) return -1; // can not be a mini row
		if (cc <2) return -1; // can not be a mini row
		register int mask = 7;
		for (int i = 0; i < 9; i++, mask <<= 3){
			if ((bf&mask) == bf){
				int ibox = i % 3, irow = i / 3;
				int * mybox = &gang27[9 * ibox],
					*tcol1 = mybox, *tcol2 = mybox + 3, *tcol3 = mybox + 6;
				int *bb = tband3[i3], v1, v2, v3, cell1;
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

	void InitGangsterBand3();

	void BuilPairinBox(int ib);
	void BuilTripletinBox(int ib);
	void BuildGUA4_6_();
	void Build_GUA4_6s_Band(int iband);

	//void GenUABands12(int limsize);// collect the base set fof UAs for ba,ds 1+2
	//void CollectUA2s();// collect GUA2s
	//void CollectUA2sBox(int ibox);
	//void CollectUA3s();//collect GUA3s
	//void CollectUA3sBox(int ibox);

	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
	// debugging code special call
	//int Afterb1b2(int option = 0);

}genb12;
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
void GEN_BANDES_12::Start(int mode) {
	modeb12 = mode;
	myband1.Initstd();
	zsol[81] = 0;
	nb12 = 0;
}
void GEN_BANDES_12::NewBand1(int iw) {
	i1t16 = iw;
	it16 = tn6_to_416[iw];
	myband1.InitG12(it16);
	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	strcpy(zsol, myband1.band);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	for (int i = 0; i < 9; i++) // init columns status
		cold[i] = 0x1ff ^ myband1.gangster[i];
	zsol[27] = 0;
	myband1.DoExpandBand( 0);// expand and set index
	cout << "i1t16=" << i1t16 << " it16=" << it16 
		<< " n auto morphs=" << n_auto_b1 << endl;
	ntc = 0;
	BitsInTable32(tc, ntc, cold[0]);// first col 6 digits in table
	for (int ip = 0; ip < 20; ip++) {//0;ip<20 setup initial values for rows columns
		GetStartB2(ip);
		Find_band2B();
	}
}
void GEN_BANDES_12::Find_band2B() {
	int * zs0= &grid0[27];
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
		if (crcb[4])if (__popcnt(fr0) < 3) goto back; // 3 clues needed here
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
	it16_2 = pband2.i416;
	i2t16 = t416_to_n6[it16_2];
	if (i2t16 < i1t16)goto next;// not canonical
	if (i2t16 == i1t16)if (it16_2 < it16)goto next;// not canonical

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
			if (ir == 1)				goto next;
			else if (!ir) {// auto morph b1 b2 store it for later
				t_auto_b1b2[n_auto_b1b2++] = p;
			}
		}
	}

	n_auto_b2b1 = 0;// possible automorph after perm b1b2
	if (i1t16 == i2t16) {// must try perm bands 12 auto morphs
		int b23[3][9];
		for (int i = 0; i < 3; i++) {// morph band1 to band2 minlex
			register int * rrd = b23[i], *rro = &grid0[9 * i];
			for (int j = 0; j < 9; j++)
				rrd[j] = pband2.map[rro[pband2.cols[j]]];
		}
		int ir = G17ComparedOrderedBand(zs0, b23[0]);// is it same as base
		if (ir == 1) 			goto next;
		else if (!ir)// auto morph b1 b2 store it for later
			t_auto_b2b1[n_auto_b2b1++].InitBase(i2t16);
		// must also test all auto morphs b2b1
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {// same automorphs b1 b2
			BANDMINLEX::PERM &pp = t_auto_b1[imorph];
			int b23_a[3][9];
			for (int i = 0; i < 3; i++) {
				register int * rrd = b23_a[i], *rro = b23[i];
				for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
			}
			int ir = G17ComparedOrderedBand(&grid0[27], b23_a[0]);
			if (ir == 1)goto next;
			else if (!ir)// auto morph b1 b2 store it for later
				t_auto_b2b1[n_auto_b2b1++] = pp;
		}
	}
	nb12++;
	if (ValidBand2())return;
	goto next;
back:
	if (--ii >= 0) goto next;
}
void Go_c17_91_go();

int GEN_BANDES_12::ValidBand2() {
	myband2.InitBand2_3(i2t16, &zsol[27], pband2);
	//_______________________ std process
	if (modeb12 < 10) {
		if ((nb12 >> 6) < skip) return 0;// here restart value, kept untouched if no band 3 found
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		}
		memcpy(gangcols, cold, sizeof gangcols);
		M10Find_band3B();
		if (0) cout << "band12=" << nb12 << "\tnband3=" << nband3 << endl;
		{// print a restart point every 64 bands 1+2 seen
			uint64_t w = genb12.nb12, w1 = w >> 6;
			w &= 63;
			if (w == 0)cout << "next skip value to use=\t" << w1 << endl;
		}
		if ((nb12 >> 6) >= last)return 1;
		return 0;
	}
	//______________________ testing options 
	switch (modeb12) {
	case 10: {// test ua collection
		if (++p_cptg[0] > sgo.vx[1]) return 1;
		cout << "deux bandes à tester" << endl;
		Go_c17_91_go();
		return 0;
	}
	case 11: {// enumeration test
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		}
		memcpy(gangcols, cold, sizeof gangcols);
		M10Find_band3B(11);
		if (nband3) {
			p_cpt[0]++;
			p_cpt[1] += nband3;
		}
		//cout << "nband3" << nband3 << endl;
		return 0;
	}
	}
	return 0;
}
void GEN_BANDES_12::M10Find_band3B(int m10) {
	BANDMINLEX::PERM pout;
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
	int ir = bandminlex.Getmin(zs0, &pout, 0);
	if (ir < 0) {//would be bug  did not come in enumeration
		cerr << "gen band 3 invalid return Getmin" << endl;
		return;
	}
	int it16_3 = pout.i416;
	i3t16 = t416_to_n6[it16_3];
	if (i3t16 < i1t16)goto next;// not canonical
	if (i3t16 < i2t16)goto next;// not canonical (must be in this case
	//==============================  b1=b2=b3 use minlex check (simplest code, not common)
	if (i1t16 == i3t16 && i3t16 == i2t16) {// 3 bands equal use diagonal test 
		BANDMINLEX::PERM * p = minlexusingbands.pout;
		p[0].InitBase(i1t16);
		p[1] = pband2;
		p[2] = pout;
		if (minlexusingbands.IsLexMinDirect(grid0, i1t16, t_auto_b1, n_auto_b1, ndiag))
			goto next;
		//int box, rows[9], cols[9], out[81];
		//rowminlexcheck(grid0, out, box, rows, cols);
		//for (int i = 0; i < 81; i++)if (out[i] < grid0[i])goto next;
		goto exit_diag;// rowminlex includes diagonal check
	}
	//========================== morphs on b1b2 base test
	if (n_auto_b1b2) {// still direct automorphism b1b2
		for (int imorph = 0; imorph < n_auto_b1b2; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1b2[imorph];
			int b23[3][9];
			// direct
			for (int i = 0; i < 3; i++) {// band 3 only
				register int * rrd = b23[i], *rro = &grid0[54 + 9 * i];
				for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
			}
			if (G17ComparedOrderedBand(&grid0[54], b23[0]) == 1)				goto next;
		}
	}
	//=========================== perm b1b2 and base test (b1=b2)
	if (n_auto_b2b1) {// possible lower band3 with a perm band1 band2
		int b23[3][9];//first morph to band 2 min lexical
		for (int i = 0; i < 3; i++) {// rows 4 to 9 as of band 2 perm
			register int * rrd = b23[i], *rro = &grid0[54 + 9 * i];
			for (int j = 0; j < 9; j++)
				rrd[j] = pband2.map[rro[pband2.cols[j]]];
		}
		for (int imorph = 0; imorph < n_auto_b2b1; imorph++) {// then apply auto morphs 
			BANDMINLEX::PERM &pp = t_auto_b2b1[imorph];
			int b23_a[3][9];
			for (int i = 0; i < 3; i++) {
				register int * rrd = b23_a[i], *rro = b23[i];
				for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
			}
			if (G17ComparedOrderedBand(&grid0[54], b23_a[0]) == 1) goto next;
		}
	}
	//========================= (b2=b3)#b1  perm b2b3 to consider (direct done)
	if (i3t16 == i2t16) {// check b3b2 on  auto morphs b1
		if (grid0[27] - 1) goto next; // must be '2' in r4c1
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			int b23[6][9];
			for (int i = 0; i < 6; i++) {// rows 4 to 9 from band 2
				register int * rrd = b23[i], *rro = &grid0[27 + 9 * i];
				for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
			}
			int ir = G17ComparedOrderedBand(&grid0[27], b23[3]);
			if (ir == 1)goto next;
			if (ir < 1 && G17ComparedOrderedBand(&grid0[54], b23[0]) == 1)goto next;
		}
	}
	//============================= b1=b3 #b2 
	if (minlexusingbands.IsLexMinDiagB(grid0, i1t16, i2t16, i3t16, t_auto_b1, n_auto_b1, ndiag))goto next;

exit_diag:
	//genb12.B3add(pout.i416);
	bands3[nband3++].InitBand2_3(i3t16,&zs[54],pout);
	//valid_bands.bands3[nband3].Init(zs0);
	//memcpy(tband3[nband3++], zs0, 27 * sizeof zs0[0]);
	goto next;
back:
	if (--ii >= 0) goto next;
	if (m10 != 1)return;
	if (nband3)		g17b.GoM10();// call the process for that entry
}


void GEN_BANDES_12::InitGangsterBand3() {//depending on the list of valid band 3
	BuildGang9x3();
	bands_pairs.SetAll_0();  bands_triplets.SetAll_0();
	memset(tbands_pairs, 0, sizeof tbands_pairs);// nb3 * sizeof bm128);
	memset(tbands_triplets, 0, sizeof tbands_triplets);//, nb3 * sizeof bm128);
	for (int ib = 0; ib < 3; ib++) {
		BuilPairinBox(ib); BuilTripletinBox(ib);
	}
}
void GEN_BANDES_12::BuilPairinBox(int ibox) {
	int cols[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
	int * mybox = &gang27[9 * ibox];
	for (int ipcol = 0; ipcol < 3; ipcol++) {
		int col1 = cols[ipcol][0], col2 = cols[ipcol][1];
		int *tcol1 = &mybox[3 * col1], *tcol2 = &mybox[3 * col2];
		for (int d1 = 0; d1 < 3; d1++) {
			for (int d2 = 0; d2 < 3; d2++) {
				int dig1 = tcol1[d1], dig2 = tcol2[d2],
					i_27 = 9 * ipcol + 3 * d1 + d2,
					i_81c = i_27 + 27 * ibox;
				pairs_cols_digits[i_81c][0] = col1 + 3 * ibox;
				pairs_cols_digits[i_81c][1] = col2 + 3 * ibox;
				pairs_cols_digits[i_81c][2] = dig1;
				pairs_cols_digits[i_81c][3] = dig2;
				uint32_t	&boxpairs = bands_pairs.bf.u32[ibox],
					bit = 1 << i_27;
				for (int iband3 = 0; iband3 < nband3; iband3++) {
					int *bb = &tband3[iband3][3 * ibox], r1, r2;
					if (bb[col1] == dig1)r1 = 0;
					else if (bb[col1 + 9] == dig1) r1 = 1;
					else r1 = 2;
					if (bb[col2] == dig2)r2 = 0;
					else if (bb[col2 + 9] == dig2) r2 = 1;
					else r2 = 2;
					if (r1 == r2) {//valid pair for that minirow
						boxpairs |= bit;
						tbands_pairs[iband3].bf.u32[ibox] |= bit;
						int i_81 = i_27 + 32 * ibox,
							start_minirow = 9 * r1 + 3 * ibox,
							bit1 = 1 << (start_minirow + 2 - ipcol);
						tipairs[iband3][i_81].u32[0] = bit1;
						tipairs[iband3][i_81].u32[1] = (7 << start_minirow) ^ bit1;
					}
				}
			}
		}
	}


}
void GEN_BANDES_12::BuilTripletinBox(int ibox) {
	int * mybox = &gang27[9 * ibox],
		*tcol1 = mybox, *tcol2 = mybox + 3, *tcol3 = mybox + 6;
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			for (int d3 = 0; d3 < 3; d3++) {
				int dig1 = tcol1[d1], dig2 = tcol2[d2], dig3 = tcol3[d3],
					digs = (1 << dig1) | (1 << dig2) | (1 << dig3),
					i_27 = 9 * d1 + 3 * d2 + d3;
				uint32_t	&boxtriplets = bands_triplets.bf.u32[ibox],
					bit = 1 << i_27;
				ncolok = 0;
				for (int iband3 = 0; iband3 < nband3; iband3++) {
					int *bb = &tband3[iband3][3 * ibox], r1, r2, r3;
					if (bb[0] == dig1)r1 = 0;
					else if (bb[9] == dig1) r1 = 1;
					else r1 = 2;
					if (bb[1] == dig2)r2 = 0;
					else if (bb[10] == dig2) r2 = 1;
					else r2 = 2;
					if (bb[2] == dig3)r3 = 0;
					else if (bb[11] == dig3) r3 = 1;
					else r3 = 2;
					if (r1 == r2 && r1 == r3) {//valid pair for that minirow
						boxtriplets |= bit;
						tbands_triplets[iband3].bf.u32[ibox] |= bit;
						tindextriplets[iband3][i_27 + 32 * ibox] =
							1 << (3 * r1 + ibox);// directly the minirow
					}
				}
			}
		}
	}
}

void GEN_BANDES_12::BuildGUA4_6_() {
	for (int iband3 = 0; iband3 < nband3; iband3++) {
		//for (int i = 0; i < 27; i++) cout << tband3[iband3][i] + 1;
		//cout << " band 3 studied" << endl;
		Build_GUA4_6s_Band(iband3);
		bands_pairs |= tbands_UA4_6s[iband3];
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
void GEN_BANDES_12::Build_GUA4_6s_Band(int iband) {
	int * tuas46 = tbands_UA4_6s_pat[iband];
	BF128 & ti81 = tbands_UA4_6s[iband];
	ti81.SetAll_0();
	int tdigits[9][3][2]; // temp status for a band, 3 times a digit in three boxes
	int *bb = tband3[iband];
	for (int i = 0; i < 27; i++) {
		int istack = C_box[i], irow = C_row[i], icol = C_col[i], digit = bb[i];
		tdigits[digit][istack][0] = irow;
		tdigits[digit][istack][1] = icol;
	}
	// test the 81 i_81
	int b2[3] = { 1, 2, 0 }, b3[3] = { 2, 0, 1 }; // other boxes
	for (int d1 = 0; d1 < 8; d1++)for (int d2 = d1 + 1; d2 < 9; d2++) {// a pair of digits
		for (int ibox = 0; ibox < 3; ibox++) {
			int * rc1 = tdigits[d1][ibox], *rc2 = tdigits[d2][ibox], i_81;

			if (rc1[0] == rc2[0] || rc1[1] == rc2[1]) continue;
			int ua = 0;
			if (rc1[1] < rc2[1])    i_81 = Get_I_81(rc1[1], d1, rc2[1], d2);
			else                    i_81 = Get_I_81(rc2[1], d2, rc1[1], d1);

			int * rc1_b = tdigits[d1][b2[ibox]], *rc1_c = tdigits[d1][b3[ibox]],
				*rc2_b = tdigits[d2][b2[ibox]], *rc2_c = tdigits[d2][b3[ibox]];
			if (rc1_b[1] == rc2_b[1]) {// can be GUA4 or GUA6 first type
				if (rc1[0] == rc2_b[0] && rc2[0] == rc1_b[0]) {// it is a GUA4
					ua |= 1 << (9 * rc1_b[0] + rc1_b[1]);
					ua |= 1 << (9 * rc2_b[0] + rc2_b[1]);
					goto guaok;
				}
				else if (rc1_c[1] == rc2_c[1]) {// it is a GUA6 first type
					ua |= 1 << (9 * rc1_b[0] + rc1_b[1]);
					ua |= 1 << (9 * rc2_b[0] + rc2_b[1]);
					ua |= 1 << (9 * rc1_c[0] + rc1_c[1]);
					ua |= 1 << (9 * rc2_c[0] + rc2_c[1]);
					goto guaok;

				}
				continue;
			}
			else if (rc1_c[1] == rc2_c[1]) {// can be GUA4
				if (rc1[0] == rc2_c[0] && rc2[0] == rc1_c[0]) {// it is a GUA4
					ua |= 1 << (9 * rc1_c[0] + rc1_c[1]);
					ua |= 1 << (9 * rc2_c[0] + rc2_c[1]);
					goto guaok;

				}
				continue;
			}
			// check now GUA6 second type  first type excluded
			int * cell1, *cell2;
			if (rc1[0] == rc2_b[0]) {// select cells in same rows
				cell1 = rc2_b;
				cell2 = rc1_c;
			}
			else {
				cell1 = rc2_c;
				cell2 = rc1_b;
			}
			if (rc2[0] != cell2[0]) continue; //must be same row
			// and must be same digit filling missing 2 cells
			int diga = bb[9 * cell1[0] + cell2[1]], digb = bb[9 * cell2[0] + cell1[1]];
			if (diga == digb) {//this is a GUA6 type 2
				ua |= 1 << (9 * cell1[0] + cell1[1]);
				ua |= 1 << (9 * cell2[0] + cell2[1]);
				ua |= 1 << (9 * cell1[0] + cell2[1]);
				ua |= 1 << (9 * cell2[0] + cell1[1]);
				goto guaok;
			}
			continue;
		guaok:// a ua has been seen, finish the task
			{
				ua |= 1 << (9 * rc1[0] + rc1[1]);
				ua |= 1 << (9 * rc2[0] + rc2[1]);
				ti81.setBit(C_To128[i_81]);
				tuas46[i_81] = ua;
			}
		}
	}
}
/*
void GEN_BANDES_12::GenUABands12(int limsize) {
	if (limsize > 2000)limsize = 2000;
	mtua.Init(tua, limsize);
	mtua.nua = 0;
	zh2b[0].InitGenUas(grid0);
	zh2b[0].GenUas2();
	zh2b[0].GenUas3();
	zh2b[0].GenUas4();
	zh2b[0].GenUas5();
	zh2b[0].CollectFinal(mtua);
}

void GEN_BANDES_12::CollectUA2s() {//try now to collect UA2s UA4s
	ntua2 = 0;
	zh2b[0].InitUas2();// raz nftemp and set limit size
	for (int ib = 0; ib < 3; ib++)  CollectUA2sBox(ib);
}
void GEN_BANDES_12::CollectUA2sBox(int ibox) {
	int cols[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
	int * mybox = &gang27[9 * ibox],
		my27ii = bands_pairs.bf.u32[ibox];

	for (int ipcol = 0; ipcol < 3; ipcol++) {
		int col1 = cols[ipcol][0], col2 = cols[ipcol][1];
		int *tcol1 = &mybox[3 * col1], *tcol2 = &mybox[3 * col2];
		for (int d1 = 0; d1 < 3; d1++) {
			for (int d2 = 0; d2 < 3; d2++) {
				int dig1 = tcol1[d1], dig2 = tcol2[d2],
					digs = (1 << dig1) | (1 << dig2),
					i_27 = 9 * ipcol + 3 * d1 + d2;
				//======= limit search to active GUA2s + GUA4s +  GUA6s
				if (!(my27ii&(1 << i_27))) continue; //only active 81
				zh2b[0].GenUAs2_minirowpair(col1 + 3 * ibox, col2 + 3 * ibox,
					dig1, dig2, i_27 + 27 * ibox);
				ntua2 = zh2b[0].CollectFinalua2s(tua2,
					1000, ntua2);
			}
		}
	}
}
void GEN_BANDES_12::CollectUA3s() {//try now to collect UA2s UA4s
	ntua3 = 0;
	zh2b[0].InitUas2();// raz nftemp and set limit size
	for (int ib = 0; ib < 3; ib++)  CollectUA3sBox(ib);
}
void GEN_BANDES_12::CollectUA3sBox(int ibox) {
	int cols[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
	int * mybox = &gang27[9 * ibox],
		my27ii = bands_triplets.bf.u32[ibox];
	int col1 = 3 * ibox, col2 = col1 + 1, col3 = col2 + 1;
	int *tcol1 = gang[col1], *tcol2 = gang[col2], *tcol3 = gang[col3];
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			for (int d3 = 0; d3 < 3; d3++) {
				int i_27 = 9 * d1 + 3 * d2 + d3,
					i_81 = 27 * ibox + i_27;
				// add later a filter for activve 81
				if (!(my27ii&(1 << i_27))) continue; //only active 81
				zh2b[0].GenUAs2_minirowtriplet(col1, col2, col3,
					tcol1[d1], tcol2[d2], tcol3[d3], i_81);
				ntua3 = zh2b[0].CollectFinalua2s(tua3, 1000, ntua3);
			}
		}
	}

}
void GEN_BANDES_12::CollectMore() {
	zhone[0].InitUasMore();
	//============= first create base gangsters bands 1 and 2
	memset(bcols, 0, sizeof bcols);
	for (iband = 0; iband < 2; iband++) {// band1 or band 2
		//============= first create base gangsters bands 1 and 2
		int * bx = &grid0[27 * iband], *bc = bcols[iband];
		for (int i = 0; i < 27; i++)	bc[i % 9] |= 1 << bx[i];
	}
	// collect first one minirow
	for (iband = 0; iband < 2; iband++) {// band1 or band 2
		for (ibox = 0; ibox < 3; ibox++) {
			for (iminirow = 0; iminirow < 3; iminirow++) {
				// first minirow pair or triplet
				int dcols = 3 * ibox, dmini = 9 * iminirow + dcols;
				for (int i = 0; i < 4; i++) {
					ncells = 2;
					switch (i) {
					case 3:
						ncells = 3;
						tcells[2] = dmini + 2; tcols[2] = dcols + 2;
					case 0:
						tcells[0] = dmini; tcols[0] = dcols;
						tcells[1] = dmini + 1; tcols[1] = dcols + 1;
						break;
					case 1:
						tcells[0] = dmini; tcols[0] = dcols;
						tcells[1] = dmini + 2; tcols[1] = dcols + 2;
						break;
					case 2:
						tcells[0] = dmini + 1; tcols[0] = dcols + 1;
						tcells[1] = dmini + 2; tcols[1] = dcols + 2;
						break;
					}
					CollectMoreOneMinirow();// try this  minirow pair or triplet
				}
			}
		}
	}
	// and later 2 minirows (have many subsets with one minirow)
	for (iband = 0; iband < 2; iband++) {// band1 or band 2
		for (ibox = 0; ibox < 3; ibox++) {
			for (iminirow = 0; iminirow < 3; iminirow++) {
				// first minirow pair or triplet
				int dcols = 3 * ibox, dmini = 9 * iminirow + dcols;
				for (int i = 0; i < 4; i++) {
					ncells = 2;
					switch (i) {
					case 3:
						ncells = 3;
						tcells[2] = dmini + 2; tcols[2] = dcols + 2;
					case 0:
						tcells[0] = dmini; tcols[0] = dcols;
						tcells[1] = dmini + 1; tcols[1] = dcols + 1;
						break;
					case 1:
						tcells[0] = dmini; tcols[0] = dcols;
						tcells[1] = dmini + 2; tcols[1] = dcols + 2;
						break;
					case 2:
						tcells[0] = dmini + 1; tcols[0] = dcols + 1;
						tcells[1] = dmini + 2; tcols[1] = dcols + 2;
						break;
					}
					int nc0 = ncells;
					for (ibox2 = ibox + 1; ibox2 < 3; ibox2++) {
						for (iminirow2 = 0; iminirow2 < 3; iminirow2++) {
							int dcols2 = 3 * ibox2, dmini2 = 9 * iminirow2 + dcols2;
							for (int j = 0; j < 4; j++) {
								ncells = nc0 + 2;
								switch (j) {
								case 3:
									ncells++;
									tcells[2 + nc0] = dmini2 + 2; tcols[2 + nc0] = dcols2 + 2;
								case 0:
									tcells[nc0] = dmini2; tcols[nc0] = dcols2;
									tcells[1 + nc0] = dmini2 + 1; tcols[1 + nc0] = dcols2 + 1;
									break;
								case 1:
									tcells[nc0] = dmini2; tcols[+nc0] = dcols2;
									tcells[1 + nc0] = dmini2 + 2; tcols[1 + nc0] = dcols2 + 2;
									break;
								case 2:
									tcells[nc0] = dmini2 + 1; tcols[+nc0] = dcols2 + 1;
									tcells[1 + nc0] = dmini2 + 2; tcols[1 + nc0] = dcols2 + 2;
									break;
								}
								CollectMore2Minirows();
							}
						}
					}
				}
			}
		}
	}
}
void GEN_BANDES_12::CollectMoreOneMinirow() {
	int db0 = 27 * iband, db1 = 27 * (1 - iband),
		*b0 = &grid0[db0], *b1 = &grid0[27 * (1 - iband)];
	myfloors = 0;
	mybf = 0;
	// prepare the columns status to go
	memcpy(mycols, bcols[1 - iband], sizeof mycols);
	for (int i = 0; i < ncells; i++) {
		int cell = tcells[i], col = tcols[i], box = C_box[col], dcolbox = 3 * box, colr = col % 3;
		int bitdigit = 1 << b0[cell];
		myfloors |= bitdigit;
		mybf |= (uint64_t)1 << C_To128[cell + db0];
		for (int j = 0; j < 3; j++) {
			int & jcol = mycols[dcolbox + j];
			if (j == colr)	jcol |= bitdigit;
			else jcol &= ~bitdigit;
		}
	}
	zhone[0].GenMoreUas(myfloors, mycols, b1, ncells);
	// and collect the results
	int oldua = mtua.nua;
	zhone[0].CollectFinal(mtua, mybf, iband);
	int newua = mtua.nua;
}
void GEN_BANDES_12::CollectMore2Minirows() {
	int db0 = 27 * iband, *b0 = &grid0[db0], *b1 = &grid0[27 * (1 - iband)];
	myfloors = 0;
	mybf = 0;
	// prepare the columns status to go
	memcpy(mycols, bcols[1 - iband], sizeof mycols);
	for (int i = 0; i < ncells; i++) {
		int cell = tcells[i], col = tcols[i], box = C_box[col], dcolbox = 3 * box, colr = col % 3;
		int bitdigit = 1 << b0[cell];
		myfloors |= bitdigit;
		mybf |= (uint64_t)1 << C_To128[cell + db0];
		for (int j = 0; j < 3; j++) {
			int & jcol = mycols[dcolbox + j];
			if (j == colr)	jcol |= bitdigit;
			else jcol &= ~bitdigit;
		}
	}
	zhone[0].GenMoreUas(myfloors, mycols, b1, ncells);
	zhone[0].CollectFinal(mtua, mybf, iband);
}
*/
