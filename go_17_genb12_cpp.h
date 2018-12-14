void GEN_BANDES_12::SGUA2::Debug(const char * lib) {
	cout << lib << " gua2 \t" << i_81 << " cols " << col1 + 1 << col2 + 1
		<< " digs " << dig1 + 1 << dig2 + 1 << endl;
}
void GEN_BANDES_12::InitialSockets2Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 & w= tsgua2[i];
		w.i_81 = i;
		register int i9 = i/9;
		int tpcol[3][2] = { {1,2},{0,2},{0,1} };
		int rdcol = i9 % 3,dcol=3*(i9/3),*p=tpcol[rdcol];
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
void GEN_BANDES_12::SecondSockets2Setup() {
	ntua2 = 0; nactive2 = 0;
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 & w = tsgua2[i];
		w.dig1 = gang27[w.id1];
		w.dig2 = gang27[w.id2];
		w.digs = (1 << w.dig1) | (1 << w.dig2);
		w.valid = 0;
		//skip if no band3 uses it gua2 gua4 gua6 / build guas 
		for (int iband3 = 0; iband3 < nband3; iband3++) {
			if(bands3[iband3].IsGua(i))// setup guas in band
				w.valid=1;
		}
		if (!w.valid)continue;
		// build revised gangster
		memcpy(w.gangcols, gangcols, sizeof gangcols);
		w.gangcols[w.col1] ^= w.digs;
		w.gangcols[w.col2] ^= w.digs;
		// find guas of the socket
		zh2b_g.InitGangster(gangcols, w.gangcols);
		zh2b5_g.sizef5 = UALIMSIZE - 2;
		zh2b5_g.modevalid = 0;
		w.nua_start = ntua2;
		//if (p_cpt2g[0] ++>2)	continue; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<first test 
		ptua2 = &tua2[ntua2];
		nua2 = nua3 = 0;

		//================== GUA collector 2 bands 
		for (int i = 0; i < 36; i++) {// find UAs 2 digits
			GuaCollect(floors_2d[i]);
		}
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			GuaCollect(floors_3d[i]);
		}
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			GuaCollect(floors_4d[i]);
		}
		for (int i = 0; i < 126; i++) {// find UAs 5digits
			//GuaCollect(0x1f^floors_4d[i]);
		}
		if (nua2) {
			tactive2[nactive2++]=i;
			ntua2 += nua2;
		}
		w.nua_end = ntua2;// store the final count
	}
	cout << "endSecondSockets2Setup ntua2=" << ntua2
		<< " nactive i81=" << nactive2 << endl;
}
void GEN_BANDES_12::GuaCollect(int fl,int diag) {//use revised gangster
	uint64_t solved_cells = zh2b5_g.FindUAsInit(fl, 1);
	if (!solved_cells) return;// one digit solved true
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	if (!zh2b5_g.nuaf5) return;
	if (diag)cout << "collect floors=0" << oct << fl << dec
		<< " nuaf5=" << zh2b5_g.nuaf5 << endl	;
	if (diag > 1) {
		for (int i = 0; i < zh2b5_g.nuaf5; i++) {
			cout<<Char2Xout(zh2b5_g.tuaf5->bf.u64)<<endl;
		}
	}
	// now collect UAs not hit by solved cells  
	genuasb12.nuaold = 0;
	{	register uint64_t R = solved_cells;
	for (uint32_t i = 0; i < genuasb12.nuab1b2; i++)
		if (!(R & genuasb12.tuab1b2[i]))
			genuasb12.tuaold[genuasb12.nuaold++] = genuasb12.tuab1b2[i];
	for (uint32_t i = 0; i < genuasb12.nua; i++)
		if (!(R & genuasb12.tua[i]))
			genuasb12.tuaold[genuasb12.nuaold++] = genuasb12.tua[i];
	}
	if (diag) {
		cout << " nuaold=" << genuasb12.nuaold << endl;
		for (uint32_t i = 0; i < genuasb12.nuaold; i++) {
			cout << Char2Xout(genuasb12.tuaold[i]) << "old" << endl;
		}
	}
	// check subsets and add to main table
	for (int i = 0; i < zh2b5_g.nuaf5; i++) {
		genuasb12.ua = zh2b5_g.tuaf5->bf.u64;
		if (genuasb12.CheckOld()) continue;// superset of a previous ua
		uint64_t cc = _popcnt64(genuasb12.ua);
		genuasb12.ua |= cc << 59;
		if (genuasb12.AddUA64(ptua2, nua2)) {
			//cout << Char2Xout(genuasb12.ua) << "added nua2="<<nua2 << endl;
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

void GEN_BANDES_12::InitialSockets3Setup() {//load permanent data
	for (int i = 0; i < 81; i++) {// initial socket 3
	}
}
void GEN_BANDES_12::SecondSockets3Setup() {
	for (int i = 0; i < 81; i++) {// initial socket 2
		SGUA2 w = tsgua2[i];
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
		if (p_cpt2g[10]++)return 1;//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
		if ((nb12 >> 6) < skip) return 0;// here restart value, kept untouched if no band 3 found
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		}
		memcpy(gangcols, cold, sizeof gangcols);
		Find_band3B();
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
		Find_band3B(11);
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
void GEN_BANDES_12::Find_band3B(int m10) {
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
	bands3[nband3++].InitBand3(i3t16,&zs[54],pout);
	//valid_bands.bands3[nband3].Init(zs0);
	goto next;
back:
	if (--ii >= 0) goto next;
	if (m10 != 1)return;
	if (nband3)		g17b.GoM10();// call the process for that entry
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
