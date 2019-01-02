

// standard first band (or unique band)

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

void STD_B416::InitBand2_3(int i16, char * ze, BANDMINLEX::PERM & p) {
	i416 = i16;
	dband = 27;
	GetUAs();
	strncpy(band, ze, 27);
	for (int i = 0; i < 27; i++) band0[i] = band[i] - '1';
	// create the cell map in output
	for (int i = 0; i < 3; i++) {
		int vr = 9 * p.rows[i], vr0 = 9 * i;
		for (int j = 0; j < 9; j++)
		map[vr0 + j] = vr + p.cols[j];
	}
	// morph all uas
	for (uint32_t i = 0; i < nua; i++) {
		register int uao = tua[i], ua = 0;
		uint32_t cc;
		while (bitscanforward(cc, uao)) {
			uint32_t bit = 1 << cc;
			uao ^= bit;
			ua |= bit;
		}
		tua[i] = ua;
	}
	SetGangster();
}
void STD_B416::PrintStatus() {
	cout << "band status i=" << i416 << "\tstart=" << dband<<endl<<"map ";
	for (int i = 0; i < 27; i++)cout << map[i] << " ";
	cout <<endl;
	cout << band << endl<<"gangster status"<<endl;;
	zh1b_g.GetBand(band0, tua);// set zhone_i
	zhone[0].InitOne_std_band();
	zh1b_g.ndigits = 9;
	zhone[0].ImageCandidats(); // gangster status
	cout << "UAs table" << endl;
	for (uint32_t i = 0; i < nua; i++)
		cout << Char27out(tua[i]) << endl;
}
void STD_B1_2::FillMiniDigsMiniPairs() {
	for (int i = 0, j = 0; i < 9; i++, j += 3) {
		int a = (1 << band0[j]), b = (1 << band0[j + 1]), c = (1 << band0[j + 2]);
		mini_digs[i] = a | b | c;
		mini_pairs[j] = b | c;// missing a  relative columns 2,3
		mini_pairs[j+1] = a | c;// missing b
		mini_pairs[j+2] = a | b;// missing c 
	}
	valid_pairs = 0;
}
int STD_B1_2::ReviseG_triplet(int imini, int ip, STD_B1_2 * bb) {
	int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
	int dcell=3*imini,dcol= dcell % 9,
		*cols = &bb->gangster[dcol],
		*myp= tp3f[ip];
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
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	int dcol = index % 9, bit0 = 1 << dcol,*pcol= tpcol[index%3];
	bcells = (7 << dcol) ^ bit0;
	uint32_t digs= mini_pairs[index];
	bb->revised_g[dcol + pcol[0]] ^= digs;
	bb->revised_g[dcol + pcol[1]] ^= digs;
	return digs;
}
void STD_B1_2::DoExpandBand(int dband) {// find all 5 and 6 clues solutions
	struct SPB3 {// spots to find band 3 minimum valid solutions
		// ====================== constant after initialization
		int  possible_cells, all_previous_cells, active_cells, iuab3;
		//uint64_t cursol;
	}spb3[7], *s3, *sn3;
	s3 = spb3;
	s3->all_previous_cells = 0;
	s3->active_cells = BIT_SET_27;// all cells active
	s3->iuab3 = 0; // copy the start table
	s3->possible_cells = tua[0];
	register XY_EXPAND * t5 = xye5, *t6 = xye6;
	n5 = 0; n6 = 0;
	nind[1] = nind[0] = 0;
	int tcells[10];
	//____________________  here start the search
next:
	uint64_t ispot = s3 - spb3;
	//cout <<ispot<<" ispot"<<endl;
	uint32_t cell;
	if (!bitscanforward(cell, s3->possible_cells))goto back;
	{// apply cell in bitfields
		register int bit = (uint64_t)1 << cell;
		tcells[ispot] = cell + dband;
		s3->possible_cells ^= bit;// clear bit
		register int filter = s3->all_previous_cells | bit,
			ac = s3->active_cells ^ bit;
		sn3 = s3 + 1;
		sn3->all_previous_cells = filter;
		// apply index 1/2
		if (ispot < 2) {
			int ind = nind[ispot]++;
			//cout << "ispot=" << ispot << "\t" << nind[ispot] << endl;
			int * t = index1[ind];
			if (ispot)t = index2[ind];
			t[0] = filter;
			t[1] = n5;
			t[2] = n6;
		}

		sn3->active_cells = s3->active_cells = ac;
		// nextspot:take the next available ua to loop		
		for (int i = s3->iuab3 + 1; i <(int) nua; i++) {
			if (tua[i] & filter)continue;
			if (ispot >= 5) 	goto next;//passing the limit, right place to loop	
			sn3->iuab3 = i;
			sn3->possible_cells = tua[i] & ac;
			s3 = sn3; // switch to next spot
			goto next;
		}

	}
	// no more ua
	if (ispot == 5)		t6[n6++].Create6(tcells);// a 6 clues puzzle
	else {// if below 6 loop  for redundant clues
		int t[32], nt = 0;
		register int ac = s3->active_cells;
		while (bitscanforward(cell, ac)) {// put active cells in table
			register int bit = (uint64_t)1 << cell;
			ac ^= bit;
			t[nt++] = cell + dband;
		}
		if (ispot == 4) {// 5 clues, look for 6
			t5[n5++].Create5(tcells);
			for (int i6 = 0; i6 < nt; i6++) {
				tcells[5] = t[i6];
				t6[n6++].Create6(tcells);
			}
		}
		else if (ispot == 3) { // valid 4 clues  
			for (int i5 = 0; i5 < nt; i5++) {
				tcells[4] = t[i5];
				t5[n5++].Create5(tcells);
				for (int i6 = i5 + 1; i6 < nt; i6++) {
					tcells[5] = t[i6];
					t6[n6++].Create6(tcells);
				}
			}
		}
		else if (ispot == 2) { // valid 3 clues  
			for (int i4 = 0; i4 < nt; i4++) {
				tcells[3] = t[i4];
				for (int i5 = i4 + 1; i5 < nt; i5++) {
					tcells[4] = t[i5];
					t5[n5++].Create5(tcells);
					for (int i6 = i5 + 1; i6 < nt; i6++) {
						tcells[5] = t[i6];
						t6[n6++].Create6(tcells);
					}
				}
			}
		}
		else if (ispot == 1) { // valid 2 clues  
			for (int i3 = 0; i3 < nt - 2; i3++) {
				tcells[2] = t[i3];
				for (int i4 = i3 + 1; i4 < nt - 1; i4++) {
					tcells[3] = t[i4];
					for (int i5 = i4 + 1; i5 < nt; i5++) {
						tcells[4] = t[i5];
						t5[n5++].Create5(tcells);
						for (int i6 = i5 + 1; i6 < nt; i6++) {
							tcells[5] = t[i6];
							t6[n6++].Create6(tcells);
						}
					}
				}
			}
		}
	}

	goto next;
	// going back, for a non empty index, count it back
back:
	if (--s3 >= spb3)goto next;
	// and set last index
	int * t = index1[nind[0]];
	t[1] = n5;
	t[2] = n6;
	t = index2[nind[1]];
	t[1] = n5;
	t[2] = n6;
}
void STD_B1_2::DebugIndex(int ind6) {
	int nn = nind[ind6];
	cout <<"debugindex ind6="<<ind6 <<" nindex="<<nn<<endl;
	for (int i = 0; i <= nn; i++) {
		int *w = (ind6) ? index2[i] : index1[i];
		cout << oct << w[0] << dec << "\t" << w[1] << "\t" << w[2] << endl;
	}
}
void STD_B1_2::PrintShortStatus() {
	cout  << band << "\t\tband main data i0-415="<<i416 << endl;
	cout << "nua    \t" << nua << endl;
	cout << "n5     \t" << n5 << endl;
	cout << "n6     \t" << n6 << endl;
	cout << "n5 ind \t" << nind[0] << endl;
	cout << "n6 ind \t" << nind[1] << endl;

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
	//cout << band << " " << oct << minirows_bf[0]
	//	<< " " << minirows_bf[1] << dec << endl;

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
	int d1 = w81.dig1, d2 = w81.dig2,r1,r2,r3,ua,cell1,cell2;
	for (int irow = 0; irow < 3; irow++) {
		int c1 = d10[9 * irow], c2 = d20[9 * irow];
		if (c1 == d1 && c2 == d2) {// gua2
			r1 = r2 = irow;		cell1 = 9 * irow+ w81.col1;
			cell2 = 9 * irow + w81.col2;			break;
		}
		if (c1 == d1) {	r1 = irow; cell1 = 9 * irow + w81.col1;		}
		else if (c2 == d2) {r2 = irow; cell2 = 9 * irow + w81.col2;	}
		else r3 = irow;
	}
	ua = (1 << cell1) | (1 << cell2);
	if (r1 == r2) {// gua2
		guas.ua_pair[i81] = ua;
		int i27 = 9 * r1 + w81.i9;// index 0-26 of the pair
		guas.pairs[i27] = i81;
		guas.isguasocket2.Set_c(i81);
		guas.ua2_imini[i81] = 3 * r1 + w81.i9 / 3;
		return 1;
	}
	// is it a gua4 gua6 catch data to use
	int tc[2], ntc = 0,digs=w81.digs;//colums with the 2 digits
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
		d1= g[w81.id1],d2= g[w81.id2],d3= g[w81.id3];
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

void STD_B3::PrintB3Status() {
	cout << "band3 status" << endl;
	for (int i = 0, ij = 0; i < 3; i++) {
		for (int j = 0; j < 9; j++, ij++) cout << band0[ij] + 1;
		cout << endl;
	}
	cout << endl<<"gua2 gua4 gua6s" << endl;
	for (int i = 0; i < 81; i++){
		int w = guas.ua_pair[i];
		if (w) cout << Char27out(w) << " i81=" << i << endl;
	}
	cout  << "gua3s" << endl;
	for (int i = 0; i < 81; i++) {
		int w = guas.ua_triplet[i];
		if (w) cout << Char27out(w) << " i81=" << i << endl;
	}
	char ws[129];
	const char* w3 = "123456789...---...123456789";
	cout << w3 << w3 << w3 << endl;;
	cout << guas.isguasocket2.String128(ws) << " sock2" << endl;
	cout << guas.isguasocket4.String128(ws) << " sock4" << endl;
	cout << guas.isguasocket3.String128(ws) << " sock3" << endl;
}

//==================== sockets UA2s UA3s control

void GENUAS_B12::Initgen() {// buil start myband1 myband2
	limsize = UALIMSIZE;
	zh2b5_g.sizef5= UALIMSIZE;
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
		tuab1b2[nuab1b2++]= myband1.tua[i]&BIT_SET_27;
	for (uint32_t i = 0; i < myband2.nua; i++) {// collect band 2
		register uint64_t  R= myband2.tua[i] & BIT_SET_27;
		R <<= 32;
		tuab1b2[nuab1b2++] = R;
	}
	//___________________________ Start collection of uas
	zh2b_g.nua = 0;// new uas 
	for (int i = 0; i < 36; i++) BuildFloorsAndCollectOlds(floors_2d[i]);
	for (int i = 0; i < 84; i++) BuildFloorsAndCollectOlds(floors_3d[i]);
	for (int i = 0; i < 126; i++)BuildFloorsAndCollectOlds(floors_4d[i]);
	for (int i = 0; i < 126; i++) BuildFloorsAndCollectOlds(0x1ff ^ floors_4d[i]);
	//==================== collect more uas 6/7 digit 
	CollectMore();
	CollectTriplets();
	CollectMore2minirows();
}
void GENUAS_B12::BuildFloorsAndCollectOlds(int fl) {
	int diag = 0;
	if (0 && fl == 0352) {
		diag = 1;
		zh2b5_g.diag = 1;
	}
	else zh2b5_g.diag = 0;
	if (diag == 1) {
		cout << "entry GENUAS_B12::BuildFloorsAndCollectOlds(int fl)" << endl;
		cout << "floors 0" << oct << fl << dec << endl;
	}
	floors = fl;// for debugging only
	uint64_t solved_cells= zh2b5_g.FindUAsInit(fl, 0); 
	//cout << Char2Xout(solved_cells) << " solved cells" << endl;
	if (!solved_cells) return;// one digit solved true
	// now collect UAs not hit by solved cells  
	nuaold = 0;
	{	register uint64_t R = solved_cells ;
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	for (uint32_t i = 0; i < nua; i++)
		if (!(R & tua[i])) tuaold[nuaold++] = tua[i];
	}
	if(diag)cout <<  " nuaold="<<nuaold << endl;
	zh2b5_g.CollectUas5();// collect uas for this set of floors
	//if (1) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< first test
	// check subsets and add to main table
	for (int i = 0; i < zh2b5_g.nuaf5; i++) {
		ua = zh2b5_g.tuaf5[i].bf.u64;
		uint64_t cc = _popcnt64(ua);
		if (cc > limsize) continue;
		if (diag) cout << Char2Xout(ua) << "verif cc=" << cc << endl;
		if (CheckOld()) continue;// superset of a previous ua
		if (diag)cout << "goadd nua="<<nua << endl;
		ua |= cc << 59;
		AddUA64(tua, nua);
		if (diag)cout << "retour add nua=" << nua << endl;
	}
	if (diag)cout << " nua =" << nua << endl;

}
void GENUAS_B12::BuilOldUAs( uint32_t r0) {
	//====extract uas not hit in the band where is the mini row
	nuaold = 0;
	register uint64_t R = BIT_SET_27 ^ r0;
	if (ib)R <<= 32;
	for (uint32_t i = 0; i < nuab1b2; i++)
		if (!(R & tuab1b2[i])) tuaold[nuaold++] = tuab1b2[i];
	for (uint32_t i = 0; i < nua; i++)
		if (!(R & tua[i])) tuaold[nuaold++] = tua[i];

}
int GENUAS_B12::CheckOld(  ) {// ua 2x27  + 5 bit length
	uint64_t * t = tuaold;
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nuaold; iua++) {
		register uint64_t R = t[iua];
		R &= BIT_SET_2X;
		if ((R&ua2x) == R)		return 1;// we have a subset
	}
	return 0;
}
int GENUAS_B12::AddUA64(uint64_t * t, uint32_t & nt) {// ua 2x27  + 5 bit length
	register uint64_t ua2x = ua & BIT_SET_2X;
	for (uint32_t iua = 0; iua < nt; iua++) {
		register uint64_t R = t[iua];
		if (R < ua) {// is it subset
			R &= BIT_SET_2X;
			if ((R&ua2x) == R) return 0;// we have a subset
		}
		else if (R == ua) return 0;
		else {
			for (uint32_t jua = nt; jua > iua; jua--)t[jua] = t[jua - 1];
			t[iua] = ua;// new inserted
			nt++;
			// is it a subset of a previous entry
			for (iua++; iua < nt; iua++) {
				if ((t[iua] & ua2x) == ua2x) {// we have a subset
					for (uint32_t k = iua + 1; k < nt; k++)t[k - 1] = t[k];
					nt--;
					iua--; //continue same position
				}
			}
			return 2;
		}
	}
	t[nt++] = ua;// added
	return 1;
}
void GENUAS_B12::CollectMore() {// special 6 7 digits minirow
	zh1b_g.tua = tuamore;
	myband1.FillMiniDigsMiniPairs();
	myband2.FillMiniDigsMiniPairs();
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (ib = 0; ib < 2; ib++) {
		//if (ib) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//check possible pairs and triplets  can be 2 or three possible digits
		ba = mybx[ib]; 
		bb = mybx[1 - ib];
		for (int imini = 0, cell = 0; imini < 9; imini++, cell += 3) {
			//if (imini != 7) continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			int dcol = cell % 9, *cols =& bb->gangster[dcol];
			int * minip = &ba->mini_pairs[cell];
			int pcol[3][2] = { {1,2},{0,2},{0,1} };// cols for each relative pair
			for (int ip = 0; ip < 3; ip++) {// 3 pairs in mini row
				//if (ip != 2) continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

				int col1= pcol[ip][0],col2= pcol[ip][1],
					g1= cols[col1],g2= cols[col2];
				digp = minip[ip];
				if ((!(g1&digp))|| (!(g2&digp)))continue;

				// must not be a 4 cells UA (all subsets)check pairs in band bb
				int * bpairs = &bb->mini_pairs[dcol + ip]; //first pair
				if ((bpairs[0] == digp) || (bpairs[9] == digp) || (bpairs[18] == digp)) continue;

				// this is a pair in minirow to test for UAs in band bb
				ba->valid_pairs |= 1 << (cell+ip); //mark it for later

				// ______________________need the 2 cells to assign in band bb
				int cell1 = cell + pcol[ip][0], cell2 = cell + pcol[ip][1];
				int dig1 = ba->band0[cell1], dig2 = ba->band0[cell2];
				uint32_t  R0 = ((1 << cell1) | (1 << cell2));
				w0 = R0;
				if (ib) w0 <<= 32; // 2 cells each ua
				BuilOldUAs(R0);
				
				// init the brute force 
				zh1b_g.GetBand(bb->band0, tuamore);
				zhone[0].InitOne_std_band();
				//zh1b_g.diag = 2;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
				// change the column for dig1 and dig2 in the revised ganster bb
				bb->InitRevisedg();
				bb->revised_g[dcol+col1] ^= digp; 
				bb->revised_g[dcol+col2] ^= digp;
				if (0) {//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					cout << "gangster status octal" << oct << endl;
					for (int i = 0; i < 9; i++)
						cout << bb->gangster[i] << "\t" << bb->revised_g[i] << endl;

					cout << dec << endl;
				}
				CollectMoreTry6_7();// then go 6/7
			}
		}
	}
}
void GENUAS_B12::CollectMoreTry6_7() {
	//____________ try 6 digits unsolved in band a
	for (int i6 = 0; i6 < 84; i6++) {
		int fl3 = floors_3d[i6], fl6 = 0x1ff ^ fl3;
		//if (fl3 != 0250)continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		if (fl3&digp) continue;// digits must be in the multi floors
		zhone[0].InitOne_std_band();
		zhone[0].Start_Uas_Mini(fl6, digp);
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		//zhone[0].ImageCandidats();//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		zhone[0].Guess6();
		if (zh1b_g.nua) EndCollectMoreStep();
	}
	//if (1) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//____________ try now 7 digits unsolved in band a
	for (int i7 = 0; i7 < 36; i7++) {
		int fl2 = floors_2d[i7], fl7 = 0x1ff ^ fl2;
		if (fl2&digp) continue;// digits must be in the multi floors
		zhone[0].InitOne_std_band();
		zhone[0].Start_Uas_Mini(fl7, digp);
		zhone[0].ApplyGangsterChanges(bb->gangster, bb->revised_g);
		zh1b_g.nua = 0;
		zhone[0].InitGuess();
		zhone[0].Guess7();
		if (zh1b_g.nua) EndCollectMoreStep();
	}
}
void GENUAS_B12::EndCollectMoreStep() {
	int diag = 0;
	//if (zh1b_g.diag) 
	if(diag)	cout << "end collect step nua=" << zh1b_g.nua << endl;
	for (uint32_t i = 0; i < zh1b_g.nua; i++) {
		ua = zh1b_g.tua[i] &= BIT_SET_27;
		if (!ib) ua <<= 32;
		ua |= w0;
		uint64_t cc = _popcnt64(ua);
		if (diag)cout << Char2Xout(ua) << "\t " << cc << "  ua to check" << endl;
		if (cc > limsize) continue;
		if (CheckOld()) continue;
		cc <<= 59;
		ua |= cc;
		if (diag)cout << "try add" << endl;
		//if (AddUA64(tua, nua))
		//		cout << Char2Xout(ua) << "\t " << (ua >> 59) << " final more ua added cycle" << endl;
		AddUA64(tua, nua);
	}
}
void GENUAS_B12::CollectTriplets() {// special 6 7 digits full minirow
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for ( ib = 0; ib < 2; ib++) {
		//check possible  triplets  must be three possible digits
		ba = mybx[ib];
		bb = mybx[1 - ib];
		// init the brute force 
		zh1b_g.GetBand(bb->band0, zh1b_g.tua);
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
				CollectMoreTry6_7();
			}
		}
	}
}
void GENUAS_B12::CollectMore2minirows() {
	STD_B1_2 * mybx[2] = { &myband1 ,&myband2 };
	for (int ib = 0; ib < 2; ib++) {
		STD_B1_2 *ba = mybx[ib], *bb = mybx[1 - ib];
		int vpairs = ba->valid_pairs,rrgang[9];
		//cout << Char27out(vpairs) << " vpairs status" << endl;
		uint32_t bcells1, bcells2;
		//_______________________ get 2  pairs
		for (int ibox = 0; ibox <3; ibox++) {
			for (int ic = 0; ic < 9; ic++) {// cell in box 1 index to valid_pairs
				int cell1 = cellsInGroup[18 + ibox][ic];
				if (!(vpairs & (1 << cell1))) continue;
				for (int ibox2 = ibox + 1; ibox2 < 3; ibox2++) {
					for (int ic2 = 0; ic2 < 9; ic2++) {// cell in box 1 index to valid_pairs
						int cell2 = cellsInGroup[18 + ibox2][ic2];
						if (!(vpairs & (1 << cell2))) continue;
						// now 2 pairs index cell1 cell2 to try
						if(0)cout << "2 pairs go ib=" << ib << " "
							<< cellsFixedData[cell1].pt << " "
							<< cellsFixedData[cell2].pt << endl;
						bb->InitRevisedg();
						digp = ba->GetMiniData(cell1, bcells1, bb);
						digp |= ba->GetMiniData(cell2, bcells2, bb);
						uint32_t  R0 = bcells1 | bcells2;
						w0 = R0;
						if (ib) w0 <<= 32; // 2 cells each ua
						BuilOldUAs(R0);
						// init the brute force
						zh1b_g.GetBand(bb->band0, tuamore);
						zhone[0].InitOne_std_band();
						CollectMoreTry6_7();// then go 6/7
					}
				}
				//if (1) continue;

				//_______________________ get 1  pair + 1 triplet
				for (int ibox2 = 0; ibox2 < 3; ibox2++) {
					if (ibox == ibox2) continue;
					for (int iminirow = 0; iminirow < 3; iminirow++) {
						for (int ip = 0; ip < 2; ip++) {// 2 false triplets in mini row
							int imini = 3 * iminirow + ibox2;
							bb->InitRevisedg();
							digp = ba->ReviseG_triplet(imini, ip, bb);
							if (!digp)continue;// not a possible perm
							digp |= ba->GetMiniData(cell1, bcells1, bb);
							if(0)cout << "possible pair + triplet imini=" << imini
								<< " ib=" << ib << " cell1=" << cell1
								<< " digp=0" << oct << digp << dec << endl;
							// ______________________need the cells to assign in band bb
							uint32_t R0 = 7 << (3 * imini) | bcells1;
							w0 = R0;
							if (ib) w0 <<= 32; // 2 cells each ua
							BuilOldUAs(R0);
							if (0) {
								cout << "temp ua table nua=" << nuaold << endl;
								for (uint32_t i = 0; i < nuaold; i++)
									cout << Char2Xout(tuaold[i]) << endl;
								cout << "end temp ua table" << endl;
							}
							//continue;
							// init the brute force
							zh1b_g.GetBand(bb->band0, tuamore);
							zhone[0].InitOne_std_band();
							CollectMoreTry6_7();
						}
					}
				}
			}
		}
		if (1) continue;// not sure it is worth to do it
		//_______________________ get 2 triplets
		for (int ibox = 0; ibox < 3; ibox++) {
			for (int iminirow = 0; iminirow < 3; iminirow++) {// first  triplet
				for (int ip = 0; ip < 2; ip++) {// 2  triplets in mini row
					bb->InitRevisedg();
					digp = ba->ReviseG_triplet(iminirow, ip, bb);
					if (!digp)continue;// not a possible perm
					//store revised g
					memcpy(rrgang, bb->revised_g, sizeof rrgang);
					for (int ibox2 = 0; ibox2 < 3; ibox2++) {
						if (ibox2 == ibox)continue;
						for (int iminirow2 = 0; iminirow2 < 3; iminirow2++) {
							for (int ip2 = 0; ip2 < 2; ip2++) {
								memcpy(bb->revised_g, rrgang,  sizeof rrgang);
								int digp2 = ba->ReviseG_triplet(iminirow, ip, bb);
								if (!digp2)continue;// not a possible perm
								digp |= digp2;// defining the compulsory digits
								uint32_t R0 = (7 << (3 * iminirow)) | (7 << (3 * iminirow2));
								w0 = R0;
								if (ib) w0 <<= 32; // 2 cells each ua
								BuilOldUAs(R0);
								// init the brute force
								zh1b_g.GetBand(bb->band0, tuamore);
								zhone[0].InitOne_std_band();
								CollectMoreTry6_7();
							}
						}
					}
				}
			}
		}
	}
}
//=============== ua cillector socket 2

void GENUAS_B12::ProcessSocket2(int i81) {
	GEN_BANDES_12::SGUA2 &wi81 = genb12.tsgua2[i81];
	zh2b_g.InitGangster(genb12.gangcols, wi81.gangcols);

}
