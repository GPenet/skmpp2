//#define VALc15 13438
//#define VALc17 1032
//#define TESTXY 0
//2536435903110145  i1=2 i2=35
//#define TESTXY2 0
//22521159232789761
//536435903110145
#define LIM3Y 2000000
//#define DEBUGLEVEL 10
void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	const char * diagband = "265738914348591267791264538";
	const char * diagpuz = ".23.........18...6..9............9.4......2..7...6.......9.2.7.8.4...........3...";
	//diag = 2; p17diag.SetAll_0();
	uint64_t diagval = 795388;
	p_cpt2g[0] ++;
	p_cpt2g[7] +=genb12.nband3;
	if (diag) {
		if (genb12.nb12 == diagval) {
			for (int i = 0; i < 54; i++)
				cout << genb12.grid0[i] + 1;
			cout << "entry m10 nb12=" << genb12.nb12 << endl;
			//if (strcmp(diagband, myband2.band)) return;
			cout << "this is the band in diag" << endl;
			if (diag == 2) {
				cout << diagpuz << " puz known" << endl;
				for (int i = 0; i < 81; i++) if (diagpuz[i] != '.')
					p17diag.Set_c(i);
			}
		}
		else return;
	}
	if (0) {
		for (int i = 0; i < 54; i++)
			cout << genb12.grid0[i] + 1;
		cout << "entry m10 "  << endl;
	}
	g17xy.g17hh0.diagh = 0;
	g17xy.go_back = 0;
	zh2b_g.test = 0;
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	myband2.DoExpandBand(27);// expand band2
	if (!(myband1.n5| myband2.n5)) return; // no 656 no 566
	int nb3 = genb12.nband3;
	int ni3 = myband2.nind[2];
	if (diag ) {
		cout << "first check" << endl;
		//myband1.PrintStatus();
		myband1.PrintShortStatus();
		//myband2.PrintStatus(); 
		myband2.PrintShortStatus();
	}
	//if (diag) 		myband2.Debug_2_3();
	if (diag) {
		myband1.DebugIndex2();
		myband2.DebugIndex2();
		//return;
	}

	//=========================== collect UAs (still room for improvement)
	zh1b_g.modegua = 0;//must be to activate filter in UAs b12 more
	if(genuasb12.Initgen()) return;
	if (diag) {
		genuasb12.DebugUas();
		cout << diagpuz << " puz known" << endl;
		//return;
	}
	genb12.BuildGang9x3();
	zh1b_g.modegua = 1;//must be to kill  filter in GUAs 6_7 more
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 
	p_cpt2g[18] += genuasb12.nua;
	p_cpt2g[19] += genb12.ntua2;
	p_cpt2g[20] += genb12.ntua3;
	p_cpt2g[21] += genb12.nactive2;
	p_cpt2g[22] += genb12.nactive3;
	Go();// standard entry point for all 
}

void G17B::Go(){// search process external loop 2X2Y
	indexstep.StartDead();
	int n1 = myband1.nind[1],
		n2 = myband2.nind[1];
	if (diag == 2)	GodebugInit(1);//  1 base 2 all UAs 
	for (int i1 = 0; i1 < n1; i1++) {
		if (g17b.debug17) {
			int * t1= myband1.index2[i1];
			if ((t1[0] & band1_17) == t1[0]) {
				if (g17b.debug17>1)cout << Char27out(t1[0]) << " step 1  BF i1=" << i1
				<< "this is the right i1" << endl;
			}
			else continue;
		}
		if (diag == 2) {
				int * t1 = myband1.index2[i1];
				if ((t1[0] & p17diag.bf.u32[0]) == t1[0]) {
					cout << Char27out(t1[0]) << " step 1  BF i1=" << i1
						<< "this is the right i1" << endl;
				}
				else continue;

		}
		for (int i2 = 0; i2 < n2; i2++) {
			//if (i2) continue;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			indexstep.Init(i1, i2);
			g17morebig.Init();
			if (!(indexstep.n51 | indexstep.n52)) continue;
			if (indexstep.ShrinkUas())continue;// dead branch
			if (g17b.debug17) {
				int * t2 = myband2.index2[i2];
				if ((t2[0] & band2_17) == t2[0]) {
					indexstep.IndexStepDebugKnown17(i1, i2);
				}
				else continue;
			}
			if (diag == 2) {
				int * t2 = myband2.index2[i2];
				if ((t2[0] & p17diag.bf.u32[1]) == t2[0]) {
					cout << Char27out(t2[0]) << " step 1  BF i2=" << i2
						<< "this is a right i2" << endl;
					//			indexstep.IndexStepDebugKnown17(i1, i2);
				}
				else continue;
			}
			// select Y mode and XY order depending on the count
			indexstep.ShrinkGuas();
			uint64_t n56 = indexstep.n51*indexstep.n62,
				n65 = indexstep.n52*indexstep.n61;
			if (sgo.vx[4]) {// mode p2b only 656
				if (!indexstep.n52)continue;
				if (n65 < LIM3Y)  	indexstep.Do65();
				else indexstep.Do65_3y();
				if (g17xy.go_back) return;
			}
			else {// p2a both 566 and 656
				if (n65 < LIM3Y) {
					if (indexstep.n52) 	indexstep.Do65();
					if (g17xy.go_back) return;
					if (n56 < LIM3Y) {
						if (indexstep.n51) {
							indexstep.ShrinkUas();
							indexstep.Do56();
						}
						if (g17xy.go_back) return;
					}
					else {
						indexstep.Do56_3y();
						if (g17xy.go_back) return;
					}
				}
				else {
					if (n56 < LIM3Y) {// do 56 first
						if (n56) {
							indexstep.Do56();
							if (g17xy.go_back) return;
						}
					}
					else {// both in 3y mode
						indexstep.Do56_3y();
						if (g17xy.go_back) return;
					}
					indexstep.Do65_3y();
					if (g17xy.go_back) return;
				}
			}
		}
	}
}
void G17INDEXSTEP::Init(int i1, int i2) {
	g17chunk.i1 = i1; g17chunk.i2 = i2;
	//cur_step_buffer_index = 0;
	// index is bitfiel;current index 5, current index 6 
	t1 = myband1.index2[i1];
	t2 = myband2.index2[i2];
	n51 = t1[4] - t1[1];	n61 = t1[5] - t1[2];
	n52 = t2[4] - t2[1];	n62 = t2[5] - t2[2];
	//update the dead situation
	XY_EXPAND xe1 = myband1.xye6[t1[2]],
		xe2 = myband2.xye6[t2[2]];// here cell is 27-53
	int cell1 = xe1.cells.u8[0], cell2 = xe2.cells.u8[0] - 27;
	int cell1bf = 1 << cell1, cell2bf = 1 << cell2;
	if (!(cell1bf&oldx1)) {// new cell 1
		oldx1 = cell1bf;
		oldx2 = t1[0];
		xfirstdead |= cell1bf;
		x2dead = xfirstdead | t1[0];
		yfirstdead = cell2bf;
		oldy1 = cell2bf;
		y2dead = t2[0];
	}
	else if ((t1[0] & (~oldx2))) {// new xcell 2
		oldx2 = t1[0];
		x2dead = x2dead | t1[0];
		yfirstdead = cell2bf;
		oldy1 = cell2bf;
		y2dead = t2[0];
	}
	else if (!(cell2bf& oldy1)) {// new first y
		oldy1 = cell2bf;
		yfirstdead |= cell2bf;
		y2dead = yfirstdead | t2[0];
	}
	else y2dead |= t2[0];
	dead = ((uint64_t)y2dead << 32) | x2dead;
	start_dead = dead;
	b1b2_2Xbf = ((uint64_t)t2[0] << 32) | t1[0];
	start_b1b2bf = b1b2_2Xbf;
}

int G17INDEXSTEP::ShrinkUas() {// shrink table of UAs for bands1+2
	// could erase new supersets to consider doing it in 2 steps
	g17more.Init();
	uint64_t * otua = genuasb12.tua;
	int  ontua = genuasb12.nua;
	register uint64_t Rw = b1b2_2Xbf, Rn = ~dead;
	ntua = 0;
	for (int iua = 0; iua < ontua; iua++) {
		register uint64_t Ru = otua[iua];
		if (Rw & Ru) continue;
		Ru &= Rn;// erase dead cells
		if (Ru) {
			//if (TESTDEBUG) fout1 << "ua " << ntua << " " << Ru << endl;
			tua[ntua++] = Ru;
		}
		else return 1;// dead branch
	}
	ntua_raw = ntua;
	if (0) {
		cout << "shrinkua nua=" << ontua << " new count " << ntua << endl;
		cout << Char2Xout(Rw) << " b1b2_2xbf ox" << hex<< Rw <<dec<< endl;
		cout << Char2Xout(Rn) << " not dead ox" << hex << Rn << dec << endl;
	}
	if (ntua > 128* NVUAS128) ntua = 128 * NVUAS128; //limit in the chunk loop

	n128vua = (ntua + 127) >> 7;

								//cout << "rawuas=" << ntua_raw << " ntua=" << ntua <<" ontua="<<ontua<< endl;
	{ // initial vectors to appropriate value and actives vector
		memset(v256a.v, 0, sizeof v256a.v);
		if (ntua <= 128)	v256a.v[0] = maskLSB[ntua];
		else {
			v256a.v[0].SetAll_1();
			if (ntua <= 256)				v256a.v[1] = maskLSB[ntua - 128];
			else {
				v256a.v[1].SetAll_1();
				if (ntua <= 384)	v256a.v[2] = maskLSB[ntua - 256];
				else {
					v256a.v[2].SetAll_1();
					v256a.v[3] = maskLSB[ntua - 384];
				}
			}
		}
	}
	for (int i = 0; i < 54; i++) v256uas[i] = v256a;
	//if (TESTDEBUG) v256a.Fout(" v256a ");

	uint32_t cc;
	for (int i = 0; i < ntua; i++) {// set uas
		register int bloc = i >> 7, ir = i & 127;
		register uint64_t Rw = tua[i] & BIT_SET_2X;
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			Rw ^= bit2;// clear bit
			v256uas[From_128_To_81[cc]].v[bloc].clearBit(ir);
		}
	}
	/*
	if (TESTDEBUG) {
		fout1 << "table des vecteurs cells_ua" << endl;
		for (int i = 0; i < 54; i++) {
			fout1 << i << " ";
			v256uas[i].Fout("<-cell ");
		}
	}
	*/
	return 0;
}

int G17INDEXSTEP::ShrinkGuasLoad(uint64_t *tu,  int itu1) {
	if (ntgua > 128 * G17BLOCGSUA) return 0;
	int n = 0;
	ntgua_raw += itu1;
	register uint64_t Rw = b1b2_2Xbf, Rn = ~dead;
	for (int iua = 0; iua < itu1; iua++) {
		register uint64_t Ru = tu[iua];
		if (!(Rw & Ru)) {
			Ru &= Rn;
			if (Ru) {
				tgua[ntgua++] = Ru;
				n++;
				if (n >= 30)return n;
			}
		}
	}
	return n;
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
void G17INDEXSTEP::Do65() {
	mode3y = 0;
	g17chunk.i56 = 0;
	g17xy.b3lim = 6;
	{//=================== move 5 b2 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n52;
		uint64_t *Ro = (uint64_t *)(&myband2.xye5[t2[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] * nw2);
	}
	

	{//=================== move 6 b1 fix place w1
		g17xy.nby = 6;
		nw1 = n61;
		uint64_t *Ro = (uint64_t *)(&myband1.xye6[t1[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof  xyew1[0] * nw1);

	}
	Do_Common();
}
void G17INDEXSTEP::Do65_3y() {
	p_cpt2g[6] ++;
	mode3y = 1;
	g17chunk.i56 = 0;
	g17xy.b3lim = 6;
	{//=================== move 5 b2 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n52;
		uint64_t *Ro = (uint64_t *)(&myband2.xye5[t2[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] * nw2);
	}

	int i1 = g17chunk.i1, i3lim = myband1.nind[2],
		i = 0, nx = 0, i3d, i3f;
	int * ind2cur = myband1.index2[i1],
		curbf = ind2cur[0], lim1 = ind2cur[2],
		lim2 = myband1.index2[i1 + 1][2];
	for (; i < i3lim; i++) {
		int n = myband1.index3[i][1];
		if (n < lim1) continue;
		if ((myband1.index3[i][0] & curbf) != curbf)continue;//empty triplet
		break;
	}
	i3d = i;
	for (; i < i3lim; i++) {
		int n = myband1.index3[i][1];
		if (n >= lim2) break;
		nx++;
	}
	i3f = i;


	uint64_t y3_b1b2_2Xbf = start_b1b2bf;
	uint64_t y3_dead = start_dead;
	for (int i3 = i3d; i3 < i3f; i3++) {
		int * t3 = myband1.index3[i3], filter = t3[0];
		uint64_t f2x = (uint64_t)filter;
		y3_dead |= f2x;		dead = y3_dead;
		b1b2_2Xbf = y3_b1b2_2Xbf | f2x;
		if(ShrinkUas()) continue;
		ShrinkGuas();

		if (g17b.debug17) {
			if ((t3[0] & g17b.band1_17) == t3[0]) {
				if (g17b.debug17>1)cout << Char27out(t3[0]) << " step 3y valid for i3=" << i3 << endl;
				//break;
			}
			else continue;
		}
		{// ===================== move 6 b2 fix place w1
			g17xy.nby = 6;
			nw1 = t3[3] - t3[1];
			uint64_t *Ro = (uint64_t *)(&myband1.xye6[t3[1]]), *Rd = (uint64_t *)(xyew1);
			memcpy(Rd, Ro, sizeof  xyew1[0] * nw1);
		}
		Do_Common();
	}
}

void G17INDEXSTEP::Do56(){
	mode3y = 0;
	g17xy.b3lim = 6;
	g17chunk.i56 = 1;
	{// ===================== move 5 b1 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n51;
		uint64_t *Ro = (uint64_t *)(&myband1.xye5[t1[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] * nw2);
	}
		 
	{// ===================== move 6 b2 fix place w1
		g17xy.nby = 6;
		nw1 = n62;
		uint64_t *Ro = (uint64_t *)(&myband2.xye6[t2[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof  xyew1[0] *nw1);
	}
	Do_Common();
}
void G17INDEXSTEP::Do56_3y() {
	p_cpt2g[6] ++;
	mode3y = 1;
	g17xy.b3lim = 6;
	g17chunk.i56 = 1;
	{// ===================== move 5 b1 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n51;
		uint64_t *Ro = (uint64_t *)(&myband1.xye5[t1[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] * nw2);
	}
	int i2 = g17chunk.i2, i3lim = myband2.nind[2],
		i = 0,nx = 0,i3d,i3f;
	int * ind2cur = myband2.index2[i2],
		curbf = ind2cur[0],lim1= ind2cur[2],
		lim2= myband2.index2[i2 + 1][2];
	for (; i <i3lim; i++) {	
		int n = myband2.index3[i][1];
		if (n < lim1) continue;
		if ((myband2.index3[i][0] & curbf) != curbf)continue;//empty triplet
		break;
	}
	i3d = i;
	for (; i < i3lim; i++) {
		int n = myband2.index3[i][1];
		if (n >= lim2) break;
		nx++;
	}
	i3f = i;

	uint64_t y3_b1b2_2Xbf =start_b1b2bf;
	uint64_t y3_dead = start_dead;
	for (int i3 = i3d; i3 < i3f; i3++) {
		int * t3 = myband2.index3[i3],filter=t3[0];
		uint64_t f2x = (uint64_t)filter << 32;
		if (0) {
			cout << "y3=" << i3 << " filter 0" <<oct<< t3[0]<<dec<<endl;
			XY_EXPAND w = myband2.xye6[t3[1]];
			cout << "xye " << (int)w.cells.u8[0] << " " << (int)w.cells.u8[1]
				<< " " << (int)w.cells.u8[2] << endl;
		}

		y3_dead |= f2x;		dead = y3_dead;
		b1b2_2Xbf= y3_b1b2_2Xbf| f2x;
		if(ShrinkUas()) continue;
		ShrinkGuas();

		if (g17b.debug17) {
			if ((t3[0] & g17b.band2_17) == t3[0]) {
				if (g17b.debug17 > 1)cout <<Char27out(t3[0] )<< " step 3y valid for i3=" << i3 << endl;
				//break;
			}
			else continue;
		}
		{// ===================== move 6 b2 fix place w1
			g17xy.nby = 6;
			nw1 = t3[3]- t3[1];
			uint64_t *Ro = (uint64_t *)(&myband2.xye6[t3[1]]), *Rd = (uint64_t *)(xyew1);
			memcpy(Rd, Ro, sizeof  xyew1[0] * nw1);
		}
		//PrintUasShrinked(1);
		Do_Common();
	}
}

void G17INDEXSTEP::Do_Common(){// after xye w1('Y')  and w2('X') have been loaded
	if (!nw1)return;
	if(diag_on)
	cout <<"common i56="<<g17chunk.i56
		<< " n51=" << n51 << " n61=" << n61
		<< " n52=" << n52 << " n62=" << n62
		<< " ntua=" << ntua << " ntgua=" << ntgua
		<< "\tnua" << genuasb12.nua 
		<< "\tnua_raw" << ntua_raw << endl
		<<endl;
	if (g17b.debug17) if (!DebugIsKnown())return;
	Do_Common_2();// build  Y tables of vectors
	Do_Common_3();// launch chunks
}
void G17INDEXSTEP::Do_Common_2(){// build  Y tables of vectors
	for (int iy = 0; iy < nw1; iy++){ //Y is bloc xyew1
		uint8_t * ty = xyew1[iy].cells.u8;
		V256_UAS & d256 = yv256uas[iy];
		d256=v256a;
		V256_GUAS & dg256 = yv256guas[iy];
		dg256=v256ga;
		for (int i = 2; i < 6; i++){
			int cell = ty[i];
			d256 &= v256uas[cell];
			dg256 &= v256guas[cell];
		}
	}
	/*
	if (TESTDEBUG>1) {
		fout1 << "table des vecteurs yua" << endl;
		for (int i = 0; i < nw1; i++) {
			fout1 << i << " ";
			yv256uas[i].Fout("<-yindex ");
		}
	}
	*/
}
void G17INDEXSTEP::Do_Common_3_BuildXvectors(){// in fact, in the chunk for 256 X maximum
	// re using the code for 'y', iy is in fact ix
	for (int ix = 0; ix <g17chunk.nxc; ix++){ //X is in the chunk
		uint8_t * ty = g17chunk.txec[ix].cells.u8;
		V256_UAS & d256 = g17chunk.xv[ix];
		d256 = v256a;
		V256_GUAS & dg256 = g17chunk.xvg[ix];
		dg256=v256ga;
		for (int i = 2; i < 5; i++){
			int cell = ty[i];
			d256 &= v256uas[cell];
			dg256 &= v256guas[cell];
		}
	}
	/*
	if (TESTDEBUG) {
		fout1 << "check cell vector status building x vectors" << endl;
		for (int i = 0; i < 54; i++) {
			fout1 << i << " ";
			v256uas[i].Fout("<-cell ");
		}
	}
	if (TESTDEBUG>1) {
		fout1 << "table des vecteurs xua" << endl;
		for (int i = 0; i < g17chunk.nxc; i++) {
			fout1 << i << " ";
			g17chunk.xv[i].Fout("<-xindex ");
		}
	}
	*/
}
void G17INDEXSTEP::Do_Common_3(){// launch chunks 256 x256 
	// now combining all 'x' to all 'y' in chunks size ..
	//if (1) return;
	int nxch = (nw2-1) / G17CHKX , nxc;
	int nych = (nw1-1) / G17CHKY, nyc;
	//cout << "Do_Common_3() nxch=" << nxch << "  nych=" << nych 
	//	<<" nw2="<<nw2<<" nw1="<<nw1<< endl;
	for (int ix = 0, iex = 0; ix <= nxch; ix++, iex += G17CHKX){
		nxc = G17CHKX;
		if (ix == nxch) nxc = nw2 - G17CHKX * ix;
		g17chunk.nxc = nxc;
		memcpy(g17chunk.txec, &xyew2[iex], sizeof g17chunk.txec);
		Do_Common_3_BuildXvectors();// build 'X' tables in the chunk
		for (int iy = 0, iey = 0; iy <= nych; iy++, iey += G17CHKY){
			nyc = G17CHKY;
			if (iy == nych) nyc = nw1 - G17CHKY * iy;
			g17chunk.nyc = nyc;
			memcpy(g17chunk.tyec, &xyew1[iey], sizeof g17chunk.tyec);
			memcpy(g17chunk.yv, &yv256uas[iey], sizeof g17chunk.yv);
			g17chunk.yvg = &yv256guas[iey];
			if (g17b.debug17) {
				if (g17chunk.DebugIsKnown()) {
					if (g17b.debug17 > 1)cout << "valid chunk ix=" << ix << " iy=" << iy
						<< "\tnxc=" << nxc << " nyc=" << nyc << endl;
				}
				else continue;
			}
			g17chunk.GoChunk();
			if (g17xy.go_back) return;
		}
	}
}
void G17INDEXSTEP::PrintUasShrinked(int opt) {
	cout << "uas shrinked table ntua="<<ntua << endl;
	for (int i = 0; i < ntua; i++)
		cout << Char2Xout(tua[i])<<" "<<i << endl;
	if (!opt) return;
	cout << "table des vecteurs  (base)" << endl;
	v256a.Cout();
	for (int i = 0; i < 54; i++) {
		cout << "cell " << i << " ";
		v256uas[i].Cout();
	}
}

void G17CHUNK::GoChunk() {// elementary 'X' x 'Y' ychunk is 256x256
	// this is the critical code in "part 1"
	//if (1) return;
	p_cpt2g[1] ++;
	p_cpt2g[2] += nxc * nyc;
	//if (g17b.debug17) {
		//cout << "go chunk nxc=" <<nxc<<" nyc="<<nyc << endl;
		//if (!DebugIsKnown())return;
	//}
	uint16_t * limit_store;
	uint16_t tstore[G17CHKX * G17CHKY]; // storing XY passing filter 1 will never reach the limit
	{//=========Find all XY passing first filter for the chunk
		register uint16_t * pstore = tstore;
		register int iy, store, ix;
		register BF128 * tvx = xv[0].v;
		for (ix = 0; ix < nxc; ix++, tvx += NVUAS128) {
			store = ix << 8; //8 low bits for iy
			register uint64_t Rx = txec[ix].stacks.u64;
			BF128 vx = tvx[0], vx2 = tvx[1], vx3 = tvx[2], vx4 = tvx[3];
			register BF128 * tvy = yv[0].v;
			for (iy = 0; iy < nyc; iy++,tvy+= NVUAS128) {
				BF128 vy = tvy[0]; 
				if ((vx & vy).isNotEmpty())continue;
				if (indexstep.n128vua > 1) {
					vy = tvy[1];
					if ((vx2 & vy).isNotEmpty()) continue;
					if (indexstep.n128vua > 2) {
						vy = tvy[2];
						if ((vx3 & vy).isNotEmpty()) continue;
						if (indexstep.n128vua > 3) {
							vy = tvy[3];
							if ((vx4 & vy).isNotEmpty()) continue;
						}
					}
				}
				register uint64_t Ry = Rx + tyec[iy].stacks.u64;
				if ((Ry & 0xff) > 6) continue;
				if ((Ry & 0xff0000) > 0x60000) continue;
				if ((Ry & 0xff00000000) > 0x600000000) continue;
				*pstore++ = store | iy;// (ix<<8) | iy
			}
		}

		limit_store = pstore;
	}
	{//=========== send such XY to next step in the process
		uint16_t * wstore = tstore;
		while (wstore < limit_store) {
			register int w = *wstore++;// catch the value
			{
				register int ix = w >> 8;//8 bits  x
				register int iy = w & 0xff;// 8 bits y
				g17xy.xxe = txec[ix];
				g17xy.yye = tyec[iy];
				g17xy.xvg = xvg[ix];
				g17xy.yvg = yvg[iy];
				g17xy.ix = ix;
				g17xy.iy = iy;
			}
			g17xy.Go_0();
			if (g17xy.go_back) return;
		}
	}
}

void G17XY::Go_0(){// check stacks, check more uas, collect Guas status
	indexstep.diag_on = 0;
	p_cpt2g[3] ++;
	Init();
	if (g17b.debug17){
		if ((cellsbf&g17b.band1_17) != g17b.band1_17) return;
		if (((cellsbf >> 32) &g17b.band2_17) != g17b.band2_17) return;
		cout << Char2Xout(cellsbf) << " expected XY" << endl;
	}

	if (indexstep.diag_on) {
		indexstep.diag_on = 1;
		if ((cellsbf&g17b.band1_17) == g17b.band1_17) 
			if (((cellsbf >> 32) &g17b.band2_17) == g17b.band2_17) {
				cout << Char2Xout(cellsbf) << " expected XY" << endl;
				indexstep.diag_on = 2;
		}

	}
	if ( g17b.diag == 2) {
		//if (cellsbf != g17b.p17diag.bf.u64[0]) return;
		if (cellsbf == g17b.p17diag.bf.u64[0]) {
			cout << Char2Xout(cellsbf) << " expected XY" << endl;
			indexstep.diag_on = 2;
		}
	}

	stacksxy.u64 = xxe.stacks.u64 + yye.stacks.u64;
	if (g17b.debug17) cout << "stacks " << stacksxy.u16[0] << stacksxy.u16[1] << stacksxy.u16[2] << endl;
	if (g17more.Check(cellsbf))return;		// check here more uas
	if (g17morebig.Check(cellsbf))return;
	Go_Guas_collect();	// collect guas still active
	if (g17b.debug17 > 1) cout << "call FirstCheckActiveBands()" << endl;
	if (indexstep.diag_on > 1)cout << "call FirstCheckActiveBands()" << endl;
	if (FirstCheckActiveBands()) return;
	p_cpt2g[4]++;
	if (g17b.debug17 > 1) cout << "call CheckValidBand12()" << endl;
	if (indexstep.diag_on > 1)cout << "call CheckValidBand12()" << endl;
	if (CheckValidBand12()) return;
	p_cpt2g[5]++;
	//if (TESTDEBUG)fout1 << "XY" << cellsbf << endl;
	//if (DEBUGLEVEL == 1) return;
	if (g17b.debug17) {
		cout << "entry Go_0 active nclues=" << nclues << " n=" << p_cptg[8]
			<< " stacks " << stacksxy.u16[0] << stacksxy.u16[1] << stacksxy.u16[2]
			<< "\t" << stacksxy.u16[3] << endl;
	}
	//genb12.Sockets2SetupForB12(cellsbf);// seems very expensive
	if (indexstep.diag_on > 1)cout << "call BuildActiveBands()" << endl;
	BuildActiveBands();// apply pairs triplets to bands 
	if (g17b.debug17)cout << "exit buidactivebands nb3=" << ntb3 << endl;
	if (indexstep.diag_on > 1)cout << "exit buidactivebands nb3=" << ntb3 << endl;
	if (!ntb3) return; // no active band b3lim and stack filters
	// at least one band is ok, time to check validity of the band 1+2
	p_cpt2g[8] ++;
	p_cpt2g[9] += ntb3;
	if (indexstep.diag_on > 1)cout << "go valide" << endl;
	Go_ValideB12();
}

void G17XY::Init() {
	cellsbf= xxe.cellsbf | yye.cellsbf;
	int n = 0;
	{
		register uint64_t R = xxe.cells.u64;
		for (int icell = 0; icell < nbx; icell++, R >>= 8) {
			tclues[n++] = (uint32_t)(R & 0xff);
		}
		R = yye.cells.u64;
		for (int icell = 0; icell < nby; icell++, R >>= 8) {
			tclues[n++] = (uint32_t)(R & 0xff);
		}
	}
	nclues = n;
	// setup gangster of clues XYgang
	int * b = genb12.grid0;
	memset(xygang, 0, sizeof xygang);
	for (int iclue = 0; iclue < n; iclue++) {
		int cell = tclues[iclue], col = cellsFixedData[cell].pl ;
		xygang[col] |= 1 << b[cell];
	}

}
void G17XY::Go_Guas_collect(){// 
	bands_active_pairs.SetAll_0();
	more_active_pairs.SetAll_0();
	more_active_triplets.SetAll_0();
	bands_active_triplets.SetAll_0();
	//more3 = 0;
	xyvg = xvg; xyvg &= yvg;
	// setup XY guas status
	for (int i = 0; i < indexstep.nactive2; i++) {
		int i81= indexstep.tactive2[i],
			b1 = indexstep.tactive2_start[i] >> 7,
			b2 = indexstep.tactive2_end[i] >> 7;

		if (b1 > NVGUAS128) continue; // lost uas 
		if (b2 > NVGUAS128)b2 = b1; // partially lost uas 
		BF128 vw = xyvg.v[b1] & indexstep.vid81s2[i81].v[b1];
		if (vw.isEmpty()) {// check second 
			if (b1 == b2) continue;
			vw = xyvg.v[b2] & indexstep.vid81s2[i81].v[b2];
			if (vw.isEmpty()) continue;// all hit
		}
		bands_active_pairs.Set_c(i81);
	}
	for (int i = 0; i < indexstep.nactive3; i++) {
		int i81 = indexstep.tactive3[i],
			b1 = indexstep.tactive3_start[i] >> 7,
			b2 = indexstep.tactive3_end[i] >> 7;
		if (b1 > NVGUAS128) continue; // lost uas 
		if (b2 > NVGUAS128)b2=b1; // partially lost uas 
		BF128 vw = xyvg.v[b1] & indexstep.vid81s3[i81].v[b1];
		if (vw.isEmpty() ) {// check second 
			if (b1 == b2) continue;
			vw = xyvg.v[b2] & indexstep.vid81s3[i81].v[b2];
			if (vw.isEmpty()) continue;// all hit
		}
		bands_active_triplets.Set_c(i81);
	}
	//more_checked_triplet= bands_active_triplets;
}
int G17XY::FirstCheckActiveBands() {// quick check of critical situations
	//if (genb12.nband3 > 10) return 0;// small chances to exit false<<<<<<<<<<<<<<<<<<<<<<
	// consider critical stack 
	ibfirst = 0;
	int maxcount = 0, stack;
	for (int i = 0; i < 3; i++) {
		int st_count = stacksxy.u16[i];
		if (st_count > maxcount) { maxcount = st_count; stack = i; }
	}
	if (maxcount < 5)return 0;// small chances to exit false
	for (ibfirst = 0; ibfirst < genb12.nband3; ibfirst++) {
		STD_B3::GUAs & myb = genb12.bands3[ibfirst].guas;
		BF128 ws2 = bands_active_pairs & myb.isguasocket2;
		BF128 ws3 = bands_active_triplets & myb.isguasocket3;
		uint32_t ws2s= ws2.bf.u32[stack],ws3s= ws3.bf.u32[stack];
		if (g17b.debug17) 
			cout << Char27out(ws2s | ws3s) << "FirstCheckActiveBands() ws2s | ws3s= " 
			<< "stack="<<stack<<" maxcount="<<maxcount<< endl;
		if(! (ws2s | ws3s)) return 0;// one non critical band found
		if (maxcount == 6)continue; // band to exclude
		// need 2 minimum clues to exclude the band
		int * tm2 = &myb.ua2_imini[27 * stack],* tm3=&myb.ua3_imini[27 * stack];
		int tix[81], ntix = 0, mini_bf1=0,mini_bf2=0, mini_bf3 = 0;
		BitsInTable32(tix, ntix, ws2s);
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int i81_rel = tix[i] ;
			int imini = tm2[i81_rel],bit=1<<imini; 
			//cout << "bit N°" << tix[i] << "\ti81=" << i81_rel << "\timini\t" << imini << endl;
			if (mini_bf2&bit) mini_bf3 |= bit;
			if (mini_bf1&bit) mini_bf2 |= bit;
			mini_bf1 |= bit;
		}
		int ntix3 = 0,mini_triplet=0;
		BitsInTable32(tix, ntix3, ws3s);
		for (int i = 0; i < ntix3; i++) {// switch from i_81 to mini rows
			int i81_rel = tix[i] ;
			int imini = tm3[i81_rel], bit = 1 << imini;
			//cout << "triplet bit N°" << tix[i]<< "\ti81="<<i81_rel << "\timini\t" << imini << endl;
			mini_triplet |= bit;
		}
		mini_triplet &= ~mini_bf1; // no double count
		uint32_t mincount = _popcnt32(mini_bf1) + _popcnt32(mini_bf3) + _popcnt32(mini_triplet);
		if (g17b.debug17)
			cout << "FirstCheckActiveBands() mincount= "<<mincount 
			<<" min 1 3 t "<<oct<< mini_bf1 <<" "<< mini_bf3 << " "<< mini_triplet << " "<<dec<< endl;
		if (mincount < 2) return 0; // non critical band found
	}
	return 1;// all bands are in critical mode
}
int G17XY::CheckValidBand12(){
	// passing more filter now check validity or add a "more" ua
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues);
	//if (indexstep.diag_on > 1)
	if (myua){
		uint64_t cc = _popcnt64(myua);
		if (cc < 12) {
			cout<<Char2Xout(cellsbf) << " invalid ua CheckValidBand12() "
				<<cellsbf<< endl; 
			cout << Char2Xout(myua) << "  ua produced " << endl;
			cout << myband1.band << myband2.band << " bands 1+2 studied num "
				<<genb12.nb12<< endl;
			genuasb12.DebugUas();
			cout << endl << endl;
			cout << Char2Xout(indexstep.b1b2_2Xbf) << " b1b2step" << endl;
			cout << Char2Xout(indexstep.dead) << " dead" << endl << endl;
			cout << Char2Xout(cellsbf) << " expected XY" << endl;
			indexstep.PrintUasShrinked(0);
			BF128 * xv = g17chunk.xv[ix].v,
				*yv = g17chunk.yv[iy].v;
			cout << Char64out(xv[0].bf.u64[0]) << " debut xv" << endl;
			cout << Char64out(yv[0].bf.u64[0]) << " debut yv" << endl;
			zh2b[0].ValidXY(tclues, nclues,1);
			go_back = 1;
			return 1;
		}
		//if (TESTDEBUG) fout1 << "0]addua " << myua << " cc="<<cc
		//	<<"\t"<< genuasb12.nua << "\t" << g17more.nt << "\t" << g17morebig.nt << endl;
		if (cc <= UALIMSIZE) {
			g17more.Add(myua);
			genuasb12.AddUACheck(myua | (cc << 59));// and update the UA table
		}
		else g17morebig.Add(myua);
		return 1; // not a solution unique for bands 1+2
	}
	return 0;
}

void G17XY::BuildActiveBands() {
	ntb3 = 0;
	//if (ibfirst)ibfirst--;// to check later why not done
	if (indexstep.diag_on > 1)cout << "BuildActiveBands()  start ibfirst=" <<ibfirst<< endl;
	for (int ib3 = ibfirst; ib3 < genb12.nband3; ib3++) {
		if (indexstep.diag_on > 1)cout << "ib3=" << ib3 << endl;
		STD_B3::GUAs & myb = genb12.bands3[ib3].guas;
		BF128 ws2 = bands_active_pairs & myb.isguasocket2;
		BF128 ws3 = bands_active_triplets & myb.isguasocket3;
		// switch to mini rows patterns
		int tix[81], ntix = ws2.Table3X27(tix);
		uint32_t mini_bf1 = 0, mini_bf2 = 0, mini_bf3 = 0,
			pairsbf = 0,
			pairs27 = 0,
			mini_triplet = 0;
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int i81= tix[i],
				imini = myb.ua2_imini[i81],
				bit = 1 << imini;
			//cout << "i81=" << i81 << " imini=" << imini << " i=" << i << endl;
			if (mini_bf2&bit) mini_bf3 |= bit;
			if (mini_bf1&bit) mini_bf2 |= bit;
			mini_bf1 |= bit;
			pairsbf |= myb.ua_pair[i81];
			pairs27 |= 1<< myb.ua2_i27[i81];
		}
		ntix = ws3.Table3X27(tix);// now triplets to mini rows
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int imini = myb.ua3_imini[tix[i]], bit = 1 << imini;
			//cout << "triplet bit N°" << tix[i] << "\timini\t" << imini << endl;
			mini_triplet |= bit;
		}

		//___________________________ prepare a new band to process later
		uint32_t all_used_minis = mini_bf1 | mini_triplet;
		memset(&wg3, 0, sizeof wg3);
		wg3.ib3 = ib3;
		wg3.count.u16[3]= mini_bf1;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		wg3.minirows_triplets = mini_triplet;
		uint32_t mincount = _popcnt32(mini_bf1) + _popcnt32(mini_bf3) + _popcnt32(mini_triplet);
		if (indexstep.diag_on > 1)cout << "mincount=" << mincount << endl;

		if (mincount > 6) continue; // too many clues needed for 6 5 6

		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs
		
		wg3.count.u16[0] = mini_bf1;	wg3.count.u16[1] = mini_bf2;
		wg3.count.u16[2] = mini_bf3;
		
		wg3.countstack.u64 = stacksxy.u64;
		wg3.pairs.u32[0] = pairs27;
		// min count per stack
		wg3.countsum.u16[3] = mincount;
		{	
			register int s = 0111;
			for(int i=0;i<3;i++,s<<=1){
				wg3.countsum.u16[i] = _popcnt32(mini_bf1&s) 
					+ _popcnt32(mini_bf2 & s)
					+ 2 * _popcnt32(mini_bf3&s) 
					+ _popcnt32(mini_triplet&s);
			}
		}

		{ // setup sum and check stacks
			register uint64_t R = stacksxy.u64 + wg3.countsum.u64;
			if (indexstep.diag_on > 1)cout << "R=ox" << hex<<R<<dec << endl;
			wg3.countstack.u64 = R;
			if ((R & 0xff) > 6) continue;
			R >>= 16; if ((R & 0xf) > 6) continue;
			R >>= 16; if ((R & 0xf) > 6) continue;
		}
		// set up pair + triplet bitfield
		wg3.pairs.u32[1] = pairsbf;
		if (wg3.minirows_triplets) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (wg3.minirows_triplets&bit)
					wg3.pairs.u32[1] |= field;
		}
		
		if (0) {// have doubts that it pays
		// try to add triplets in empty minirows
			all_used_minis ^= 0777;// now free minis
			if (all_used_minis && (p_cpt2g[8] < 30)) {
				Find_More_Sockets3(all_used_minis);
			}
		}
		g17tb3go[ntb3++] = wg3;
	}
}
int G17XY::Find_More_Sockets3(uint32_t free_minis) {
	int rmini = 0;
	//cout << "test more sockets 3" << endl;
	STD_B3::GUAs & myb = genb12.bands3[wg3.ib3].guas;
	for (int imini = 0,bit=1; imini < 9; imini++,bit<<=1) {
		if (!(free_minis & bit))continue;
		int i81 = myb.triplet[imini];
		// is there a socket3 from a previous check
		if (more_checked_triplet.On_c(i81)) {// already cheked
			if (more_active_triplets.Off_c(i81)) {// already cheked
				more_active_triplets.Set_c(i81);
			}
			continue;
		}

		GEN_BANDES_12::SGUA3 &w = genb12.tsgua3[i81];
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = bita | bitb | bitc;
		int triplet_perms[2][3];
		if (g17xy.xygang[w.col1] & digs)continue;// not valid here
		if (g17xy.xygang[w.col1+1] & digs)continue;// not valid here
		if (g17xy.xygang[w.col1 + 2] & digs)continue;// not valid here

		int * p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, genb12.gangb12, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			// this is a possible  socket3 try to find a gua
			//cout << "call MoreSocket3 i81=" << i81 << endl;
			zh2b_g.test = 0;
			uint64_t new_ua = zh2b_g.MoreSocket2(genb12.gangb12, rgangcols,
				g17xy.tclues, g17xy.nclues, digs);
			zh2b_g.test = 0;
			//cout << Char2Xout(new_ua) << " retour zhou" << endl;
		}
	}
	return rmini;
}

void G17XY::Go_ValideB12(){// UA2 and UA3 known not yet dead with min clues in band 3
	//cout <<Char2Xout(g17xy.cellsbf)<< "entry Go_ValideB12() ntb3=" << ntb3 << endl;
	/*
	if (g17xy.cellsbf == TESTXY) {
		cout << Char2Xout(g17xy.cellsbf) << " enter XY test mode" << endl;
		g17hh0.diagh = 1;
	}
	if (g17xy.cellsbf == TESTXY2) {
		cout << Char2Xout(g17xy.cellsbf) << " leave XY test mode" << endl;
		g17hh0.diagh = 0;
	}	
	*/
	int zinitdone = 0;
	for (int i3 = 0; i3 < ntb3; i3++){
		wg3 = g17tb3go[i3];
		int nmiss= 17 - wg3.countstack.u16[3];
		if(nmiss<2)p_cpt2g[10+nmiss] ++;
		else if (nmiss < 5)p_cpt2g[12] ++;
		else p_cpt2g[13] ++;
		if(wg3.minirows_triplets)p_cpt2g[14] ++;
		if (g17b.debug17>1) 	wg3.Debug();
		if (indexstep.diag_on > 1)wg3.Debug();

		//============= collect Gua46 and uas b3 for the band split them "in-field" "out-field"
		nuasb3_1 = nuasb3_2 = 0;
		register int  Rfilt = g17tb3go[i3].pairs.u32[1];
		// first GUA46 usually shorter than UAs band3
		BF128  socket4 = genb12.bands3[wg3.ib3].guas.isguasocket4;// i81 3X
		socket4 &= bands_active_pairs;
		int * ua_46 = genb12.bands3[wg3.ib3].guas.ua_pair; // ua pattern
		if (0) {
			char ws[82];
			cout << socket4.String3X(ws) << "socket4" << endl;
		}
		int i81;
		while ((i81 = socket4.getFirstCell()) >= 0) {
			socket4.Clear_c(i81);// clear bit
			register int Ru = ua_46[i81];
			if (Ru & Rfilt)	uasb3_1[nuasb3_1++] = Ru;
			else		uasb3_2[nuasb3_2++] = Ru;
		}
		if ((!nmiss) && nuasb3_2) continue; // critical + outfield uas
		uint32_t * to = genb12.bands3[wg3.ib3].tua;
		for (uint32_t i = 0; i < genb12.bands3[wg3.ib3].nua; i++) {
			register int Ru = to[i];
			if (Ru & Rfilt)			uasb3_1[nuasb3_1++] = Ru;
			else uasb3_2[nuasb3_2++] = Ru;
		}
		if ((!nmiss) && nuasb3_2) continue; // critical + outfield uas
		p_cpt2g[15] ++;
		if (!zinitdone) {
			zinitdone = 1;
			if (zhou[0].PartialInitSearch17(tclues, nclues))
				return;// would be  bug 
		}
		memcpy(&genb12.grid0[54], genb12.bands3[wg3.ib3].band0, 4*27);
		//if (DEBUGLEVEL == 2) continue;
		//if (++p_cpt2g[19] > 10) return;
		//wg3.Debug();
		//for (int i = 54; i < 81; i++)
			//cout << zh_g2.grid0[i]+1;
		//cout << " band3 en zh_h2" << endl;


		g17hh0.Go();
//<<<<<<<<<<<<<<<<< a revoir
//		zhou[0].InitBand3PerDigit(genb12.bands3[wg3.ib3].band0);

		//if (g17b.debug17 > 1)		cout << " valide b12 go"  << endl;
	}
}

void G17XY::FoutValid17(int bf3, int ib3){
	char zs[82];
	int *g = genb12.grid0;
	strcpy(zs, empty_puzzle);
	for (int i = 0; i < nclues; i++) {
		int cell = tclues[i];
		zs[cell] = g[cell] + '1';
	}
	int bit = 1;
	g = genb12.bands3[ib3].band0;
	for (int i = 0; i < 27; i++, bit <<= 1)if (bf3&bit)
		zs[i + 54] = g[i] + '1';
	if (g17b.debug17) {
		fout1 << zs << ";" << g17b.npuz << endl;
		cout << zs << " valid found" << endl;
	}
	else fout1 << zs << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
}

//================ part 2  band procvessing
void ZHOU::InitBand3PerDigit(int * grid0b3) {
}
/*
void ZHOU::InitBand3PerDigit(int * grid0b3){
	memset(glb.band3digits, 0, sizeof glb.band3digits);
	register int * t = glb.band3digits;
	for (int i = 0; i < 27; i++){
		t[grid0b3[i]] |= 1 << i;
	}
	glb.digsols = zhxy27[0].glb.digsols;// catch pointer to solution per digit
}*/
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
	if (_popcnt32(digitsbf < 8)) return 1;// can not be one solution
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (diag) ImageCandidats();
	//__________end assign last lot start solver
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	int ir = Full17Update();
	if (diag) {
		cout << "after update" << endl;
		ImageCandidats();
	}
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0,diag);

	return zh_g.nsol;  
}

int ZHOU::Apply17SingleOrEmptyCells() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched 
	BF128 R1 = FD[0][0], R2 = R1 & FD[1][0]; 	R1 |= FD[1][0];
	BF128 R3 = R2 & FD[2][0]; R2 |= R1 & FD[2][0]; R1 |= FD[2][0];
	BF128 R4= R3 & FD[3][0]; 
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
	if (diag) {
		char ws[82];
		cout << zh_g.pairs.String3X(ws)<< " pairs "<< endl;
		cout << zh_g2.triplets.String3X(ws) << " triplets " << endl;
		cout << zh_g2.quads.String3X(ws) << " quads " << endl;
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
	if (diag) {
		cout << "guess17 index=" << index << " cell " << cellsFixedData[cell].pt << endl;
	}
	// if first step try first false
	if(!index)	ClearCandidate_c(digit, cell);// force false
	for (int idig = 0; idig < 9; idig++)
		if (FD[idig][0].On(xcell))tdig[ndig++] = idig;
	for (int idig = 0; idig < ndig; idig++) {
		ZHOU * mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(tdig[idig], cell, xcell);
		if (diag)cout <<"guess index="<<index<<" "
			<< tdig[idig]+1<< cellsFixedData[cell].pt << endl;
		mynext->Compute17Next(index+1,diag);
		if (zh_g.go_back) return;

	}
	if (!index) {
		FD[digit]->Set_c(cell);// restore the candidate
		SetaCom(digit, cell, xcell);
		if (diag)cout << "guess last index=" << index << " "
			<< digit + 1 << cellsFixedData[cell].pt << endl;
		Compute17Next(index, diag);

	}
}

void ZHOU::Compute17Next(int index, int diag) {
	int ir = Full17Update();
	if (!ir) return;// locked 
	if (diag>1) { cout << "index=" << index << endl; ImageCandidats(); }
	if (ir == 2) {//solved
		if (index)zh_g.nsol++;
		zh_g.go_back = 1;// closed anyway
		return;
	}
	Guess17(index , diag);// continue the process
}

int G17B3HANDLER::IsMultiple(int bf){
	if (diagh) 
		cout<<Char27out(bf)  << "call is multple diag="<<diagh << endl;
	if (_popcnt32(bf) > 25) return 0;
	// check first if all tuab3 is hit
	for (int i = 0; i < g17xy.ntuab3; i++)
		if (!(bf&g17xy.tuab3[i])) return 1;// previous UA not hit
	int ir = zhou[1].CallMultipleB3(zhou[0], bf , diagh);
	if (ir) {//consider store the fresh ua b3
	}
	return ir;
}

void G17B3HANDLER::Critical2pairs() {// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (wg3.count.u16[1]) {// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++) {
			int shrink = TblShrinkMask[wg3.count.u16[1] & (0111 << ist)];
			if (shrink) {// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~wg3.pairs.u32[0]);// and set the common cell as assigned
			}
		}
		wg3.count.u16[1] = 0;
	}
}
void  G17B3HANDLER::Not_Critical_wactive(){// find active cells
	// if the band is full all cells outfield are dead
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	register int Rst = 07007007,// stack 0 pattern
		Ractive_out = BIT_SET_27,//  active cells
		Rpairs = wg3.pairs.u32[1];//  pairs pattern
	//Rdead &= Rn;
	for (int ist = 0; ist < 3; ist++){
		register int Rw = Rst << (3 * ist); // current stack pattern
		if (wg3.countstack.u16[ist] < 6){
			Ractive_out ^= (Rw & Rpairs);// outfield if not pair
			continue;
		}
		Ractive_out ^= Rw;// no cell outfield if stack full
		int shrink = TblShrinkMask[wg3.count.u16[1] & (0111 << ist)];
		if (shrink){// minirows 2 pairs in that stack
			register int Mask = tbitsrows[shrink] << (3 * ist);
			//adjust count and known
			known_b3 |= Mask & (~wg3.pairs.u32[0]);// and set the common cell as assigned
			wg3.count.u16[1] &= (~(0111 << ist)); // clear the 2pairs bit(s) in stack
			wg3.pairs.u32[1] &= (~Mask);// clear the pairs in the in field bf
			wg3.pairs.u32[0] &= (~Mask);// clear the pairs in the pair bf
			ncritical -= _popcnt32(shrink);
		}
	}
	wactive0 = Ractive_out&ndead& BIT_SET_27;// be sure to have no extra bit
}
void  G17B3HANDLER::Go(){
	wg3 = g17xy.wg3;
	nmiss = 17 - wg3.countstack.u16[3];
	if (diagh && nmiss==1 ){
		cout << "start a new band3 ==========================" << endl;
		wg3.Debug();
	}
	known_b3 = 0;
	ncritical = wg3.countsum.u16[3];// to update assigning critical 2 pairs
	ndead = BIT_SET_27;
	if (g17b.debug17)	cout << " call handler go   nmiss="<<nmiss << endl;
	if (indexstep.diag_on > 1)	cout << " call handler go   nmiss=" << nmiss << endl;
	rknown_b3 = 0;
	g17xy.ntuab3 = 0;
	//if (DEBUGLEVEL == 3) return;
	if (nmiss){
		//if (DEBUGLEVEL == 3 && nmiss > 0) return;
		Not_Critical_wactive();// if the band is full all cells outfield are dead
		if (nmiss == g17xy.b3lim){// no ua2_3 expand the uas table within the limits
			SPOT17_NOUAS2_3 s;
			memset(&s, 0, sizeof s);
			s.active = wactive0;
			s.stacks = wg3.countstack;
			s.stacks.u16[3] = nmiss;
			s.newspot(g17xy.uasb3_2, g17xy.nuasb3_2, 0, 0);
		}
		else switch (nmiss){// add in priority cells outside the critical field
		case 1: Go_Not_Critical_miss1(); return;
		case 2: Go_Not_Critical_miss2(); return;
		case 3: Go_Not_Critical_miss3(); return;
		case 4: Go_Not_Critical_miss4(); return;
		case 5: Go_Not_Critical_miss5(); return;
		}
	}
	else {
		//if(!diagh )
		Go_Critical();
	}
}
//================= critical process
void G17B3HANDLER::CriticalAssignCell(int Ru){// assign a cell within the critical cells
	// Ru is usually a regidster containing a 27 bits field with one bit on
	// 2 pairs in a miniriow have already been applied
	known_b3 |= Ru;
	uint32_t cell;
	bitscanforward(cell, Ru); // catch the cell
	if (g17b.debug17>1) {
		cout << Char27out(known_b3) << " new known status after  CriticalAssignCell  "<<dec << cell << endl;
	}
	register int mini = C_minirow[cell],// minirow to clear
		bit = 1 << mini,
		Mask = 7 << (3 * mini);
	if (bit & wg3.count.u16[2]){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		wg3.count.u16[2] ^= bit; // now only a pairto hit
		wg3.count.u16[0] |= bit;
	}
	else{// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_b3 &= (~Mask); // kill the minirow as active
		wg3.count.u16[0] &= ~bit;
		wg3.minirows_triplets &= ~bit;
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
void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	p_cpt2g[16]++;
	active_b3 = wg3.pairs.u32[1];
	Critical2pairs();// assign 2 pairs in minirow to common cell
	rknown_b3 = known_b3| active_b3;
	if (!active_b3){
		if (diagh)cout << "final check1" << endl;
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (diagh) {
		cout << "critical status after assign 2 pairs" << endl;
		cout << Char27out(active_b3) << " active b3"<< endl;
		cout << Char27out(known_b3) << " known b3" << endl;
		wg3.Debug();
	}
	if(IsMultiple(rknown_b3))	return;

	if (g17b.debug17) {
		if ((rknown_b3 & g17b.band3_17) != g17b.band3_17) return;
		else 	cout << Char27out(rknown_b3) << " known b3 start critical" << endl;

	}
	if (diagh)cout << "back from is multiple)" << endl;
	if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
	int wknown = known_b3 | active_b3;
	if (!active_b3) {
		if (diagh)cout << "final check2" << endl;
		CriticalFinalCheck();
		return; // should be seen earlier if possible
	}
	if (rknown_b3 != wknown)
	if (IsMultiple(wknown))	return;// not valid using all cells
	rknown_b3 = wknown;
	if (diagh || g17b.debug17)cout << "start  irloop="<<irloop << endl;
	if (irloop)		CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalLoop(int *to, int no){
	if (g17b.debug17>1)cout << "critical loop" << endl;
	if (ShrinkUas1(to, no)) return;
	if (irloop)CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalExitLoop(int *uasb3, int nuasb3){
	int nmissb = g17xy.b3lim - _popcnt32(known_b3);// missing clues
	if (diagh || g17b.debug17>1 ){
		cout << Char27out(known_b3) << "known exit loop nuasb3= " << nuasb3 << endl;
		cout << Char27out(active_b3) << "active exit loop nmissb="<<nmissb << endl;
	}
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
	if (g17b.debug17>1) {
		cout<<Char9out(wg3.count.u16[0]) <<   "\tlook for next uab3 nuab3="<<nuasb3   << endl
			<<" priority to a smallest nuab3 if nmissb=1 else priority to pair if exists"<<endl;
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
		else if (wg3.count.u16[0]){	// use in priority an unsolved pair it is a smallest
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[0]);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3&mask;// catch the minirow
		}
		else{
			for (int i = 0; i < nuasb3; i++){
				register int ua = uasb3[i], cc = _popcnt32(ua);
				if (cc < sizeua){ wua = ua; sizeua = cc; }
				if (cc < 3)break; // this is the minimum 
			}
			if (sizeua >= 2 && wg3.minirows_triplets){// use the triplet in priority
				uint32_t mini;
				bitscanforward(mini, wg3.minirows_triplets);
				int  shift = 3 * mini, mask = 7 << shift;
				wua = active_b3 &mask;// catch the minirow

			}
		}
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			if (g17b.debug17>1)cout << Char27out(wua) << " wua to assign in g17b.debug17 mode" << endl;
		}

		while (bitscanforward(cell, wua)){
			register int bit = 1 << cell;
			wua ^= bit;// clear bit
			// clean the bit in active_b3, this is now a dead cell downstream
			active_b3 ^= bit;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bit);
			hn.CriticalLoop(uasb3, nuasb3);
		}
	}
	else Critical_0_UA(); // no more ua, some non assigned pairs or triplets
}
void G17B3HANDLER::Critical_0_UA(){
	int nmissb = g17xy.b3lim - _popcnt32(known_b3);// missing clues
	if (g17b.debug17>1 ){
		cout << Char27out(known_b3) << "known 0 ua" << endl;
		cout << Char27out(active_b3) << "active 0 ua " << endl;
	}
	if (nmissb < 0)return;
	if (!nmissb){// nothing more to assign (granted at first call in a branch)
		CriticalFinalCheck();
		return;
	}
	if (wg3.count.u16[2])	{// in active minirows with 3 pairs, assign 2
		while (wg3.count.u16[2]){
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[2]);
			int shift = 3 * mini, bit = 1 << shift;
			wg3.count.u16[2] ^= 1 << mini; //clear bit the mini row is always killed
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
	if (wg3.count.u16[0]){// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, wg3.count.u16[0]);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		//wg3.count.u16[0] ^= 1 << mini; //clear bit the mini row is always killed
		int x = active_b3&mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		wg3.count.u16[0] ^= 1 << mini;// and clear the minirow bit as active
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
	if (wg3.minirows_triplets){// safety control should always be
		uint32_t mini;
		bitscanforward(mini, wg3.minirows_triplets);
		int shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		//wg3.minirows_triplets ^= 1 << mini; //clear bit the mini row is always killed
		active_b3 &= ~mask;// clear the minirow
		wg3.minirows_triplets ^= 1 << mini;// and clear the minirow bit as active
		for (int i = 0; i < 3; i++){
			int bb = bit << i;
			G17B3HANDLER hn = *this;
			hn.CriticalAssignCell(bb);
			hn.Critical_0_UA();
		}
	}
}
void G17B3HANDLER::CriticalFinalCheck(){// no more ua is it a valid solution 
	//if (p_cpt2g[17] == 2252)
	int ncl = _popcnt32(known_b3);
	//if (g17b.debug17) cout << "final check test ncl=" << ncl  << endl;
	if (ncl != g17xy.b3lim) return; // should be seen earlier if possible
	if (diagh)cout <<Char27out(known_b3 )<< "critical final check" << p_cpt2g[17] << endl;
	register int ir = IsMultiple(known_b3);// , g17b.debug17);
	if (diagh)cout << "retour" << ir << endl;
	if (ir){
		if (g17b.debug17&& nmiss != g17xy.b3lim) {
			cout << "final check retour faux ir=" <<ir<< endl;
			//diagh = 1;
			//cout << "retour check debug=" << IsMultiple(known_b3) << endl;;
		}
		return;// mode debug pour voir
	}
	g17xy.FoutValid17(known_b3, wg3.ib3);
	g17b.a_17_found_here = 1;
}
//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow(){
	if (g17b.debug17)cout << "entry Go_SubcriticalMiniRow() ndead="<<ndead << endl;
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	G17B3HANDLER hn;
	int bit = 1<<ndead, mask = 7<<(3*ndead), stack = ndead%3;
	for (int i = ndead; i < 9; i++, stack++, bit <<= 1, mask <<= 3){
		if (stack >2) stack = 0;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (g17b.debug17)cout << Char27out(M) << " mini row to process i="<<i << endl;
		if (bit & wg3.count.u16[0])// it was a gua2 pair assign both 
			hn.SubMini(this, M, mask, stack, 1, bit);
		else if (bit & wg3.count.u16[1])// it was 2 gua2 pair assign 2 out of 3 
			for (int j = 0; j < 3; j++){
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(this, M, mask, stack, 2, bit);
			}
		else if (bit & wg3.count.u16[2])// it was 3 gua2 pair assign 3 out of 3 
			hn.SubMini(this, M, mask, stack, 3, bit);
		else if (bit & wg3.minirows_triplets)// it was a gua3 triplet assign 2 out of 3
			for (int j = 0; j < 3; j++){
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.wg3.minirows_triplets ^= bit;
				hn.SubMini(this, M, mask, stack, 4, bit);
			}
		else // second add in the mini row one residual cell take it
			hn.SubMini(this, M, mask, stack, 0, 0);
	}
}
void G17B3HANDLER::SubMini(G17B3HANDLER * o, int M, int mask, int stack, int imode, int bit){
	*this = *o;
	if (imode){
		if (imode<4)wg3.count.u16[imode-1] ^= bit;
		else wg3.minirows_triplets ^= bit;
	}
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added 
	active_b3 &= ~mask;
	active_sub ^= M;
	if (g17b.debug17) {
		if (g17b.debug17>1)cout << Char27out(known_b3) << " known sub mini  imode=" << imode << endl;
		if ((known_b3 & g17b.band3_17) != known_b3) {
			cout << "not the right one, goback" << endl;
			return;
		}
		if (g17b.debug17 > 1)cout << Char27out(active_b3) << " active_b3 nmiss=" << nmiss << endl;
	}
	// now adjust the stack count
	wg3.countstack.u16[stack]++;
	wg3.countstack.u16[3]++;
	if (wg3.countstack.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else{	// leave sub critical mode and enter the critical mode 
		if (g17b.debug17>1) wg3.Debug();
		Critical2pairs();// assign 2 pairs in minirow to common cell
		if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
		rknown_b3 = known_b3 | active_b3;
		if (IsMultiple(rknown_b3))return;// not valid using all cells
		if (g17b.debug17>1)cout << "still valide" << endl;
		if (irloop)		CriticalLoop(uasb3, nuasb3);
		else CriticalExitLoop(uasb3, nuasb3);
	}
}
void G17B3HANDLER::Go_Subcritical(int docheck){// nmiss to select in the critical field
	p_cpt2g[16]++;
	active_b3 = active_sub = wg3.pairs.u32[1];
	// check first if a global solution  is still possible
	if(docheck)if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	int cct = _popcnt32(wg3.pairs.u32[1]) - ncritical;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	for (int ist = 0; ist < 3; ist++){// check stacks 
		if (wg3.countstack.u16[ist]>5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
	uasb3 = g17xy.uasb3_1;
	nuasb3 = g17xy.nuasb3_1;
	if (g17b.debug17)cout <<Char27out(active_sub )
		<< " active sub initial call to Go_SubcriticalMiniRow() nmiss="<<nmiss << endl;
	Go_SubcriticalMiniRow();// find the first miss
}
//======================================================================= not critical sequence


#define FINDWUA {register int Ra = wactive0, Rfilt =  known_b3;\
for (int iua = 0; iua < g17xy.nuasb3_2; iua++){\
register int Ru = g17xy.uasb3_2[iua];\
if (Ru & Rfilt) continue;\
Ru &= Ra;\
register int cc = _popcnt32(Ru);\
if (cc < ncells)	{ ncells = cc; wua = Ru; rawua = g17xy.uasb3_2[iua]; }\
if (!cc)return;\
if (cc > 6 && ncells < 5) break;}}


#define APPLYWUA(G)uint32_t res;\
int x = wua;\
while (bitscanforward(res, x)){\
int bit = 1 << res;x ^= bit;wactive0 ^= bit;\
G17B3HANDLER hn = *this;hn.nmiss--;hn.known_b3 |= bit;\
int stack = C_box[res];\
hn.wg3.countstack.u16[stack]++;\
hn.wg3.countstack.u16[3]++;\
if (hn.wg3.countstack.u16[stack] > 5){\
hn.wactive0 &= ~(07007007 << (3 * stack));}\
G;}
void G17B3HANDLER::Go_Not_Critical_miss5(){ // select and assign a small ua outfield then miis4
	// or subfield 2 cells
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			cout << Char27out(wua) << " wua nmiss5  ncells=" << ncells << endl;
		}
	APPLYWUA(hn.Go_Not_Critical_miss4())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 4 cells
}
void G17B3HANDLER::Go_Not_Critical_miss4(){ // select and assign a small ua outfield then miis2
	// or subfield 2 cells
	if (g17b.debug17) {
		if (known_b3 & (~g17b.band3_17)) return;
		cout << Char27out(known_b3) << " known entry nmiss4" << endl;
	}
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			cout << Char27out(wua) << " wua nmiss 4    ncells=" << ncells << endl;
		}
	APPLYWUA(hn.Go_Not_Critical_miss3())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 3 cells
}
void G17B3HANDLER::Go_Not_Critical_miss3(){ // select and assign a small ua outfield then miis2
	// or subfield 2 cells
	if (g17b.debug17) {
		cout << Char27out(known_b3) << " known entry nmiss3" << endl;
		if (known_b3 & (~g17b.band3_17)) return;
	}
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			cout << Char27out(wua) << " wua nmiss 3  ncells=" << ncells << endl;
		}
	APPLYWUA(hn.Go_Not_Critical_miss2())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 3 cells
}
void G17B3HANDLER::Go_Not_Critical_miss2(){ // select and assign a small ua outfield then miis1
	// or subfield 2 cells
	if (g17b.debug17)	{
		cout << Char27out(known_b3) << "entry nmiss2" << endl
			<< Char27out(g17b.band3_17) << endl;
		if (known_b3 & (~g17b.band3_17)) return;
		cout << "start critical miss2" << endl;
		//diagh = 1;
	}

	if (IsMultiple(wg3.pairs.u32[1] | wactive0 | known_b3))return;
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
		if (g17b.debug17) {
			wua &= g17b.band3_17;
			cout << Char27out(wua) << " wua nmiss 2   ncells=" << ncells << endl;
		}
	APPLYWUA(hn.Go_Not_Critical_miss1())
		if (ncells == 27)Go_Subcritical();// and finish in Subcritical mode here miss 2 cells
}
void G17B3HANDLER::Go_Not_Critical_miss1(){
	//fout1 << Char27out(known_b3) << "entry nmiss1" << endl;
	//cout << Char27out(known_b3) << "entry nmiss1" << endl;	
	if (diagh || g17b.debug17) {
		cout << Char27out(known_b3) << "entry nmiss1 g17xy.nuasb3_2=" << g17xy.nuasb3_2 << endl;
		if (g17b.debug17) {
			if (known_b3 & (~g17b.band3_17)) return;
		}
	}
	int nn = 0,wact_miss1=wactive0;
	{//================= first collect UAs and GUAs46 100% outfield 
		register int Ra = wactive0, Rfilt =  known_b3;
		for (int iua = 0; iua < g17xy.nuasb3_2; iua++){
			register int Ru = g17xy.uasb3_2[iua];
			if (Ru & Rfilt) continue;// hit or cells in the critical field
			Ra &= Ru; // one of these cells must be hit
			// next test will say no for any UA having no cell in the studied field
			nn++;
			if (!Ra) break;// need more than one cell to hit all uas located 
		}
		wact_miss1 = Ra;// now cells hitting all UAs with no link to the partial critical field
	}
	if (diagh|| g17b.debug17) cout << Char27out(wact_miss1) << " wact_miss1 nn= " << nn << endl;
	if (!nn){
		if (IsMultiple(wg3.pairs.u32[1] | known_b3)){// must be one cell out
			//wact_miss1 &= zhou[0].glb.go_back;
			if (!wact_miss1) return;
			goto one_out_forced;
		}
		if (diagh ||g17b.debug17) {
			int field = wg3.pairs.u32[1] | known_b3;
			if (g17b.debug17) {
				if((field & g17b.band3_17)!= g17b.band3_17)
					goto wactive_to_test;
			}
			cout << "miss1 call  hn.Go_Subcritical" << endl;
			cout << Char27out(known_b3) << " known b3 "  << endl;
			cout << Char27out(field) << " known b3 + pairs" << endl;
			//return;
		}
		G17B3HANDLER hn = *this;
		hn.Go_Subcritical(0);// try Subcritical mode in a new handler
		// must try each cell in wactiveas redundant
		goto wactive_to_test;
	}
	//else return;
one_out_forced:
	{
		uint32_t bf = (wg3.pairs.u32[1] | wact_miss1 | known_b3);
		if (diagh ||g17b.debug17) {
			cout << " miss 1 entry one_out_forced" << endl;
			cout << Char27out(bf) << " bf " << endl;
		}
		if (wact_miss1&&IsMultiple(bf))	return;
		if (0) {// try to reduce the count using stacks to be revised, not ok here
			//if (g17b.debug17 > 1)cout << " one_out_forced:" << endl;
			int st = 07007007, rmiss1 = wact_miss1;
			for (int ist = 0; ist < 3; ist++) {
				int wastack = wact_miss1 & st;
				if (wastack && (wastack != rmiss1)) {
					int bfstack = wg3.pairs.u32[1] | wastack | known_b3;
					if (IsMultiple(bfstack)) {
						//wact_miss1 ^= wastack;
						//wact_miss1 &= zhou[0].glb.go_back;// use received ua to reduce the count
					}
				}
			}

		}

	}
wactive_to_test:
	if (g17b.debug17) {
		cout << Char27out(wact_miss1) << " wactive_to_test :" << endl;
		wact_miss1 &= g17b.band3_17;
		cout << Char27out(wact_miss1) << " wactive_to_test valid for 17" << endl;
	}
	if (!wact_miss1) return;
	p_cpt2g[17]++;
	if (diagh) {
		cout << Char27out(wact_miss1) << " wactive_to_test p_cpt2g[17]:"
			<< p_cpt2g[17] << endl;
		cout << Char27out(known_b3) << " known b3" << endl;
		wg3.Debug();

	}
	//if (1) return;
	uint32_t res;
	while (bitscanforward(res, wact_miss1)){
		int bit = 1 << res;
		wact_miss1 ^= bit;// clear bit
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		int stack = C_box[res];
		hn.wg3.countstack.u16[stack]++;
		hn.wg3.countstack.u16[3]++;
		if (diagh) cout << "go critical res=" << res << endl;
		//if(diagh && res==23)
		hn.Go_Critical();
	}
	return;
}

//==================================== handler direct 17 to test

void SPOT17_NOUAS2_3::newspot(int * oua, int onua, SPOT17_NOUAS2_3 * sp, int cell, int bit){
	if (sp){
		*this = *sp;// copy previous except if first
		// apply cell to stack count and check for the limit in stack
		known |= bit;
		int stack = C_box[cell];
		stacks.u16[stack]++;
		if (stacks.u16[stack] > 5)	active &= ~(07007007 << (3 * stack));
	}
	stacks.u16[3]--;// reduce the "missing cells" number
	// early check if 2 stacks filled 
	if (_popcnt32(active)<10 && g17xy.g17hh0.IsMultiple(known | active)) return;
	tua = &oua[onua];
	nua = 0;
	int wua = active, ncells = 27;
	// if last step, must hit all remaining uas
	if (stacks.u16[3]){//================= first find a small ua
		register int Ra = active, Rfilt = known;
		for (int iua = 0; iua < onua; iua++){//and now UAs band3
			register int Ru = oua[iua];
			if (Ru & Rfilt) continue;// alreadyhit 
			Ru &= Ra;// no dead cells
			register int cc = _popcnt32(Ru);
			if (!cc)return; // dead branch
			if (cc < ncells)	{ ncells = cc; wua = Ru; }
			tua[nua++] = Ru;
		}
	}
	else{
		register int Ra = active, Rfilt = known;
		for (int iua = 0; iua < onua; iua++){//and now UAs band3
			register int Ru = oua[iua];
			if (Ru & Rfilt) continue;// alreadyhit 
			Ra &= Ru;// common cells
			if (!Ra)return; // can not hit all uas
		}
		if (g17xy.g17hh0.IsMultiple(known | Ra)) return;// early check for last
		wua = Ra;
	}
	// apply the found wua
	uint32_t res;
	int x = wua;
	while (bitscanforward(res, x)){
		int bit = 1 << res;
		x ^= bit;// clear bit
		active ^= bit;// clear cell in that branch (dead cell)
		if (stacks.u16[3]){
			SPOT17_NOUAS2_3 sn;
			sn.newspot(tua, nua, this, res, bit);
		}
		else{// last step
			g17xy.g17hh0.known_b3 = known | bit;
			g17xy.g17hh0.CriticalFinalCheck();
		}
	}
}


//================ debugging code
void G17Debug_CoutEntry(int * g) {
for (int i = 0; i < 54; i++) cout << g[i] + 1;
cout << "; entry canonical" << endl;
}
const char * g17tl1[10] = { "entry", "steps1", "steps2", "steps3", "steps4",
"common", "f1", "dead", "n51", "n62" };


const char * g17tl[40] = { "0__", "falluas", "fstack1", "UasGuas", "fb12v", "  b3v\t\t", "fnov", //0-6
"h3go", "critb", "ncritb", "m1go", "subgo", "subvv ",   //12
"gcrit", "gcvv", "gc_exit", "m2", "m2vv", "m1_ss", "19_", "fin check",//20
"21_", "22_", "23_", "24_", "25_", "26__", "gocheck", "check0", "nocheck0",//29
"30nua", "nua2", "nua3", "nb3", "__34", "vt1", "vt2", "", "T1", "T2",//39
};

//=============================== debugging sequences
void G17B::GodebugInit(int mode) {
	if (!mode) return;
	int n1 = myband1.nind[1],
		n2 = myband2.nind[1];
	//zhou_i.ImageCandidats();
	//myband1.PrintShortStatus(); myband2.PrintShortStatus();
	cout << "n bands3      \t" << genb12.nband3 << endl;
	cout << "ua bands1+2   \t" << genuasb12.nua << endl;
	cout << "guas socket2  \t" << genb12.ntua2 << endl;
	cout << "guas socket3  \t" << genb12.ntua3 << endl;
	cout << "active socket2\t" << genb12.nactive2 << endl;
	cout << "active socket3\t" << genb12.nactive3 << endl;
	cout << "n1=" << n1 << " n2=" << n2 << " n5 "<< myband1.n5<<";"<< myband2.n5
		<< "n6  " << myband1.n6 << ";" << myband2.n6 << endl;
	uint64_t t1 = (uint64_t)myband1.n5 *(uint64_t)myband2.n6;
	uint64_t t2 = (uint64_t)myband2.n5 *(uint64_t)myband1.n6;
	cout << t1 << endl << t2 << endl<<t1+t2<< endl;
	if (mode == 2) {
		cout << "sockets 2 table" << endl;
		int n2 = 0;
		for (int i = 0; i < genb12.nactive2; i++) {
			int i81 = genb12.tactive2[i];
			GEN_BANDES_12::SGUA2 & w = genb12.tsgua2[i81];
			cout << i81 << " " << w.nua << endl;
			n2 += w.nua;
		}
		cout << "cumul=" << n2 << endl;
		cout << "sockets 3 table" << endl;
		int n3 = 0;
		for (int i = 0; i < genb12.nactive3; i++) {
			int i81 = genb12.tactive3[i];
			GEN_BANDES_12::SGUA3 & w = genb12.tsgua3[i81];
			cout << i81 << " " << w.nua << endl;
			n3 += w.nua;
		}
		cout << "cumul=" << n3 << " total=" << n2 + n3 << endl;
	}

}
int G17B::GodebugFindKnown17() {// locate the known 17 
	if (0) {
		XY_EXPAND * tx1 = myband1.xye6;
		cout << "table 6 bande 1 index2" << endl;
		for (int i = 0; i < myband1.nind[1]; i++) {
			int * t = myband1.index2[i];
			cout << Char27out(t[0]) << " i=" << i << " start=" << t[2]
				<< " end=" << t[5] << endl;
		}
	}
	if (0) {
		XY_EXPAND * tx1 = myband1.xye6;
		cout << "table 5 bande 22 index2" << endl;
		for (int i = 0; i < myband2.nind[1]; i++) {
			int * t = myband2.index2[i];
			cout << Char27out(t[0]) << " i=" << i << " start=" << t[1]
				<< " end=" << t[4] << endl;
		}

	}
	if (1) {
		XY_EXPAND * tx1 = myband1.xye6;
		cout << "table 6 bande 1" << endl;
		for (int i = 1872; i < 2097 /*myband1.n6*/; i++)
			cout << Char2Xout(tx1[i].cellsbf) << " i=" << i << endl;

	}
	if (0) {
		XY_EXPAND * tx2 = myband2.xye5;
		cout << "table 5 bande 2" << endl;
		for (int i = 3778; i < 4073 /*myband2.n5*/; i++)
			cout << Char2Xout(tx2[i].cellsbf) << "i=" << i << endl;

	}
	return 0;
}
int G17B::GodebugCheckUas(const char * lib) {
	uint32_t nua = genuasb12.nua;
	uint64_t * tua= genuasb12.tua;
	for (uint32_t i = 0; i < nua; i++) {
		if (tua[i] & g17b.band12_17) continue;
		cout << lib << "check ua failed" << endl;
		cout << Char2Xout(tua[i]) << " not hit by the known 17" << endl;
		return 1;
	}
	return 0;
}
int G17INDEXSTEP::DebugIsKnown() {
	if (_popcnt32(g17b.band1_17) == 6) {
		if (g17b.debug17>1)cout << "check known in step  mode 6_5 " << endl;
		XY_EXPAND * tx1 = xyew1, *tx2 = xyew2;
		// w1 => Y  w2 => X
		for (int i = 0; i < nw1; i++) {
			if ((tx1[i].cellsbf &g17b.band1_17) == g17b.band1_17){
				if (g17b.debug17 > 1)cout << "seen band1 i=" << i << endl;;
				goto check2;
			}
		}
		return 0; // first band not seen
		check2:;
		for (int i = 0; i < nw2; i++) {
			if (((tx2[i].cellsbf>>32) &g17b.band2_17) == g17b.band2_17) {
				if (g17b.debug17 > 1)cout << "seen band2 i=" << i << endl;
				return 1;
			}
		}
		return 0;//  band2 not seen
	}
	else {
		if (g17b.debug17 > 1)cout << "check known in step  mode 5_6 " << endl;
		XY_EXPAND * tx1 = xyew2, *tx2 = xyew1;
		// w1 => Y  w2 => X
		for (int i = 0; i < nw2; i++) {
			if ((tx1[i].cellsbf &g17b.band1_17) == g17b.band1_17) {
				if (g17b.debug17 > 1)cout << "seen band1 i=" << i << endl;;
				goto check2a;
			}
		}
		return 0; // first band not seen
	check2a:
		for (int i = 0; i < nw1; i++) {
			if (((tx2[i].cellsbf >> 32) &g17b.band2_17) == g17b.band2_17) {
				if (g17b.debug17 > 1)cout << "seen band2 i=" << i << endl;
				return 1;
			}
		}
		return 0;//  band2 not seen

	}
	return 0;
}

void V256_UAS::Debug(const char * lib, int mirror) {
	cout << " v256_ua for " << lib << endl;
	uint32_t * t32=v[0].bf.u32;
	for (int i = 0; i < 8; i++){
		int w = t32[i];
		if (mirror) w =~w;
		if(w)
		cout << Char32out(w) << " " << 32 * i << "-" << 32 * (i + 1) - 1 << endl;
	}
}
void V256_UAS::Fout(const char * lib) {
	fout1 << lib << v[0].bf.u64[0] << " " << v[0].bf.u64[1] << " "
		<< v[1].bf.u64[0] << " " << v[1].bf.u64[1] << endl;
}
void V256_UAS::Cout() {
	cout << Char64out(v[0].bf.u64[0]);
	if(indexstep.ntua>64)	cout << Char64out(v[0].bf.u64[1]);
	if (indexstep.ntua > 128)cout << Char64out(v[1].bf.u64[0]);
	if (indexstep.ntua > 192)cout << Char64out(v[1].bf.u64[1]);
	cout<< endl;
}

void G17INDEXSTEP::DebugIndex(int i1, int i2) {
	int t0 = n51 * n62 + n52 * n61;
	fout1<<i1<<" "<<i2 << " debug index"  
		<<" n51="<<n51<<" n61="<<n61
		<< " n52=" << n52 << " n62=" << n62
		<<" ntua="<<ntua << " ntgua=" << ntgua
		<<"\tnua"<<genuasb12.nua<<endl
		<<t0 << " attendus  faits" << p_cpt2g[2]<< " futur"<< p_cpt2g[2] +t0<< endl;
	//PrintUasShrinked();
	if (0) {
		v256a.Debug("v256a");
		for (int ic = 0; ic < 54; ic++) {
			cout << "cell=" << ic << endl;
			v256uas[ic].Debug(" mirror cell", 1);
		}

	}
}
int G17CHUNK::DebugIsKnown() {
	if (_popcnt32(g17b.band1_17) == 6) {
		if (g17b.debug17 > 1)cout << "chunk check known in chunk  mode 6_5 " << endl;
		XY_EXPAND * tx1 = tyec, *tx2 = txec;
		// w1 => Y  w2 => X
		for (int i = 0; i < nyc; i++) {
			if ((tx1[i].cellsbf &g17b.band1_17) == g17b.band1_17) {
				if (g17b.debug17 > 1)cout << "chunk seen band1 i=" << i << endl;;
				iy17 = i;
				goto check2;
			}
		}
		return 0; // first band not seen
	check2:;
		for (int i = 0; i < nxc; i++) {
			if (((tx2[i].cellsbf >> 32) &g17b.band2_17) == g17b.band2_17) {
				ix17 = i;
				cout << "chunk seen band2 ix17=" << ix17<<" iy17="<<iy17 
					<<"\ti1="<<i1<<"i2="<<i2<< endl;
				goto checkok;
			}
		}
		return 0;//  band2 not seen
	}
	else {
		if (g17b.debug17 > 1)cout << "chunk check known in chunk  mode 5_6 " << endl;
		XY_EXPAND * tx1 = txec, *tx2 = tyec;
		// w1 => Y  w2 => X
		for (int i = 0; i < nxc; i++) {
			if ((tx1[i].cellsbf &g17b.band1_17) == g17b.band1_17) {
				if (g17b.debug17 > 1)cout << "chunk seen band1 i=" << i << endl;;
				ix17 = i;
				goto check2a;
			}
		}
		return 0; // first band not seen
	check2a:;
		for (int i = 0; i < nyc; i++) {
			if (((tx2[i].cellsbf >> 32) &g17b.band2_17) == g17b.band2_17) {
				iy17 = i;
				cout << "chunk seen band2 ix17=" << ix17 << " iy17=" << iy17 << endl;
				goto checkok;
			}
		}
		return 0;//  band2 not seen

	}
	return 0;

checkok: 
	if (g17b.debug17 > 1) { // print the X vector U vector
		xv[ix17].Cout();
		yv[iy17].Cout();

	}
	return 1;
}

void G17INDEXSTEP::IndexStepDebugKnown17(int i1, int i2) {
	if (g17b.debug17 > 1) {

		cout << Char27out(t2[0]) << " step 2 BF i2=" << i2 << "this is the right i2" << endl;
		cout << " n51=" << n51 << " n61=" << n61 << " n52=" << n52 << " n62=" << n62
			<< " n51*n62=" << n51 * n62
			<< "\tnua step=" << ntua << "\tp_cpt2g[3]=" << p_cpt2g[3] << endl;
		cout << "band1 indexstep" << endl;
		int * t1 = myband1.index2[i1];
		cout << Char27out(t1[0]) << " i1=" << i1
			<< " start6=" << t1[2] << " end6=" << t1[5]
			<< " start5=" << t1[1] << " end5=" << t1[4] << endl;
		cout << "band2 indexstep" << endl;
		int * t2 = myband2.index2[i2];
		cout << Char27out(t2[0]) << " i2=" << i2
			<< " start6=" << t2[2] << " end6=" << t2[5]
			<< " start5=" << t2[1] << " end5=" << t2[4] << endl;
	}
}
void G17B::PrintEndPuz() {
	cout << endl;
	for (int i = 0; i < 10; i++) if (p_cpt1[i])
		cout << g17tl1[i] << "\t" << p_cpt1[i] << endl;
	cout << endl;
	for (int i = 0; i < 32; i++) if (p_cpt[i])
		cout << g17tl[i] << "\t" << p_cpt[i] << endl;
}

void G17TB3GO::Debug() {
	cout << "status G17TB3GO ib3=" << ib3 << endl;
	//char ws[82];
	//cout << g17xy.bands_active_triplets.String3X(ws)
	//	<< "g17xy active triplets" << endl;
	cout << genb12.bands3[ib3].band << " \tbande3" << endl;
	cout << Char27out(pairs.u32[1]) << " in field bf" << endl << endl;
	cout << Char27out(pairs.u32[0]) << " pairs bf" << endl<<endl;
	cout << Char9out(count.u16[3]) << " all minis2" << endl;
	cout << Char9out(count.u16[0]) << "     minis2/1" << endl;
	cout << Char9out(count.u16[1]) << "     minis2/2" << endl;
	cout << Char9out(count.u16[2]) << "     minis2/3" << endl<<endl;
	cout << Char9out(minirows_triplets) << " mini triplets" << endl << endl;

	cout << countsum.u16[0] << "\t" << countsum.u16[1] << "\t" << countsum.u16[2]
		<< "\t total " << countsum.u16[3] << endl;
	cout << countstack.u16[0] << "\t" << countstack.u16[1] << "\t" << countstack.u16[2]
		<< "\t total " << countstack.u16[3] << endl;
}