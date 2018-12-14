G17B::G17B() {
	register uint64_t R = (uint64_t)&xye6[0][1];
	R -= (uint64_t)&xye6[0][0];
	xyexpand_size = (int)(R >> 3);

	cout << "G17 constructor chunksize (64 bits) for XYexpand one size=" << xyexpand_size << endl;
}

void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	if (p_cpt1g[0] ++)	return; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<first test 
	int * zs0 = genb12.grid0;
	myband2.DoExpandBand(27);// expand band2
	if (!(myband1.n5| myband2.n5)) return; // no 656 no 566
	int nb3 = genb12.nband3;
	//=========================== collect UAs (still room for improvement)
	genuasb12.Initgen();
	genb12.BuildGang9x3();
	if (0) {
		cout << "end gen UAs bands 1+2 nb3="<<nb3<<" nua="<< genuasb12.nua << endl;
		zh2b[0].Init_std_bands();
		zh2b[0].ImageCandidats();
		cout << "controle de solution bande 1 2 par digit" << endl;
		zh2b[0].DebugSol();
		cout << "gangster band 3" << endl;
		for (int i = 0; i < 9; i++) {
			int * w = genb12.gang[i];
			cout << w[0] + 1 << w[1] + 1 <<w[2] + 1<<" ";
		}
		cout << endl;
		cout << "start collect gua2s gua3s" << endl;
	}
	if (1) {// print first uas in the global table for debugging
		cout << "first uas in the global table" << endl;
		for (uint32_t i = 0; i < genuasb12.nua; i++){
			if (i == 6) break;
			cout << Char2Xout(genuasb12.tua[i]) << endl;
		}

	}


	genb12.SecondSockets2Setup();// collect GUA2s 
	if (0) {// check bands 
		for (int i = 0; i < 1; i++)
			genb12.bands3[i].PrintB3Status();
		return;
	}
	if (1) return;
	//	genb12.CollectUA2s();// collect GUA2s
	//	genb12.CollectUA3s();//collect GUA3s
	//	genb12.CollectMore();// Add more medium size UAs bands 1+2 with special patterns 
	Go();// standard entry point for all 
	if (p_cpt1g[0])g17b.PrintEndPuz();
}

void G17B::Go(){// search process external loop 2X2Y
	g17more.Init();
	g17morebig.Init();
	for (int i = 0; i < 3; i++)g17moreg[i].Init();
	indexstep.StartDead();
	zh_g.grid0 = genb12.grid0;
	int n1 = myband1.nind[0],
		n2 = myband2.nind[1];
	for (int i1 = 0; i1 < n1; i1++){
		if (GTEST17_ON)if (Godebug17_1(i1)) continue;
		for (int i2 = 0; i2 < n2; i2++){
			if (GTEST17_ON)if(Godebug17_2(i2)) continue;

			indexstep.added_more = 0;
			indexstep.Init(i1, i2);
			g17more.Init();
			g17morebig.Init();
			if (indexstep.ShrinkUas())continue;// dead branch
			indexstep.ShrinkGuas();

			if (GTEST17_ON)Godebug17_3();// dead ntua ntgua

			if (indexstep.n52){
				indexstep.Do65();
				if (genuasb12.nua <4000 && indexstep.n64vua< G17BLOCSUA-2){
					indexstep.GoAddNewUas_sub_step();
				}
				indexstep.GoAddNewGUas_sub_step();
			}
			if (indexstep.n51)				indexstep.Do56();

			//indexstep.Do66(); // keep it for a 18 222222222 search
			GoAddNewUas();
			GoAddNewGUas();
			g17more.Init();// restart with empty more tables for the next step
			g17morebig.Init();
			for (int i = 0; i < 3; i++)g17moreg[i].Init();
		}
		if (GTEST17_ON || g17b.debug17>1)
			Godebug17_4(i1, n2);// i1 i1 nua nua2 nua3 p_cpt[4]
	}
}

void G17B::GoAddNewUas(){// end step include small more in the the main uas table
	// tua is supposed to have a capacity of TUA64_12SIZE uas
	register int N = g17more.nt, ndeb = genuasb12.nua;
	if (N){
		uint32_t & nua = genuasb12.nua;
		uint64_t * tt = genuasb12.tua, *tt0 = g17more.t;
		while (N-- && nua < TUA64_12SIZE){
			genuasb12.AddUA(*tt0);
			tt0++;
		}
	}
	if (0){
		cout << "more end step" << g17more.nt << " ndeb=" << ndeb
			<< " final count=" << genuasb12.nua << " ecart " << genuasb12.nua - ndeb << endl;
		g17more.PrintUasDirect();
	}


}

void G17B::GoAddNewGUas(){// end step include small more guas in the guas tables
	if (genb12.ntua2 >= 3000) return;
	register int nt = g17moreg[0].nt;
	register uint64_t * t = g17moreg[0].t;
	while (nt--){
		register  uint64_t v = *t, vid = v >> 56;
		//cout << Char64out(v) << "v" << endl;
		if (vid>80){
			if (genb12.ntua3 < 3000) InsertGua(genb12.tua3, genb12.ntua3, v);
		}
		else if (genb12.ntua2 < 3000)	InsertGua(genb12.tua2, genb12.ntua2, v);
		t++;
	}

}

void G17B::InsertGua(uint64_t * tu, uint32_t & ntu, uint64_t gu){
	GINT64 guw; guw.u64 = gu;
	if (guw.u8[7] >= 81)guw.u8[7] -= 81;
	register uint64_t w = guw.u64;
	for (uint32_t iua = 0; iua < ntu; iua++){
		register uint64_t r = tu[iua];
		if (w > r) continue;
		if (w == r) return;
		// insert the new here
		uint64_t temp[8000]; // temporary for copy forward
		__movsq(temp, &tu[iua], (ntu - iua));
		__movsq(&tu[iua + 1], temp, (ntu - iua));
		tu[iua] = w;
		ntu++;
		return;
	}
	tu[ntu++] = w; // add to the end
}


void G17INDEXSTEP::GoAddNewUas_sub_step(){// partial step build in vector 
	if (n64vua >= G17BLOCSUA) return; // do nothing if no room
	// open a new vector with available small uas and store them
	register int N = g17more.nt, ndeb = genuasb12.nua;
	if (N <= 15)goto morebig;// do nothing if small count
	// open a new vector
	register int bloc = n64vua;
	indexstep.vauas[bloc] = maskLSB[N].u64[0];
	//uint32_t & nua = genuasb12.nua;
	//uint64_t * tt = genb12.tua;

	uint32_t cc;
	for (int i = 0; i <N; i++){// set uas
		register uint64_t Rw = g17more.t[i], biti = 1;
		genuasb12.AddUA(Rw);
		biti <<= i;
		while (bitscanforward64(cc, Rw)){// look for  possible cells
			register uint64_t bit2 = 1; bit2 <<= cc;
			Rw ^= bit2;// clear bit
			vcellsuas[From_128_To_81[cc]][bloc] ^= biti;
		}
	}
	if (0){
		cout << "more step" << g17more.nt << " ndeb=" << ndeb
			<< " final count=" << genuasb12.nua << " ecart " << genuasb12.nua - ndeb
			<< " n64vua avant =" << n64vua << endl;
		g17more.PrintUasDirect();
	}

	g17more.Init();
	n64vua++;
	g17chunk.n64vua = n64vua;
morebig:// just put big in vectors 
	//if (1) return;
	N = g17morebig.nt;
	if (N < 20) return;
	if (n64vua >= G17BLOCSUA) return; // do nothing if no room
	// open a new vector
	bloc = n64vua;
	vauas[bloc] = maskLSB[N].u64[0];
	for (int i = 0; i <N; i++){// set uas
		register uint64_t Rw = g17morebig.t[i], biti = 1;
		biti <<= i;
		while (bitscanforward64(cc, Rw)){// look for  possible cells
			register uint64_t bit2 = 1; bit2 <<= cc;
			Rw ^= bit2;// clear bit
			vcellsuas[From_128_To_81[cc]][bloc] ^= biti;
		}
	}
	g17morebig.Init();
	n64vua++;
	g17chunk.n64vua = n64vua;
}


void G17INDEXSTEP::GoAddNewGUas_sub_step(){// end step include small more guas in the guas tables
	// copy mores in the appropriate table
	int imore = added_more++;
	G17TMORE * tm = g17moreg0;
	memcpy(tm, g17moreg, sizeof g17moreg);
	// store small in the main tables
	register int nt = tm[0].nt;
	register uint64_t * t = tm[0].t;
	while (nt--){
		register  uint64_t v = *t, vid = v >> 56;
		if (vid>80){
			if (genb12.ntua3 < 3000) g17b.InsertGua(genb12.tua3, genb12.ntua3, v);
		}
		else if (genb12.ntua2 < 3000)	g17b.InsertGua(genb12.tua2, genb12.ntua2, v);
		t++;
	}
	for (int i = 0; i < 3; i++)g17moreg[i].Init();// reset more tables
	// build vectors for tm 
	//initial vectors to appropriate value and actives vector
	//register uint64_t In = ~0;
	memset(vmcellsguas, 255, sizeof vmcellsguas);
	memset(vmid81, 255, sizeof vmid81);
	for (int i = 0; i < 3; i++)vmuas[i] = maskLSB[tm[i].nt].u64[0];

	uint32_t cc;
	// and setup vectors for each table
	for (int itab = 0; itab < 3; itab++){
		int nu = tm[itab].nt;
		register uint64_t * tu = tm[itab].t;
		for (int i = 0; i < nu; i++){// set uas
			register uint64_t Rw = tu[i], biti = (uint64_t)1<<i, Rid = Rw >> 56;
			Rw &= 0xffffffffffffff;// forget the i81 id
			vmid81[Rid][itab] ^= biti;// flag it to clear later the id
			while (bitscanforward64(cc, Rw)){// look for  possible cells
				register uint64_t bit2 = (uint64_t)1 << cc;
				Rw ^= bit2;// clear bit
				vmcellsguas[cc][itab] ^= biti;
			}
		}
	}
}

void G17INDEXSTEP::Init(int i1, int i2){
	t1 = myband1.index2[i1];
	t2 = myband2.index2[i2];
	n51 = t1[4] - t1[1];
	n61 = t1[5] - t1[2];
	n52 = t2[4] - t2[1];
	n62 = t2[5] - t2[2];

	//update the dead situation
	XY_EXPAND xe1 = xye6[0][t1[2]],
		xe2 = xye6[1][t2[2]];// here cell is 27-53
	int cell1 = xe1.cells.u8[0], cell2 = xe2.cells.u8[0] - 27;
	int cell1bf = 1 << cell1, cell2bf = 1 << cell2;
	if (!(cell1bf&oldx1)){// new cell 1
		oldx1 = cell1bf;
		oldx2 = t1[0];
		xfirstdead |= cell1bf;
		x2dead = xfirstdead | t1[0];
		yfirstdead = cell2bf;
		oldy1 = cell2bf;
		y2dead = t2[0];
	}
	else if ((t1[0] & (~oldx2))){// new xcell 2
		oldx2 = t1[0];
		x2dead = x2dead | t1[0];
		yfirstdead = cell2bf;
		oldy1 = cell2bf;
		y2dead = t2[0];
	}
	else if (!(cell2bf& oldy1)){// new first y
		oldy1 = cell2bf;
		yfirstdead |= cell2bf;
		y2dead = yfirstdead | t2[0];
	}
	else y2dead |= t2[0];

	dead = dead54 = y2dead;
	dead <<= 32;
	dead |= x2dead;
	dead54 <<= 27;
	dead54 |= x2dead;
	b1b2_2Xbf = b1b2_54bf = t2[0];
	b1b2_2Xbf <<= 32;
	b1b2_54bf <<= 27;
	b1b2_2Xbf |= t1[0];
	b1b2_54bf |= t1[0];
}
int G17INDEXSTEP::ShrinkUas(){// shrink table of UAs for bands1+2
	// could erase new supersets to consider doing it in 2 steps
	uint64_t * otua = genuasb12.tua;
	int  ontua = genuasb12.nua;
	register uint64_t Rw = b1b2_2Xbf, Rn = ~dead;
	ntua = 0;
	for (int iua = 0; iua < ontua; iua++){
		register uint64_t Ru = otua[iua];
		if (Rw & Ru) continue;
		Ru &= Rn;// erase dead cells
		if (Ru)		tua[ntua++] = Ru;
		else return 1;// dead branch
	}
	if (ntua > 64 * G17BLOCSUA) ntua = 64 * G17BLOCSUA; //working with a mawimum of 64 * G17BLOCSUA UAs for an index chunk
	n64vua = (ntua + 63) >> 6;
	g17chunk.n64vua = n64vua;
	{ // initial vectors to appropriate value and actives vector
		memset(vcellsuas, 255, sizeof vcellsuas);
		register uint64_t * Ra = vauas;
		register int Rn = ntua;
		for (int i = 0; i < G17BLOCSUA; i++, Rn -= 64){
			if (Rn > 64)vauas[i] = BIT_SET_64;
			else if (Rn > 0)vauas[i] = maskLSB[Rn].u64[0];
			else vauas[i] = 0;
		}
	}
	uint32_t cc;
	for (int i = 0; i < ntua; i++){// set uas
		register int bloc = i >> 6, ir = i & 63;
		register uint64_t Rw = tua[i], biti = (uint64_t)1<< ir;
		while (bitscanforward64(cc, Rw)){// look for  possible cells
			register uint64_t bit2 = (uint64_t)1<< cc;
			Rw ^= bit2;// clear bit
			vcellsuas[From_128_To_81[cc]][bloc] ^= biti;
		}
	}
	return 0;
}
void G17INDEXSTEP::ShrinkGuas(){// shrink table for gangster uas2 ua3s
	register uint64_t Rn;
	{// put dead in 54 bit mode GUAs are in the 54 first bits
		register uint64_t w = dead;
		w <<= 32; w >>= 32; Rn = w;// low 27 bits
		w = dead; w >>= 32; w <<= 27; Rn |= w; // high 27 bits
		Rn = ~Rn; // complement to "and" guas
	}
	uint64_t * otua = genb12.tua2;
	int  ontua = genb12.ntua2;
	pairsok.SetAll_0(); // if empty ua found
	tripletsok.SetAll_0();
	register uint64_t Rw = b1b2_54bf;
	ntgua = 0;
	for (int iua = 0; iua < ontua; iua++){
		register uint64_t Ru = otua[iua];
		if (!(Rw & Ru)){
			Ru &= Rn;
			tgua[ntgua++] = Ru;
		}
	}
	// now Gua3s
	otua = genb12.tua3;
	ontua = genb12.ntua3;
	for (int iua = 0; iua < ontua; iua++){
		register uint64_t Ru = otua[iua];
		Ru += 0x5100000000000000; //socket switched from 0_80 to 81_161
			if (!(Rw & Ru)){
				Ru &= Rn;
				tgua[ntgua++] = Ru;
			}
	}
	if (ntgua > 64 * G17BLOCGSUA) ntgua = 64 * G17BLOCGSUA; //working with a mawimum of 640 GUAs for an index chunk
	n64vgua = (ntgua + 63) >> 6;
	g17chunk.n64vgua = n64vgua;
	{ // initial vectors to appropriate value and actives vector
		memset(vmcellsguas, 255, sizeof vmcellsguas);
		memset(vmid81, 255, sizeof vmid81);
		register uint64_t * Ra = vaguas;
		register int Rn = ntgua;
		for (int i = 0; i < G17BLOCGSUA; i++, Rn -= 64){
			if (Rn > 64)vaguas[i] = BIT_SET_64;
			else if (Rn > 0)vaguas[i] = maskLSB[Rn].u64[0];
			else vaguas[i] = 0;
		}
	}

	uint32_t cc;
	for (int i = 0; i < ntgua; i++){// set uas
		register int bloc = i >> 6, ir = i & 63;
		register uint64_t Rw = tgua[i], biti = 1, Rid = Rw >> 56;
		Rw &= 0xffffffffffffff;// forget the i81 id
		biti <<= ir;
		vid81[Rid][bloc] ^= biti;// flag it to clear later the id
		while (bitscanforward64(cc, Rw)){// look for  possible cells
			register uint64_t bit2 = 1; bit2 <<= cc;
			Rw ^= bit2;// clear bit
			vcellsguas[cc][bloc] ^= biti;
		}
	}
}
void G17INDEXSTEP::Do65(){
	{//=================== move 6 b1 fix place w1
		g17xy.nby = 6;
		nw1 = n61;
		uint64_t *Ro = (uint64_t *)(&xye6[0][t1[2]]), *Rd = (uint64_t *)(xyew1);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw1);
	}
	{//=================== move 5 b2 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n52;
		uint64_t *Ro = (uint64_t *)(&xye5[1][t2[1]]), *Rd = (uint64_t *)(xyew2);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw2);
	}

	g17xy.b3lim = 6;
	Do_Common();
}
void G17INDEXSTEP::Do56(){
	{// ===================== move 6 b2 fix place w1
		g17xy.nby = 6;
		nw1 = n62;
		uint64_t *Ro = (uint64_t *)(&xye6[1][t2[2]]), *Rd = (uint64_t *)(xyew1);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw1);
	}
	{// ===================== move 5 b1 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n51;
		uint64_t *Ro = (uint64_t *)(&xye5[0][t1[1]]), *Rd = (uint64_t *)(xyew2);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw2);
	}
	g17xy.b3lim = 6;
	Do_Common();
}
void G17INDEXSTEP::Do66(){// assumed after do 56 x is band 1 y is band2
	if (!indexstep.n51){// ====== move 6 b2 fix place w1 if not yet done
		g17xy.nby = 6;
		nw1 = n62;
		uint64_t *Ro = (uint64_t *)(&xye6[1][t2[2]]), *Rd = (uint64_t *)(xyew1);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw1);
	}

	{// ===================== move 6 b1 fix place w2
		g17xy.nbx = 6;
		ncluesw2 = 6;
		nw2 = n61;
		uint64_t *Ro = (uint64_t *)(&xye6[0][t1[2]]), *Rd = (uint64_t *)(xyew2);
		__movsq(Rd, Ro, g17b.xyexpand_size*nw2);
	}
	g17xy.b3lim = 5;
	Do_Common();
}


void G17INDEXSTEP::Do_Common(){// after xye w1('Y')  and w2('X') have been loaded
	// find vector x start and vector y start (2 cells in each)
	if (!nw1)return;
	p_cpt1[5]++;
	// setup vectors for the 2 first cells X and Y
	Do_Common_2();// build  Y tables of vectors
	Do_Common_3();// launch chunks

}
void G17INDEXSTEP::Do_Common_2(){// build  Y tables of vectors
	for (int iy = 0; iy < nw1; iy++){ //Y is bloc xyew1
		uint8_t * ty = xyew1[iy].cells.u8;
		uint64_t *d = &tvuas[iy*n64vua],
			*dg = &tvguas[iy*n64vgua];// vectors to build in vuas vguas
		for (int j = 0; j < n64vua; j++)d[j] = vauas[j];
		for (int j = 0; j < n64vgua; j++)dg[j] = vaguas[j];
		for (int i = 2; i < 6; i++){
			int cell = ty[i];
			uint64_t * v = vcellsuas[cell], *vg = vcellsguas[cell];
			for (int j = 0; j < n64vua; j++)d[j] &= v[j];
			for (int j = 0; j < n64vgua; j++)dg[j] &= vg[j];
		}
		if (added_more){
			uint64_t *dg = &tvm1guas[iy*3];// vectors to build in vuas vguas
			for (int j = 0; j < 3; j++)dg[j] = vmuas[j];
			for (int i = 2; i < 6; i++){
				int cell = ty[i];
				uint64_t  *vg = vmcellsguas[cell];
				for (int j = 0; j <3; j++)dg[j] &= vg[j];
			}
		}
	}

}
void G17INDEXSTEP::Do_Common_3_BuildXvectors(){// in fact, in the chunk for 256 X maximum
	// re using the code for 'y', iy is in fact ix
	for (int iy = 0; iy <g17chunk.nxc; iy++){ //X is in the chunk
		uint8_t * ty = g17chunk.txec[iy].cells.u8;
		uint64_t *d = &g17chunk.vxc[iy*n64vua],
			*dg = &g17chunk.vxgc[iy*n64vgua];// vectors to build in vuas vguas				
		for (int j = 0; j < n64vua; j++)d[j] = vauas[j];
		for (int j = 0; j < n64vgua; j++)dg[j] = vaguas[j];
		for (int i = 2; i <g17xy.nbx; i++){
			int cell = ty[i];
			uint64_t * v = vcellsuas[cell], *vg = vcellsguas[cell];
			for (int j = 0; j < n64vua; j++)d[j] &= v[j];
			for (int j = 0; j < n64vgua; j++)dg[j] &= vg[j];
		}
		if (added_more){
			uint64_t *dg = &g17chunk.vm1xgc[iy * 3];// vectors to build in vuas vguas
			for (int j = 0; j < 3; j++)dg[j] = vmuas[j];
			for (int i = 2; i < 6; i++){
				int cell = ty[i];
				uint64_t  *vg = vmcellsguas[cell];
				for (int j = 0; j <3; j++)dg[j] &= vg[j];
			}
		}
	}
}
void G17INDEXSTEP::Do_Common_3(){// launch chunks 256 x256 
	// now combining all 'x' to all 'y' in chunks size ..
	int nxch = (nw2-1) / G17CHKX , nxc;
	int nych = (nw1-1) / G17CHKY, nyc;

	for (int ix = 0, iex = 0; ix <= nxch; ix++, iex += G17CHKX){
		nxc = G17CHKX;
		if (ix == nxch) nxc = nw2 - G17CHKX * ix;
		g17chunk.nxc = nxc;
		{// move XY_EXPAND in the chunk from w2
			uint64_t *Ro = (uint64_t *)(&xyew2[iex]), *Rd = (uint64_t *)(g17chunk.txec);
			__movsq(Rd, Ro, g17b.xyexpand_size*nxc);
		}
		//cout << "before buildx" << endl;
		Do_Common_3_BuildXvectors();// build 'X' tables in the chunk
		//cout << "after buildx" << endl;
		for (int iy = 0, iey = 0; iy <= nych; iy++, iey += G17CHKY){
			nyc = G17CHKY;
			if (iy == nych) nyc = nw1 - G17CHKY * iy;
			g17chunk.nyc = nyc;
			{// move Y_expand in the chunk from w1
				uint64_t *Ro = (uint64_t *)(&xyew1[iey]), *Rd = (uint64_t *)(g17chunk.tyec);
				__movsq(Rd, Ro, g17b.xyexpand_size * nyc);
			}
			// move buas and vguas of the chunk to the chunk
			__movsq(g17chunk.vyc, &tvuas[n64vua*iey], n64vua * nyc);
			__movsq(g17chunk.vygc, &tvguas[n64vgua*iey], n64vgua * nyc);
			
			if (added_more)// and more moreguas as well
				__movsq(g17chunk.vm1ygc, &tvm1guas[3*iey], 3 * nyc);
			//cout << "gochunk" << endl;
			g17chunk.GoChunk();
		}
	}
}

void G17CHUNK::GoChunk(){// elementary 'X' x 'Y' ychunk is 256x256
	// this is the critical code in "part 1"
	p_cpt1[6]++;
	uint16_t * limit_store;
	{//=========Find all XY passing first filter for the chunk
		register uint16_t * pstore = tstore;
		register int iy, store, ix, ipx = 0, Rvua = n64vua;
		register uint64_t *px = vxc;
		switch (Rvua){
		case 1:{
			for (ix = 0; ix < nxc; ix++, px ++){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++, py++){
					if (px[0] & py[0]) continue;
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		case 2:{
			for (ix = 0; ix < nxc; ix++, px = &px[2]){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++, py=&py[2]){
					if (px[0] & py[0]|| px[1] & py[1]) continue;
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		case 3:{
			for (ix = 0; ix < nxc; ix++, px = &px[3]){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++,  py=&py[3]){
					if (px[0] & py[0]|| px[1] & py[1]|| px[2] & py[2]) continue;
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		case 4:{
			for (ix = 0; ix < nxc; ix++, px = &px[4]){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++,  py=&py[4]){
					if (px[0] & py[0] || px[1] & py[1] || px[2] & py[2] ||
						px[3] & py[3]) continue;	
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		default:{
			for (ix = 0; ix < nxc; ix++, px = &px[Rvua]){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++,  py=&py[Rvua]){
					if (px[0] & py[0] || px[1] & py[1] || px[2] & py[2] ||
						px[3] & py[3]) continue;	
					for(int j=4;j<Rvua;j++)if (px[j] & py[j]) goto exitj;
					*pstore++ = store | iy;// (ix<<8) | iy
				exitj:{}
				}
			}
			break;
		}
		}// end switch

		limit_store = pstore;
	}
	{//=========== send such XY to next step in the process
		uint16_t * wstore = tstore;
		int ix, iy;
		while (wstore < limit_store){
			register int w = *wstore++;// catch the value
			ix = w >> 8;//8 bits  x
			iy = w & 0xff;// 8 bits y
			g17xy.xxe = txec[ix];
			g17xy.yye = tyec[iy];
			__movsq(g17xy.vxg, &vxgc[indexstep.n64vgua * ix], G17BLOCGSUA);
			__movsq(g17xy.vyg, &vygc[indexstep.n64vgua * iy], G17BLOCGSUA);
			__movsq(g17xy.vm1xg, &vm1xgc[3 * ix], 3);
			__movsq(g17xy.vm1yg, &vm1ygc[3 * iy], 3);
			g17xy.Go_0();
			//p_cpt1[7]++;
		}
	}
}



void G17XY::Go_0(){// check stacks, check more uas, collect Guas status
	
	if (g17b.debug17){
		register uint64_t R = indexstep.searched17_54,
			R1 = xxe.cellsbf | yye.cellsbf;
		if (GTEST17_ON){ if (R1 != R) return; }
		if (R1 == R)cout << Char54out(R1) << "Rxy" << endl;
		if (a_17_found_here == 1) return;
	}

	{ // compute and check stack count limit 6
		register uint64_t R = xxe.stacks.u64 + yye.stacks.u64;
		stacksxy.u64 = R;
		if ((R & 0xff) > 6) return;
		R >>= 16; if ((R & 0xf) > 6) return;
		R >>= 16; if ((R & 0xf) > 6) return;
	}
	{	// check here more uas
		register uint64_t Rxy = xxe.cellsbf;
		Rxy |= yye.cellsbf;
		// reshape Rxy to be in 2x mode
		Rxy = (Rxy & 0xffffffff) | ((Rxy >> 27) << 32);
		if (g17more.Check(Rxy))return;
		if (g17morebig.Check(Rxy))return;
	}
	Go_Guas_collect();// collect guas still active
	Go_Guas_collectMore();// collect guas still active
	BuildActiveBands();// apply pairs triplets to bands 

	if (!ntb3) return; // no active band b3lim and stack filters
	// at least one band is ok, time to check validity of the band 1+2
	p_cpt[3]++;
	if (CheckValidBand12()) return;
//	if (g17b.debug17 > 1)		cout << " valide b 12" << endl;

	Go_ValideB12();
}
void G17XY::Go_Guas_collect(){// collect guas still active
	register int imore = indexstep.added_more;
	uint64_t vw1[3];
	uint32_t iua;
	// setup XY guas status
	for (int i = 0; i < indexstep.n64vgua; i++)		vw[i] = vxg[i] & vyg[i];
	if (imore)for (int i = 0; i < 3; i++)		vw1[i] = vm1xg[i] & vm1yg[i];
	bands_active_pairs.SetAll_0();
	bands_active_triplets.SetAll_0();
	more_active_pairs.SetAll_0();
	more_active_triplets.SetAll_0();
	more3 = 0;
	for (int i = 0; i < indexstep.n64vgua; i++){//scan used  blocks 
		register uint64_t x = vw[i];
		//x &= indexstep.vaguas[i]; //already done
		while (bitscanforward64(iua, x)){
			register uint64_t uaid = (indexstep.tgua[iua + 64 * i]) >> 56;
			uint64_t *tv = indexstep.vid81[uaid];// clear all same id iv vectors
			x &= tv[i];// clear bit and other bits same id
			vw[i + 1] &= tv[i + 1];// clear id in next block 
			if (imore){
				uint64_t *tmv = indexstep.vmid81[uaid];
				for (int j = 0; j < 3; j++)vw1[j] &= tmv[j];
			}
			if (uaid <81)bands_active_pairs.setBit(C_To128[uaid]);
			else bands_active_triplets.setBit(C_To128[uaid - 81]);
		}
	}
	if (!imore) return;	//  look for missing ids add them in "or" mode 
	uint64_t * tx = vw1; 
	G17TMORE * tm = g17moreg0;		
	for (int i = 0; i < 3; i++){//scan used  blocks 
		register uint64_t x = tx[i];
		//x &= indexstep.vaguas[i];
		while (bitscanforward64(iua, x)){
			register uint64_t uaid = (tm[i].t[iua]) >> 56;
			uint64_t *tv = indexstep.vmid81[uaid];// clear all same id iv vectors
			x &= tv[i];// clear bit and other bits same id
			if (uaid <81)bands_active_pairs.setBit(C_To128[uaid]);
			else bands_active_triplets.setBit(C_To128[uaid - 81]);
		}
	}
}

void G17XY::Go_Guas_collectMore(){// collect guas from more tables
	register uint64_t R = xxe.cellsbf;	R |= yye.cellsbf;
	for (int it = 0; it < 3; it++){
		register int nt = g17moreg[it].nt;
		register uint64_t * t = g17moreg[it].t;
		while (nt--){
			register  uint64_t v = *t;
			if (!(v&R)){
				v >>= 56; //now socket 0-161
				if (v <81)bands_active_pairs.setBit(C_To128[v]);
				else bands_active_triplets.setBit(C_To128[v - 81]);
			}
			t++;
		}
	}
}


uint64_t G17XY::GetPairs(int ib3){
	uint32_t iuaid;
	register uint64_t pairbfa = 0;
	register uint64_t * ti = &genb12.tipairs[ib3][0].u64;
	{ // process first box
		register uint32_t X = genb12.tbands_pairs[ib3].bf.u32[0];
		X &= bands_active_pairs.bf.u32[0];
		while (bitscanforward(iuaid, X)){
			X &= ~(1 << iuaid);
			pairbfa |= ti[iuaid];
		}
	}
		{ // process second box
			register uint32_t X = genb12.tbands_pairs[ib3].bf.u32[1];
			X &= bands_active_pairs.bf.u32[1];
			while (bitscanforward(iuaid, X)){
				X &= ~(1 << iuaid);
				pairbfa |= ti[iuaid + 32];
			}
		}
		{ // process last box
			register uint32_t X = genb12.tbands_pairs[ib3].bf.u32[2];
			X &= bands_active_pairs.bf.u32[2];
			while (bitscanforward(iuaid, X)){
				X &= ~(1 << iuaid);
				pairbfa |= ti[iuaid + 64];
			}
		}
		return pairbfa;
}
int G17XY::GetTriplets(int ib3){
	register int minirows_triplets = 0;
	uint32_t iuaid;
	register int * ti = genb12.tindextriplets[ib3];
	{ // process first box
		register uint32_t X = genb12.tbands_triplets[ib3].bf.u32[0];
		X &= bands_active_triplets.bf.u32[0];
		while (bitscanforward(iuaid, X)){
			X &= ~(1 << iuaid);
			minirows_triplets |= ti[iuaid];
		}
	}
	{ // process second box
		register uint32_t X = genb12.tbands_triplets[ib3].bf.u32[1];
		X &= bands_active_triplets.bf.u32[1];
		while (bitscanforward(iuaid, X)){
			X &= ~(1 << iuaid);
			minirows_triplets |= ti[iuaid + 32];
		}
	}
	{ // process third box
		register uint32_t X = genb12.tbands_triplets[ib3].bf.u32[2];
		X &= bands_active_triplets.bf.u32[2];
		while (bitscanforward(iuaid, X)){
			X &= ~(1 << iuaid);
			minirows_triplets |= ti[iuaid + 64];
		}
	}
	return minirows_triplets;
}
void G17XY::BuildActiveBands(){
	ntb3 = 0;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++){
		GINT64 & ccx = wg3.count, &ccxs = wg3.countsum;
		memset(&wg3, 0, sizeof wg3);
		wg3.ib3 = ib3;
		register uint64_t pairbfa = GetPairs(ib3);
		wg3.pairs.u64 = pairbfa;
		if (pairbfa){// parse the pair status 
			register int bit = 1, Ror = (int)pairbfa;
			for (int i = 0; i < 9; i++){// 9 minirows to explore
				register int n = __popcnt(Ror & 7);
				if (n){
					ccx.u16[0] |= bit;
					ccx.u16[n] |= bit;
				}
				bit <<= 1;
				Ror >>= 3;
			}
		}
		register int A = GetTriplets(ib3);
		{// add possible triplets in empty mini rows (not very common)
			int tbr[8] = { 0, 07, 070, 077, 0700, 0707, 0770, 0777 };
			A &= (~wg3.count.u16[0]);// keep only minirows with no pair
			wg3.minirows_triplets = A;
			if ((!pairbfa) && (!A)){
				wg3.countstack.u64 = stacksxy.u64;// already seen as no stack filter
				goto band_to_keep;
			}
			ccxs.u16[3] = __popcnt(ccx.u16[0]) + __popcnt(ccx.u16[3]) + __popcnt(A);
			if (ccxs.u16[3] > b3lim)				goto dont_keep;
			// fill active triplets in wg3.pairs.u32[1] 
			wg3.pairs.u32[1] |= tbr[A & 7];
			wg3.pairs.u32[1] |= (tbr[(A >> 3) & 7] << 9);
			wg3.pairs.u32[1] |= (tbr[(A >> 6) & 7] << 18);

		}


		{// set up minimum count per stack (pairs + 3 pairs+triplet ) in each mini row
			A |= wg3.count.u16[0];//now all used mini rows
			register int B = ccx.u16[3]; // mini rows with 3 pairs are counted as 2 minimum clues
			ccxs.u16[0] = __popcnt(A & 0111) + __popcnt(B & 0111);
			ccxs.u16[1] = __popcnt(A & 0222) + __popcnt(B & 0222);
			ccxs.u16[2] = __popcnt(A & 0444) + __popcnt(B & 0444);
			register uint64_t R = stacksxy.u64 + ccxs.u64;
			wg3.countstack.u64 = R;
			if ((R & 0xff) > 6) goto dont_keep;
			R >>= 16; if ((R & 0xf) > 6) goto dont_keep;
			R >>= 16; if ((R & 0xf) > 6) goto dont_keep;
		}
	band_to_keep:
		g17tb3go[ntb3++] = wg3;
		continue;
	dont_keep:{}// forced continue in any place
	}
}
int G17XY::CheckValidBand12(){
	// passing more filter now check validity or add a "more" ua
	zh2b_i1.Init_x_(xxe.cells, nbx);
	register uint64_t myua = zh2b[0].Init_y_(yye.cells, nby);
	if (myua){
		if (_popcnt64(myua) <= 20)		g17more.Add(myua);
		else g17morebig.Add(myua);
		return 1; // not a solution unique for bands 1+2
	}
	return 0;
}
void G17XY::Go_ValideB12(){// UA2 and UA3 known not yet dead with min clues in band 3
	p_cpt[4]++;
	int tcluesb3[20], ntb3a = 0;
	for (int i = 0; i < nbx; i++)tcluesb3[ntb3a++] = xxe.cells.u8[i];
	for (int i = 0; i < nby; i++)tcluesb3[ntb3a++] = yye.cells.u8[i];
	if (zhou[0].PartialInitSearch17(tcluesb3, ntb3a))return;// would be  bug 
	for (int i3 = 0; i3 < ntb3; i3++){
		wg3 = g17tb3go[i3];
		__movsd((unsigned long *)&genb12.grid0[54], 
			(unsigned long*)genb12.bands3[wg3.ib3].band0, 27);
		zhou[0].InitBand3PerDigit(genb12.bands3[wg3.ib3].band0);
		if (more3 && Rebuild())continue;
		g17hh0.Init(wg3);

		//if (g17b.debug17 > 1)		cout << " valide b12 go"  << endl;

		//============= collect Gua46 and uasb3 for the band split them "in-field" "out-field"
		nuasb3_1 = nuasb3_2 = 0;
		// first GUA46 usually shorter than UAs band3
		register int  Rfilt = g17tb3go[i3].pairs.u32[1];
		int * tpat = genb12.tbands_UA4_6s_pat[wg3.ib3];
		for (int i = 0; i < 3; i++){
			register unsigned long dres = 27 * i;
			register uint32_t res;
			register int x = bands_active_pairs.bf.u32[i];				
			x |= more_active_pairs.bf.u32[i];
			x &= genb12.tbands_UA4_6s[wg3.ib3].bf.u32[i];
			while (bitscanforward(res, x)){
				x ^= 1 << res;// clear bit
				res += dres;
				register int Ru = tpat[res];
				if (Ru & Rfilt)	uasb3_1[nuasb3_1++] = Ru;
				else		uasb3_2[nuasb3_2++] = Ru;
			}
		}
		uint32_t * to = genb12.bands3[wg3.ib3].tua;
		for (uint32_t i= 0; i < genb12.bands3[wg3.ib3].nua; i++){
			register int Ru = to[i];
			if (Ru & Rfilt)			uasb3_1[nuasb3_1++] = Ru;
			else uasb3_2[nuasb3_2++] = Ru;
		}
		g17hh0.Go();
	}
}

int G17XY::Rebuild(){
	int rib3 = wg3.ib3;// reset 
	GINT64 & ccx = wg3.count , &ccxs = wg3.countsum;
	memset(&wg3, 0, sizeof wg3);
	wg3.ib3 = rib3;
	uint32_t iuaid;
	register uint64_t * ti = &genb12.tipairs[wg3.ib3][0].u64;
	register uint64_t pairbfa = 0;
	for (int ibox = 0; ibox < 3; ibox++){
		register uint32_t X = bands_active_pairs.bf.u32[ibox];
		X |= more_active_pairs.bf.u32[ibox];
		X &= genb12.tbands_pairs[wg3.ib3].bf.u32[ibox];
		while (bitscanforward(iuaid, X)){
			X ^= 1 << iuaid;
			iuaid += 32 * ibox;
			pairbfa |= ti[iuaid];
		}
	}
	wg3.pairs.u64 |= pairbfa;// add new pairs to olds 
	if (pairbfa){// parse the more pair status
		register int bit = 1, Ror = (int)pairbfa;
		for (int i = 0; i < 9; i++, bit <<= 1, Ror >>= 3){// 9 minirows to explore
			register int n = __popcnt(Ror & 7);
			if (n){ ccx.u16[0] |= bit;	ccx.u16[n] |= bit; }
		}
	}
	int tbr[8] = { 0, 07, 070, 077, 0700, 0707, 0770, 0777 };
	register int A = wg3.minirows_triplets;// collect new triplets + active
	// kill olds mini rows triplets in wg3.pairs.u32[1]
	wg3.pairs.u32[1] &= ~tbr[A & 7];
	wg3.pairs.u32[1] &= ~(tbr[(A >> 3) & 7] << 9);
	wg3.pairs.u32[1] &= ~(tbr[(A >> 6) & 7] << 18);

	register int * ti3 = genb12.tindextriplets[wg3.ib3];
	for (int ibox = 0; ibox < 3; ibox++){
		register uint32_t  X = bands_active_triplets.bf.u32[ibox];
		X |= more_active_triplets.bf.u32[ibox];
		X &= genb12.tbands_triplets[wg3.ib3].bf.u32[ibox];
		while (bitscanforward(iuaid, X)){
			X ^= 1 << iuaid;
			iuaid += 32 * ibox;
			A |= ti3[iuaid];
		}
	}
	{// add possible triplets in empty mini rows (not very common)
		A &= (~wg3.count.u16[0]);// keep only minirows with no pair
		wg3.minirows_triplets = A;
		if ((!wg3.pairs.u64) && (!A)){
			wg3.countstack.u64 = stacksxy.u64;// already seen as no stack filter
			return 0; // continue process
		}
		ccxs.u16[3] = __popcnt(ccx.u16[0]) + __popcnt(ccx.u16[3]) + __popcnt(A);
		if (ccxs.u16[3] >b3lim) return 1; // cancelled
		// fill  new active triplets in wg3.pairs.u32[1]
		wg3.pairs.u32[1] |= tbr[A & 7];
		wg3.pairs.u32[1] |= (tbr[(A >> 3) & 7] << 9);
		wg3.pairs.u32[1] |= (tbr[(A >> 6) & 7] << 18);
	}

	{// reset  minimum count per stack (pairs + 3 pairs+triplet ) in each mini row
		A |= wg3.count.u16[0];//now all used mini rows
		register int B = ccx.u16[3]; // mini rows with 3 pairs are counted as 2 minimum clues
		ccxs.u16[0] = __popcnt(A & 0111) + __popcnt(B & 0111);
		ccxs.u16[1] = __popcnt(A & 0222) + __popcnt(B & 0222);
		ccxs.u16[2] = __popcnt(A & 0444) + __popcnt(B & 0444);
		register uint64_t R = stacksxy.u64 + ccxs.u64;
		wg3.countstack.u64 = R;
		if ((R & 0xff) > 6) return 1;
		R >>= 16; if ((R & 0xf) > 6) return 1;
		R >>= 16; if ((R & 0xf) > 6) return 1;
	}
	return 0;
}

void G17XY::FoutValid17(int bf3, int ib3){
	char zs[82];
	int *g = genb12.grid0;
	strcpy(zs, empty_puzzle);
	register uint64_t w = xxe.cellsbf, bit = 1;
	w |= yye.cellsbf;
	for (int i = 0; i < 54; i++, bit <<= 1)if (w&bit)
		zs[i] = g[i] + '1';
	bit = 1;
	g = genb12.bands3[ib3].band0;
	for (int i = 0; i < 27; i++, bit <<= 1)if (bf3&bit)
		zs[i + 54] = g[i] + '1';
	fout1 << zs << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
}

//================ part 2  band procvessing
void G17B3HANDLER::Init(G17TB3GO & wgo){
	wg3 = wgo;
	nmiss = g17xy.b3lim - wg3.countsum.u16[3];
	known_b3 = 0;
	ncritical = wg3.countsum.u16[3];// to update assigning critical 2 pairs
	ndead = BIT_SET_27;
}
int G17B3HANDLER::IsMultiple(int bf){
	//if (g17b.debug17 > 1) cout << "is multiple ntuab3=" << g17xy.ntuab3 << endl;
	// check first if all tuab3 is hit
	for (int i = 0; i < g17xy.ntuab3; i++)
		if (!(bf&g17xy.tuab3[i])) return 1;// previous UA not hit
	int rua = zhou[1].CallMultipleB3(zhou[0], bf);
	//if (g17b.debug17 > 1) cout<<Char27out(rua) << "===========rua  is multiple all hit"  << endl;
	if (rua){
		if (g17xy.ntuab3 < 40)g17xy.tuab3[g17xy.ntuab3++] = rua;
		int nrua = __popcnt(rua);
		//cout << Char2Xout(zhou[0].glb.b12nextua) << " ";
		if (nrua < 4){
			int socket = genb12.GetSocket(rua, wg3.ib3),
				nrua = __popcnt(rua);
			// pack rua 54 bits
			GINT64 w;
			{
				uint64_t R = zh_g.b12nextua, R1 = R & 0xffffffff;
				R >>= 32; R <<= 27; R1 |= R;
				w.u64 = R1;
			}
			if (nrua == 2){
				g17xy.more_active_pairs.Set_c(socket);
				g17xy.more3 |= 1;
				w.u8[7] = socket;
				//cout << " socket b3_2=" << socket << endl;
			}
			else {
				g17xy.more_active_triplets.Set_c(socket);
				g17xy.more3 |= 2;
				w.u8[7] = socket + 81;
				//cout << " socket b3_3=" << socket+81 << endl;
			}
			uint64_t nb12 = _popcnt64(zh_g.b12nextua);
			if (nb12 <= 14)g17moreg[0].Add_If_New(w.u64);
			else if (nb12 <= 18)g17moreg[1].Add_If_New(w.u64);
			else g17moreg[2].Add_If_New(w.u64);
		}
		//else 		cout << Char27out(rua) << "ua de sortie"<< endl;
	}

	return rua;
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
		int shrink = TblShrinkMask[wg3.count.u16[2] & (0111 << ist)];
		if (shrink){// minirows 2 pairs in that stack
			register int Mask = tbitsrows[shrink] << (3 * ist);
			//adjust count and known
			known_b3 |= Mask & (~wg3.pairs.u32[0]);// and set the common cell as assigned
			wg3.count.u16[2] &= (~(0111 << ist)); // clear the 2pairs bit(s) in stack
			wg3.pairs.u32[1] &= (~Mask);// clear the pairs in the bit field
			ncritical -= __popcnt(shrink);
		}
	}
	wactive0 = Ractive_out&ndead;
}
void  G17B3HANDLER::Go(){
	if (0 &&g17b.debug17 > 1){
		cout << " call handler go" << endl;
		Print_Ok3();
	}
	p_cpt[7]++;
	rknown_b3 = 0;
	g17xy.ntuab3 = 0;
	if (nmiss){
		p_cpt[20 + nmiss]++;
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
		if (g17xy.nuasb3_2)return;// some uas outfield 
		p_cpt[8]++;
		Go_Critical();
	}
}
//================= critical process
void G17B3HANDLER::Critical2pairs(){// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (wg3.count.u16[2]){// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++){
			int shrink = TblShrinkMask[wg3.count.u16[2] & (0111 << ist)];
			if (shrink){// minirows 2 pairs in that stack
				register int Mask = tbitsrows[shrink] << (3 * ist);
				active_b3 &= (~Mask); // clear the minirow
				known_b3 |= Mask & (~wg3.pairs.u32[0]);// and set the common cell as assigned
			}
		}
	}
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
	if (bit & wg3.count.u16[3]){// the cell is in a minirow with 3 pairs active
		active_b3 &= ~Ru; //clear the cell
		wg3.count.u16[3] ^= bit; // now only a pairto hit
		wg3.count.u16[1] |= bit;
	}
	else{// either one pair or a triplet in the minirow
		//active_b3 &= ~Ru; //clear the cell
		active_b3 &= (~Mask); // kill the minirow as active
		wg3.count.u16[1] &= ~bit;
		wg3.minirows_triplets &= ~bit;
	}
}

void G17B3HANDLER::Go_Critical(){// critical situation all clues in pairs tripl:ets
	//if (g17b.debug17 > 1 && known_b3)cout << Char27out(known_b3) << " entry critical" << endl;
	active_b3 = wg3.pairs.u32[1];
	rknown_b3 = known_b3| active_b3;
	if(IsMultiple(rknown_b3))	return;	
	Critical2pairs();// assign 2 pairs in minirow to common cell
	if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
	int wknown = known_b3 | active_b3;
	if (rknown_b3 != wknown)
	if (IsMultiple(wknown))	return;// not valid using all cells
	rknown_b3 = wknown;
	if (irloop)		CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalLoop(int *to, int no){
	if (ShrinkUas1(to, no)) return;
	if (irloop)CriticalLoop(uasb3, nuasb3);
	else CriticalExitLoop(uasb3, nuasb3);
}
void G17B3HANDLER::CriticalExitLoop(int *uasb3, int nuasb3){
	int nmissb = g17xy.b3lim - __popcnt(known_b3);// missing clues
	if (0 &&g17b.debug17 > 1){
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
	if (rknown_b3 != wknown)
	if (IsMultiple(wknown))	return;// not valid using all cells
	rknown_b3 = wknown;

	if (nuasb3){		// find the smallest ua and apply it
		// use in priority an unsolved pair it is a smallest
		int wua = 0, sizeua = 27;
		unsigned long cell;
		if (wg3.count.u16[1]){
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[1]);
			int  shift = 3 * mini, mask = 7 << shift;
			wua = active_b3&mask;// catch the minirow
		}
		else{
			for (int i = 0; i < nuasb3; i++){
				register int ua = uasb3[i], cc = __popcnt(ua);
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
	int nmissb = g17xy.b3lim - __popcnt(known_b3);// missing clues
	if (0 &&g17b.debug17 > 1){
		cout << Char27out(known_b3) << "known 0 ua" << endl;
		cout << Char27out(active_b3) << "active 0 ua " << endl;
	}
	if (nmissb < 0)return;
	if (!nmissb){// nothing more to assign (granted at first call in a branch)
		CriticalFinalCheck();
		return;
	}
	if (wg3.count.u16[3])	{// in active minirows with 3 pairs, assign 2
		while (wg3.count.u16[3]){
			uint32_t mini;
			bitscanforward(mini, wg3.count.u16[3]);
			int shift = 3 * mini, bit = 1 << shift;
			wg3.count.u16[3] ^= 1 << mini; //clear bit the mini row is always killed
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
	if (wg3.count.u16[1]){// active pair in minirow
		uint32_t mini;
		bitscanforward(mini, wg3.count.u16[1]);
		int  shift = 3 * mini, bit = 1 << shift, mask = 7 << shift;
		wg3.count.u16[1] ^= 1 << mini; //clear bit the mini row is always killed
		int x = active_b3&mask;// catch the minirow
		active_b3 &= ~mask;// clear the minirow
		wg3.count.u16[1] ^= 1 << mini;// and clear the minirow bit as active
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
		wg3.minirows_triplets ^= 1 << mini; //clear bit the mini row is always killed
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
	int ncl = __popcnt(known_b3);
	//if (g17b.debug17 > 1) cout << "final check test=" << (ncl < g17xy.b3lim) << endl;
	if (ncl < g17xy.b3lim) return; // should be seen earlier if possible
	register int ir = IsMultiple(known_b3);
	if (ir){
		//if (g17b.debug17 > 1)cout << Char27out(ir) << "final check retour faux" << endl;
		return;// mode debug pour voir
	}
	g17xy.FoutValid17(known_b3, wg3.ib3);
	a_17_found_here = 1;
}
//=============== sub critical process   missing(s)  in the critical area
void G17B3HANDLER::Go_SubcriticalMiniRow(){
	int c2[3] = { 3, 5, 6 };// 2 cells in a mini row
	G17B3HANDLER hn;
	int bit = 1<<ndead, mask = 7<<(3*ndead), stack = ndead%3;
	for (int i = ndead; i < 9; i++, stack++, bit <<= 1, mask <<= 3){
		if (stack >2) stack = 0;
		register int M = active_sub & mask;
		if (!M)continue;
		ndead = i;
		if (bit & wg3.count.u16[1])// it was a gua2 pair assign both 
			hn.SubMini(this, M, mask, stack, 1, bit);
		else if (bit & wg3.count.u16[2])// it was 2 gua2 pair assign 2 out of 3 
			for (int j = 0; j < 3; j++){
				M = c2[j] << (3 * i);// 2cell bit field in the mini row
				hn.SubMini(this, M, mask, stack, 2, bit);
			}
		else if (bit & wg3.count.u16[3])// it was 3 gua2 pair assign 3 out of 3 
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
		if (imode<4)wg3.count.u16[imode] ^= bit;
		else wg3.minirows_triplets ^= bit;
	}
	known_b3 |= M;// assign 1 or 2
	nmiss--;// one added 
	active_b3 &= ~mask;
	active_sub ^= M;
	// now adjust the stack count
	wg3.countstack.u16[stack]++;
	wg3.countstack.u16[3]++;
	if (wg3.countstack.u16[stack] > 5)active_sub &= ~(07007007 << (3 * stack));
	if (nmiss) Go_SubcriticalMiniRow();// continue till a"no missing clue condition"
	else{	// leave sub critical mode and enter the critical mode 
		Critical2pairs();// assign 2 pairs in minirow to common cell
		if (ShrinkUas1(g17xy.uasb3_1, g17xy.nuasb3_1)) return;// dead branch
		rknown_b3 = known_b3 | active_b3;
		if (IsMultiple(rknown_b3))return;// not valid using all cells
		if (irloop)		CriticalLoop(uasb3, nuasb3);
		else CriticalExitLoop(uasb3, nuasb3);
	}
}
void G17B3HANDLER::Go_Subcritical(int docheck){// nmiss to select in the critical field
	p_cpt[11]++;
	active_b3 = active_sub = wg3.pairs.u32[1];
	// check first if a global solution  is still possible
	if(docheck)if (IsMultiple(known_b3 | active_b3))return;// not valid using all cells
	p_cpt[12]++;
	int cct = __popcnt(wg3.pairs.u32[1]) - ncritical;
	if (cct < nmiss)return;// not enough remaining cells in GUA2s GUA3s to reach the count
	for (int ist = 0; ist < 3; ist++){// check stacks 
		if (wg3.countstack.u16[ist]>5)active_sub &= ~(07007007 << (3 * ist));// kill the stack for more clues
	}
	ndead = 0;
	uasb3 = g17xy.uasb3_1;
	nuasb3 = g17xy.nuasb3_1;
	Go_SubcriticalMiniRow();// find the first miss
}
//======================================================================= not critical sequence


#define FINDWUA {register int Ra = wactive0, Rfilt =  known_b3;\
for (int iua = 0; iua < g17xy.nuasb3_2; iua++){\
register int Ru = g17xy.uasb3_2[iua];\
if (Ru & Rfilt) continue;\
Ru &= Ra;\
register int cc = __popcnt(Ru);\
if (cc < ncells)	{ ncells = cc; wua = Ru; rawua = g17xy.uasb3_2[iua]; }\
if (!cc)return;\
if (cc > 6 && ncells < 5) break;}}


#define APPLYWUA(G)unsigned long res;\
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
	APPLYWUA(hn.Go_Not_Critical_miss4())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 4 cells
}
void G17B3HANDLER::Go_Not_Critical_miss4(){ // select and assign a small ua outfield then miis2
	// or subfield 2 cells
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
	APPLYWUA(hn.Go_Not_Critical_miss3())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 3 cells
}
void G17B3HANDLER::Go_Not_Critical_miss3(){ // select and assign a small ua outfield then miis2
	// or subfield 2 cells
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
	APPLYWUA(hn.Go_Not_Critical_miss2())
		if (ncells == 27)	Go_Subcritical();// and finish in Subcritical mode here miss 3 cells
}
void G17B3HANDLER::Go_Not_Critical_miss2(){ // select and assign a small ua outfield then miis1
	// or subfield 2 cells
	if (IsMultiple(wg3.pairs.u32[1] | wactive0 | known_b3))return;
	int wua = wactive0, ncells = 27, rawua = 0;
	FINDWUA
	//if (g17b.debug17 > 1)	cout<<Char27out(wua) << " wua nmiss2 unique" << endl;
	APPLYWUA(hn.Go_Not_Critical_miss1())
		if (ncells == 27)Go_Subcritical();// and finish in Subcritical mode here miss 2 cells
}
void G17B3HANDLER::Go_Not_Critical_miss1(){
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
	if (0 &&g17b.debug17 > 1)	
		cout << Char27out(known_b3) << " known b3 nmiss1 entry"<<endl 
		<< Char27out(wact_miss1) << " wact_miss1 nn= " <<nn<< endl;
	if (!nn){
		if (IsMultiple(wg3.pairs.u32[1] | known_b3)){// must be one cell out
			//wact_miss1 &= zhou[0].glb.go_back;
			if (!wact_miss1) return;
			goto one_out_forced;
		}
		G17B3HANDLER hn = *this;
		hn.Go_Subcritical(0);// try Subcritical mode in a new handler
		// must try each cell in wactiveas redundant
		goto wactive_to_test;
	}
	if (wact_miss1){
		if (IsMultiple(wg3.pairs.u32[1] | wact_miss1 | known_b3))	return;
	}
	//else return;
	one_out_forced:
	//if (g17b.debug17 > 1)cout << " one_out_forced:" << endl;
	p_cpt[10]++;
	int st = 07007007;// try to reduce the count using stacks
	for (int ist = 0; ist < 3; ist++){
		int wastack = wact_miss1 & st;
		if (wastack){
			int bfstack = wg3.pairs.u32[1] | wastack | known_b3;
			if (IsMultiple(bfstack)){
				wact_miss1 ^= wastack;
				//wact_miss1 &= zhou[0].glb.go_back;// use received ua to reduce the count
			}
		}
	}
	wactive_to_test:
	//if (g17b.debug17 > 1)cout << Char27out(wact_miss1) << " wactive_to_test:" << endl;
	p_cpt[18]++;
	uint32_t res;
	while (bitscanforward(res, wact_miss1)){
		int bit = 1 << res;
		wact_miss1 ^= bit;// clear bit
		G17B3HANDLER hn = *this;
		hn.known_b3 |= bit;
		int stack = C_box[res];
		hn.wg3.countstack.u16[stack]++;
		hn.wg3.countstack.u16[3]++;
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
	if (__popcnt(active)<10 && g17xy.g17hh0.IsMultiple(known | active)) return;
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
			register int cc = __popcnt(Ru);
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
char * g17tl1[10] = { "entry", "steps1", "steps2", "steps3", "steps4",
"common", "f1", "dead", "n51", "n62" };


char * g17tl[40] = { "0__", "falluas", "fstack1", "UasGuas", "fb12v", "  b3v\t\t", "fnov", //0-6
"h3go", "critb", "ncritb", "m1go", "subgo", "subvv ",   //12
"gcrit", "gcvv", "gc_exit", "m2", "m2vv", "m1_ss", "19_", "fin check",//20
"21_", "22_", "23_", "24_", "25_", "26__", "gocheck", "check0", "nocheck0",//29
"30nua", "nua2", "nua3", "nb3", "__34", "vt1", "vt2", "", "T1", "T2",//39
};

//=============================== debugging sequences
void G17B::PrintEndPuz() {
	cout << endl;
	for (int i = 0; i < 10; i++) if (p_cpt1[i])
		cout << g17tl1[i] << "\t" << p_cpt1[i] << endl;
	cout << endl;
	for (int i = 0; i < 32; i++) if (p_cpt[i])
		cout << g17tl[i] << "\t" << p_cpt[i] << endl;
}
int G17B::Godebug17_1(int i1) { // only if g17test_on 
	if (g17b.debug17) {// skip false band1
		int * t1 = myband1.index2[i1];
		int n51 = t1[4] - t1[1];
		int n61 = t1[5] - t1[2];
		XY_EXPAND xe1 = xye6[0][t1[2]];// here cell is 27-53
		int cell1 = xe1.cells.u8[0], cell2 = xe1.cells.u8[1];
		int cellsbf = (1 << cell1) | (1 << cell2);
		register uint64_t R = cellsbf;
		R &= (~indexstep.searched17_54);
		if (R) return 1;
		cout << Char27out(cellsbf) << " band 1 2cells to search" << endl;
	}
	return 0;
}
int G17B::Godebug17_2(int i2) {// only if g17test_on 
	if (g17b.debug17) {// skip false band1
		int * t1 = myband2.index2[i2];
		int n51 = t1[4] - t1[1];
		int n61 = t1[5] - t1[2];
		XY_EXPAND xe1 = xye6[1][t1[2]];// here cell is 27-53
		int cell1 = xe1.cells.u8[0], cell2 = xe1.cells.u8[1];
		uint64_t a = 1, b = 1;
		a <<= cell1; b <<= cell2;
		uint64_t cellsbf = a | b;
		register uint64_t R = cellsbf;
		R &= (~indexstep.searched17_54);
		//if (g17b.debug17>2)
		if (R) return 1;
		cout << Char54out(cellsbf) << " band 2 2cells  to search" << endl;
	}
	return 0;
}
void G17B::Godebug17_3() {// only if g17test_on 
	if (g17b.debug17) {
		cout << Char54out(indexstep.dead54) << " dead in active branch " << endl;
		cout << "indexstep.ntua=" << indexstep.ntua
			<< " indexstep.ntgua= " << indexstep.ntgua << endl;
	}
}
void G17B::Godebug17_4(int i1, int i2) {// only if g17test_on 
	cout << i1 << "\t" << i2 << "\t" << genuasb12.nua
		<< "\t" << genb12.ntua2 << "\t" << genb12.ntua3
		<< "\t" << p_cpt[4] << endl;

}

void G17B3HANDLER::Print_Ok3() {
	for (int i = 54; i < 81; i++) cout << genb12.grid0[i] + 1;
	cout << " band 3 studied wg3.ib3=" << wg3.ib3 << endl;
	cout << Char27out(wg3.pairs.u32[0]) << " pair mode pairs" << endl;
	cout << Char27out(wg3.pairs.u32[1]) << " pair mode bit field" << endl;
	cout << "minirows status  all 1 pair 2 pairs  3 pairs  triplet b3lim=" << g17xy.b3lim << endl;
	for (int j = 0; j < 4; j++) cout << Char9out(wg3.count.u16[j]) << endl;
	cout << Char9out(wg3.minirows_triplets) << endl;
	for (int j = 0; j < 4; j++)cout << wg3.countstack.u16[j] << "\t";
	cout << " compte status full grid nmiss=" << nmiss << endl;
	for (int j = 0; j < 4; j++)cout << wg3.countsum.u16[j] << "\t";
	cout << " compte status bande 3 total=" << wg3.countsum.u16[3] << endl;
	cout << Char27out(known_b3) << " assigned  cells" << endl;

}
