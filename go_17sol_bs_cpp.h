

void G17B::GoM10(){// processing an entry 656 566 with the relevant table of ba,ds3
	zh2b_g.test = 0;
	p_cpt2g[0] ++;
	if (p_cpt2g[0]>1 )	return; //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<first test 
	myband2.MorphUas();
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	int * zs0 = genb12.grid0;
	myband2.DoExpandBand(27);// expand band2
	if (!(myband1.n5| myband2.n5)) return; // no 656 no 566
	int nb3 = genb12.nband3;
	//=========================== collect UAs (still room for improvement)
	genuasb12.Initgen();
	genb12.BuildGang9x3();
	genb12.SecondSockets2Setup();// collect GUA2s 
	genb12.SecondSockets3Setup();// collect GUA3s 
	Go();// standard entry point for all 
}

void G17B::Go(){// search process external loop 2X2Y
	// temporary code to test the external loop
	//g17more.Init();
	//g17morebig.Init();
	//for (int i = 0; i < 3; i++)g17moreg[i].Init();
	//zh_g.grid0 = genb12.grid0;
	zh_g.grid0= genb12.grid0;
	indexstep.StartDead();
	int n1 = myband1.nind[1],
		n2 = myband2.nind[1];
	if (1) {
		//myband1.PrintShortStatus(); myband2.PrintShortStatus();
		cout << "n bands3      \t" << genb12.nband3 << endl;
		cout << "ua bands1+2   \t" << genuasb12.nua << endl;
		cout << "guas socket2  \t" << genb12.ntua2 << endl;
		cout << "guas socket3  \t" << genb12.ntua3 << endl;
		cout << "active socket2\t" << genb12.nactive2 << endl;
		cout << "active socket3\t" << genb12.nactive3 << endl;
		cout << "n1=" << n1 << " n2=" << n2 << " could be dynamic<<<<<<<<" << endl;
	}
	for (int i1 = 0; i1 < n1; i1++) {
		if (0) {
			cout << "i1=" << i1 << "\tn51=" << indexstep.n51 
				<< " n61=" << indexstep.n61		<< endl;
		}
		//if (i1) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		//if (Godebug17_1(i1)) continue;
		for (int i2 = 0; i2 < n2; i2++) {
			//if (i2>6) return;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			//if (Godebug17_2(i2)) continue;
			indexstep.added_more = 0;
			indexstep.Init(i1, i2);
			if (!(indexstep.n51 | indexstep.n52)) continue;
			if (indexstep.ShrinkUas())continue;// dead branch
			if (0) {
				cout << " i2=" << i2<< " n52=" << indexstep.n52	 << " n62=" << indexstep.n62  
					<< "\tnua step" << indexstep.ntua <<"\tnua"<< genuasb12.nua
					<< "\tgnua step" << indexstep.ntgua << "\tngua" << indexstep.ntgua_raw
					<<"\tp_cpt2g[3]="<< p_cpt2g[3] << endl;
			}
			indexstep.ShrinkGuas();
			g17more.Init();
			g17morebig.Init();
			genb12.Init_tmore_sockets2();
			if (indexstep.n52) {
				indexstep.Do65();
				//if (genuasb12.nua < TUA64_12SIZE && indexstep.n64vua < G17BLOCSUA - 2) {
					//indexstep.GoAddNewUas_sub_step();
				//}
				//indexstep.GoAddNewGUas_sub_step();
			}
			if (indexstep.n51)				indexstep.Do56();
			//indexstep.Do66(); // keep it for a 18 222222222 search
		}
	}

/* 
	for (int i1 = 0; i1 < n1; i1++){
		for (int i2 = 0; i2 < n2; i2++){
			g17more.Init();// restart with empty more tables for the next step
			g17morebig.Init();
			for (int i = 0; i < 3; i++)g17moreg[i].Init();
		}
		if (GTEST17_ON || g17b.debug17>1)
			Godebug17_4(i1, n2);// i1 i1 nua nua2 nua3 p_cpt[4]
	}
*/
}
void G17INDEXSTEP::Init(int i1, int i2) {
	cur_step_buffer_index = 0;
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
	b1b2_2Xbf = ((uint64_t)t2[0] << 32) | t1[0];
}
int G17INDEXSTEP::ShrinkUas() {// shrink table of UAs for bands1+2
	// could erase new supersets to consider doing it in 2 steps
	uint64_t * otua = genuasb12.tua;
	int  ontua = genuasb12.nua;
	register uint64_t Rw = b1b2_2Xbf, Rn = ~dead;
	ntua = 0;
	for (int iua = 0; iua < ontua; iua++) {
		register uint64_t Ru = otua[iua];
		if (Rw & Ru) continue;
		Ru &= Rn;// erase dead cells
		if (Ru)		tua[ntua++] = Ru;
		else return 1;// dead branch
	}
	if (ntua > 256) ntua = 256; //limit in the chunk loop
	n64vua = (ntua + 63) >> 6;
	vcellsuas_size = 54 * n64vua;
	vcellsuas = AllocLockStepBuffer(vcellsuas_size);
	g17chunk.n64vua = n64vua;
	{ // initial vectors to appropriate value and actives vector
		memset(vcellsuas, 255, 8*vcellsuas_size);
		register uint64_t * Ra = vauas;
		register int Rn = ntua;
		for (int i = 0; i < G17BLOCSUA; i++, Rn -= 64) {
			if (Rn > 64)vauas[i] = BIT_SET_64;
			else if (Rn > 0)vauas[i] = maskLSB[Rn].u64[0];
			else vauas[i] = 0;
		}
	}
	uint32_t cc;
	for (int i = 0; i < ntua; i++) {// set uas
		register int bloc = i >> 6, ir = i & 63;
		register uint64_t Rw = tua[i], biti = (uint64_t)1 << ir;
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			Rw ^= bit2;// clear bit
			vcellsuas[From_128_To_81[cc] * n64vua + bloc] ^= biti;
		}
	}
	return 0;
}
int G17INDEXSTEP::ShrinkGuasLoad(uint64_t *tu,  int itu1) {
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
			}
		}
	}
	return n;
}
void G17INDEXSTEP::ShrinkGuas() {// shrink table for gangster uas2 ua3s
// group all gua2 gua2 soskets in one table
	int tactive2_start[81], tactive2_end[81],
		tactive3_start[81], tactive3_end[81];
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
	if (ntgua > 64 * G17BLOCGSUA) ntgua = 64 * G17BLOCGSUA;
	//working with a limited number of GUAs for an index chunk
	n64vgua = (ntgua + 63) >> 6;
	g17chunk.n64vgua = n64vgua;

	vcellsguas_size = 54 * n64vgua;
	vcellsguas = AllocLockStepBuffer(vcellsguas_size);
	{ // initial vectors to appropriate value and actives vector
		memset(vcellsguas, 255, 8 *  vcellsguas_size);
		register uint64_t * Ra = vaguas;
		register int Rn = ntgua;
		for (int i = 0; i < G17BLOCGSUA; i++, Rn -= 64) {
			if (Rn > 64)vaguas[i] = BIT_SET_64;
			else if (Rn > 0)vaguas[i] = maskLSB[Rn].u64[0];
			else vaguas[i] = 0;
		}
	}
	uint32_t cc;
	for (int i = 0; i < ntgua; i++) {// set uas
		register int bloc = i >> 6, ir = i & 63;
		register uint64_t Rw = tgua[i], biti = (uint64_t)1 << ir;
		Rw &= BIT_SET_2X;
		while (bitscanforward64(cc, Rw)) {// look for  possible cells
			register uint64_t bit2 = (uint64_t)1 << cc;
			Rw ^= bit2;// clear bit
			vcellsguas[From_128_To_81[cc] * n64vgua + bloc] ^= biti;
		}
	}

	// ____________________________________create id81 vectors of corresponding guas
	vid81_size = n64vgua * (nactive2 + nactive3);
	vid81 = AllocLockStepBuffer(vid81_size);
	memset(vid81, 0, 8 * vid81_size);
	uint64_t * vx = vid81;
	for (int i = 0; i < nactive2; i++) {
		SetV(vx, tactive2_start[i], tactive2_end[i]);
		vx += n64vgua;
	}
	for (int i = 0; i < nactive3; i++) {
		SetV(vx, tactive3_start[i], tactive3_end[i]);
		vx += n64vgua;
	}
}
void G17INDEXSTEP::SetV(uint64_t * v, int i1, int i2) {
	while (1) {// *v initial value set to 0 outside the call
		if (i1 >= 64) {
			v++;
			i1 -= 64;
			i2 -= 64;
			continue;
		}
		*v = ~0; // active bloc initial value to all '1'
		if (i1) {// clear the corresponding 
			register uint64_t R = (uint64_t)1 << i1;
			R -= 1;
			*v ^= R; // bits cleared
			i1 = 0;
		}
		if (i2 > 64) {
			v++;
			i2 -= 64;
		}
		else{// last active bloc, clear corresponding bits
			if (i2 == 64)return;
			register uint64_t R = (uint64_t)1 << i2;
			R -= 1;
			*v &= R; // bits to keep
			return;
		}
	}
}
void G17INDEXSTEP::GoAddNewUas_sub_step(){// partial step build in vector 
	/*
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
			vcellsuas[From_128_To_81[cc]*n64vua+bloc] ^= biti;
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
			vcellsuas[From_128_To_81[cc]*n64vua+bloc] ^= biti;
		}
	}
	g17morebig.Init();
	n64vua++;
	g17chunk.n64vua = n64vua;
	*/
}
void G17INDEXSTEP::GoAddNewGUas_sub_step(){// end step include small more guas in the guas tables
	// copy mores in the appropriate table
	/*
	int imore = added_more++;
	G17TMORE * tm = g17moreg0;
	memcpy(tm, g17moreg, sizeof g17moreg);
	// store small in the main tables
	register int nt = tm[0].nt;
	register uint64_t * t = tm[0].t;
	while (nt--){
		register  uint64_t v = *t, vid = v >> 56;
		if (vid>80){
			//if (genb12.ntua3 < 3000) g17b.InsertGua(genb12.tua3, genb12.ntua3, v);
		}
		//else if (genb12.ntua2 < 3000)	g17b.InsertGua(genb12.tua2, genb12.ntua2, v);
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
	*/
}

void G17INDEXSTEP::Do65(){
	{//=================== move 6 b1 fix place w1
		g17xy.nby = 6;
		nw1 = n61;
		uint64_t *Ro = (uint64_t *)(&myband1.xye6[t1[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof  xyew1[0] * nw1);

	}
	{//=================== move 5 b2 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n52;
		uint64_t *Ro = (uint64_t *)(&myband2.xye5[t2[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] *nw2);
	}

	g17xy.b3lim = 6;
	Do_Common();
}
void G17INDEXSTEP::Do56(){
	{// ===================== move 6 b2 fix place w1
		g17xy.nby = 6;
		nw1 = n62;
		uint64_t *Ro = (uint64_t *)(&myband2.xye6[t2[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof  xyew1[0] *nw1);
	}
	{// ===================== move 5 b1 fix place w2
		g17xy.nbx = 5;
		ncluesw2 = 5;
		nw2 = n51;
		uint64_t *Ro = (uint64_t *)(&myband1.xye5[t1[1]]), *Rd = (uint64_t *)(xyew2);
		memcpy(Rd, Ro, sizeof  xyew2[0] *nw2);
	}
	g17xy.b3lim = 6;
	Do_Common();
}
void G17INDEXSTEP::Do66(){// assumed after do 56 x is band 1 y is band2
	/*
	if (!indexstep.n51){// ====== move 6 b2 fix place w1 if not yet done
		g17xy.nby = 6;
		nw1 = n62;
		uint64_t *Ro = (uint64_t *)(&myband2.xye6[t2[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof  xyew1[0] *nw1);
	}

	{// ===================== move 6 b1 fix place w2
		g17xy.nbx = 6;
		ncluesw2 = 6;
		nw1 = n61;
		uint64_t *Ro = (uint64_t *)(&myband1.xye6[t1[2]]), *Rd = (uint64_t *)(xyew1);
		memcpy(Rd, Ro, sizeof xyew2[0] * nw1);
	}
	g17xy.b3lim = 5;
	Do_Common();
	*/
}


void G17INDEXSTEP::Do_Common(){// after xye w1('Y')  and w2('X') have been loaded
	if (!nw1)return;
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
			uint64_t * v =& vcellsuas[cell*n64vua], 
				*vg = &vcellsguas[cell*n64vgua];
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
			uint64_t * v = &vcellsuas[cell*n64vua], 
				*vg = &vcellsguas[cell*n64vgua];
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
			memcpy(Rd, Ro, sizeof  g17chunk.txec[0] *nxc);
		}
		Do_Common_3_BuildXvectors();// build 'X' tables in the chunk
		for (int iy = 0, iey = 0; iy <= nych; iy++, iey += G17CHKY){
			nyc = G17CHKY;
			if (iy == nych) nyc = nw1 - G17CHKY * iy;
			g17chunk.nyc = nyc;
			{// move Y_expand in the chunk from w1
				uint64_t *Ro = (uint64_t *)(&xyew1[iey]), *Rd = (uint64_t *)(g17chunk.tyec);
				memcpy(Rd, Ro, sizeof   g17chunk.tyec[0] * nyc);
			}
			// dynamic move buas and vguas  to the chunk
			g17chunk.vyc = &g17chunk.vxc[n64vua * nxc];
			memcpy(g17chunk.vyc, &tvuas[n64vua*iey], 8*n64vua * nyc);
			memcpy(g17chunk.vygc, &tvguas[n64vgua*iey], 8*n64vgua * nyc);
			
			if (added_more)// and more moreguas as well
				memcpy(g17chunk.vm1ygc, &tvm1guas[3*iey], 8*3 * nyc);
			g17chunk.GoChunk();
		}
	}
}
void G17INDEXSTEP::PrintUasShrinked() {
	cout << "uas shrinked table" << endl;
	for (int i = 0; i < ntua; i++)
		cout << Char2Xout(tua[i]) << endl;
}

void G17CHUNK::GoChunk(){// elementary 'X' x 'Y' ychunk is 256x256
	// this is the critical code in "part 1"
	p_cpt2g[1] ++;
	p_cpt2g[2] += nxc * nyc;
	uint16_t * limit_store;
	uint16_t tstore[256 * 256]; // storing XY passing filter 1 will never reach the limit
	{//=========Find all XY passing first filter for the chunk
		register uint16_t * pstore = tstore;
		register int iy, store, ix; 
		register uint64_t *px = vxc;
		switch (n64vua){
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
			for (ix = 0; ix < nxc; ix++, px += 2){
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
			for (ix = 0; ix < nxc; ix++, px +=3){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++,  py=&py[3]){
					if (px[0] & py[0]|| (px[1] & py[1])| (px[2] & py[2])) continue;
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		case 4:{
			for (ix = 0; ix < nxc; ix++, px +=4){
				store = ix << 8; //8 low bits for iy
				register uint64_t *py = vyc;
				for (iy = 0; iy < nyc; iy++,  py=&py[4]){
					if (px[0] & py[0] || 	(px[1] & py[1]) | 
						(px[2] & py[2]) |(px[3] & py[3])  ) continue;	
					*pstore++ = store | iy;// (ix<<8) | iy
				}
			}
			break;
		}
		}// end switch

		limit_store = pstore;
	}

	{//=========== send such XY to next step in the process
		uint16_t * wstore = tstore;
		while (wstore < limit_store){
			register int w = *wstore++;// catch the value
			{
				register int ix = w >> 8;//8 bits  x
				register int iy = w & 0xff;// 8 bits y
				g17xy.xxe = txec[ix];
				g17xy.yye = tyec[iy];
				memcpy(g17xy.vxg, &vxgc[indexstep.n64vgua * ix], sizeof g17xy.vxg);
				memcpy(g17xy.vyg, &vygc[indexstep.n64vgua * iy], sizeof g17xy.vyg);
				memcpy(g17xy.vm1xg, &vm1xgc[3 * ix], sizeof g17xy.vm1xg);
				memcpy(g17xy.vm1yg, &vm1ygc[3 * iy], sizeof g17xy.vm1yg);
			}
			g17xy.Go_0();
		}
	}
}

void G17XY::Go_0(){// check stacks, check more uas, collect Guas status
	p_cpt2g[3] ++;
	Init();
	/*
	if (g17b.debug17){
		register uint64_t R = indexstep.searched17_54,
			R1 = xxe.cellsbf | yye.cellsbf;
		if (GTEST17_ON){ if (R1 != R) return; }
		if (R1 == R)cout << Char54out(R1) << "Rxy" << endl;
		if (a_17_found_here == 1) return;
	}
	*/

	{ // compute and check stack count limit 6
		register uint64_t R = xxe.stacks.u64 + yye.stacks.u64;
		stacksxy.u64 = R;
		stacksxy.u16[3] = stacksxy.u16[0] + stacksxy.u16[1] + stacksxy.u16[2];
		if ((R & 0xff) > 6) return;
		R >>= 16; if ((R & 0xff) > 6) return;
		R >>= 16; if ((R & 0xff) > 6) return;
	}

	if (g17more.Check(cellsbf))return;		// check here more uas
	if (g17morebig.Check(cellsbf))return;
	Go_Guas_collect();						// collect guas still active
	if (FirstCheckActiveBands()) return;
	p_cpt2g[4]++;
	if (CheckValidBand12()) return;
	p_cpt2g[5]++;
	if (0) {
		cout << Char2Xout(cellsbf) << " entry Go_0 active nclues=" << nclues << " n=" << p_cptg[8]
			<< " stacks " << stacksxy.u16[0] << stacksxy.u16[1] << stacksxy.u16[2]
			<< "\t" << stacksxy.u16[3] << endl;
	}
	genb12.Sockets2SetupForB12(cellsbf);

	BuildActiveBands();// apply pairs triplets to bands 
	if (!ntb3) return; // no active band b3lim and stack filters
	// at least one band is ok, time to check validity of the band 1+2
	p_cpt2g[8] ++;
	p_cpt2g[9] += ntb3;
//	if (g17b.debug17 > 1)		cout << " valide b 12" << endl;
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
	more3 = 0;
	int n = indexstep.n64vgua;
	uint64_t vw[G17BLOCGSUA];
	for (int i = 0; i < n; i++) vw[i] = vxg[i] & vyg[i];
	// setup XY guas status
	for (int i = 0; i < indexstep.nactive2; i++) {
		register uint64_t * pw = &indexstep.vid81[n * i];
		register uint64_t vwi = (vw[0] & pw[0]);
		for (int j = 1; j < n; j++) vwi |= (vw[j] & pw[j]);
		int i81 = indexstep.tactive2[i];
		if (vwi)bands_active_pairs.Set_c(i81);
	}
	uint64_t * myvid = &indexstep.vid81[indexstep.nactive2 * n];
	for (int i = 0; i < indexstep.nactive3; i++) {
		register uint64_t * pw = &myvid[n * i];
		register uint64_t vwi = (vw[0] & pw[0]);
		for (int j = 1; j < n; j++) vwi |= (vw[j] & pw[j]);
		int i81 = indexstep.tactive3[i];
		if (vwi)bands_active_triplets.Set_c(i81);
	}
	more_checked_triplet= bands_active_triplets;
}
int G17XY::FirstCheckActiveBands() {// quick check of critical situations
	//if (genb12.nband3 > 10) return 0;// small chances to exit false<<<<<<<<<<<<<<<<<<<<<<
	// consider critical stack 
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
		if(! (ws2s | ws3s)) return 1;// one non critical band found
		if (maxcount == 6)continue; // band to exclude
		// need 2 minimum clues to exclude the band
		int * tm2 = &myb.ua2_imini[27 * stack],* tm3=&myb.ua3_imini[27 * stack];
		int tix[81], ntix = 0, mini_bf1=0,mini_bf2=0, mini_bf3 = 0;
		BitsInTable32(tix, ntix, ws2s);
		for (int i = 0; i < ntix; i++) {// switch from i_81 to mini rows
			int i81 = tix[i] + 27 * stack;
			int imini = tm2[i81],bit=1<<imini; 
			//cout << "bit N°" << tix[i] << "\timini\t" << imini << endl;
			if (mini_bf2&bit) mini_bf3 |= bit;
			if (mini_bf1&bit) mini_bf2 |= bit;
			mini_bf1 |= bit;
		}
		int ntix3 = 0,mini_triplet=0;
		BitsInTable32(tix, ntix3, ws3s);
		for (int i = 0; i < ntix3; i++) {// switch from i_81 to mini rows
			int i81 = tix[i] + 27 * stack;
			int imini = tm3[i81], bit = 1 << imini;
			//cout << "triplet bit N°" << tix[i] << "\timini\t" << imini << endl;
			mini_triplet |= bit;
		}
		uint32_t mincount = _popcnt32(mini_bf1) + _popcnt32(mini_bf3) + _popcnt32(mini_triplet);
		if (mincount > 2) return 1; // non critical band found
	}
	return 0;
}
int G17XY::CheckValidBand12(){
	// passing more filter now check validity or add a "more" ua
	register uint64_t myua = zh2b[0].ValidXY(tclues,nclues);
	if (myua){
		uint64_t cc = _popcnt64(myua);
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
	for (int ib3 = ibfirst; ib3 < genb12.nband3; ib3++) {
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

		uint32_t mincount = _popcnt32(mini_bf1) + _popcnt32(mini_bf3) + _popcnt32(mini_triplet);
		if (mincount > 6) continue; // too many clues needed for 6 5 6

		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs
		
		wg3.count.u16[0] = mini_bf1;
		wg3.count.u16[1] = mini_bf2;
		wg3.count.u16[2] = mini_bf3;
		
		wg3.countstack.u64 = stacksxy.u64;
		wg3.pairs.u32[0] = pairs27;
		// min count per stack
		wg3.countsum.u16[3] = mincount;
		register int s = 0111;
		wg3.countsum.u16[0] = _popcnt32(mini_bf1&s) + _popcnt32(mini_bf2 & s)
			+ _popcnt32(mini_bf3&s)		+ _popcnt32(mini_triplet&s);
		s = 0222;
		wg3.countsum.u16[1] = _popcnt32(mini_bf1 & s) + _popcnt32(mini_bf2 & s)
			+ _popcnt32(mini_bf3 & s) + _popcnt32(mini_triplet & s);
		s = 0444;
		wg3.countsum.u16[2] = _popcnt32(mini_bf1 & s) + _popcnt32(mini_bf2 & s)
			+ _popcnt32(mini_bf3 & s) + _popcnt32(mini_triplet & s);
		{ // setup sum and check stacks
			register uint64_t R = stacksxy.u64 + wg3.countsum.u64;
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
		
		// try to add triplets in empty minirows
		all_used_minis ^= 0777;// now free minis
		if (all_used_minis &&(p_cpt2g[8]<30)) {
			Find_More_Sockets3(all_used_minis);
		}

		g17tb3go[ntb3++] = wg3;
	}
}
int G17XY::Find_More_Sockets3(uint32_t free_minis) {
	int rmini = 0;
	cout << "test more sockets 3" << endl;
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
			cout << "call MoreSocket3 i81=" << i81 << endl;
			zh2b_g.test = 1;
			uint64_t new_ua = zh2b_g.MoreSocket2(genb12.gangb12, rgangcols,
				g17xy.tclues, g17xy.nclues, digs);
			zh2b_g.test = 0;
			cout << Char2Xout(new_ua) << " retour zhou" << endl;
		}
	}
	return rmini;
}

void G17XY::Go_ValideB12(){// UA2 and UA3 known not yet dead with min clues in band 3
	//cout <<Char2Xout(g17xy.cellsbf)<< "entry Go_ValideB12() ntb3=" << ntb3 << endl;
	int zinitdone = 0;
	for (int i3 = 0; i3 < ntb3; i3++){
		wg3 = g17tb3go[i3];
		int nmiss= 17 - wg3.countstack.u16[3];
		if(nmiss<2)p_cpt2g[10+nmiss] ++;
		else if (nmiss < 5)p_cpt2g[12] ++;
		else p_cpt2g[13] ++;
		if(wg3.minirows_triplets)p_cpt2g[14] ++;
		p_cpt2g[19] ++;
		if (p_cpt2g[19] > 10) continue;
		if (1) {
			wg3.Debug();

		}
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
		while ((i81 = socket4.getFirsCell()) >= 0) {
			socket4.Clear_c(i81);// clear bit
			continue;
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
		if (p_cpt2g[15]!=2) continue;
		if (!zinitdone) {
			zinitdone = 1;
			if (zhou[0].PartialInitSearch17(g17xy.tclues, g17xy.nclues))
				return;// would be  bug 
			//zhou[0].ImageCandidats();
		}
		memcpy(&genb12.grid0[54], genb12.bands3[wg3.ib3].band0, 4*27);
		zhou[0].InitBand3PerDigit(genb12.bands3[wg3.ib3].band0);

		//if (g17b.debug17 > 1)		cout << " valide b12 go"  << endl;
		g17hh0.Go();
	}
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

int G17B3HANDLER::IsMultiple(int bf){
	//if (g17b.debug17 > 1) 
	cout << "is multiple ntuab3=" << g17xy.ntuab3 << endl;

	// check first if all tuab3 is hit
	for (int i = 0; i < g17xy.ntuab3; i++)
		if (!(bf&g17xy.tuab3[i])) return 1;// previous UA not hit
	int rua = zhou[1].CallMultipleB3(zhou[0], bf);
	//if (g17b.debug17 > 1) 
	cout<<Char27out(rua) << "===========rua  is multiple all hit"  << endl;
	if (0 &&rua){
		/*
		if (g17xy.ntuab3 < 40)g17xy.tuab3[g17xy.ntuab3++] = rua;
		int nrua = _popcnt32(rua);
		//cout << Char2Xout(zhou[0].glb.b12nextua) << " ";
		if (nrua < 4){
			int socket = genb12.GetSocket(rua, wg3.ib3),
				nrua = _popcnt32(rua);
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
		*/
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
			ncritical -= _popcnt32(shrink);
		}
	}
	wactive0 = Ractive_out&ndead;
}
void  G17B3HANDLER::Go(){
	wg3 = g17xy.wg3;
	nmiss = 17 - wg3.countstack.u16[3];
	known_b3 = 0;
	ncritical = wg3.countsum.u16[3];// to update assigning critical 2 pairs
	ndead = BIT_SET_27;
	cout << " call handler go p_cpt2g[8]=" << p_cpt2g[8] << endl;
	if (0 &&g17b.debug17 > 1){
		//Print_Ok3();
	}
	rknown_b3 = 0;
	g17xy.ntuab3 = 0;
	if (nmiss){
		if (1) return;
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
		Go_Critical();
	}
}
//================= critical process
void G17B3HANDLER::Critical2pairs(){// assign 2 pairs in minirow to common cell
	int tbitsrows[8] = { 0, 07, 07000, 07007, 07000000, 07000007, 07007000, 07007007 };
	if (wg3.count.u16[1]){// and 2 pairs in minirow forced to common cell
		register int Rst = 07007007;// stack 0 pattern 
		for (int ist = 0; ist < 3; ist++){
			int shrink = TblShrinkMask[wg3.count.u16[1] & (0111 << ist)];
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
	Critical2pairs();// assign 2 pairs in minirow to common cell
	rknown_b3 = known_b3| active_b3;
	cout << Char27out(known_b3) << " known b2" << endl;
	cout << Char27out(active_b3) << " active_b3" << endl;
	cout<<Char27out(rknown_b3) << "entry critical appel ismultiple" << endl;
	if(IsMultiple(rknown_b3))	return;
	if (1)return;
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
	int nmissb = g17xy.b3lim - _popcnt32(known_b3);// missing clues
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
	int ncl = _popcnt32(known_b3);
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
	int cct = _popcnt32(wg3.pairs.u32[1]) - ncritical;
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
register int cc = _popcnt32(Ru);\
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
	{
		//if (g17b.debug17 > 1)cout << " one_out_forced:" << endl;
		p_cpt[10]++;
		int st = 07007007;// try to reduce the count using stacks
		for (int ist = 0; ist < 3; ist++) {
			int wastack = wact_miss1 & st;
			if (wastack) {
				int bfstack = wg3.pairs.u32[1] | wastack | known_b3;
				if (IsMultiple(bfstack)) {
					wact_miss1 ^= wastack;
					//wact_miss1 &= zhou[0].glb.go_back;// use received ua to reduce the count
				}
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
void G17B::PrintEndPuz() {
	cout << endl;
	for (int i = 0; i < 10; i++) if (p_cpt1[i])
		cout << g17tl1[i] << "\t" << p_cpt1[i] << endl;
	cout << endl;
	for (int i = 0; i < 32; i++) if (p_cpt[i])
		cout << g17tl[i] << "\t" << p_cpt[i] << endl;
}
int G17B::Godebug17_1(int i1) { // only if g17test_on 
	if (!GTEST17_ON) return 0;
/*
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
*/
	return 0;
}
int G17B::Godebug17_2(int i2) {// only if g17test_on 
	if (!GTEST17_ON) return 0;
	/*
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
	*/
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


void G17TB3GO::Debug() {
	cout << "status G17TB3GO ib3=" << ib3 << endl;
	char ws[82];
	cout << g17xy.bands_active_triplets.String3X(ws)
		<< "g17xy active triplets" << endl;
	cout << genb12.bands3[ib3].band << " \tbande3" << endl;
	cout << Char27out(pairs.u32[1]) << " in field bf" << endl << endl;
	cout << Char27out(pairs.u32[0]) << " pairs bf" << endl<<endl;
	cout << Char9out(count.u16[3]) << " all minis2" << endl;
	cout << Char9out(count.u16[0]) << "     minis2/1" << endl;
	cout << Char9out(count.u16[1]) << "     minis2/2" << endl;
	cout << Char9out(count.u16[2]) << "     minis2/3" << endl<<endl;
	cout << Char9out(minirows_triplets) << "mini triplets" << endl << endl;

	cout << countsum.u16[0] << "\t" << countsum.u16[1] << "\t" << countsum.u16[2]
		<< "\t total " << countsum.u16[3] << endl;
	cout << countstack.u16[0] << "\t" << countstack.u16[1] << "\t" << countstack.u16[2]
		<< "\t total " << countstack.u16[3] << endl;
}