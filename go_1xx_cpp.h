
// working in solver mode

/*
DynamicForcingChain=85,
DynamicForcingChainPlus=90,
// 90 is DynamicForcingChainPlus
// all "events"  claiming, pointing, pair, hidden pair, XWing
// follows an empty step 85
// consider each false candidate as start
// search for new bi values. If none, skip it
// look for new false thru basic sets
NestedForcingChain=95,
Nestedmultiple chains=100,
NestedForcingChain dynamic subchains dynamic=105,
NesttedLevel5=110
*/


extern PM_GO pm_go;
void Go_c110(){// template serate mode
    if (!sgo.finput_name) return;
	long long cptg[20];
	memset(cptg, 0, sizeof cptg);
	cout << "Go_110 entry " << sgo.finput_name << " input" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		if(sgo.bfx[8])cout << finput.ze << "to process npuz=" << npuz << endl;
		zh_g2.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << endl;
			continue;
		}
		pm_go.SolveSerate110();
		if (sgo.vx[1] && npuz >= sgo.vx[1]) break;
	}
}
//fout_diam, fout_pearl, fout_l45, fout_l65, fout_solved, fout_unsolved;
//fout1=l45;  fout2 solved;   fout3 unsolved
void Go_c111(){// fast serate mode 
	if (!sgo.finput_name) return;
	uint64_t cptg[20];
	memset(cptg, 0, sizeof cptg);
	cout << "Go_111 entry " << sgo.finput_name << " input" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	if (!sgo.foutput_name)return;
    {// create output files to split rated puzzles
		char zn[200];
		strcpy(zn, sgo.foutput_name);
		int ll = (int)strlen(zn);
		strcpy(&zn[ll], "_0diam.txt");		fout_diam.open(zn);// quasi diam 45_,,
		strcpy(&zn[ll], "_1pearl.txt");		fout_pearl.open(zn);// quasi pearl 45_,,
		strcpy(&zn[ll], "_2l65.txt");		fout_l65.open(zn);// rating 45 62
  	}

	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (0)cout << finput.ze << "to process npuz=" << npuz << endl;
		zh_g2.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << endl;
			continue;
		}
		pm_go.SolveSerate111();
	}
	

}
int  ParsePm(char * ze, int * cells) {
	int n = 0, nze = 0;
	while (ze[nze]) {
		register int c = ze[nze++];
		if (c<'1' || c>'9') continue;
		int & mycell = cells[n++];
		if (n > 9)return 0;
		mycell = 1 << (c - '1');
		while (ze[nze]) {
			register int c2 = ze[nze++];
			if (c2<'1' || c2>'9') break;
			mycell |= 1 << (c2 - '1');
		}
	}
	return n;
}
void Go_c118() { // study a sub grid
	if (!sgo.finput_name) return;
	cout << "Go_118 entry " << sgo.finput_name << " input" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char zp[82]; zp[81] = 0;
	if (!finput.GetPuzzle(zp)) {
		cout << "missing puzzle to study="   << endl;
		return;
	}
	if (zh_g.Go_InitSolve(zp)) {
		cout << finput.ze << "invalid or multiple solutions "  << endl;
		return;
	}
	cout << zp << "puzzle to study" << endl;
	cout << zh_g2.zsol << endl;
	// now the pm in 9 lines
	ZHOU myz;
	memset(&myz, 0, sizeof myz);
	myz.cells_unsolved = zhoustart[18];
	for (int i = 0; i < 9; i++)myz.FD[i][0].bf.u32[3] = 0x1ff;

	char *ze = finput.ze;
	for (int i = 0; i < 9; i++) {
		if (!finput.GetLigne()) {
			cout << "missing pm i="<<i << endl;
			return;
		}
		int tc[10], ntc = ParsePm(ze, tc);
		if (ntc != 9) {
			cout << ze << "illegal line n=" << i << endl;
			return;
		}
		for (int ic = 0; ic < 9; ic++) {
			int cell = 9 * i + ic,digs=tc[ic];
			for (int id = 0, bit = 1; id < 9; id++, bit <<= 1)if (digs&bit)
				myz.FD[id][0].Set_c(cell);
		}
	}
	myz.ImageCandidats();
	int  invalid_subpuzzle = 0;
	for (int i = 0; i < 81; i++) {
		int d = zh_g2.zsol[i] - '1';
		if (myz.FD[d][0].Off_c(i)) {
			invalid_subpuzzle = 1;
			pm_go.is_valid_puzzle = 0;
			break;
		}
	}
	if (invalid_subpuzzle) cout << "invalid_subpuzzle" << endl;
	else cout << "valid_subpuzzle" << endl;
	int ir = myz.FullUpdate();
	cout << "after full update ir ="<<ir << endl;
	myz.ImageCandidats();
	if (!ir)  return;
	zhou_solve = myz;
	pm_go.SolveDet(44,7,1);
	zhou_solve.ImageCandidats();
	pm_go.SolveDet(85, 7, 1);
	zhou_solve.ImageCandidats();
}
void Go_c199(){// test on demand
	if (!sgo.finput_name) return;
	cout << "Go_199 entry " << sgo.finput_name << " input" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)){
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		//cout << finput.ze << "to process npuz=" << npuz << endl;
		zh_g2.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << endl;
			continue;
		}
		pm_go.Solve199test();
		if (sgo.vx[1] && npuz >= sgo.vx[1]) break;
	}
}
void Go_c1xx(int subcase) {// test on demand
	if (!sgo.finput_name) return;
	cout << "Go_1xx entry  subcase "<< subcase
		<<" input "<< sgo.finput_name   << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = 0;
	uint32_t npuz = 0;
	while (finput.GetPuzzle(ze)) {
		npuz++;
		if (npuz < sgo.vx[0]) continue;
		//cout << finput.ze << "to process npuz=" << npuz << endl;
		zh_g2.npuz = npuz;
		if (zh_g.Go_InitSolve(ze)) {
			cout << finput.ze << "invalid or multiple solutions npuz=" << npuz << endl;
			continue;
		}
		switch (subcase) {
		case 120:pm_go.Solve120_MultiAnalysis(); break;
		case 125:pm_go.Solve125_Find_Vloop(); break;
		case 130:pm_go.Solve130_Find_JExocet(); break;
		}
		if (sgo.vx[1] && npuz >= sgo.vx[1]) break;
	}
}