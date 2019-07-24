
/*
442 dat file to txt file

..\release\skmpp    -i"extractclues"  -c"446"   -v24         "//split =idpos clues"
 #  19    1.2/1.2/1.2 - gsf
 1.2/1.2/1.2 - gsf
*/

//____________________________________________________________

/*	400 small tasks
		0x add 0,sequence 1,string0 2 nclues 
			9 std puzzle
			1x replace 0,empty'.' 1,erase'"'  2 cut to vx[1]
	401 .dat to .txt
	402 morph a puzzle

*/
void Go_c400_40p(char *ze, int task) {// work on bitfields std puz
	BF128 pattern; pattern.SetAll_0();
	int d_bands[6], d_units[27],digits=0;
	int c_bands[6], c_units[27];
	memset(d_bands, 0, sizeof d_bands);
	memset(d_units, 0, sizeof d_units);
	memset(c_bands, 0, sizeof c_bands);
	memset(c_units, 0, sizeof c_units);
	const char * rcb  = "RCB";
	for (int i = 0; i < 81; i++)if (ze[i] != '.') {
		pattern.Set_c(i);
		int d = ze[i] - '1',bit=1<<d;
		if (d < 0 || d>8) return;// invalid puzzle
		CELL_FIX & cf = cellsFixedData[i];
		int band = cf.el / 3, stack = cf.pl / 3 + 3;
		d_bands[band] |= bit; d_bands[stack] |= bit; digits |= bit;
		d_units[cf.el] |= bit; d_units[cf.plu] |= bit; d_units[cf.ebu] |= bit;
		c_bands[band]++; c_bands[stack]++;
		c_units[cf.el]++; c_units[cf.plu]++; c_units[cf.ebu]++;
	}
	ze[81] = 0;
	fout1 << ze;
	switch (task) {
	case 40: // count digits
		fout1 << ";" << _popcnt32(digits);
		break;
	case 41: // count given per band
		for (int i = 0; i < 6; i++) fout1 << ";" << c_bands[i];
		break;
	case 42: // count digits per band
		for (int i = 0; i < 6; i++) fout1 << ";" << _popcnt32(d_bands[i]);
		break;
	case 43: // count given/digits per unit
		for (int i1 = 0, iu = 0; i1 < 3; i1++) {
			fout1 << ";" << rcb[i1];
			for (int i = 0; i < 9; i++,iu++) fout1 << ";" << c_units[iu];
		}
		break;
	}
	fout1 << endl;
}
void Go_c400() {// small tasks on entry -v0- is the task
	int task = sgo.vx[0];// parameter defining the small task
	cout << "Go_400 entry " << sgo.finput_name 
		<< " small tasks on entry subtask "<<task << endl;
	uint64_t npuz = 0;
	char * ze = finput.ze;
	char zout[200];
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()) {
		int ll = (int)strlen(ze);
		if (ll < 81 || ll > 250) continue;
		if (task >= 40 && task<50) {	Go_c400_40p(ze, task); continue;}
		switch (task) {

		//=========== add to ze

		case 0:// add sequence
			fout1 << ze << ";" << npuz++ << endl;			break;
		case 1:// add string 0
			fout1 << ze << ";" << sgo.s_strings[0] << endl; break;
		case 2:// add nclues 
			{	int nc = 0;
				for (int i = 0; i < 81; i++)
					if (ze[i]>= '1' && ze[i]<= '9') nc++;
				fout1 << ze << ";" << nc << endl; break;
			}
		case 9:// add std puzzle to entry out puz + std puz
			{
				char zs[82];	zs[81] = 0;
				if (ll < 81) break;
				ze[81] = 0;
				for (int i = 0; i < 81; i++) {
					char c = ze[i];
					if (c > '0' && c <= '9') zs[i] = c;
					else zs[i] = '.';
				}
				fout1 << ze << ";" << zs << endl;
				break;
		}

		//===================== replace 

		case 10://'.' if not '0' '9' in puzzle area 
			{	if (ll < 81) break;
				for (int i = 0; i < 81; i++) {
					int c = ze[i];
					if (c < '1' || c > '9')ze[i] = '.';
				}
				fout1 << ze  << endl; break;
			}
		case 11://erase " in the entry
			for (int i = 0; i < ll; i++) if (ze[i] - '"') fout1 << ze[i];
			fout1 << endl; break;
		case 12:// cut entry to sgo.vx[1]
			if (ll < (int)sgo.vx[1]) break;
			ze[sgo.vx[1]] = 0;
			fout1 << ze << endl; break;
		case 15:// maxtext in output
			{	int known = 0;
				char vec[9], vi = '9';
				for (int i = 0; i < 81; i++)if (ze[i] - '.') {
					int c = ze[i] - '1', bit = 1 << c;;
					if (! (known&bit)) {
						vec[c] = vi--;
						known|=bit;
					}
					ze[i] = vec[c];
				}
				fout1 << ze << endl; 
				break;
			}
		case 16:// mintext in output
		{	int known = 0;
			char vec[9], vi = '1';
			for (int i = 0; i < 81; i++)if (ze[i] - '.') {
				int c = ze[i] - '1', bit = 1 << c;;
				if (!(known&bit)) {
					vec[c] = vi++;
					known |= bit;
				}
				ze[i] = vec[c];
			}
			fout1 << ze << endl;  
			break;
		}
		case 19:// int ratings to floating ratings pot hardest
		{	 // default location parameters ,2,3,4
			strncpy(zout, ze, 84);// puz;nn
			if(!sgo.bfx[0])sgo.bfx[0] = 2 + 4 + 8; //.111
			sgo.ParseInt(ze,';');
			int er = sgo.tparse[1], ep = sgo.tparse[2], ed = sgo.tparse[3];
			ze[81] = 0;
			fout1 << ze<<";"<<er/10<<"."<<er%10 
				<< ";" << ep / 10 << "." << ep % 10
				<< ";" << ed / 10 << "." << ed % 10 << endl;
			break;
		}

		//================= extract  
		case 21:// extract 81 start sgo.vx[1]
			if (ll >= (int)sgo.vx[1]) {
				ze[sgo.vx[1] + 81] = 0;
				fout1 << &ze [sgo.vx[1]] << endl;
			}
			break;
		case 22:// extract first sgo.vx[1] puzzles
			if (npuz++ < sgo.vx[1]) fout1 << ze << endl;
			else fout2 << ze << endl;
			break;
		case 23://sampling skip sgo.vx[1] one every sgo.vx[2]
			if (npuz++ < sgo.vx[1]) break;
			{	int rn = (npuz - sgo.vx[1]) % sgo.vx[2];
				if (!rn)fout1 << ze << endl;
			}
			break;
		//====================== split
		case 32:// split on  nclues 
			if (sgo.vx[1] < 18 || sgo.vx[1] >30) {
				cerr << "invalid split nclues value -v1- " << sgo.vx[1] << endl;
				return;
			}
			else {
				int nc = 0, lim = sgo.vx[1];
				for (int i = 0; i < 81; i++)if (ze[i] >= '1' && ze[i] <= '9') nc++;
				if (nc < lim)	fout1 << ze << endl;
				else if (nc == lim)fout2 << ze << endl;
				else fout3 << ze << endl;
				break;
			}
		}


	}
}
void Go_c401() {// .dat to .txt
	cout << "Go_401 entry " << sgo.finput_name << "dat file to txt file" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char zs[200];  int n = 0;
	while (!finput.eof()) {
		char c;
		finput.get(c); if (finput.eof()) break; // fin ligne normale
		if (c < 20) {//fin ligne 0a
			if (n) {// no line length 0
				zs[n] = 0;
				fout1 << zs << endl;
				n = 0;
			}
		}
		else {
			if (n < 190) zs[n++] = c; // cut line more than 190
		}
	}
	if (n) {// last line if exist
		zs[n] = 0;
		fout1 << zs << endl;
	}
}
void Go_c402() {// morph row s[1] cols s[2] v[1] diag
	cout << "Go_c402 entry " << sgo.finput_name << " morph a puzzle" << endl;
	char ze[82],zout[82]; ze[81] = zout[81] = 0;
	char * sr = sgo.s_strings[1], *sc = sgo.s_strings[2];
	if (!sr || !sc) {
		cout << "check missing -s1- for rows and/or -s2- for cols " << endl;
		return;
	}
	int rows = 0, cols = 0, tr[9], tc[9];
	int lrows = (int)strlen(sr), lcols = (int)strlen(sc);
	if (lrows != 9 || lcols != 9) {
		cout << "check -s1- for rows -s2- for cols must digits 1-9" << endl;
		return;
	}
	for (int i = 0; i < 9; i++) {
		int r = sr[i], c = sc[i];
		if (r<'1' || r>'9' || c<'1' || c>'9') {
			cout << "check -s1- -s2- for  digits 1-9" << endl;
			return;
		}
		r -= '1'; c -= '1';
		rows |= 1 << r; cols |= 1 << c;
		tr[i] = r; tc[i] = c;
	}
	if(rows!=0x1ff || cols!=0x1ff) {
		cout << "check -s1- -s2- for  digits 1-9 must all be there" << endl;
		return;
	}

	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	// build the morphing table
	int tmorph[81];
	for (int ir = 0; ir < 9; ir++) for(int ic=0;ic<9;ic++){
		tmorph[9 * ir + ic] = 9 * tr[ir] + tc[ic];
	}
	if (sgo.vx[1]) for(int i=0;i<81;i++)// then make it diagonal 
		tmorph[i]= C_transpose_d[tmorph[i]];
	while (finput.GetPuzzle(ze)) {
		for (int i = 0; i < 81; i++) zout[i] = ze[tmorph[i]];
		fout1 << zout << endl;
	}
}


void Go_cxx() {
	cout << "Go_xx entry " << sgo.finput_name << " xxxxx" << endl;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	//char * ze = finput.ze;
	while (finput.GetLigne()) {
	}

	char zep[82]; zep[81] = 0;
	while (finput.GetPuzzle(zep)) {
	}
}

//_________________________________________________________________________

void Go_c440(){
	cout << "Go_440 entry " << sgo.finput_name  <<" game results to parse"<< endl;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()){
		int ll = (int)strlen(ze);
		if (ll<96 || ll > 200) continue;
		char zout[200];
		for (int i = 0; i < 81; i++)
			if (ze[i] - '0')zout[i] = ze[i];
			else zout[i] = '.';
		zout[81] = ';';
		int n = 82;
		for (int i = 88; i < 102; i++){
			if (ze[i] == '/') {
				zout[n++] = ';'; continue;
			}
			if (ze[i] == ' ')continue;
			zout[n++] = ze[i];
		}
		zout[n++] = ';';
		if (ze[103] == '+')zout[n++] = '+'; else zout[n++] = ' ';
		zout[n++] = ';';
		strcpy(&zout[n], &ze[104]);
		fout1 << zout << endl;
	}

}
void Go_c445(){// filter on one integer parameter
	cout << "Go_445 entry " << sgo.finput_name << " param=" << sgo.bfx[0] << endl;
	if (_popcnt32(sgo.bfx[0]) != 1) return;//pointer to the  parameter to consider
	uint32_t ipar; bitscanforward(ipar, sgo.bfx[0]);
	cout << "split o, parameter rank=" << ipar << " file 1 <=" << sgo.vx[0] << endl;
	int v = sgo.vx[0];
	//ipar--;// switch to index;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()){
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()){
		int ll = (int)strlen(ze);
		if (ll<82 || ll > 200) continue;
		sgo.ParseInt(ze, ';');
		//cout << ze << " check v=" << sgo.tparse[ipar] << endl;
		if (sgo.tparse[ipar] <= v) fout1 << ze << endl;
		else fout2 << ze << endl;
		//break;
	}

}

void Go_c470() {//extract not equal in regressive test
	cout << "Go_470 entry " << sgo.finput_name << " first entry" << endl;
	FINPUT fin2;
	if (!sgo.foutput_name) {
		cerr << "missing output root" << endl;
		return;
	}
	if (!sgo.s_strings[0]) {
		cerr << "missing  file 2 in" << endl;
		return;
	}
	cout << "Go_470 entry " << sgo.s_strings[0] << " second entry" << endl;

	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	fin2.open(sgo.s_strings[0]);
	if (!fin2.is_open()) {
		cerr << "error open second entry " << sgo.s_strings[0] << endl;
		return;
	}
	uint64_t cpt = 0,cpt2=0;
	//char zc1[100], zc2[100];
	while (1) {
		if (!finput.GetLigne()) break;
		if (!fin2.GetLigne()) break;
		cpt++;
		if (strncmp(finput.ze, fin2.ze, 81)) {
			cout << "synchro broken" << endl;
			break;
		}
		switch (sgo.vx[0]) {
			case 0:// 2 outputs of skmppv2
				if (strcmp(finput.ze, fin2.ze)) {
					cout << finput.ze << endl << fin2.ze << endl;
					fout1 << finput.ze  << endl;
					fout2 << fin2.ze << endl;
					fout3 << finput.ze << &fin2.ze[81] << endl;
					cpt2++;
				}
				break; 
		}
	}
	cerr << cpt<<" records processed "<<cpt2 <<" not equal"<<endl;
}


void Go_c446() {// split on ER EP ED potential hardest
	cout << "Go_446 entry " << sgo.finput_name << " file1 pot hardest"  << endl;
	int er = sgo.vx[0],ep = sgo.vx[1], ed = sgo.vx[2];
	cout << "er " << er << " ep " << ep << " ed " << ed << endl;
	//ipar--;// switch to index;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()) {
		int ll = (int)strlen(ze);
		if (ll < 82 || ll > 200) continue;
		sgo.ParseInt(ze, ';');
		//cout <<ze <<  " check v=" << sgo.tparse[1] 
		//	<< " " << sgo.tparse[2] << " " << sgo.tparse[3] << endl;
		if (sgo.tparse[1] < 100) continue;
		if (sgo.tparse[1] >=er) goto isok;
		if (sgo.tparse[2] >= ep) goto isok;
		if (sgo.tparse[3] >= ed) goto isok;
		fout2 << ze << endl;
		continue;
	isok:
		fout1 << ze << endl;
		//break;
	}

}
void Go_c480() {//add  compressed clues to entry
	cout << "Go_481 entry " << sgo.finput_name << " base check" << endl;
	if (!sgo.foutput_name) {
		cerr << "missing output root" << endl;
		return;
	}
	if (!sgo.s_strings[0]) {
		cerr << "missing base in" << endl;
		return;
	}
	if (!sgo.s_strings[1]) {
		cerr << "missing base out" << endl;
		return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char ze[200]; ze[81] = ';';
	int lcmp = 0;
	while (finput.GetPuzzle(ze)) {
		int  p = 82;
		for (int i = 9; i < 81; i++)if (ze[i] - '.')
			ze[p++] = ze[i];

		int ll = p-82;
		if (lcmp && lcmp != ll) {
			cerr << "stop wrong input file" << endl;
		}
		else lcmp = ll;
		ze[p] = 0;
		fout1 << ze << endl;
	}

}
void Go_c481() {//base check -i ads -s1- base -o add root
	FINPUT fin2;
	cout << "Go_481 entry " << sgo.finput_name << " base check" << endl;
	if (!sgo.foutput_name) {
		cerr << "missing output root" << endl;
		return;
	}
	if (!sgo.s_strings[0]) {
		cerr << "missing base in" << endl;
		return;
	}
	int update  = sgo.vx[0];
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	fin2.open(sgo.s_strings[0]);
	if (!fin2.is_open()) {
		cerr << "error open base " << sgo.s_strings[0] << endl;
		return;
	}
	

	//__________________________________
	int lcmp = 0;
	char scomp[10]; scomp[0] = 0;
	char * zcompress = &ze[81],*zeb= fin2.ze;;
	zeb[0] = 0;
	while (finput.GetLigne()) {
		int ll = (int)strlen(ze);
		if (ll < 91 || ll > 120) {
			cerr << "wrong add file cancelled " << endl;
			return;
		}
		if(lcmp){
			if (ll != lcmp) {
				cerr << "stop wrong add file" << endl;
				return;
			}
		}
		else lcmp = ll;
		//int lbase = (int)strlen(zcompress);
		// read/ out if update base as long as below
	loop_base:
		int icomp = strcmp(zcompress, zeb);
		if (!icomp)goto next_add;
		if (icomp > 0) {// next in base
			if (zeb[0] && update) fout2 << zeb << endl; //write old if update
			if(fin2.GetLigne())	goto loop_base;
		}
		// new or end of file in old base
		if (update) {
			fout2 << zcompress << endl;// new in base
			ze[81] = 0;
			fout1 << ze << endl;// just the puzzle in output
		}
		else fout1 << ze << endl;// just recopy in output
			//char zout[200];
	next_add:;
	}
	if (update && zeb[0] !=-1) {// copy the rest of the base
		fout2 << zeb << endl;
		while(fin2.GetLigne())fout2 << zeb << endl;
	}
}


void Go_c484() {
	cout << "Go_484 entry " << sgo.finput_name << " data base to restore" << endl;
	if (!sgo.s_strings[0]) {
		cerr << "missing -s1- pattern " << endl;
		return;
	}
	if (strlen(sgo.s_strings[0]) != 81) {
		cerr << "-s1- not length 81 " << endl;
		return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze;
	char zout[82];
	strcpy(zout, empty_puzzle);
	int tclues[50], nclues = 0,puz_int[81];
	memset(puz_int, 0, sizeof puz_int);
	char *w = sgo.s_strings[0];
	for (int i = 0; i < 81; i++)
		if (w[i] >= '1' && w[i] <= '9')puz_int[i] = w[i];
	for (int i = 0; i < 9; i++) if (puz_int[i]) zout[i] = (char)puz_int[i];
	for (int i = 9; i < 81; i++) 
		if (puz_int[i]&& nclues<40) tclues[nclues++]=i;
	if (nclues > 30) {
		cerr << "-s1- too many clues cancel " << endl;
		return;
	}
	uint64_t ne = 0,nd= sgo.vx[1],nf= sgo.vx[2];
	while (finput.GetLigne()) {
		if (++ne < nd)continue;
		if (nf && ne >= nf)break;
		ze[1] = 0;// temp code to clear the warning
		cout << ze << tclues[0] << endl;
	}

}


void Go_c485() {// check close to first entry in data base canonical
	FINPUT fin2; //data base
	cout << "Go_485 entry " << sgo.finput_name << " puzzle to check" << endl;
	if (!sgo.s_strings[0]) {
		cerr << "missing -s1- database full name " << endl;
		return;
	}
	fin2.open(sgo.s_strings[0]);
	if (!fin2.is_open()) {
		cerr << "error open file " << sgo.s_strings[0] << endl;
		return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	char * ze = finput.ze,*ze2= fin2.ze;
	char zout[82];
	int tclues[50], ticlues[50],nclues = 0;
	while (finput.GetLigne()) {
		strncpy(zout, ze, 9);// first row un touched
		strcpy(&zout[9], &empty_puzzle[9]);

		for (int i = 9; i < 81; i++) {
			int c = ze[i];
			if (c > '0' && c <= '9') {
				ticlues[nclues] = i;
				tclues[nclues++] = c;
			}
			if (nclues > 30) {
				cerr << " too many clues cancel " << endl;
				return;
			}
		}
		while (fin2.GetLigne()) {
			if ((int)strlen(ze2) != nclues) {
				cerr << " wrong clues number data base "
					<< strlen(ze2) << " expected "<<nclues<< endl;
				return;
			}
			int nm = 0;
			for (int i = 0; i < nclues; i++) {
				int iout = ticlues[i], c = tclues[i];
				zout[iout] = ze2[i];
				if ((int)ze2[i] != c)nm++;
			}
			if (nm <= 5)fout1 << zout << ";" << nm << endl;
		}
		break;// only one puzzle per run
	}

}


/*


void GO_MISC::Do_70(){  // start a game from gremlin pattenr
	char * ze=myin->ze,zs[82];
	zs[81]=0;
	for(int i=0;i<9;i++){
		myin->GetLigne();
		int length=strlen(ze);
		(*myout1)<<ze<<endl;
		for(int j=0;j<9;j++)
			zs[9*i+j]=ze[2*j];
	}
	(*myout1)<<endl<<endl<<zs<<endl;
	for(int i=0;i<81;i++)
		if(zs[i]-'.') zs[i]='1';
	(*myout1)<<endl<<zs<<endl<<endl;
	char *a="1__2__3__",*b="4__5__6__",*c="7__8__9__";
	(*myout1)<<a<<a<<a<<b<<b<<b<<c<<c<<c<<endl;

}

void Go_c484(){  // restore game database for a pattern
	char * ze=myin->ze,zs[82];
	zs[81]=0;
	PUZ0 puzd;
	puzd.Empty();
	if((!options.first) ||	(strlen(options.first)-81) ){
		(*myout2)<< " pattern missing or length not correct "<<endl;
		return;
	}
	(*myout2)<<options.first<< " pattern "<<endl;
	int nd=0;
	for(int i=9;i<81;i++)
		if(options.first[i]-'.') nd++;

	while(myin->GetLigne()){
		if(strlen(ze)-nd) continue;
		int nx=0,lig1='1';
		for(int i=0;i<81;i++)
			if(options.first[i]=='.')
				zs[i]='.';
			else
				if(i<9)
					zs[i]=lig1++;
				else
					zs[i]=ze[nx++];

		(*myout1)<<zs<<endl;
	}

}

void GO_MISC::Do_72(){  // find puzzles close to a pattern
	char * ze=myin->ze,ztarget[82],maxf=11;
	PUZ0 puzd;
	puzd.Empty();
	if((!options.first) ||	(strlen(options.first)-81) ){
		(*myout2)<< " pattern missing or length not correct "<<endl;
		return;
	}
	(*myout2)<<options.first<< " target "<<endl;
	int nd=0;
	for(int i=9;i<81;i++)
		if(options.first[i]-'.') ztarget[nd++]=options.first[i];

	while(myin->GetLigne()){
		if(strlen(ze)-nd) {
			(*myout2)<< " length not correct in the file"<<endl;
			return;
		}
		int nx=0;
		for(int i=0;i<nd;i++)
			if(ztarget[i]-ze[i]) nx++;
		if(nx>=maxf) continue;
		if(nx>=6)maxf=nx;
		(*myout2)<<ze<<";"<<nx<<endl;
		if(nx<5)
			(*myout1)<<ze<<endl;

	}

}

*/

//================= symmetry analysis

int t_perms_boxes[9][9] = {// 9 boxes perms to have each box central
	{1,0,2,4,3,5,7,6,8}, // 3 central
	{0,1,2,3,4,5,6,7,8}, // 4 central and start
	{0,2,1,3,5,4,6,8,7}, // 5 central
	{1,0,2,7,6,8,4,3,5}, // 6 central
	{0,1,2,6,7,8,3,4,5}, // 7 central
	{0,2,1,6,8,7,3,5,4}, // 8 central
	{4,3,5,1,0,2,7,6,8}, // 0 central
	{3,4,5,0,1,2,6,7,8}, // 1 central  
	{3,5,4,0,2,1,6,8,7}, // 2 central
};
int rc_perms[9][2][9] = {
	0,1,2,3,4,5,6,7,8,	3,4,5,0,1,2,6,7,8,// 3 central
	0,1,2,3,4,5,6,7,8,	0,1,2,3,4,5,6,7,8,// 4 central
	0,1,2,3,4,5,6,7,8,	0,1,2,6,7,8,3,4,5,// 5 central
	0,1,2,6,7,8,3,4,5,	3,4,5,0,1,2,6,7,8,// 6 central
	0,1,2,6,7,8,3,4,5,	0,1,2,3,4,5,6,7,8,// 7 central
	0,1,2,6,7,8,3,4,5,	0,1,2,6,7,8,3,4,5,// 8 central
	3,4,5,0,1,2,6,7,8,	3,4,5,0,1,2,6,7,8,// 0 central
	3,4,5,0,1,2,6,7,8,	0,1,2,3,4,5,6,7,8,// 1 central
	3,4,5,0,1,2,6,7,8,	0,1,2,6,7,8,3,4,5,// 2 central
};
int t_perms_boxes_rc[9][6] = {// 9 boxes perms to have each cell central
	{0,1,2,1,0,2}, //3 cols 0 1
	{0,1,2,0,1,2}, //4 no morph
	{0,1,2,0,2,1}, //5 cols 12
	{0,2,1,1,0,2}, //6 cols 01
	{0,2,1,0,1,2}, //7 rows 12
	{0,2,1,0,2,1}, //8 rows 12 cols 12
	{1,0,2,1,0,2}, //0  rows 01 cols 01
	{1,0,2,0,1,2}, //1 rows 0 1 
	{1,0,2,0,2,1}, //2 rows 0 1 cols 12
};
int tpermsdiag[6][6] = {// row 1 "3 cases" x cols 2perms
	{0,1,2,0,1,2},
	{0,1,2,0,2,1},
	{1,0,2,0,1,2},
	{1,0,2,0,2,1},
	{2,0,1,0,1,2},
	{2,0,1,0,2,1},
};
int t_perms_boxes_diagonal[6][9] = {// for main diagonal symmetry
	{0,1,2,3,4,5,6,7,8}, // box 0 base
	{0,2,1,3,5,4,6,8,7}, //		perm stack
	{1,0,2,4,3,5,7,6,8}, // box 1 base
	{1,2,0,4,5,3,7,8,6}, //		perm satck
	{2,0,1,5,3,4,8,6,7}, // box 2 base
	{2,1,0,5,4,3,8,7,6}, //		perm stack

};
int tpdiag_c[6][9] = {// col order of each perm (rows unchanged)
	0,1,2,3,4,5,6,7,8,// box 0 base
	0,1,2,6,7,8,3,4,5,
	3,4,5,0,1,2,6,7,8,// box 1 base
	3,4,5,6,7,8,0,1,2,
	6,7,8,0,1,2,3,4,5,// box 2 base
	6,7,8,3,4,5,0,1,2,
};

struct PUZC_SYM {
	char puz[82],puz2[82], puz3[82]; // normalized puzzle
	uint32_t boxes[9],
		ctb[9], wctb[9], // count per box  
		ctbs[6],//band stack count
		nclues,cb_parity,sym_boxes;
	BF128 pat;
	void InitPat(char* ze) {
		strncpy(puz, ze, 81);
		pat.SetAll_1();
		UpdatePat(ze);
	}
	void Init(char * ze) {
		nclues = 0;
		memset(boxes, 0, sizeof boxes);
		memset(ctb, 0, sizeof ctb);
		strncpy(puz, ze, 81);
		puz[81] = puz2[81] = puz3[81] = 0;
		for (int i = 0; i < 81; i++) {
			if (puz[i] < '1' || puz[i] > '9') {
				puz[i] = '.';
				continue;
			}
			boxes[C_box[i]] |= 1 << C_box_rel[i];
			ctb[C_box[i]]++;
			nclues++;
		}
		cb_parity = nclues & 1;
	}
	void UpdatePat(char *ze) {
		BF128 wpat; wpat.SetAll_0();
		for (int i = 0; i < 81; i++)if (ze[i] != '.')wpat.Set_c(i);
		if (wpat.Compare(pat) < 0) pat = wpat;
	}
	int GetSym();
	int MorphPatB5(PUZC_SYM & pat,int diag=0);
	void PatCentral() {
		int tp[12][5] = {
			0,1,2,3,4,	0,2,1,3,4,
			1,0,2,3,4,	1,2,0,3,4,
			2,0,1,3,4,	2,1,0,3,4,
			0,1,2,5,4,	0,2,1,5,4,
			1,0,2,5,4,	1,2,0,5,4,
			2,0,1,5,4,	2,1,0,5,4,
		};
		char ws[82], puz2d[82];;
		for (int i1 = 0; i1 < 12; i1++)for (int i2 = 0; i2 < 12; i2++) {
			int *pr = tp[i1], *pc = tp[i2];
			for (int ir = 0; ir < 5; ir++)for (int ic = 0; ic < 5; ic++) {
				int d1 = 9 * ir + ic, d2 = 9 * ir + 8 - ic, d3 = 80 - d2, d4 = 80 - d1;
				int o1=9*pr[ir]+pc[ic],o2= 9 * pr[ir] + 8-pc[ic], o3 = 80 - o2, o4 = 80 - o1;
				puz2[d1] = puz[o1]; puz2[d2] = puz[o2]; puz2[d3] = puz[o3]; puz2[d4] = puz[o4];
				puz2d[C_transpose_d[d1]] = puz[o1];
				puz2d[C_transpose_d[d2]] = puz[o2];
				puz2d[C_transpose_d[d3]] = puz[o3];
				puz2d[C_transpose_d[d4]] = puz[o4];
			}
			UpdatePat(puz2);
			UpdatePat(puz2d);
		}
		puz[81] = 0;
		fout1 << puz << ";" << pat.String3X(ws) << endl;
	}

	void PatDiagonal() {// try all perms on the diagonal
		char ws[82];
		for (int ipr = 0; ipr < 1296; ipr++) {// perms on diagonal
			int * rc = tpermorder[ipr]; // same perm for both
			for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
				int i0 = 9 * ir + ic, i1 = 9 * rc[ir] + rc[ic];
				puz2[i0] = puz[i1];
			}
			UpdatePat(puz2);
		}
		puz[81] = 0;
		fout1 << puz<<";"<<pat.String3X(ws)<<endl;
	}

};


int PUZC_SYM::GetSym() {
	for (int iperm = 0; iperm < 9; iperm++) {
		for (int i = 0; i < 9; i++) 			wctb[i] = ctb[t_perms_boxes[iperm][i]];
		// central box must have the right parity to be central
		if (cb_parity != (wctb[4] & 1)) continue;
		for (int i = 0; i < 5; i++)	if (wctb[i] != wctb[8 - i])continue;
		// morph puz2 to morph
		int * rr = rc_perms[iperm][0], *cc = rc_perms[iperm][1];
		for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
			int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic];
			puz2[i0] = puz[i1];
		}
		// and test all morphs
		for (int ipr = 0; ipr < 216; ipr++) {// perms on rows
			for (int ipc = 0; ipc < 216; ipc++) {// perms on columns
				int * rr = tpermorder[ipr], *cc = tpermorder[ipc];
				for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
					int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic];
					puz3[i0] = puz2[i1];
				}
				for (int i0 = 0; i0 < 40; i0++) {
					int  i1 = 80-i0, ch1 = puz3[i0], ch2 = puz3[i1];
					if (ch1 == '.' && ch2 != '.') goto nextperm_rc_central;
					if (ch1 != '.' && ch2 == '.') goto nextperm_rc_central;
				}
				return 1;// valid puz3 found available in puz3
			nextperm_rc_central:;
			}
		}

	}
	//try now diagonal 6 perm to have box1 b0 or b3 or b6 top left
	for (int iperm = 0; iperm < 6; iperm++) {// six perms box1 from stack1
		for (int i = 0; i < 9; i++)wctb[i] = ctb[t_perms_boxes_diagonal[iperm][i]];
		if( (wctb[1] != wctb[3])|| (wctb[2] != wctb[6])|| (wctb[5] != wctb[7]))continue;
		// morph puz2 to morph
		int  *cc = tpdiag_c[iperm];
		for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
			int i0 = 9 * ir + ic, i1 = 9 * ir + cc[ic];
			puz2[i0] = puz[i1];
		}
		// and test all morphs
		for (int ipr = 0; ipr < 216; ipr++) {// perms on rows
			for (int ipc = 0; ipc < 216; ipc++) {// perms on columns
				int * rr = tpermorder[ipr], *cc = tpermorder[ipc];
				for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
					int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic];
					puz3[i0] = puz2[i1];
				}
				for (int ir = 0; ir < 8; ir++)for (int ic = ir + 1; ic < 9; ic++) {
					int i0 = 9 * ir + ic, i1 = 9 * ic + ir, ch1 = puz3[i0], ch2 = puz3[i1];
					if (ch1 == '.' && ch2 != '.') goto nextperm_rc;
					if (ch1 != '.' && ch2 == '.') goto nextperm_rc;
				}
				//cout << puz3 << "morph ok" << endl;
				return 2;// valid puz3 found available in puz3
			nextperm_rc:;
			}
		}
	}
	return 0;
}
void Go_c490() {
	cerr << "Go_490 entry " << sgo.finput_name << " extract potential symmetry" << endl;
	PUZC_SYM psym;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()) {
		if (strlen(ze) < 81) continue;
		psym.Init(ze);
		int smode = psym.GetSym();
		if (smode) {
			//cout << ze << "ok" << endl;
			fout1<< psym.puz3 <<";" << ze << ";" << smode << endl;
		}
		else {}//			cout << ze << "\t pas ok" << endl;
	}

}

void Go_c492() {//puzzle;sym(1 central 2 diagonal
	cerr << "Go_492 entry " << sgo.finput_name << " find minimal pattern" << endl;
	PUZC_SYM psym;
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	while (finput.GetLigne()) {
		if (strlen(ze) < 83) continue;
		int smode = ze[82] - '1';
		psym.InitPat(ze);
		//cout << ze << " process smode="<<smode << "ze[82]"<< ze[82] << endl;
		switch (smode) {
		case 0: {//central symmetry
			psym.PatCentral();
		}
				break;
		case 1: {// diagonal symetry
			//cout <<ze<< "diag go" << endl;
			//psym.PatDiagonal();
		}
				break;
		}

	}

}
int PUZC_SYM::MorphPatB5(PUZC_SYM & pat,int diag) {
	char  puzcorner[82], puzcornerdiag[82]; // normalized puzzle
	puzcorner[81] = puz2[81] = 0;
	for (int iperm = 0; iperm < 9; iperm++) {// start with central symmetry
	// reorder boxes and build pattern
		for (int i = 0; i < 9; i++) wctb[i] = ctb[t_perms_boxes[iperm][i]];
		// central box count must be ==
		if (pat.ctb[4] != wctb[4]) continue;
		// do the box perm for the given entry
		int * rr = rc_perms[iperm][0], *cc = rc_perms[iperm][1];
		for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
			int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic];
			puz2[i0] = puz[i1];
		}
		//===========================================================
		// 8 perms four boxes in top left + diagonal symmetry
		int  wctb2[9];
		int p2[8][9] = {
			{0,1,2,3,4,5,6,7,8},
			{2,1,0,5,4,3,8,7,6},
			{6,7,8,3,4,5,0,1,2},
			{8,7,6,5,4,3,2,1,0},
		};
		int rc_perms2[4][2][9] = {
			0,1,2,3,4,5,6,7,8,	0,1,2,3,4,5,6,7,8,// box 0 top left
			0,1,2,3,4,5,6,7,8,  6,7,8,3,4,5,0,1,2,// box 2 top left
			6,7,8,3,4,5,0,1,2,	0,1,2,3,4,5,6,7,8,// box 6 top left
			6,7,8,3,4,5,0,1,2,  6,7,8,3,4,5,0,1,2,// box 8 top left
		};
		for (int icorner = 0; icorner < 4; icorner++) {
			int * pp = p2[icorner];
			for (int i = 0; i < 9; i++) {
				wctb2[i] = wctb[pp[i]];
				if (wctb2[i] != pat.ctb[i]) goto nextcorner;
			}
			// morph puz2 to puzcorner
			int * rr = rc_perms2[icorner][0], *cc = rc_perms2[icorner][1];
			for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
				int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic];
				puzcorner[i0] = puz2[i1];
			}
			// and test all morphs 
			for (int ipr = 0; ipr < 216; ipr++) {// perms on rows
				for (int ipc = 0; ipc < 216; ipc++) {// perms on columns
					int * rr = tpermorder[ipr], *cc = tpermorder[ipc];
					for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
						int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic],
							ch= puzcorner[i1];
						puz3[i0] = ch;
						if (ch == '.' && pat.puz[i0] != '.') goto nextperm_rc;
						if (ch != '.' && pat.puz[i0] == '.') goto nextperm_rc;
					}
					return 1;// valid puz3 found available in puz3
				nextperm_rc:;
				}
			}
			// test also diagonal (if central symmetry, can be compulsory)
			for (int i = 0; i < 81; i++)puzcornerdiag[i] = puzcorner[C_transpose_d[i]];
			for (int ipr = 0; ipr < 216; ipr++) {// perms on rows
				for (int ipc = 0; ipc < 216; ipc++) {// perms on columns
					int * rr = tpermorder[ipr], *cc = tpermorder[ipc];
					for (int ir = 0; ir < 9; ir++)for (int ic = 0; ic < 9; ic++) {
						int i0 = 9 * ir + ic, i1 = 9 * rr[ir] + cc[ic],
							ch = puzcornerdiag[i1];
						puz3[i0] = ch;
						if (ch == '.' && pat.puz[i0] != '.') goto nextperm_rc_diag;
						if (ch != '.' && pat.puz[i0] == '.') goto nextperm_rc_diag;
					}
					return 1;// valid puz3 found available in puz3
				nextperm_rc_diag:;
				}
			}
		nextcorner:;
		}
	}
	return 0;
}

void Go_c491() {// morph entry to a given pattern
	cout << "Go_490 entry " << sgo.finput_name << " morph to pattern"	<< endl;
	if (!sgo.s_strings[0]) {
		cerr << "missing pattern -s0-"  << endl;
		return;
	}
	if (strlen(sgo.s_strings[0]) != 81) {
		cerr<< sgo.s_strings[0] << "invalid pattern -s0-" << endl;
		return;
	}
	cout << sgo.s_strings[0] << " pattern to use" << endl;
	//entry must fit with the pattern
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	PUZC_SYM psym,psym_pat; 
	psym_pat.Init(sgo.s_strings[0]);
	for (int i = 0; i < 9; i++) cout << psym_pat.ctb[i];
	cout << "  compte par boite" << endl;
	char  *ze= finput.ze;  
	while (finput.GetLigne( )) {
		if (strlen(ze) < 81)continue;
		//cout << ze << "to process" << endl;
		psym.Init(ze);// create boxes 
		if (psym.nclues != psym_pat.nclues)continue; // invalid entry
		// start with b5 (as in c490 central) count = b5 pat
		//fout2 << ze << endl;
		if (psym.MorphPatB5(psym_pat)) {
			//cout<< psym.puz3 << "retour ok" << endl;
			fout1 << psym.puz3 << endl;
			fout2 << psym.puz3 << ";" << ze<<endl;
		}
		//if (1) break;
	}
}

//================== file compare


void Go_c495() {//compare log files
	FINPUT fin2;
	cout << "Go_495 entry " << sgo.finput_name << " base check" << endl;
	if (!sgo.s_strings[0]) {
		cerr << "missing file in 2" << endl;
		return;
	}
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file 1 " << sgo.finput_name << endl;
		return;
	}
	fin2.open(sgo.s_strings[0]);
	if (!fin2.is_open()) {
		cerr << "error open file2 " << sgo.s_strings[0] << endl;
		return;
	}

	int option = sgo.vx[0];
	uint64_t cpt = 0;
	while (1) {
	read1:
		cpt++;
		if (!finput.GetLigne()) break;
		if (0) goto read1;
	read2:
		if (!fin2.GetLigne()) break;
		if (0) goto read2;
		int ir = strcmp(finput.ze, fin2.ze);
		if (!ir)goto read1;
		cout << "files don't match cpt=" << cpt << endl;
		cout << finput.ze << endl;
		cout << fin2.ze << endl;
		return;
	}


}