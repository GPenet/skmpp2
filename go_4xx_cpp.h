
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
	{0,1,2,3,4,5,6,7,8}, // 4 central and start
	{0,2,1,3,5,4,6,8,7}, // 5 central
	{1,0,2,7,6,8,4,3,5}, // 6 central
	{0,1,2,6,7,8,3,4,5}, // 7 central
	{0,2,1,6,8,7,3,5,4}, // 8 central
	{4,3,5,1,0,2,7,6,8}, // 0 central
	{3,4,5,0,1,2,6,7,8}, // 1 central  
	{3,5,4,0,2,1,6,8,7}, // 2 central
	{1,0,2,4,3,5,7,6,8}, // 3 central
};
int t_perms_boxes_rc[9][6] = {// 9 boxes perms to have each box central
	{0,1,2,0,1,2}, //4 no morph
	{0,1,2,0,2,1}, //5 cols 12
	{0,1,2,1,0,2}, //6 cols 0 1
	{0,2,1,0,1,2}, //7 rows 12
	{0,2,1,0,2,1}, //8 rows 12 cols 12
	{1,0,2,1,0,2}, //0  rows 01 cols 01
	{1,0,2,0,1,2}, //1 rows 0 1 
	{1,0,2,0,2,1}, //2 rows 0 1 cols 12
	{0,1,2,1,0,2}, // 3 cols 01
};
int tpermsdiag[6][6] = {// row 1 "3 cases" x cols 2perms
	{0,1,2,0,1,2},
	{0,1,2,0,2,1},
	{1,0,2,0,1,2},
	{1,0,2,0,2,1},
	{2,0,1,0,1,2},
	{2,0,1,0,2,1},
};
int diag_remapping_index[6][9];
struct SYM_SPOT {
	uint32_t boxes[9],// current boxes bit fields
		rel_band_stack[6][3],//relatice order (row/col) inside the band/stack
		band_stack_order[6]; // bands/stack/reordering
	void Init(uint32_t * bb) {
		memcpy(band_stack_order, t_perms_boxes_rc, sizeof band_stack_order);
		memcpy(rel_band_stack[0], t_perms_boxes_rc, sizeof band_stack_order);
		memcpy(rel_band_stack[2], t_perms_boxes_rc, sizeof band_stack_order);
		memcpy(rel_band_stack[4], t_perms_boxes_rc, sizeof band_stack_order);
		memcpy(boxes, bb, sizeof boxes);
	}
 	void MovBandsStacks(SYM_SPOT & e, int ip) {
		int * perm = t_perms_boxes[ip], *permbs = t_perms_boxes_rc[ip];
		for (int i = 0; i < 9; i++)boxes[i] = e.boxes[perm[i]];
		for (int i = 0; i < 6; i++) {
			band_stack_order[i] = e.band_stack_order[permbs[i]];
			memcpy(rel_band_stack[i], e.rel_band_stack[permbs[i]],
				sizeof rel_band_stack[0]);
		}
	}
	void MoveMinirows( int ib, int ip) {
		int *permbs = t_perms_boxes_rc[ip];
		register int b = boxes[ib];
		boxes[ib] = 0;
		switch (permbs[0]) {
		case 0:boxes[ib] |= b & 7; break;
		case 1:boxes[ib] |= (b>>3) & 7; break;
		case 2: boxes[ib] |= (b >> 6) & 7; break;
		}
		switch (permbs[1]) {
		case 0:boxes[ib] |= (b<<3) & 070; break;
		case 1:boxes[ib] |= b & 070; break;
		case 2: boxes[ib] |= (b >> 3) & 070; break;
		}
		switch (permbs[2]) {
		case 0:boxes[ib] |= (b << 6) & 0700; break;
		case 1:boxes[ib] |= (b << 3) & 0700;  break;
		case 2: boxes[ib] |= b & 0700; break;
		}
	}
	void MoveMinicols(int ib, int ip) {
		int *permbs = &t_perms_boxes_rc[ip][3];
		register int b = boxes[ib];
		boxes[ib] = 0;
		switch (permbs[0]) {
		case 0:boxes[ib] |= b & 0111; break;
		case 1:boxes[ib] |= (b >> 1) & 0111; break;
		case 2: boxes[ib] |= (b >> 2) & 0111; break;
		}
		switch (permbs[1]) {
		case 0:boxes[ib] |= (b << 1) & 0222; break;
		case 1:boxes[ib] |= b & 0222; break;
		case 2: boxes[ib] |= (b >> 1) & 0222; break;
		}
		switch (permbs[2]) {
		case 0:boxes[ib] |= (b << 2) & 0444; break;
		case 1:boxes[ib] |= (b << 1) & 0444;  break;
		case 2: boxes[ib] |= b & 0444; break;
		}

	}

	void MorphStackToPerm( int stack, int * perm){
		//cout << "stack to perm " << perm[0] << perm[1] << perm[2] << endl;
		uint32_t * rstack = rel_band_stack[3 + stack], rrstack[3];
		memcpy(rrstack, rstack, sizeof rrstack);
		for (int ib = stack; ib < 9; ib += 3) {
			int rb = box_col_to_row[boxes[ib]],b=0;
			for (int ir = 0; ir < 3; ir++) { //rows after diag transpose
				int rr = perm[ir];
				rstack[ir] = rrstack[rr];
				b |= ((rb>>(3*rr)) & 7)<<(3*ir);
			}
			boxes[ib]= box_col_to_row[b];
		}
	}
	void MorphBandToPerm(int band, int * perm) {
		int boxd = band * 3, boxf = boxd + 3;
		uint32_t * rband = rel_band_stack[band], rrband[3];
		memcpy(rrband, rband, sizeof rrband);
		//cout << band << " band to perm " << perm[0] << perm[1] << perm[2]
		//	<<"\trel rows\t"<< rrband[0] << rrband[1] << rrband[2] << endl;
		for (int ib = boxd; ib < boxf; ib ++) {
			int rb = boxes[ib], b = 0;
			for (int ir = 0; ir < 3; ir++) { //mini rows  
				int rr = perm[ir];
				rband[ir] = rrband[rr];
				b |= ((rb >> (3 * rr)) & 7) << (3 * ir);
			}
			boxes[ib] = b;
		}
		//cout << band<< " band final rel rows\t" << rband[0] << rband[1] << rband[2] << endl;

	}

	void BandRelativeOrder(int band, int ip) {
		int old[3];
		register uint32_t *o = rel_band_stack[band];
		memcpy(old, o, sizeof old);
		register int *p = t_perms_boxes_rc[ip];
		for (int i = 0;  i < 3; i++) {
			o[i] = old[p[i]];
		}
	}
	void StackRelativeOrder(int band, int ip) {
		int old[3];
		register uint32_t *o = rel_band_stack[band+3];
		memcpy(old, o, sizeof old);
		register int *p = &t_perms_boxes_rc[ip][3];
		for (int i = 0; i < 3; i++) {
			o[i] = old[p[i]];
		}
	}	

	void MovBox5(SYM_SPOT * e, int ip) {// reorder rows 345 columns 345 
		int * perm = t_perms_boxes[ip], *permbs = t_perms_boxes_rc[ip];
		*this = *e; // start with previous status
		MoveMinirows(3, ip);		MoveMinirows(5, ip);
		MoveMinicols(1, ip);		MoveMinicols(7, ip);
		// adjust central band stack relative order
		BandRelativeOrder(1, ip); StackRelativeOrder(1, ip);
		{// and morph central box
			register int b = boxes[4];
			boxes[4] = 0;
			for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
				if (b&bit) boxes[4] |= 1 << perm[i]	;
			}
		}
		//if (1)  Debug("morph b5 after mov");			 
	}
	int GoB5Central(SYM_SPOT * e, int ip) {
		int diag = 0;
		//if (ip == 7) diag = 1;
		if (diag) cout << "entry gob5central" << endl;
		// apply the perm to central band/stack
		MovBox5(e, ip);
		{	//check fit minirows 3 5
			register int b1 = boxes[3], b2 = boxes[5];
			if (_popcnt32(b1 & 7) != _popcnt32(b2 & 0700))return 0;
			if (_popcnt32(b1 & 070) != _popcnt32(b2 & 070))return 0;
			if (_popcnt32(b1 & 0700) != _popcnt32(b2 & 7))return 0;
		}		 
		{	//check fit 1 7
			register int b1 = boxes[1], b2 = boxes[7];
			if (_popcnt32(b1 & 0111) != _popcnt32(b2 & 0444))return 0;
			if (_popcnt32(b1 & 0222) != _popcnt32(b2 & 0222))return 0;
			if (_popcnt32(b1 & 0444) != _popcnt32(b2 & 0111))return 0;
		}		 
 		if (diag) cout << "try align " << endl;

		// Align minicol 1 to minicol 9 in band 3
		int db3 = box_col_to_row[boxes[3]],
			db5 = box_col_to_row[boxes[5]],
			db3_mini1 = db3 & 7,
			db3_mini2 = (db3 >> 3) & 7;// mini col 1;2 in box 3
		for (int ip5 = 0; ip5 < 6; ip5++) {// box5 matches to box3
			int * tp = tperm6[ip5], //column order in box 5
				cr2 = tp[1], cr3 = tp[2];
			int mini3 = treverse_mini[(db5 >> (3 * cr3)) & 7];
			if (mini3 != db3_mini1) continue;
			int mini2 = treverse_mini[(db5 >> (3 * cr2)) & 7];
			if (mini2 != db3_mini2) continue;
			int db1_mini1 = boxes[1] & 7,
				db1_mini2 = (boxes[1] >> 3) & 7;// mini row 1;2 in box 1
			if (diag) cout << "band match ip="<<ip5 << endl;

			for (int ip5s = 0; ip5s < 6; ip5s++) {//box7 matches to box1
				int * tps = tperm6[ip5s], //row order in box 7
					rr2 = tps[1], rr3 = tps[2];
				int minis3 = treverse_mini[(boxes[7] >> (3 * rr3)) & 7];
				if (minis3 != db1_mini1) continue;
				int minis2 = treverse_mini[(boxes[7] >> (3 * rr2)) & 7];
				if (minis2 != db1_mini2) continue;
				if (diag) cout << "stack match ip=" << ip5s << endl;
				// morph stack3 to perm b5 band 3 to perm b7
				SYM_SPOT * sn = this; sn++;
				*sn = *this;
				if (sn->GoCentralBand3Stack3(tps, tp))// tps for band tp for stack
					return 1;
			}
		}
		return 0;
	}
	int GoCentralBand3Stack3(int *tpb, int *tps) {// morph band3 stack3 final check
		//if (1)Debug("entry GoCentralBand3Stack3");
		MorphBandToPerm(2, tpb);
		MorphStackToPerm(2, tps);
		//if(1)Debug("morph b5 band/stack morphed");
		// the four corners must be ok (central symmetry)
		for (int ib = 0; ib < 3; ib += 2) {// pairing 0 8 and 2 6
			int b1 = boxes[ib], b2 = boxes[8 - ib];
			if((b1 & 7)!= treverse_mini[(b2 >> 6) & 7]) return 0;
			if (((b1>>3) & 7) != treverse_mini[(b2 >> 3) & 7]) return 0;
			if (((b1 >> 6) & 7) != treverse_mini[(b2 ) & 7]) return 0;
		}
		MoveFinal();
		return 1;
	}

	int IsNotDiagonal(int ib, int ip) {
		int * p = diag_remapping_index[ip], v = 0, b = boxes[ib];
		for (int i = 0, bit = 1; i < 9; i++, bit <<= 1) {
			int bit2 = 1 << p[i];
			if (bit2&b) v |= bit;
		}
		return (v != box_col_to_row[v]);
	}
	int GoMainDiagDet(SYM_SPOT * e, int ip0, int ip4, int ip8) {//morph bans and stacks
		*this = *e;
		// morph bands and stacks to ip0 ip4 ip8
		MorphBandToPerm(0, tpermsdiag[ip0]);
		MorphBandToPerm(1, tpermsdiag[ip4]);
		MorphBandToPerm(2, tpermsdiag[ip8]);
		MorphStackToPerm(0, tpermsdiag[ip0]);
		MorphStackToPerm(1, tpermsdiag[ip4]);
		MorphStackToPerm(2, tpermsdiag[ip8]);
		// return 0 if now not diagonal		
		if (boxes[1] != box_col_to_row[boxes[3]])return 0;
		if (boxes[2] != box_col_to_row[boxes[6]])return 0;
		if (boxes[5] != box_col_to_row[boxes[7]])return 0;
		MoveFinal();
		return 1;
	}
	int GoMainDiag() {
		int diag = 0;
		//if (ip == 7) diag = 1;
		if (diag)this->Debug("entry gob5 main diag" );
		// box 0;8 in diagonal mode
		for (int ip0 = 0; ip0 < 6; ip0++) {// box0 must be diagonal 
			if (IsNotDiagonal(0, ip0)) continue;
			for (int ip4 = 0; ip4 < 6; ip4++) {// box0 must be diagonal 
				if (IsNotDiagonal(4, ip4)) continue;
				for (int ip8 = 0; ip8 < 6; ip8++) {// box0 must be diagonal 
					if (IsNotDiagonal(8, ip8)) continue;
					if (diag) cout << "try diagonal " << endl;
					SYM_SPOT * sn = this; sn++;
					if (sn->GoMainDiagDet(this, ip0, ip4, ip8))
						return 1;
				}
			}
		}
		return 0;
	}

	void MoveFinal();
	BF128 BuilBF128() {
		BF128 w;	w.bf.u32[3]=0;
		for (int ib = 0; ib < 3; ib++) {
			uint32_t b[3],band=0; 
			memcpy(b, &boxes[3 * ib], sizeof b);
			for (int ir = 0, shift = 0; ir < 3; ir++) for (int ibx = 0; ibx < 3; ibx++){
				uint32_t x = b[ibx] & 7;
				b[ibx] >>= 3;// clear the mini row moved
				band |= x << shift;
				shift += 3;
			}
			w.bf.u32[ib] = band;
		}
		return w;
	}

	void BuilSwappingTable(int * t) {
		for (int ib = 0; ib < 3; ib++) {
			int band = band_stack_order[ib];
			for (int irow = 0; irow < 3; irow++) {
				int row = 3 * band + rel_band_stack[ib][irow];
				for (int is = 0; is < 3; is++) {
					int stack = band_stack_order[is+3];
					for (int icol = 0; icol < 3; icol++) {
						int col = 3 * stack + rel_band_stack[is+3][icol];
						t[27 * ib + 9 * irow + 3 * is + icol] = 9 * row + col;
					}
				}
			}
		}
	}

	void Debug(const char * lib) {
		char ws[82];
		BF128 w = BuilBF128();
		cout << w.String3X(ws) << " "<<lib << endl;
	}
}sym_spot_final;
void SYM_SPOT::MoveFinal() { sym_spot_final = *this; }



struct PUZC_SYM {
	SYM_SPOT tsym_spot[5]; 
	char puz[82],puz2[82], puz3[82]; // normalized puzzle
	uint32_t boxes[9], wboxes[9],// box pattern //start, after perm
		ctb[9], wctb[9], // count per box  
		ctbs[6],//band stack count
		nclues,cb_parity,sym_boxes;
	void BuildRemappingIndex();
	void Init(char * ze);
	int GetSym();
	int MorphPatB5(PUZC_SYM & pat,int diag=0);
	int MorphPatB1S1B2S2(PUZC_SYM & pat,int * pb0,int * ps0);

};
void PUZC_SYM::BuildRemappingIndex() {
	// build diag remapping index
	for (int iperm = 0; iperm < 6; iperm++) {// six perms box1 from stack1
		// reorder boxes and build pattern
		int *bs_order = tpermsdiag[iperm], *bs_cells = diag_remapping_index[iperm];
		for (int ib = 0, i = 0; ib < 3; ib++) {
			int ibo = bs_order[ib];
			for (int is = 0; is < 3; is++, i++) {
				int iso = bs_order[3 + is],
					io = 3 * ibo + iso;
				bs_cells[i] = io;
			}
		}
	}
}
void PUZC_SYM::Init(char * ze) {
	nclues = 0;
	memset(boxes, 0, sizeof boxes);
	memset(ctb, 0, sizeof ctb);
	strncpy(puz, ze, 81);
	for (int i = 0; i < 81; i++) {
		if (puz[i] < '1' || puz[i] > '9') {
			puz[i] = '.';
			continue;
		}
		boxes[C_box[i]]|= 1<< C_box_rel[i];
		ctb[C_box[i]]++;
		nclues++;
	}
	cb_parity = nclues & 1;
	tsym_spot[0].Init(boxes); 	
}
int PUZC_SYM::GetSym() {
	tsym_spot[1] = tsym_spot[0];
	int irs = 0, irsym = 0;
	SYM_SPOT  *so = &tsym_spot[1], *sn = so + 1;
	for (int iperm = 0; iperm < 9; iperm++) {// start with central symmetry
	// reorder boxes and build pattern
		for (int i = 0; i < 9; i++) {
			wboxes[i] = boxes[t_perms_boxes[iperm][i]];
			wctb[i] = ctb[t_perms_boxes[iperm][i]];
		}
		// central box must have the right parity to be central
		if (cb_parity != (wctb[4] & 1)) continue;
		for (int i = 0; i < 5; i++)	if (wctb[i] != wctb[8 - i])continue;
		// build spot1 to do more checks
		memcpy(so->band_stack_order, t_perms_boxes_rc[iperm],
			sizeof so->band_stack_order);
		memcpy(so->boxes, wboxes, sizeof wboxes);
		int v0[9], v[9];
		{	register int B = wboxes[4];// central box at start
			for (int i = 0; i < 9; i++, B >>= 1) v0[i] = B & 1;
		}
		for (int iperm = 0; iperm < 9; iperm++) {// each cell central
			int *p = t_perms_boxes[iperm];
			for (int i = 0; i < 9; i++)v[i] = v0[p[i]];
			if (v[4] != cb_parity)continue; // must be if central
			for (int i = 0; i < 4; i++)	if (v[i] != v[8 - i])continue;
			if (sn->GoB5Central(so, iperm))// apply the perm to central band/stack
				return 1;
		}
	}
	//try now diagonal 6 perm to have box1 b0 or b3 or b6 top left
	//cout << "trydiagonal" << endl;
	for (int iperm = 0; iperm < 6; iperm++) {// six perms box1 from stack1
		// reorder boxes and build pattern
		int *bs_order = tpermsdiag[iperm],
			*bs_cells = diag_remapping_index[iperm];
		for (int i = 0; i < 9; i++) {
			int io = bs_cells[i];
			wboxes[i] = boxes[io];
			wctb[i] = ctb[io];
		}
		if (wctb[1] != wctb[3])continue;
		if (wctb[2] != wctb[6])continue;
		if (wctb[5] != wctb[7])continue;
		*sn = *so;
		memcpy(sn->boxes, wboxes, sizeof wboxes);
		memcpy(sn->band_stack_order, bs_order, sizeof sn->band_stack_order);
		if (sn->GoMainDiag()) return 2;
	}
	return 0;
}
void Go_c490() {
	cerr << "Go_490 entry " << sgo.finput_name << " extract potential symmetry" << endl;
	PUZC_SYM psym;
	psym.BuildRemappingIndex();
	char * ze = finput.ze;
	finput.open(sgo.finput_name);
	if (!finput.is_open()) {
		cerr << "error open file " << sgo.finput_name << endl;
		return;
	}
	int diag = 0;
	while (finput.GetLigne()) {
		if(diag)cout << "go" << endl;
		if (strlen(ze) < 81) continue;
		if (diag)cout << ze << endl;
		psym.Init(ze);
		psym.sym_boxes = psym.GetSym();
		if (psym.sym_boxes) {
			char ws[82];
			BF128 w = sym_spot_final.BuilBF128();
			if (diag)cout << w.String3X(ws) << "morphed final to" << endl;
			int rbs[6][3],//relatice order (row/col) inside the band/stack
				bso[6]; // bands/stack/reordering
			memcpy(rbs, sym_spot_final.rel_band_stack, sizeof rbs);
			memcpy(bso, sym_spot_final.band_stack_order, sizeof bso);
			// reorder the source
			if (diag)cout << "band_stack_order " << bso[0] << bso[1] << bso[2]
				<< "\t" << bso[3] << bso[4] << bso[5] << endl;
			if (diag) {
				cout << " relative band stack order" << endl;
				for (int i = 0; i < 6; i++) {
					int *w = rbs[i];
					cout << "\t" << w[0] << w[1] << w[2] << "\t";
				}
				cout << endl << "reshaped entry" << endl;
			}
			for (int ib = 0; ib < 3; ib++) {
				int band = bso[ib];
				int * wib = rbs[ib];
				for (int ir = 0; ir < 3; ir++) {
					int row = 3 * band + wib[ir];
					for (int is = 0; is < 3; is++) {
						int stack = bso[3 + is];
						int * wis = rbs[3+is];
						for (int ic = 0; ic < 3; ic++) {
							int col = 3 * stack + wis[ic];
							if (diag)cout << ze[9 * row + col];
							fout1 << ze[9 * row + col];
						}
					}
					if (diag)cout << endl;
				}
			}
			ze[81] = 0;
			fout1<<";" << ze << ";" << psym.sym_boxes << endl;
		}

		//strcpy(&zout[n], &ze[104]);
		//fout1 << zout << endl;
	}

}

int PUZC_SYM::MorphPatB1S1B2S2(PUZC_SYM & pat, int * pb0, int * ps0) {
	SYM_SPOT  &s2 = tsym_spot[2], &s3 = sym_spot_final;
	for (int ipb1 = 0; ipb1 < 6; ipb1++)
		for (int ips1 = 0; ips1 < 6; ips1++) {
			int * pb1 = tperm6[ipb1], *ps1 = tperm6[ips1];
			s2 = tsym_spot[1];
			s2.MorphBandToPerm(1, pb1);
			s2.MorphStackToPerm(1, ps1);
			uint32_t * bx = s2.boxes, *bpx = pat.boxes;
			if ((bpx[4]!= bx[4])||( bpx[1]!= bx[1])||(bpx[3]!= bx[3])) continue;
			for (int ipb2 = 0; ipb2 < 6; ipb2++)
				for (int ips2 = 0; ips2 < 6; ips2++) {
					int * pb2 = tperm6[ipb2], *ps2 = tperm6[ips2];
					s3 = s2;
					s3.MorphBandToPerm(2, pb2);
					s3.MorphStackToPerm(2, ps2);
					uint32_t * bx2 = s3.boxes;
					if (bpx[8]!= bx2[8]) continue;
					if (bpx[2]!= bx2[2]||bpx[5]!= bx2[5]) continue;
					if (bpx[6]!= bx2[6]||bpx[7]!= bx2[7]) continue;
					return 1;
				}
		}

	return 0;
}


int PUZC_SYM::MorphPatB5(PUZC_SYM & pat,int diag) {
	for (int iperm = 0; iperm < 9; iperm++) {// start with central symmetry
	// reorder boxes and build pattern
		for (int i = 0; i < 9; i++) {
			wboxes[i] = boxes[t_perms_boxes[iperm][i]];
			wctb[i] = ctb[t_perms_boxes[iperm][i]];
		}
		// central box count must be ==
		if (pat.ctb[4] != wctb[4]) continue;
		// do the box perm for the given entry
		strcpy(puz2, empty_puzzle);
		for (int i = 0; i < 9; i++)  {// box 4 not moved
			int i0 = t_perms_boxes[iperm][i];
			byte *bp = cellsInGroup[18 + i], *bp0 = cellsInGroup[18 + i0];
			for (int ip = 0; ip < 9; ip++) {
				puz2[bp[ip]] = puz[bp0[ip]];
			}
		}
		if (diag) {
			SYM_SPOT  &s = tsym_spot[1];
			tsym_spot[1] = tsym_spot[0];
			memcpy(s.boxes, wboxes, sizeof wboxes);
			s.Debug("z2 control");
			cout << puz2 <<"puz2"<< endl;
		}
		//===========================================================
		// 8 perms four boxes in top left + diag symmetru
		int wb2[9], wctb2[9];
		int p2[8][9] = {
			{0,1,2,3,4,5,6,7,8},
			{0,3,6,1,4,7,2,5,8},
			{2,1,0,5,4,3,8,7,6},
			{6,3,0,7,4,1,8,5,2},
			{8,7,6,5,4,3,2,1,0},
			{8,5,2,7,4,1,6,3,0},
			{6,7,8,3,4,5,0,1,2},
			{2,5,8,1,4,7,0,3,6}
		};
		for (int ip2 = 0; ip2 < 8; ip2++) {
			int * pp = p2[ip2];
			for (int i = 0; i < 9; i++) {
				wb2[i] = wboxes[pp[i]];
				wctb2[i] = wctb[pp[i]];
				if (wctb2[i] != pat.ctb[i]) goto nextip2;
			}
			strcpy(puz3, puz2);// copy unchanged
			for (int i = 0; i < 9; i++) if(i-4){// box 4 not moved
				int i0 = pp[i];
				byte *bp = cellsInGroup[18+i],*bp0= cellsInGroup[18 + i0];
				for (int ip = 0; ip < 9; ip++) {
					puz3[bp[ip]] = puz2[bp0[ip]];
				}
			}
			if (diag) {
				SYM_SPOT  &s = tsym_spot[1];
				tsym_spot[1] = tsym_spot[0];
				memcpy(s.boxes, wb2, sizeof wb2);
				s.Debug("z3 control");
				cout << puz3<<"puz3" << endl;
			}
			//global count match, try more 
			for (int ipb0 = 0; ipb0<6;ipb0++)
				for (int ips0 = 0; ips0 < 6; ips0++) {
					int * pb0 = tperm6[ipb0], *ps0 = tperm6[ips0];
					SYM_SPOT  &s = tsym_spot[1];
					tsym_spot[1] = tsym_spot[0];
					memcpy(s.boxes, wb2, sizeof wb2);
					s.MorphBandToPerm(0, pb0);
					s.MorphStackToPerm(0, ps0);
					if (pat.boxes[0]!= s.boxes[0]) continue;
					if (MorphPatB1S1B2S2(pat,pb0,ps0)) return 1;
				}
		nextip2:	;
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
	psym.BuildRemappingIndex();
	psym_pat.Init(sgo.s_strings[0]);
	for (int i = 0; i < 9; i++) cout << psym_pat.ctb[i];
	cout << "  compte par boite" << endl;

	char  ze[82]; ze[81] = 0;
	while (finput.GetPuzzle(ze)) {
		psym.Init(ze);// create boxes 
		if (psym.nclues != psym_pat.nclues)continue; // invalid entry
		// start with b5 (as in c490 central) count = b5 pat
		//fout2 << ze << endl;
		if (psym.MorphPatB5(psym_pat)) {
			//cout << "retour ok" << endl;
			int t[81];
			sym_spot_final.BuilSwappingTable(t);
			char wout[82]; 
			strcpy(wout, empty_puzzle);
			for (int i = 0; i < 81; i++) {
				if (t[i] < 0 || t[i]>80)continue;
				int c = psym.puz3[t[i]];
				wout[i] = c; 			}
			fout1 << wout  << endl;
		}
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