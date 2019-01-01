
//========================================
const char * zh_g_cpt[10] = { "npuz", "guess", "close_d ", "upd1 ", "upd2 ",
"fupd ", "hpair ", "htripl ", " ", " " };

const char * libs_c17_00_cpt2g[20] = {
	"bands1+2 processed entries M10",//0
	"chunks processed entries gochunk",//1
	"XY count",//2
	"XY passing UA filter",//3
	"XY brute force",//4
	"valid brute force",//5
	"more sockets2 searched ",//6
	"more sockets2 searched found",//7
	"valid b12 after build active",//8
	"band3 active after build active",//9
	"",//10
	"",//11
	"",//12
	"",//13
	"",//14
	"",//15
	"",//16
	"",//17
	"",//18
	"control valib12// band3",//19
};
void Go_c17_00( ) {// p2 process
	cout << "Go_c17_00 search batch 17 clues 656 566 " << endl;
	cout << sgo.vx[0] << " -v0- band 0_415" << endl;
	cout << sgo.vx[2] << " -v2- skip first nnn restart after batch failure" << endl;
	cout << sgo.vx[3] << " -v3- last entry number for this batch must be > vx[2]" << endl;
	int it16_start = sgo.vx[0];
	genb12.skip = sgo.vx[2];
	genb12.last = sgo.vx[3];
	if (sgo.vx[2] < 0) {
		cerr << "invalid value for skip" << endl;
		return;
	}
	if (sgo.vx[3] < sgo.vx[2]) {
		cerr << "invalid value for last to process" << endl;
		return;
	}
	memset(p_cpt, 0, sizeof p_cpt);// band2 and band 3 count
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	memset(p_cpt1g, 0, sizeof p_cpt1g);// used in debugging sequences only
	memset(p_cpt2g, 0, sizeof p_cpt2g);// used in debugging sequences only
	genb12.Start(0);
	genb12.NewBand1(sgo.vx[0]);
	cout << "print final stats" << endl;
	for (int i = 0; i < 10; i++) {
		if (!p_cpt2g[i])continue;
		cout << p_cpt2g[i] << "\t\t" << libs_c17_00_cpt2g[i] << endl;
	}
}
//=========================== regressive tests

void Go_c17_91_go() {
	if (p_cptg[0] > 1) return;
	myband1.PrintStatus();
	myband2.PrintStatus();
	genuasb12.Initgen();
}
void Go_c17_80() {// enumeration test
	cout << "Go_c17_20 phase 2a enumeration test " << endl;
	cout << sgo.vx[0] << " -v0- first id 0_415" << endl;
	cout << sgo.vx[1] << " -v1- second id 0_415" << endl;
	int it16_start = sgo.vx[0], it16_end = sgo.vx[1];
	genb12.skip = 0;
	genb12.last = 100000;

	if (it16_start > 415 || it16_end > 415 || it16_start > it16_end) {
		cerr << "invalid it16_start it16_end" << endl;
		return;
	}
	memset(p_cptg, 0, sizeof p_cptg);// used in debugging sequences only
	genb12.Start(11);
	for (int i1t16 = it16_start; i1t16 <= it16_end; i1t16++) {
		memset(p_cpt, 0, sizeof p_cpt);// band2 and band 3 count
		genb12.nb12 = 0;
		genb12.NewBand1(i1t16);
		cout << genb12.i1t16 << "\t" << genb12.it16 << "\t" << p_cpt[0]
			<< "\t" << p_cpt[1] << endl;
		p_cptg[0] += p_cpt[0];
		p_cptg[1] += p_cpt[1];
	}
	cout << "total\t\t" << p_cptg[0]
		<< "\t" << p_cptg[1] << endl;
}
void Go_c17_90() {// creating band UAs
	//the table exists, this is somehow a regressive test  
	//or to rebuild the table
	cout << "Go_10 entry  creating band uas test" << endl;
	STD_B416 wband;
	wband.Initstd();
	int maxua = 0;
	zh1b_g.diag = 0;
	for (int i416 = 0; i416 < 416; i416++) {
		//for (int i416 = 64; i416 <= 64; i416++) {	
		wband.InitC10(i416);
		int nua = t16_nua[i416]; // reference known
		zhone[0].InitOne_std_band();
		cout << wband.band << " i416=" << i416 << endl;

		for (int i = 0; i < 36; i++) {// find UAs 2 digits
			if (zhone[0].Start_nFloors(floors_2d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess2();
		}
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			if (zhone[0].Start_nFloors(floors_3d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess3();
		}
		if (t16_nua[i416] == zh1b_g.nua)goto ok;
		for (int i = 0; i < 126; i++) {// find UAs 4 digits
			if (zhone[0].Start_nFloors(floors_4d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess4();
		}
		for (int i = 0; i < 126; i++) {// find UAs 5 digits
			if (zhone[0].Start_nFloors(0x1ff ^ floors_4d[i])) continue;
			zhone[0].InitGuess();
			if (0 && floors_4d[i] == 0113) {
				zhone[0].ImageCandidats();
				zh1b_g.diag = 1;
			}
			else zh1b_g.diag = 0;
			zhone[0].Guess5();
		}
		for (int i = 0; i < 84; i++) {// find UAs 6 digits
			if (zhone[0].Start_nFloors(0x1ff ^ floors_3d[i])) continue;
			zhone[0].InitGuess();
			zhone[0].Guess6();
		}
		if (t16_nua[i416] == zh1b_g.nua)goto ok;

		//if(1)zh1b_g.PrintTua();

		zh1b_g.FindMissingUAs();

		if (t16_nua[i416] == zh1b_g.nua)goto ok;
		{
			cout << "uas final table" << endl;
			zh1b_g.PrintTua();
			cout << "uas known table" << endl;
			int id = t16_indua[i416], iff = id + t16_nua[i416];
			for (int i = id; i < iff; i++)
				cout << i - id << "\t"
				<< Char27out(t16_UAs[i]) << "\t" << _popcnt32(t16_UAs[i]) << endl;
			continue;
		}
	ok:
		cout << "ua count equal" << endl;
	}
}
void Go_c17_91() {// UAs collector 2 bands
	cout << "Go_11 test UAs collector 2 ba,ds" << endl;
	cout << sgo.vx[0] << " -v0- id bande low order 0_415" << endl;
	cout << sgo.vx[0] << " -v1- number of band2" << endl;
	memset(p_cpt, 0, sizeof p_cpt);// used in debugging sequences only
	memset(p_cpt1, 0, sizeof p_cpt1);// used in debugging sequences only
	genb12.Start(10);
	genb12.NewBand1(sgo.vx[0]);
}
/* sample used for Go_c17_92()
1 0     0       0       0       0       729     18225
28 27   0       38      2025    26285   157845  610977
29 28   9       567     7398    51516   237762  803574
412 411 0       0       0       324     10503   112491
413 412 0       0       0       0       729     18225
32 31   0       0       0       546     12291   101597
		*/
void Go_c17_92() {// test UAs 5 6 expand
	int tb1[6] = { 0,27,28,411,412,31 };
	int tb5[6] = { 0,26285,51516,324,0,546 };
	int tb6[6] = { 729,157845,237762,10503,729,12291 };


	cout << "Go_c17_92 test expand" << endl;
	STD_B1_2 wband;
	wband.Initstd();
	int maxua = 0;
	zh1b_g.diag = 0;
	for (int i = 0; i < 6; i++) {
		int i416 = tb1[i];
		wband.InitG12(i416);
		int nua = t16_nua[i416]; // reference known
		cout << wband.band << " i416=" << i416 
			<<"nua="<<nua<<"\t"<< t16_nua[i416]	<< endl;
		wband.DoExpandBand(0);
		cout << "ua count equal" << endl;
		cout << "n5\t" << tb5[i] << "\t" << wband.n5<< endl;
		cout << "n6\t" << tb6[i] << "\t" << wband.n6 << endl;
		//wband.DebugIndex(0);
		//wband.DebugIndex(1);
	}
}
