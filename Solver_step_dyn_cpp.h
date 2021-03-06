﻿/* Dynamic plus mode 90 + length
   same as dynamic but including several "objects" ALS patterns
   all "events"  claiming, pointing, pair, hidden pair, XWing


*/
#define SKMPPV2
ZHOU zdyn,zdynstep,zhou_solve_start;
ZH_GLOBAL zhgdyn;
void BuildHiddenBiv_Xwings(PM3X& pm, PM3X& pmdiag, HID_BIV & rbiv, HID_BIV & xbiv,int diag=0) {
	// phase 1 set up bi values per digit
	uint32_t dig_bivs[9];
	memset(dig_bivs, 0, sizeof dig_bivs);
	memset(&rbiv, 0, sizeof rbiv);// 36 pairs of digits
	memset(&xbiv, 0, sizeof xbiv);// 36 pairs of rows/cols 9 digits
	for (int iu = 0; iu < 27; iu++) {
		int bit = 1 << iu;// bit to set
		BF128 w0 = units3xBM[iu];
		for (int idig = 0; idig < 9; idig++) {
			BF128 w = w0 & pm.pmdig[idig];
			if (w.Count() == 2)dig_bivs[idig] |= bit;
		}
	}
	//cout << "pairs biv" << endl;
	// phase 2 check pairs of digits biv same unit
	for (int i1 = 0, i12 = 0; i1 < 8; i1++) {
		for (int i2 = i1 + 1; i2 < 9; i2++, i12++) {
			uint32_t sets_com = dig_bivs[i1] & dig_bivs[i2];
			if (!sets_com)continue;
			int tsets[30], ntsets=0;
			BitsInTable32(tsets, ntsets, sets_com);// potential  hidden pairs
			for (int iu = 0; iu < ntsets; iu++) {// must be same cells
				int unit = tsets[iu];
				BF128 w = units3xBM[unit];// cells of the unit
				w &= (pm.pmdig[i1] & pm.pmdig[i2]);
				if (w.Count() == 2) {// this is a hidden bi value i1;i2
					rbiv.sets_biv[i12] |= 1 << unit;
				}
			}
		}
	}
	if(diag) cout << "xwings" << endl;
	// same rows columns for xwings
	for (int i1 = 0, i12 = 0; i1 < 8; i1++) {
		int i1b = i1 / 3, i1shift = 9 * (i1 % 3);
		for (int i2 = i1 + 1; i2 < 9; i2++, i12++) {
			int i2b = i2 / 3, i2shift = 9 * (i2 % 3);
			int mask = floors_2d[i12];// rows or cols
			for (int idig = 0; idig < 9; idig++) {				
				int dbiv = dig_bivs[idig];
				if (_popcnt32(dbiv&mask) == 2) {// can be xwing row
					int r1 = (pm.pmdig[idig].bf.u32[i1b] >> i1shift) & 0x1ff;
					int r2 = (pm.pmdig[idig].bf.u32[i2b] >> i2shift) & 0x1ff;
					if (r1 == r2) xbiv.sets_biv[i12] |= (1 << idig);
				}
				dbiv >>= 9;// now columns in the low 8 bits
				if (_popcnt32(dbiv&mask) == 2) {// can be xwing column
					int r1 = (pmdiag.pmdig[idig].bf.u32[i1b] >> i1shift) & 0x1ff;
					int r2 = (pmdiag.pmdig[idig].bf.u32[i2b] >> i2shift) & 0x1ff;
					if (r1 == r2) xbiv.sets_biv[i12] |= (1 << (idig+9));
				}
			}
		}
	}
}
struct XYSEARCHDYNPLUS{
	BF128 pairs, pairs_step_old,
		cells_biv_all, cells_all, cells_biv_true, dig_b_true[9];
	BF128 loop, wb;
	BF128 cells_unsolved_live;
	int dig_cells_live[81];
	struct PATH {
		GINT64  t[400];// u8 cell,digit,sign,step u16[2] source 
		int nt, dig, cell;
	}paths[9];

	int idig, digit,  cell, sign,ddig,dcell,dsign,dind,
		nt,ntd,ntp, maxpas,maxrating, npaths,elim_done,
		nsteps, c1, c2, 
		locdiag, diagexpand,
		mode,fastmode,expandmode,opprint;
	int dig_sets[9][27];
	uint32_t d1,d2;
	HID_BIV hp_biv_0, hp_biv_old, hp_biv_step, //hidden pairs status (d1 d2 unit)
		   xw_biv_0, xw_biv_old, xw_biv_step;  // xwing status (rc1 rc2 digit)
	BF32 dig_bivsets[9],dig_sets3[9];
	GINT64 tback[300],t[400];// cell,digit,digit2,source
	PM3X //used_off_digits,  
		used_off_old,used_off_new,
		//used_on_digits, 
		cleang, cleanstart,
		active_all,active_unit,dbiv;
	GINT telims[50];	int ntelims ;
	//====== for dynamic and more
	int is_contradiction; // direct index for storing tables
	PM3X off_off_status;
	//PM3X off_status[400],contradiction,elim_stored;
	//GINT64 off_path[100][400];
	//int tex[9][81]; // index to first "on" index in the path
	void SearchInit(int fast);
	void Addt(int cell, int di, int sign, int source);
	inline void AddLastInCell(int ce, int di){	Addt(ce, di, 0, 0x1000 | cell);	}	
	inline void AddLastInUnit(int ce, int di,int unit){	Addt(ce, di, 0, 0x2000 | unit);	}
	inline int DynLock1(int band) {
		int shrink = (TblShrinkMask[band & 0x1FF]) |
			(TblShrinkMask[(band >> 9) & 0x1FF] << 3) |
			(TblShrinkMask[(band >> 18) & 0x1FF] << 6);
		return band ^ (band &TblComplexMask[shrink]);
	}
	void OffToOn_Plus(int i);
	void OnToOff_Plus(int i);
	void OffToOff_Plus();
	void OffToOff_Plus_New_Pairs_Status(BF128 *pmz, BF128 &wp_new);
	int Is_To_Clean(int rating);
	int AddElim(int d, int c, int rat);
	int AddElimFast(PM3X & elims);
	void Clean_d_bf(int d, BF128 & b,int source);
	void Clean_digs_cell(int digs, int cell,int source);
	int SearchDynPlus(int fast);
	void SearchDynPassPlus(int nmax);
	void SearchDynPassMulti(int nmax);
	void ExpandPlusInit(GINT cand,int off=0);
	void ExpandPlusGo();
	int ExpandDynamicToElim(GINT cand,GINT targ);
	int BackDynamic(GINT64 target, GINT64 * tb, int ntb);
	int BackDynamicOff(GINT target);
	void DynamicSolveContradiction(GINT cand,PM3X cont);
	void DynamicSolveContradiction(int dig1, int cell1, int dig2, int cell2, PM3X cont);
	void PrintTback();
	void PrintBackMulti(int elim_dig, int elim_cell);
	void PrintBackCom(const char * lib,GINT64 * ptback, int nback,int mode );
	void DebugT();
}xysearchplus;

//============================================= XYSEARCHDYNPLUS
void XYSEARCHDYNPLUS::SearchInit(int fast) {
	fastmode = fast;
	zhou_solve_start = zhou_solve;
	expandmode = 0;// start search no way back
	elim_done = ntelims = 0;
	maxrating = 200;
	pm_go.xysearch.elim_stored.SetAll_0();
	if(opprint)cout << "start build " << endl;
	BuildHiddenBiv_Xwings(zh_g2.pm, zh_g2.pmdiag, hp_biv_0, xw_biv_0);
	if (opprint)cout << "end build " << endl;
}

void XYSEARCHDYNPLUS::Addt(int ce , int di, int sign, int source) {
	if (nt > 380) return;
	GINT64 & tx = t[nt++];
	tx.u32[0] = ce + (di << 8) + (sign << 16) + (nsteps << 24);
	tx.u32[1] = source;
	if (diagexpand) {
		cout << "Addt ";
		if (sign) cout << "~";
		cout << di + 1 << cellsFixedData[ce].pt << endl;
	}
	if (sign) {
		pm_go.xysearch.used_off_digits.Set_c(di,ce);
		dig_cells_live[ce] &= ~(1 << di); 
		zdyn.ClearCandidate_c(di, ce);
	}
	else pm_go.xysearch.used_on_digits.Set_c(di, ce);
}
void XYSEARCHDYNPLUS::Clean_d_bf(int d, BF128 & b,int source) {
	int tc[80], ntc = b.Table3X27(tc);
	for (int i = 0; i < ntc; i++)Addt(tc[i], d, 1, source);
}
void XYSEARCHDYNPLUS::Clean_digs_cell(int digs, int cell, int source) {
	if (!digs) return;
	int td[10], ntd = 0;
	BitsInTable32(td, ntd,digs);
	for (int i = 0; i < ntd; i++)Addt(cell, td[i], 1, source);
}

void XYSEARCHDYNPLUS::OffToOn_Plus(int isource){

	if (pairs.On_c(cell)){
		int dig = zh_g2.dig_cells[cell] ^ (1 << digit);
		bitscanforward(d2, dig);
		if (pm_go.xysearch.used_on_digits.Off_c(d2, cell)){
			pm_go.xysearch.tex[d2][cell] = nt;// priority to direct
			Addt(cell, d2,0, isource);
			//used_on_digits.Set_c(d2, cell);
		}
		else{// check to offset a last on same step
			int oldi = pm_go.xysearch.tex[d2][cell];
			GINT64 w = t[oldi];
			if( (w.u16[2] & 0x3000)&& oldi>= ntd){// old was not direct same step
				t[oldi].u16[2] = (uint16_t)isource;
			}
		}
	}
	else{// and now last in cell or empty cell
		int ndc = _popcnt32(dig_cells_live[cell]);
		if(ndc>1)goto exit_last_in_cell;// nothing to do
		if (ndc==1){// true last in cell put it on if not yet
			//cout << "add last in cell " << endl;
			bitscanforward(d1, dig_cells_live[cell]);
			if (pm_go.xysearch.used_on_digits.Off_c(d1, cell)){
				pm_go.xysearch.tex[d1][cell] = nt;// first
				AddLastInCell(cell, d1);
				//used_on_digits.Set_c(d1, cell);
			}
		}
		else{// empty cell (no "on") set "on" all "off"  
			int digs = zh_g2.dig_cells[cell] ^ (1 << digit);;
			while ( digs){
				bitscanforward(d2, digs);
				digs ^= 1 << d2;
				if (pm_go.xysearch.used_on_digits.Off_c(d2, cell)){
					pm_go.xysearch.tex[d2][cell] = nt;// first
					AddLastInCell(cell, d2);
					//used_on_digits.Set_c(d2, cell);
				}
			}
		}
	exit_last_in_cell:;
	}
	CELL_FIX & cf = cellsFixedData[cell];
	int tu[3];
	cf.GetRegions(tu);
	for (int iu = 0; iu < 3; iu++){// process each set
		int unit = tu[iu];
		BF128 wu = units3xBM[unit]; wu &= zh_g2.pm.pmdig[digit];
		int  nn_start = wu.Count();
		BF128 w81 = wu - pm_go.xysearch.used_off_digits.pmdig[digit];
		if ((w81&pm_go.xysearch.used_on_digits.pmdig[digit]).isNotEmpty()) continue;// no on
		int nbc = w81.Count(),wcell2;
		if (nbc > 1)continue; // not last in region
		if (nbc == 1){
			wcell2 = w81.getFirstCell();
			if (nn_start == 2){// bi value
				if (pm_go.xysearch.used_on_digits.Off_c(digit, wcell2)){
					pm_go.xysearch.tex[digit][wcell2] = nt;// priority to direct
					Addt(wcell2, digit, 0,isource);
				}
				else{// offset not direct same step
					int oldi = pm_go.xysearch.tex[digit][wcell2];
					GINT64 w = t[oldi];
					if ((w.u16[2] & 0x3000) && oldi >= ntd){// old was not direct same step
						t[oldi].u16[2] = (uint16_t)isource;
					}
				}
			}
			else {//last in region
				if (pm_go.xysearch.used_on_digits.Off_c(digit, wcell2)){
					pm_go.xysearch.tex[digit][wcell2] = nt;// first
					AddLastInUnit(wcell2, digit,unit);
					//used_on_digits.Set_c(digit, wcell2);
				}
				else if(nn_start==3){//offset not direct if smaller size same step
					int oldi = pm_go.xysearch.tex[digit][wcell2];
					GINT64 w = t[oldi];
					if ((w.u16[2] & 0x2000) && oldi >= ntd){// old was not direct same step
						t[oldi].u16[2] = (uint16_t)(0x2000 | unit); 
					}
				}
			}
		}
		else {//empty unit set(no "on)  "on" all  
			wu.Clear_c(cell);
			while ((wcell2 = wu.getFirstCell()) >= 0){
				wu.Clear_c(wcell2);
				if (pm_go.xysearch.used_on_digits.Off_c(digit, wcell2)){
					pm_go.xysearch.tex[digit][wcell2] = nt;// first
					AddLastInUnit(wcell2, digit, unit);
					//used_on_digits.Set_c(digit, wcell2);
				}
			}
		}
	}
}
void XYSEARCHDYNPLUS::OnToOff_Plus(int isource){// no bi value filter
	int digs = zh_g2.dig_cells[cell] ^ (1 << digit);// can be more than one
	while ( digs){
		bitscanforward(d2, digs);
		digs ^= 1 << d2;
		if (pm_go.xysearch.used_off_digits.On_c(d2, cell))continue;
		if ((cell == dcell) && ((int)d2 == ddig)) continue;
		Addt(cell, d2,1, isource);
		//used_off_digits.Set_c(d2, cell);
	}
	CELL_FIX & cf = cellsFixedData[cell];
	BF128 wb = cell_z3x[cell];
	wb &= zh_g2.pm.pmdig[digit];
	wb -= pm_go.xysearch.used_off_digits.pmdig[digit];
	wb.Clear_c(dcell);// be sure not to use it
	//used_off_digits.pmdig[digit] |= wb;
	while ((c2 = wb.getFirstCell()) >= 0){
		wb.Clear_c(c2);
		Addt(c2, digit,1, isource);
	}
}
void XYSEARCHDYNPLUS::DebugT(){
	cout << "\nsummary of expand nt="<<nt << endl;
	int curstep = -1;
	for (int i = 0; i < nt; i++){
		GINT64 w = t[i];
		int wc = w.u16[0], wd = w.u16[1], ws = w.u16[3], source = w.u16[2];
		if (ws != curstep){
			curstep = ws;
			cout << "\nstep " << curstep << "->";
			cout << wd + 1 << cellsFixedData[wc].pt << ";0x" << hex << source << dec << " ";
		}

	}
	cout << endl;
}
int XYSEARCHDYNPLUS::AddElim(int d, int c, int rating){
	if (fastmode || (rating <= pm_go.rat_er)){// then do it
		ntelims = 0;
		maxrating = pm_go.rat_er;
		if (zhou_solve.IsOnCandidate_c(d, c)){
			zhou_solve.ClearCandidate_c(d, c);
			if (opprint)	cout << "immediate elim " << d + 1 << cellsFixedData[c].pt << " rating=" << rating
				<< "  pm_go.rat_er=" << pm_go.rat_er << " nt=" << nt << endl;
			elim_done = 1;
			return 1;
		}
		return 0;
	}
	if (rating > maxrating) return 0;
	if (elim_done)return 0;
	if (opprint)	cout << "addelim " << d + 1 << cellsFixedData[c].pt << " rating=" << rating
		<< "  pm_go.rat_er=" << pm_go.rat_er << " nt=" << nt << " ntelim=" << ntelims << endl;
	if (rating < maxrating){
		maxrating = rating;
		ntelims = 0;
		pm_go.xysearch.elim_stored.SetAll_0();
		if (opprint)cout<< "reset elim_count to 0" << endl;
	}
	if (ntelims < 40 && pm_go.xysearch.elim_stored.Off_c(d,c)){// must continue to look for smaller
		telims[ntelims++].u32 = c | (d << 16);
		pm_go.xysearch.elim_stored.Set_c(d, c);
		return 2;
	}
	return 0;
}
int XYSEARCHDYNPLUS::AddElimFast(PM3X & elims) {
	int iret = 0;
	for (int idig = 0; idig < 9; idig++) {
		BF128 ee = elims.pmdig[idig];
		if (ee.isEmpty())continue;
		int tc[80], ntc = ee.Table3X27(tc);
		for (int ic = 0; ic < ntc; ic++) {
			int c = tc[ic];
			if (zhou_solve.IsOnCandidate_c(idig, c)) {
				zhou_solve.ClearCandidate_c(idig, c);
				if (opprint)	cout << "immediate elim " << idig + 1
					<< cellsFixedData[c].pt << " rating=96" << endl;
				elim_done = iret = 1;
			}
		}
	}
	return iret;
}
//____________________ XY DYNPLUS  base 90
void XYSEARCHDYNPLUS::ExpandPlusInit(GINT cand,int off) {// start with cand on
	zdyn = zhou_solve_start;
	zhgdyn = zh_g;
	// initial status for cells pair hid pairs and xwings
	memcpy(dig_cells_live, zh_g2.dig_cells, sizeof dig_cells_live);
	pairs_step_old = zh_g.pairs;
	hp_biv_old = hp_biv_0;
	xw_biv_old = xw_biv_0;

	ddig = cand.u8[1];
	dcell = cand.u8[0];
	dsign = off;
	dind = cand.u16[1];
	nsteps = is_contradiction = 0;// start with 1 step 
	if (zh_g2.zerobased_sol[dcell] == idig || dsign)
		is_contradiction = 2;// skip test if valid or start off
	//nt = 1;	
	//t[0].u64 = dcell | (ddig << 8)|(dsign<<16);	// source to 0
	pm_go.xysearch.used_off_digits.SetAll_0();
	pm_go.xysearch.used_on_digits.SetAll_0();
	used_off_old.SetAll_0();
	nt = 0;
	Addt(dcell, ddig,dsign, 0);

	//if(dsign)used_off_digits.Set_c(ddig, dcell);
	//else used_on_digits.Set_c(ddig, dcell);
}
void XYSEARCHDYNPLUS::OffToOff_Plus_New_Pairs_Status(BF128 *pmz ,BF128 &wp_new) {
	BF128 wpairs; wpairs.SetAll_0();
	for (int i = 0; i < 81; i++)
		if (_popcnt32(dig_cells_live[i]) == 2) wpairs.Set_c(i);
	wp_new = wpairs - pairs_step_old;
	pairs_step_old = wpairs;

}

void XYSEARCHDYNPLUS::OffToOff_Plus() {
	if (is_contradiction==1 || nt>350) return;
	if (diagexpand) cout << "entry off to off" << endl;
	zdynstep = zdyn;// be sure to avoid cumulative action
	used_off_new = pm_go.xysearch.used_off_digits;
	used_off_new -= used_off_old;
	int tactive[10], nactive = 0,source=0;
	PM3X pm,pmdiag;
	BF128 * pmz=pm.pmdig, wp_new;
	for (int i = 0; i < 9; i++){
		pmz[i] = zdynstep.FD[i][0];
		if (used_off_new.pmdig[i].isNotEmpty())	tactive[nactive++] = i;
	}
	pmdiag.Diag3x27(pm);// need diagonald status for columns analyssi
	OffToOff_Plus_New_Pairs_Status(pmz, wp_new);// find new pairs 
	// locked in box/row/col
	for (int id = 0; id < nactive; id++) {
		int dig = tactive[id];
		BF128 w = zdynstep.FD[dig][0], wdiag;
		wdiag.Diag3x27(w);
		for (int iband = 0; iband < 3; iband++) {
			int rband = DynLock1(w.bf.u32[iband]);
			if (rband) {
				if (diagexpand) cout << "seen lock dig " << dig + 1 << " band " << iband << endl;
				int tw[27], ntw=0;
				BitsInTable32(tw, ntw, rband);
				for (int i = 0; i < ntw; i++) {
					int ce = tw[i] + 27 * iband;
					if (expandmode) { source = 0; }// to see later
					// here need a BF128 with digit + cells
					Addt(ce, dig, 1, source);
				}
			}
			int rbanddiag = DynLock1(wdiag.bf.u32[iband]);
			if (rbanddiag) {
				if (diagexpand) cout << "seen lock diag dig " << dig + 1 << " band " << iband << endl;
				int tw[27], ntw=0;
				BitsInTable32(tw, ntw, rbanddiag);
				for (int i = 0; i < ntw; i++) {
					int ce = C_transpose_d[tw[i] + 27 * iband];
					if (expandmode) { source = 0; }// to see later
					// here need a BF128 with digit + cells
					Addt(ce, dig, 1, source);
				}

			}
		}// band loop
	}// digit loop
	if(wp_new.isNotEmpty()){// look here for new  naked pairs 
		int tnak[80], ntnak = wp_new.Table3X27(tnak);
		for (int inak = 0; inak < ntnak; inak++) {// see row col box
			int celln1 = tnak[inak];
			BF128 seenpairs = cell_z3x[celln1]; 
			seenpairs &= pairs_step_old;// here new status
			int tdig[2], ndig = 0;
			for (int i = 0; i < 9; i++) {
				BF128 p = pmz[i];
				if (p.Off_c(celln1))continue;
				tdig[ndig++] = i;
				seenpairs &= p;// must be the same pair
			}
			if (seenpairs.isEmpty())continue;
			if (diagexpand) cout << "gopair with cell " << cellsFixedData[celln1].pt << endl;
			int tunit[3];
			cellsFixedData[celln1].GetRegions(tunit);
			for (int iu = 0; iu < 3; iu++) {// check the 3 units
				BF128 wu= units3xBM[tunit[iu]];
				wu &= seenpairs;// 0 or 1 cell (more???)
				if (wu.isEmpty())continue;
				int celln2 = wu.getFirstCell();
				if (diagexpand) cout << "gopair2 with cell " << cellsFixedData[celln2].pt << endl;
				BF128 clean= cell_z3x[celln1];
				clean &= cell_z3x[celln2];//cells to clean
				BF128 cleand1 = clean & pmz[tdig[0]];
				cleand1 -= pm_go.xysearch.used_off_digits.pmdig[tdig[0]];//cells to clean dig1
				BF128 cleand2 = clean & pmz[tdig[1]];
				cleand2 -= pm_go.xysearch.used_off_digits.pmdig[tdig[1]];//cells to clean dig2
				Clean_d_bf(tdig[0], cleand1, source);// see later for source
				Clean_d_bf(tdig[1], cleand2, source);

			}
		}
	}
	BuildHiddenBiv_Xwings(pm, pmdiag, hp_biv_step, xw_biv_step);
	// hidden pair 36 pairs of digits 27 units floors_2d[36]
	if (hp_biv_step.Comp(hp_biv_old)) {// a new hidden pair has been seen
		for (int id1d2 = 0; id1d2 < 36; id1d2++) {// check the 36 pairs of digits
			int units_new = hp_biv_step.sets_biv[id1d2] ^ hp_biv_old.sets_biv[id1d2];
			if (!units_new)continue;
			int d1d2 = floors_2d[id1d2], tunits[27], ntu = 0;;
			bitscanforward(d1, d1d2); bitscanreverse(d2, d1d2);
			BitsInTable32(tunits, ntu, units_new);
			for (int iu = 0; iu < ntu; iu++) {
				int unit = tunits[iu];
				BF128 wu= units3xBM[unit];
				wu &= pm.pmdig[d1];
				wu &= pm.pmdig[d2];// now wu has 2 cells see if naked
				wu-= pairs_step_old;// nothing to do in bivalue cells
				if (wu.isEmpty())continue;
				int c1 = wu.getFirstCell(), c2 = wu.getLastCell();
				if (diagexpand) {
					char ws[82];
					cout<<wu.String3X(ws) << "try hidden pairs unit "<<unit << endl;
					cout << d1 + 1 << d2 + 1 << " " << cellsFixedData[c1].pt << " " << cellsFixedData[c2].pt << endl;
				}
				int digs1 = dig_cells_live[c1]^d1d2, digs2 = dig_cells_live[c2]^d1d2;
				Clean_digs_cell(digs1, c1, source);
				Clean_digs_cell(digs2, c2, source);
			}
		}

		hp_biv_old = hp_biv_step;
	}
	//Xwing must be hidden biv 2 rows or 2 cols
	if (xw_biv_step.Comp(xw_biv_old)) {// a new xwing  has been seen
		for (int ir1r2 = 0; ir1r2 < 36; ir1r2++) {// check the 36 pairs of digits
			int digits_new = xw_biv_step.sets_biv[ir1r2] ^ xw_biv_old.sets_biv[ir1r2];
			if (!digits_new)continue;// unit is here rows or cols
			int r1r2 = floors_2d[ir1r2], tdigits[18], ntd = 0,rr1,rr2;
			bitscanforward(rr1, r1r2); bitscanreverse(rr2, r1r2);
			BitsInTable32(tdigits, ntd, digits_new);
			for (int id = 0; id < ntd; id++) {
				int digit = tdigits[id];
				if (digit < 9) {// rows Xwing
					BF128 wu = units3xBM[rr1];
					wu |= units3xBM[rr2];
					wu &= pm.pmdig[digit];
					int c1 = wu.getFirstCell(), c2 = wu.getLastCell();
					int col1 = cellsFixedData[c1].plu, col2 = cellsFixedData[c2].plu;
					BF128 wu2 = units3xBM[col1];
					wu2 |= units3xBM[col2];
					wu2 &= pm.pmdig[digit];
					wu2 -= wu;// cells to clear
					if (wu2.isEmpty())continue;
					if (diagexpand) {
						char ws[82];
						cout << wu.String3X(ws) << "xwing row digit " << digit+1 << " "
						     << rr1 + 1 << rr2 + 1 << " " << cellsFixedData[c1].pt << " " << cellsFixedData[c2].pt << endl;
					}
					Clean_d_bf(digit, wu2, source);
				}
				else {//columns Xwing
					BF128 wu = units3xBM[rr1+9];
					wu |= units3xBM[rr2+9];
					wu &= pm.pmdig[digit-9];
					int c1 = wu.getFirstCell(), c2 = wu.getLastCell();
					int row1 = cellsFixedData[c1].el, row2 = cellsFixedData[c2].el;
					BF128 wu2 = units3xBM[row1];
					wu2 |= units3xBM[row2];
					wu2 &= pm.pmdig[digit-9];
					wu2 -= wu;// cells to clear
					if (wu2.isEmpty())continue;
					if ( diagexpand) {
						char ws[82];
						cout << wu.String3X(ws) << "xwing column digit " << digit -8 << " "
						 << rr1 -8 << rr2 -8 << " " << cellsFixedData[c1].pt << " " << cellsFixedData[c2].pt << endl;
					}
					Clean_d_bf(digit-9, wu2, source);

				}
			}
		}

		xw_biv_old = xw_biv_step;
	}
	if (diagexpand) {
		cout << "exit off to off" << endl;
		//zdyn.ImageCandidats();
	}

}
void XYSEARCHDYNPLUS::ExpandPlusGo(){// start with cand on
			//cout <<  ddig + 1 << cellsFixedData[dcell].pt << endl;
	diagexpand = 0;
	//if (pm_go.cycle == 2 && ddig == 6 && dcell == 9)		diagexpand = 2;
	if (diagexpand) {
		cout << "start xyexpand " << ddig + 1 << cellsFixedData[dcell].pt
			<< " dind=" << dind << " maxpas=" << maxpas << " cont="<< is_contradiction << endl;
		zdyn.ImageCandidats();
	}
	ntd = 0;
	while (nsteps++ < maxpas){ 
		int ntpr = ntd;
		ntp = ntd;
		ntd =nt;
		if (diagexpand >1) cout << "step " << nsteps << " ntp=" << ntp << " ntd=" << ntd 
			<< " cont=" << is_contradiction << endl;
		for (int i = ntp; i < ntd; i++){// std off/on and on/off
			cell = t[i].u8[0];	 
			digit = t[i].u8[1]; 
			sign= t[i].u8[2];
			if (diagexpand > 1) {
				cout << "process ";
				if (sign)cout << "~";
				cout << digit + 1 << cellsFixedData[cell].pt << endl;
			}
			if (sign)	OffToOn_Plus(i);		
			else  OnToOff_Plus(i);
			if (nt > 350) break;
		}
		// if new off in the previous step, consider  off to off
		if (memcmp(&pm_go.xysearch.used_off_digits, &used_off_old, 9 * 16)) {
			PM3X r_off = pm_go.xysearch.used_off_digits;
			if (diagexpand > 1)cout << "start off to off step" << endl;
			OffToOff_Plus();
			used_off_old = r_off;
		}

		if (nt == ntd) 	break;
		if (!is_contradiction){
			pm_go.xysearch.contradiction = pm_go.xysearch.used_off_digits;  
			pm_go.xysearch.contradiction &= pm_go.xysearch.used_on_digits;
			if (!pm_go.xysearch.contradiction.IsEmpty()){
				if (diagexpand) {
					cout << ddig + 1 << cellsFixedData[dcell].pt << " contradiction for npas=" << nsteps << " nt=" << nt << endl;
					pm_go.xysearch.contradiction.Print("contradiction ");
					if (diagexpand > 2)DebugT();
				}
				is_contradiction = 1;
				if (fastmode && (!dsign)) {
					AddElim(ddig, dcell, 96);
					if (diagexpand)cout << "94 fastelim " << ddig + 1 << cellsFixedData[dcell].pt << endl;
					pm_go.xysearch.off_status[dind] = pm_go.xysearch.used_off_digits;
					return;
				}
				//if (maxpas > nsteps + 2)maxpas = nsteps + 2;
//<<<<<<<<<<<< //DynamicSolveContradiction(cand,contradiction);
			}
		}
		if (diagexpand) {
			cout << "end of step=" << nsteps << " nt=" << nt << endl;
			zdyn.ImageCandidats();
		}
		if ((maxpas<=12 && nt > 150)|| nt>300)break;
	}
	if (diagexpand) {
		cout << "expand closed nt=" << nt << " off count" << pm_go.xysearch.used_off_digits.Count() << endl;
		pm_go.xysearch.used_off_digits.Print("dyn off status ");
	}
	if(dsign)off_off_status= pm_go.xysearch.used_off_digits;
	else pm_go.xysearch.off_status[dind] = pm_go.xysearch.used_off_digits;// store off status 

}
/*
int XYSEARCHDYNPLUS::ExpandDynamicToElim(GINT cand,GINT target){// start with cand on
	int tdig = target.u16[1], tcell = target.u16[0];
	ddig = cand.u16[1];
	dcell = cand.u16[0];
	int locdiag = 0;
	//if (pm_go.cycle == 5 && ddig == 7 && dcell == 19 && tdig ==4 && tcell== 20)		locdiag = 1;
	if (locdiag)cout << "expand to elim diag mode target " << cand.u16[1] + 1 << cellsFixedData[cand.u16[0]].pt
		<< " to "<< tdig + 1 << cellsFixedData[tcell].pt << endl;
	//if( tcell == 33) locdiag = 1;
	nsteps = 0;
	nt = 1;	
	t[0].u64 = dcell | (ddig << 16);	// source to 0
	used_off_digits.SetAll_0();	used_on_digits.SetAll_0();
	used_on_digits.Set_c(ddig, dcell);
	ntd = 0;
	while (nsteps++ < 40){// should never pass the target
		ntp = ntd;
		ntd = nt;
		for (int i = ntp; i < ntd; i++){
			cell = t[i].u16[0];		 digit = t[i].u16[1];
			if (nsteps & 1){
				OnToOff_Dyn(i);
				if (used_off_digits.On_c(tdig, tcell)) return 1;
			}
			else OffToOn_Dyn(i);
			if (nt > 350)return  0;// TEMP SAFETY TEST
		}
		if (nt > 350) {
			if (locdiag) {
				cout << "limit reached  npas=" << nsteps << " nt=" << nt << endl;
				if (nsteps & 1)used_off_digits.Print("off status");
				else used_on_digits.Print("on status");
			}
			return 0;
		}
		if (locdiag ){
			cout << "end npas=" << nsteps << " nt=" << nt<< endl;
			if (nsteps&1)used_off_digits.Print("off status");
			else used_on_digits.Print("on status");
		}
		if (nt == ntd ) 	break;// should never be
	}
	return 0;
}



void XYSEARCHDYNPLUS::SearchDynPassMulti(int nmax){// try multi chains if nothing low
	if (elim_done) return;
	if (ntelims && maxrating <= 88) return;
	if (opprint)cout << "try cells not bi values" << endl;
	// try all bi values in mode x->~a and y->~a adding one in length
	BF128 wp = zh_g.cells_unsolved_e-pairs;
	GINT target,p;
	GINT64 tbn[9][400];
	int ntbn[9], nx,xcell;
	PM3X welims;
	while ((xcell = wp.getFirsCell()) >= 0){
		wp.Clear_c(xcell);
		int digs = zh_g.dig_cells[xcell], ndigs=_popcnt32(digs);
		if (nmax < 8 && ndigs>3)continue;
		welims.SetAll_1();
		while ( digs){
			bitscanforward(d2, digs);
			digs ^= 1 << d2;
			int i1 = ind_pm[d2][xcell], coff = tcands[i1].u16[1];
			welims &= off_status[coff];
		}
		if (welims.IsEmpty()) continue;
		if (fastmode){// do it in bloc and return
			if (opprint) {
				cout << "cell multi analysis " << cellsFixedData[xcell].pt << endl;
				welims.Print(" elims seen ");
			}
			zhou_solve.Clean(welims);
			elim_done = 1;
			continue;
		}
		for (int id = 0; id < 9; id++) if (welims.pmdig[id].isNotEmpty()){
			BF128 elimd = welims.pmdig[id];
			int elim_cell;
			while ((elim_cell = elimd.getFirsCell()) >= 0){
				elimd.Clear_c(elim_cell);
				target.u32 = elim_cell | (id << 16);
				int length = 0, n = 0;
				digs = zh_g.dig_cells[xcell];
				while (digs){
					bitscanforward(d2, digs);
					digs ^= 1 << d2;
					p.u32 = xcell | (d2 << 16);
					if (!ExpandDynamicToElim(p, target)) {// redo expansion should work
						cout<< "anomaly in cell redo "<<d2+1<<cellsFixedData[xcell].pt
							<<" to target "<< id+1<< cellsFixedData[elim_cell].pt << endl;
						fout1 << "anomaly in cell redo" << endl;
						return;
					}
					nx=ntbn[n] = BackDynamicOff(target);
					if (locdiag){
						cout << "nx=" << nx  << endl;
						PrintBackCom("locdiag on path ", tback, nx, 1);
					}
					length += nx;
					memcpy(tbn[n++], tback, nx * 8);

				}
				int rating = pm_go.hint.ChainLengthAdjusted(85, length);
				if (AddElim(id, elim_cell, rating)){
					if (opprint){
						cout << "cleaned or stored rating " << rating << endl;
						for (int ip = 0; ip < n; ip++)
							PrintBackCom("", tbn[ip], ntbn[ip], 1);
					}
				}

			}
		}

	}
	if (opprint)cout << "try units not bi values" << endl;
	for (int id1 = 0; id1 < 9; id1++){
		for (int iu = 0; iu < 27; iu++){
			if (dig_bivsets[id1].On(iu))continue;
			BF128 wu = units3xBM[iu]; wu &= zh_g.pm.pmdig[id1];
			if (wu.isEmpty())continue;
			int tcu[10], ntcu = wu.Table3X27(tcu);
			welims.SetAll_1();
			for (int icu = 0; icu < ntcu; icu++){
				int xcell = tcu[icu]	;
				int i1 = ind_pm[id1][xcell], coff = tcands[i1].u16[1];
				welims &= off_status[coff];
			}
			if (welims.IsEmpty()) continue;
			if (0 && opprint) {
				cout << "digit " << id1 + 1 << " unit " << iu  << " multi analysis"<< endl;
				welims.Print(" unit elims seen ");
			}
			if (fastmode){// do it in bloc and return
				zhou_solve.Clean(welims);
				elim_done = 1;
				continue;
			}
			for (int id = 0; id < 9; id++) if (welims.pmdig[id].isNotEmpty()){
				BF128 elimd = welims.pmdig[id];
				int elim_cell;
				while ((elim_cell = elimd.getFirsCell()) >= 0){
					elimd.Clear_c(elim_cell);
					target.u32 = elim_cell | (id << 16);
					if (locdiag)cout << "try unit elim" << id + 1 << cellsFixedData[elim_cell].pt << endl;
					int length = 0, n = 0,xcell;
					for (int icu = 0; icu < ntcu; icu++){
						xcell = tcu[icu];
						p.u32 = xcell | (id1 << 16);
						if (!ExpandDynamicToElim(p, target)) {// redo expansion should work
							cout << "trying unit elim" << id + 1 << cellsFixedData[elim_cell].pt << endl;
							cout << "anomaly in unit/cell redo start "
								<<id1+1<<cellsFixedData[xcell].pt<< endl;
							fout1 << "anomaly in unit/cell redo" << endl;
							return;
						}
						nx = ntbn[n] = BackDynamicOff(target);
						if (locdiag){
							cout << "nx=" << nx << endl;
							PrintBackCom("locdiag on path ", tback, nx, 1);
						}
						length += nx;
						memcpy(tbn[n++], tback, nx * 8);

					}
					int rating = pm_go.hint.ChainLengthAdjusted(85, length);
					if (AddElim(id, elim_cell, rating)){
						if (opprint){
							cout << "cleaned or stored rating " << rating << endl;
							for (int ip = 0; ip < n; ip++)
								PrintBackCom("", tbn[ip], ntbn[ip], 1);
						}
					}
				}
			}
		}
	}
}
*/

void XYSEARCHDYNPLUS::SearchDynPassPlus(int nmax) {	// try a  pass limited to nmax steps
	//nind = 0;
	maxpas = nmax;
	for (int icand = 0; icand < pm_go.xysearch.ntcands; icand++) {// all candidates processed 
		ExpandPlusInit(pm_go.xysearch.tcands[icand]);
		ExpandPlusGo();
		ExpandPlusInit(pm_go.xysearch.tcands[icand],1);
		ExpandPlusGo();
		// check si signe de 2 est moins
		// on peut aussi tester dans expand arrêt si pas de - -> -
	}
	if (opprint) {
		cout << "search pass end phase 1 elim_done status=" << elim_done
			<< "  maxrating=" << maxrating
			<< endl;
	}
	if (elim_done) return;
	if (1) return;
	/*
	if (opprint)cout << "try cells bi values" << endl;
	// try all bi values in mode x->~a and y->~a adding one in length
	BF128 wp = pairs;
	uint32_t  dc1, dc2;
	while ((cell = wp.getFirsCell()) >= 0) {
		if (0 && opprint) {
			cout << "cells " << cellsFixedData[cell].pt << endl;
		}
		wp.Clear_c(cell);
		int digs = zh_g.dig_cells[cell];
		bitscanforward(dc1, digs);
		bitscanreverse(dc2, digs);
		int i1 = ind_pm[dc1][cell], i2 = ind_pm[dc2][cell];// pointers to tcands
		GINT cand1 = tcands[i1], cand2 = tcands[i2];
		PM3X welims = off_status[cand1.u16[1]]; welims &= off_status[cand2.u16[1]];
		if (welims.IsEmpty()) continue;
		//<<<<<<<<<DynamicSolveContradiction(dc1, cell, dc2, cell, welims);
	}
	if (opprint)cout << "try unit bi values" << endl;
	for (int id = 0; id < 9; id++) {
		//if (nmax > 10)cout << "digit=" << id + 1 << endl;
		for (int iu = 0; iu < 27; iu++) {
			if (dig_bivsets[id].Off(iu))continue;
			BF128 wu = units3xBM[iu]; wu &= zh_g.pm.pmdig[id];
			int cell1 = wu.getFirsCell();
			wu.Clear_c(cell1);
			int cell2 = wu.getFirsCell();
			int i1 = ind_pm[id][cell1], i2 = ind_pm[id][cell2];// pointers to tcands
			GINT cand1 = tcands[i1], cand2 = tcands[i2];
			PM3X welims = off_status[cand1.u16[1]]; welims &= off_status[cand2.u16[1]];
			if (welims.IsEmpty()) continue;
			DynamicSolveContradiction(id, cell1, id, cell2, welims);
		}
	}

	//if (nmax == 6)SearchDynPassMulti(6);
	//else
	SearchDynPassMulti(nmax);
	*/
}
int XYSEARCHDYNPLUS::SearchDynPlus(int fast){
	opprint = pm_go.opprint;
	PM3X welims;
	SearchInit(fast);
	if (opprint) {// debug hp_biv
		cout << "end init status" << endl;
		hp_biv_0.Debug("digit bi_values at start");
		xw_biv_0.Debug("xwings at start");
	}
	if (fast) {
		maxpas = 35;
		if (opprint)cout << "fast 1" << endl;
		for (int icand = 0; icand < pm_go.xysearch.ntcands; icand++) {// all candidates processed 
			ExpandPlusInit(pm_go.xysearch.tcands[icand]);
			if ((int)ddig == (int)zh_g2.zerobased_sol[dcell]) continue;
			ExpandPlusGo();
		}
		if(elim_done) return 1;
		if (opprint)cout << "fast 2" << endl;
		for (int icand = 0; icand < pm_go.xysearch.ntcands; icand++) {// all candidates processed 
			ExpandPlusInit(pm_go.xysearch.tcands[icand]);
			if ((int)ddig == (int)zh_g2.zerobased_sol[dcell])
				ExpandPlusGo();// add missing candidates off status (no contradiction)
		}
		// test on and off (includes all bivalues) no contradiction
		for (int icand = 0; icand < pm_go.xysearch.ntcands; icand++) {// all candidates processed 
			ExpandPlusInit(pm_go.xysearch.tcands[icand], 1);
			ExpandPlusGo();// off_off_status
			off_off_status &= pm_go.xysearch.off_status[dind];
			if (!off_off_status.IsEmpty()) {
				if(AddElimFast(off_off_status))
					if(opprint)	cout << "elim fast on/off "<<ddig+1
						<<cellsFixedData[dcell].pt<<endl;
			}
		}
		if (elim_done) return 1;
		if (opprint)cout << "fast 3" << endl;
		// test multis cells >2 digits units >2 cells
		if (opprint)cout << "try fast multi cells" << endl;
		// try all bi values in mode x->~a and y->~a adding one in length
		BF128 wp = zh_g2.cells_unsolved_e - pairs;
		int  xcell;
		while ((xcell = wp.getFirstCell()) >= 0) {
			wp.Clear_c(xcell);
			int locdiag = 0;
			//if (pm_go.cycle == 7 && xcell == 44)		locdiag = 2;
			int digs = zh_g2.dig_cells[xcell];
			if (_popcnt32(digs) < 3)continue;
			welims.SetAll_1();
			if (locdiag)cout << "multi cell dplus r5c9" << endl;
			while (digs) {
				bitscanforward(d2, digs);
				digs ^= 1 << d2;
				int i1 = pm_go.xysearch.ind_pm[d2][xcell], 
					coff = pm_go.xysearch.tcands[i1].u16[1];
				if (locdiag) {
					cout << "for digit " << d2 + 1 <<" coff="<<coff<< endl;
					//char ws[82];
					pm_go.xysearch.off_status[coff].Print("off status");
				}
				welims &= pm_go.xysearch.off_status[coff];
			}
			if (welims.IsEmpty()) continue;
			if (AddElimFast(welims))
				if (opprint)	cout << "elim fast multicell " 
					<< cellsFixedData[xcell].pt << endl;
		}
		if (opprint)cout << "try fast multi units" << endl;
		for (int id1 = 0; id1 < 9; id1++) {
			for (int iu = 0; iu < 27; iu++) {
				if (dig_bivsets[id1].On(iu))continue;
				BF128 wu = units3xBM[iu]; wu &= zh_g2.pm.pmdig[id1];
				if (wu.isEmpty())continue;
				int tcu[10], ntcu = wu.Table3X27(tcu);
				welims.SetAll_1();
				for (int icu = 0; icu < ntcu; icu++) {
					int xcell = tcu[icu];
					int i1 = pm_go.xysearch.ind_pm[id1][xcell], 
						coff = pm_go.xysearch.tcands[i1].u16[1];
					welims &= pm_go.xysearch.off_status[coff];
				}
				if (welims.IsEmpty()) continue;
				if (AddElimFast(welims))
					if (opprint)	cout << "elim fast multi unit digit "
						<<id1+1 << "unit " <<iu<< endl;
			}
		}
		return elim_done;
	}

	SearchDynPassPlus(10);// try a first pass limited to 6 steps
	if (elim_done) return 1;
	if (1) return elim_done;
	if (ntelims && maxrating <= 93){// do it if small enough 
		pm_go.hint.rating_done = (USHORT)maxrating;
		for (int i = 0; i < ntelims; i++){
			GINT w = telims[i];
			zhou_solve.ClearCandidate_c(w.u16[1], w.u16[0]);
		}
		return 1;
	}
	if (opprint) {
		cout << "search phase 1 closed ntelims="<< ntelims	<< endl;
	}
	//if (ntelims) SearchDynPass(12);//  just secure the found rating
	//else 
	//SearchDynPass(25);// try a second  pass "no limit"
	if (elim_done) return 1;
	if (ntelims ){// do it  
		pm_go.hint.rating_done = (USHORT)maxrating;
		for (int i = 0; i < ntelims; i++){
			GINT w = telims[i];
			zhou_solve.ClearCandidate_c(w.u16[1], w.u16[0]);
		}
		return 1;
	}
	return 0;
}
/*
void XYSEARCHDYNPLUS::DynamicSolveContradiction(GINT cand,PM3X cont){// find path and elims low rating a-> x/~x
	if (fastmode){// do it in bloc and return
		zhou_solve.ClearCandidate_c(cand.u8[1], cand.u8[0]);
		elim_done = 1;
		return;
	}
	int diag = 0;
	//if (pm_go.cycle == 14 &&cand.u8[1] == 8 && cand.u8[0] == 15)diag = 1;
	if (diag) cout << "solve cont for " << cand.u8[1] + 1 << cellsFixedData[cand.u8[0]].pt << endl;
	for (int id = 0; id < 9; id++) if (cont.pmdig[id].isNotEmpty()){
		BF128 elimd = cont.pmdig[id];
		int elim_cell;
		while ((elim_cell = elimd.getFirsCell()) >= 0){
			elimd.Clear_c(elim_cell);
			GINT64 target_on, target_off;
			int n = 0;
			for (int i = 0; i < nt; i++){// locate targets 
				GINT64 w = t[i];
				if (w.u16[1] != id)continue;
				if (w.u16[0] != elim_cell)continue;
				n++;
				if (w.u16[3] & 1)target_off = w; else target_on = w;
			}
			if (n != 2) continue;// should never be
			if (diag){
				cout << "target_on" << target_on.u16[1] + 1 << cellsFixedData[target_on.u16[0]].pt
					<< " step" << target_on.u16[3] << " source=" << target_on.u16[2] << endl;
				cout << "target_off" << target_off.u16[1] + 1 << cellsFixedData[target_off.u16[0]].pt
					<< " step" << target_off.u16[3] << " source=" << target_off.u16[2] << endl;

			}
			int n1 = BackDynamic(target_on, t, nt);
			GINT64 t1b[200];
			memcpy(t1b, tback, n1 * 8);
			int n2 = BackDynamic(target_off, t, nt);
			int rating = pm_go.hint.ChainLengthAdjusted(85, n1 + n2);
			//cout << "rating " << rating<<endl<<endl;
			if (AddElim(tback[0].u16[1], tback[0].u16[0], rating)){
				if (opprint){
					cout << "cleaned or stored rating " << rating <<" ntelim="<<ntelims<< endl;
					PrintBackCom("off path ", tback, n2, 1);
					PrintBackCom("on  path ", t1b, n1, 1);
					cout << "done" << endl;
				}
			}
		}
	}
}

void XYSEARCHDYNPLUS::DynamicSolveContradiction(int dig1, int cell1, int dig2, int cell2, PM3X cont){
	if (fastmode){// do it in bloc and return
		zhou_solve.Clean(cont);
		elim_done = 1;
		return;
	}
	int locdiag = 0;
	//if (pm_go.cycle == 13 && dig1==2 &&dig2==2 && cell1>=72)		locdiag = 1;
	//if (dig1 != dig2 && cell1 == 62) locdiag = 1;
	if (locdiag){
		cout << "DynamicSolveContradiction " << dig1 + 1 << cellsFixedData[cell1].pt << " "
			<< dig2 + 1 << cellsFixedData[cell2].pt  << endl;
		cont.Print("for elims");

	}
	GINT p1, p2,target;
	p1.u32 = cell1 | (dig1 << 16); 
	p2.u32 = cell2 | (dig2 << 16);
	for (int id = 0; id < 9; id++) if (cont.pmdig[id].isNotEmpty()){
		BF128 elimd = cont.pmdig[id];
		int elim_cell;
		while ((elim_cell = elimd.getFirsCell()) >= 0){
			elimd.Clear_c(elim_cell);
			target.u32 = elim_cell | (id <<16);
			if (locdiag) cout << "to elim " << id + 1 << cellsFixedData[elim_cell].pt << endl;
			if(!ExpandDynamicToElim(p1,target)) continue;// redo expansion should work
			int n1 = BackDynamicOff(target);
			if (locdiag){
				cout << "n1=" << n1 << endl;
				PrintBackCom("locdiag on path ", tback, n1, 1);
			}
			if (!n1) continue; // should never be
			GINT64 t1b[200];
			memcpy(t1b, tback, n1 * 8);
			if (!ExpandDynamicToElim(p2, target)) continue;// redo expansion should work
			int n2 = BackDynamicOff(target);
			if (locdiag){
				cout << "n2=" << n2 << endl;
				PrintBackCom("locdiag off path ", tback, n2, 1);
				if (!n2){
					cout << "pas trouve cible nt=" << nt << endl;
					used_off_digits.Print("used off status");
				}
			}
			if (!n2) continue; // should never be
			int rating = pm_go.hint.ChainLengthAdjusted(85, n1 + n2);
			if (AddElim(id, elim_cell, rating)){
				if (opprint){
					cout << "cleaned or stored rating " << rating << " ntelim=" << ntelims << endl;
					PrintBackCom("off path ", tback, n2, 1);
					PrintBackCom("on  path ", t1b, n1, 1);
				}
			}
		}
	}
}
int XYSEARCHDYNPLUS::BackDynamic(GINT64 target, GINT64 * tb, int ntb) {
	int diag = 0;
	//if (pm_go.cycle == 14 && target.u16[1] == 5 && target.u16[0] == 31)diag = 1;

	PM3X back_bf; back_bf.SetAll_0(); // no possible conflict in the way back
	int itret = 1, itret1 = 0;
	GINT64 tretw[300];// cell;digit;source/last in;step/sign
	tretw[0] =target;
	back_bf.Set_c(target.u16[1], target.u16[0]);
	while (itret1 < itret){// && itret < 100 ) { // solve each entry back
		GINT64 x = tretw[itret1];
		if (diag){
			cout << "BackDynamic " << x.u16[1]+1 << cellsFixedData[x.u16[0]].pt
				<< " step " << x.u16[3] << " source=0x" << hex << x.u16[2] << dec << endl;
		}
		if (!x.u16[3]) { itret1++;	continue; }  // start point
		if (!(x.u16[2] & 0x3000)) { // this is direct
			GINT64 y = tb[x.u16[2]];
			if (diag){
				cout << "source direct BackDynamic " << y.u16[1] + 1 << cellsFixedData[y.u16[0]].pt
					<< " status  " << back_bf.On_c(y.u16[1], y.u16[0] )<< endl;
			}
			if (back_bf.On_c(y.u16[1], y.u16[0])){ itret1++;	continue; }
			back_bf.Set_c(y.u16[1], y.u16[0]);
			tretw[itret++] = y;
			itret1++;
			continue;
		}
		// now this comes from a set last in cell/region
		int xcell = x.u16[0], xdig = x.u16[1];
		switch (x.u16[2] >> 12){
		case 1:{// last in cell
			int digs = zh_g.dig_cells[xcell] ^ (1 << xdig);
			while (digs){
				bitscanforward(d2, digs);
				digs ^= 1 << d2;
				if (back_bf.On_c(d2, xcell))	continue;
				back_bf.Set_c(d2, xcell);
				for (int i = 0; i < ntb; i++){// must exist
					GINT64 w = tb[i];
					if (w.u16[0] == xcell && w.u16[1] == d2 && (w.u16[3] & 1))
						tretw[itret++] = w;
				}
			}
			break;
		}
		case 2:{// last in region
			int unit = x.u8[4], ucell;
			BF128 wu = units3xBM[unit]; wu &= zh_g.pm.pmdig[xdig];
			if (diag){
				char ws[82];
				cout << wu.String3X(ws) << " wu to backload for last" << endl;
			}
			while ((ucell = wu.getFirsCell()) >= 0){
				wu.Clear_c(ucell);
				if (back_bf.On_c(xdig, ucell))	continue;
				back_bf.Set_c(xdig, ucell);
				//cout << " look for cell " << cellsFixedData[ucell].pt << endl;
				for (int i = 0; i < ntb; i++){// must exist
					GINT64 w = tb[i];
					if (w.u16[0] == ucell && w.u16[1] == xdig && (w.u16[3] & 1)){
						tretw[itret++] = w;
						break;
					}
				}
			}

		}
		}//  end switch
		itret1++;
	}	
	// send it back increasing order of step
	int stepl = tretw[0].u16[3];
	int nb = 0;
	for (int ist = 0; ist <= stepl; ist++){
		for (int j = 0; j < itret; j++) if (tretw[j].u16[3] == ist)
			tback[nb++] = tretw[j];
	}
	return nb;
}
*/
/*

int XYSEARCHDYNPLUS::BackDynamicOff(GINT target){
		// locate target in the last step (must be a on to off step
	for (int i = ntd; i < nt; i++){
		if (t[i].u32[0]==target.u32)
			return BackDynamic(t[i],t,nt);
	}
	return 0;
}
*/

