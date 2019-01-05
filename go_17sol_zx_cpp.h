void ZHOU::SetUab3() {
	zh_g.band3nextua = 0;
	for (int idig = 0; idig < 9; idig++) {// check sol against solution grid
		zh_g.band3nextua |= FD[idig][0].bf.u32[2] 
			& (~zh_g.band3digits[idig]);
	}
}
void ZHOU::SetUab12() {
	uint64_t rd = 0;
	for (int idig = 0; idig < 9; idig++) {// check sol against solution grid
		register uint64_t r = ~zh_g.digsols[idig];
		r &= FD[idig][0].bf.u64[0];
		rd |= r;
	}
	zh_g.b12nextua = rd;
}

void ZHOU::InitBand3PerDigit(int * grid0b3){
	memset(zh_g.band3digits, 0, sizeof zh_g.band3digits);
	register int * t = zh_g.band3digits;
	for (int i = 0; i < 27; i++){
		t[grid0b3[i]] |= 1 << i;
	}
	zh_g.digsols = zh2b_g.digsols;// catch pointer to solution per digit
}
int ZHOU::PartialInitSearch17(uint32_t * t, int n){
	zh_g.digitsbf = 0;
	memset(zh_g.Digit_cell_Assigned, 0, sizeof zh_g.Digit_cell_Assigned);
	*this=zhou_i;
	for (int icell = 0; icell<n; icell++)   {
		int cell = t[icell], digit = zh_g.grid0[cell];
		zh_g.digitsbf |= 1 << digit;
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g.Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i<9; i++)  FD[i][0] &= w | zh_g.Digit_cell_Assigned[i];
	return 0;
}
int ZHOU::EndInitSearch17(ZHOU & o, int * t, int n) {
	*this = o;
	BF128 Digit_cell_Assigned[9];
	memcpy(Digit_cell_Assigned, zh_g.Digit_cell_Assigned, sizeof Digit_cell_Assigned);
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = zh_g.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		Digit_cell_Assigned[digit].Set(xcell);
	}
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | Digit_cell_Assigned[i];
	return 0;
}
void ZHOU::GuessB12() {// band3 solved "false"
	uint32_t res;
	if (zh_g.nsol > zh_g.lim) return;
	if (cells_unsolved.isEmpty()) {
		zh_g.nsol++;
		zh_g.go_back = zh_g.band3nextua;
		SetUab12();
		return;
	}
	int ndig;
	uint64_t R3 = zh_g.pairs.bf.u64[0];
	if (bitscanforward64(res, R3)) {// first naked pair is ok to go
		ndig = 2;
		int cell = From_128_To_81[res], tdig[2], ndig = 0;
		for (int idig = 0; idig < 9; idig++)
			if (FD[idig][0].On(res))tdig[ndig++] = idig;
		ZHOU * mynext = this + 1; // start next guess
		*mynext=*this;
		mynext->SetaCom(tdig[0], cell, res);
		mynext->ComputeNextB12();
		SetaCom(tdig[1], cell, res);
		Upd1(tdig[1]);// try to save one update cycle
		ComputeNextB12();
		return;
	}
	// can godirectly to first cell, kind of escape lane
	int xcell = cells_unsolved.getFirst128(), cell = From_128_To_81[xcell];
	for (int idig = 0; idig < 9; idig++) {
		if (FD[idig][0].On(xcell)) {// one valid digit
			ZHOU * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->ComputeNextB12();
		}
	}

}
void ZHOU::GuessB3() {//guess priority to small ua in band 3
	if (zh_g.go_back) return;
	if (cells_unsolved.isEmpty()) {//if band 3 is "false" use it as UA
		SetUab3();
		if (zh_g.band3nextua) {// stop at first ua
			zh_g.nsol++;
			zh_g.go_back = zh_g.band3nextua;
			SetUab12();
		}
		return;
	}
	if (zh_g.pairs.bf.u32[2]) {
		GuessCellB3(zh_g.pairs.bf.u32[2]);
		return;
	}
	{// default use any cell in band 3 
		register int R3 = cells_unsolved.bf.u32[2];
		if (!R3) {// if all good in band 3 game over
			SetUab3();
			if (!zh_g.band3nextua) return; // this is the good solution
			//if not,must look for a valid puzzle in bands 1+2
			GuessB12();
			return;
		}
		GuessCellB3(R3);
	}
}
void ZHOU::GuessCellB3(int field) {
	uint32_t res;
	register int R3 = field;// zh_g.pairs.bitmap128.m128i_u32[2];
	if (bitscanforward(res, R3))return; // no  cell  should  not be

	int cell = res + 54, xcell = res + 64, bit = 1 << res;
	int digit = genb12.grid0[cell];
	{// select first the "good" digit
		ZHOU * mynext = this + 1; // start next guess
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->ComputeNextB3();
	}
	if (zh_g.go_back)return; // job is done
	//if job is not yet done, select the "false"
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue; // it is the "good" value
		if (FD[idig][0].bf.u32[2] & bit) {// one valid digit	
			ZHOU * mynext = this + 1; // start next guess
			*mynext = *this;
			mynext->SetaCom(idig, cell, xcell);
			mynext->ComputeNextB3();
			if (zh_g.go_back)return;
		}
	}
}
int ZHOU::GetNextUaB3x() {
	zh_g.band3nextua = zh_g.nsol = zh_g.lim = zh_g.go_back = 0;
	ComputeNextB3();
	return zh_g.band3nextua;
}
int ZHOU::MultipleTrial() {
	GuessB3();
	return zh_g.go_back;
}
int ZHOU::CallMultipleB3(ZHOU & o, int bf) {
	*this = o;
	zh_g.go_back = 0;
	zh_g.nsol = 0; zh_g.lim = 0;
	return MultipleB3(bf);
}
int ZHOU::MultipleB3(int bf) {// bf is band3 known cells bit field
	BF128 dca[9];
	int digitsbf = zh_g.digitsbf;
	memcpy(dca, zh_g.Digit_cell_Assigned, sizeof dca);
	{	uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = zh_g.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc+64; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	if (_popcnt32(digitsbf < 8)) return 1;// can not be one solution
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	if (FullUpdateV2() == 2) return 0;// solved can not be multiple

	ImageCandidats();
	BF128 tres[1000]; // buffer for perms
	GuessV2(tres);
	if (1) return 0;

	// solve a stack in band 3 as guess take the smallest
	int R3 = cells_unsolved.bf.u32[2];
	if (!cells_unsolved.bf.u32[2]) return 0; //band 3 solved, unique
	int bit = 1, mask = 07007007, nbf = 10, bfmx, nmx = 0;
	for (int i = 0; i < 3; i++, mask <<= 3) {
		int bf = cells_unsolved.bf.u32[2] & mask;
		if (bf) {
			nmx++;
			register  int n = _popcnt32(bf);
			if (n < nbf) { nbf = n;  bfmx = bf; }
		}
	}
	if (nmx > 1) {// more than one stack, loop
		ZHOU * mynext = this + 1; // start next guess
		*mynext = *this;
		int rp = zh_g.pairs.bf.u32[2]; // store cells with pairs
		if (mynext->MultipleB3loop(bfmx))return zh_g.go_back;
		zh_g.pairs.bf.u32[2] = rp; // restore cells with pairs status
		MultipleTrial();
	}
	else {// last stack trial and error from now on
		MultipleTrial();
	}
	return zh_g.go_back;
}



int ZHOU::MultipleB3loop(int bf){// bf is band3 known cells bit field
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)){
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = zh_g.grid0[cell];
			int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
			if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
			SetaCom(digit, cell, xcell);
		}

	}
	if (FullUpdate() == 2) return 0;// solved can not be multiple
	int R3 = cells_unsolved.bf.u32[2];
	if (!cells_unsolved.bf.u32[2]) return 0; //band 3 solved, unique
	int bit = 1, mask = 07007007, nbf = 10, bfmx, nmx = 0;
	for (int i = 0; i < 3; i++, mask <<= 3){
		int bf = cells_unsolved.bf.u32[2] & mask;
		if (bf){
			nmx++;
			register  int n = _popcnt32(bf);
			if (n < nbf){ nbf = n;  bfmx = bf; }
		}
	}
	if (nmx > 1){// more than one stack, loop
		ZHOU * mynext = this + 1; // start next guess
		*mynext = *this;
		int rp = zh_g.pairs.bf.u32[2]; // store cells with pairs
		if (mynext->MultipleB3loop(bfmx))return zh_g.go_back;
		zh_g.pairs.bf.u32[2] = rp; // restore cells with pairs status
		MultipleTrial();
	}
	else{// last stack trial and error from now on
		MultipleTrial();
	}
	return zh_g.go_back;
}
