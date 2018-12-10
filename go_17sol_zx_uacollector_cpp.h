

//================== UA collector one band 
/*
void GENUAS_1B::Initgen(int floords, int * cols, int * g0, int ncells_out) {// g0 solution; cols gangster to use
	//==============set solutions bit fields
	ismore = 0;
	limsize = 18 - ncells_out;
	memset(digsols, 0, sizeof digsols);// set all bit maps to 0
	register int R = 1;
	for (int i = 0; i < 27; i++, R <<= 1)digsols[g0[i]] |= R;
	if (0) {
		cout << "check digsols" << endl;
		for (int i = 0; i < 9; i++) {
			cout << Char27out(digsols[i]) << " dig=" << i + 1 << endl;
		}
	}
	//==============set pm bit fields
	memset(digpm, 0, sizeof digpm);
	for (int icol = 0; icol < 9; icol++) {// apply the gangster
		register int bits = 01001001 << icol, cdigs = cols[icol], R = 1;
		for (int idig = 0; idig < 9; idig++, R <<= 1) {
			if (cdigs & R)digpm[idig] |= bits;
		}
	}
	if (0) {
		cout << "check digpm" << endl;
		for (int i = 0; i < 9; i++) {
			cout << Char27out(digpm[i]) << " dig=" << i + 1 << endl;
		}
	}


	//============== order initial floors
	nfloorsd = __popcnt(floords);
	nfloors = ntdmore = 0;
	//cout <<" ncells_out="<<ncells_out<< " initial set of digits ";
	for (int i = 0; i < 9; i++)
		if ((floords & (1 << i))) {
			ordered_digpm[nfloors] = digpm[i];
			ordered_digsols[nfloors++] = digsols[i];
			//cout << i + 1;
		}
		else {
			tdmorebf[ntdmore] = 1 << i;
			tdmore[ntdmore++] = i;
		}
	//cout << endl;
	memset(ntemp, 0, sizeof ntemp);
}

void GENUAS_1B::InitUasMore() {// keep base UAs for subsets
	limsize = 18;
	runsolved = 07000000000;
}

int GENUAS_1B::AddUA(int  ux, int cc) {// add if and clear subsets supersets
	if ((!cc) || cc > limsize) return 0;// should never be but !!!
	register int w = ux | ismore, nw = ~w;
	for (int i = 0; i < cc; i++) {// is it a superset of another same cycle
		register int * tt = ztemp[i];
		for (int j = 0; j < ntemp[i]; j++)
			if (!(tt[j] & nw))return 0; // subset foundforget the new one
	}
	{ // possible redundancy, not supposed to come here
		register int * tt = ztemp[cc];
		for (int j = 0; j < ntemp[cc]; j++)
			if (tt[j] == w) return 0; // should not be here but !!
		if (ntemp[cc] < 200)tt[ntemp[cc]++] = w;
		else return 0;
		//cout << Char27out(w) << " added" << endl;
	}

	{// is w superset of a previous entry
		for (int i = cc + 1; i <= limsize; i++) {
			register int * tt = ztemp[i];
			for (int j = 0; j < ntemp[i]; j++) {
				register uint64_t x = ~tt[j];
				if (!(x & w)) {//w is subset of w, kill x
					for (int k = j + 1; k < ntemp[i]; k++)tt[k - 1] = tt[k];
					ntemp[i]--;
				}
			}
		}
	}
	return 1;// new loaded
}

void GENUAS_1B::PrintTemp(int ideb) {
	cout << "Print  Temporary table " << endl;
	int n = 0;
	for (int i = ideb; i <= limsize; i++)if (ntemp[i])
		for (int j = 0; j < ntemp[i]; j++) {
			int w = ztemp[i][j];
			if (w >> 27) // it is a true more
				cout << Char27out(w) << " " << ++n << " " << __popcnt(w&BIT_SET_27) << endl;
		}
}
*/


