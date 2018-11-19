

/* GO_CAN encapsulates all data and functions speficic to canonicalisation

   class V3C works on a pattern zero based puzzle
     Blocs(USHORT * res,char * lig) processes a ligne or a column
       . finds the number of clues per triplet row/box ...
	   . send back the six permutations in a 6 bits fields
	     aa bb cc  where aa,bb,cc are numer of clues 0_3
	     the six permutations are defined in the constant iii

	class BAND_STACK describes one permutation for a specific band/stack


*/
struct GO_CAN{
	struct FIRST_IN_BAND{
		int rc;        // row or column absolute index		
		VV9    perm,      // rows or columns order
			   puz;       // digit mode	
	};

	struct BAND_STACK{
		FIRST_IN_BAND fib;
		int rows[2]; // rows values for that peermutation
		VV9 relab; //relab status after relabelling 27
		int nextdigit;
	};

	struct TWO_BANDS{
		USHORT ibs,ra,rb,rc;
		VV9 relab; // relab status after relabelling 54 
		int nextdigit;
	};
	struct WSYM_PERM{ // a valid perm
		int target_is_row,  // 1 if row_ column  0 if column row
			   perm1[9],perm2[9];
	};
	struct WSYM{ // data for collection of Study_Perms_Sym()
		WSYM_PERM tp[100],wp;
		USHORT nfinal,band2,maxrow,  // size of final collection of 
			   cptr[3],cptrg[3][3],
			   otherband2[2][2],
			   ig1,ig2,ig3, // groups for first row
			   bandb,ir[9],ir2[9];   // first band / stack

	void Prepare1(BF128 & patx);
	void Prepare2();
	int Check(int a,int b);
	void Store();
	void Image_Out2();
	};

	static int iii[6][3],bf16_starts_row1[20][2];
	int bf16_directs[64];  //same as bf16_starts_row1 in direct access
	  // options of the command line
	//USHORT cop23,   add_fam_id,	   add_seq,	   input_0,	   output_0;
	//UINT   idpos;
	//char fam[40];

	VV9 relab_init, // set to all 0 by the constructor
		relab,      // new label in relabelling
		v9pat_0based[18],
		v9puz[18];

	int	max_store1, // max value found for store on pattern 
			t_store1[110], // table for store on pattern max in theory 18 x 6 = 108
			n_store1,      // rows columns storage
			n_store2,	  // band+perm in band storage
			n_store54,
			n_fib,		  // ib storage
			n_rows,		  // usuall 9 but can be 6 or 6
			wrc_store1,   // row studied from strore1
			wp_store1;    //perm studied from store1   
	BF16    bf16_store1;  // pattern equivalent to max_store1 

	FIRST_IN_BAND fibw,  t_fib[432]; // usually small but 216 x 2 if 9 clues

	BAND_STACK   bsw,       // to prepare the move
		         bs[2000];   // all permutation leading to max_store2 usually <10
	                       // but should limit the band number of clues

	TWO_BANDS w_2b,t_2b[100]; // all max for 2 bands usually < 5

	WSYM	wsym; // working area for search of valid perms
	//BF81    patx;
	//PUZC	pze;

	char	puz_in[82],        // pointer to pze_puz
			puz_pat[81],   //
		    w_store2[28],  //27 + string end
			max_store2[28],   //27 + string end
			w_store54[55],
			max_store54[55],
	        max_store81[82],
			w_store81[82];
	int		nextdigit;   // for relabeling
    GO_CAN();             // initial data in the class
	void InitPuzStd(char *ze ); // ze in the entry 81 or motr
	void InitPattern(char *ze); // ze in the entry 81 or motr
	int Canonical();  // standard entry for full caninical form
	//inline char * Canonical(PUZC pe){pze=pe; Canonical();return max_store81;} // for internal call 
	void Canonical_Band1();  // find all possible band 1
	void Canonical_Band2();  // find all possible band 1 + band 2
	void Canonical_Band3();  // find all possible band 1 + band 2 / band 3 final step
	void FindStarts(int i);
	void Store_Max_Pattern(int ie);
	void Store_Max_Band(int istore1);
	void Store_Max_Band(int r2, int r3);
	int  Store_2_Bands(int ra,int rb,int rc);
	void Store_3_Bands(int ra,int rb,int rc);
	void New_Store1();  // prepare perms for one entry of the store1 table

	//void Do_10();		// batch process standard
	//void Do_11();		// batch add ;can to entry
	//void Do_12();		// only can in output
	//void Do_20();		// can on clues
	//void Do_21();		// add can on clues "11.11" to entry
	//void Do_30();		// max or minlex no perm
	//void Do_31();		// max or minlex on pat
	//void Do_99();		// tets in can mode

	void Study_Perms_Sym(char * pate);
	void SPS_Row1();
	void SPS_Last(USHORT ir2a,USHORT ir2b);
	void MaxMinLex(char * puz);

}go_can;

GO_CAN::GO_CAN(){
   for(int i=0;i<20;i++) {
		bf16_directs[bf16_starts_row1[i][0]]=bf16_starts_row1[i][1];
	//	cout <<bf16_starts_row1[i][0]<<" "<<bf16_starts_row1[i][1]<<endl;
   }
  // puz_in=pze.puz;
   w_store2[27]=max_store2[27]=w_store54[54]=max_store54[54]
                =max_store81[81]=w_store81[81]=puz_in[81]=0;
   for(int i=0;i<9;i++)
	    relab_init.v[i]=0;
}


int GO_CAN::iii[6][3]={0,1,2, 0,2,1, 1,2,0, 1,0,2, 2,0,1, 2,1,0};

 /* table of possible pattern for row/column 1
 */
 int GO_CAN::bf16_starts_row1[20][2] ={
	       (3<<4)			,0x1c0,	  //111 000 000
		   (3<<4)+(1<<2)	,0x1e0,   //111 100 000
		   (3<<4)+(1<<2)+1	,0x1e4,   //111 100 100
		   (3<<4)+(2<<2) 	,0x1f0,   //111 110 000 
		   (3<<4)+(2<<2)+1	,0x1f4,   //111 110 100 
		   (3<<4)+(2<<2)+2	,0x1f6,   //111 110 110 
		   (3<<4)+(3<<2)	,0x1f8,   //111 111 000 
		   (3<<4)+(3<<2)+1	,0x1fc,   //111 111 100 
		   (3<<4)+(3<<2)+2	,0x1fe,   //111 111 110 
		   (3<<4)+(3<<2)+3	,0x1ff,   //111 111 111
		   (2<<4)			,0x180,   //110 000 000
		   (2<<4)+(1<<2)	,0x1a0,   //110 100 000
		   (2<<4)+(1<<2)+1	,0x1a4,   //110 100 100
		   (2<<4)+(2<<2)	,0x1b0,   //110 110 000
		   (2<<4)+(2<<2)+1	,0x1b4,   //110 110 100
		   (2<<4)+(2<<2)+2	,0x1b6,   //110 110 110
		   (1<<4)+(1<<2)	,0x120,   //100 100 000
		   (1<<4)+(1<<2)+1	,0x124,   //100 100 100
 };
//=================================== functions called
 /* Store_Max_Pattern  processes a row or a column
	 . finds the number of clues per triplet row/box ...
	 . chexk wheter the six permutations can be a start
	 . if yes, store the entry in a table and keep the start value in max_store1
 */
  void GO_CAN::Store_Max_Pattern(int ie) {
	 char * puz = v9pat_0based[ie].v;
	 UCHAR bb[3] = { 0,0,0 }, res;
	 for (int i = 0, k = 0; i < 3; i++)
		 for (int j = 0; j < 3; j++)
			 if (puz[k++])
				 bb[i]++;
	 for (int i = 0; i < 6; i++) {
		 int * iix = iii[i];
		 res = (bb[iix[0]] << 4) + (bb[iix[1]] << 2) + bb[iix[2]];

		 if (res > max_store1) { // reset to null if new max
			 n_store1 = 0;
			 max_store1 = res;
		 }
		 if (res == max_store1) // and store if equal
			 t_store1[n_store1++] = (i << 5) + ie;/// rrrrrreeeee  stored	
	 }
 }
  /* create a table of perm for the first line
	 bf16_store1 is the pattern of digits to find
	 wrc_store1,   is the row/col studied
	 wp_store1;    is the perm index studied (perm for block)
	 in each bloc, we try the 6 perms looking for
	   match with  bf16_store1

  */
  void GO_CAN::New_Store1() {
	  char * source = v9puz[wrc_store1].v;
	  char * spat = v9pat_0based[wrc_store1].v;
	  //============= bloc1
	  int base91 = 3 * iii[wp_store1][0];
	  char * pos_bloc1 = &spat[base91];
	  for (int iperm1 = 0; iperm1 < 6; iperm1++) {// 6 perms per bloc
		  int * pw1 = iii[iperm1];
		  int aig1 = 1;
		  for (int j1 = 0, ji1 = 8; j1 < 3; j1++, ji1--) {
			  int iw1 = (bf16_store1.f & (1 << ji1)) ? 1 : 0,
				  ir1 = pw1[j1],
				  cur_pos1 = base91 + ir1;
			  if (iw1 ^ pos_bloc1[ir1]) { // stop if it does not match 
				  aig1 = 0;
				  break;
			  }
			  // if match, put it in output  
			  fibw.perm.v[j1] = cur_pos1;
			  fibw.puz.v[j1] = source[cur_pos1];
		  }
		  if (!aig1)  //not a valid perm
			  continue;
		  //==================bloc 2
		  int base92 = 3 * iii[wp_store1][1];
		  char * pos_bloc2 = &v9pat_0based[wrc_store1].v[base92];
		  for (int iperm2 = 0; iperm2 < 6; iperm2++) {// 6 perms per bloc
			  int * pw2 = iii[iperm2];
			  int aig2 = 1;
			  for (int j2 = 0, ji2 = 5; j2 < 3; j2++, ji2--) {
				  int iw2 = (bf16_store1.f & (1 << ji2)) ? 1 : 0,
					  ir2 = pw2[j2],
					  cur_pos2 = base92 + ir2;
				  if (iw2 ^ pos_bloc2[ir2]) { // stop if it does not match 
					  aig2 = 0;
					  break;
				  }
				  // if match, put it in output  
				  fibw.perm.v[3 + j2] = cur_pos2;
				  fibw.puz.v[3 + j2] = source[cur_pos2];
			  }
			  if (!aig2)  //not a valid perm
				  continue;
			  //================== bloc 3
			  int base93 = 3 * iii[wp_store1][2];
			  char * pos_bloc3 = &v9pat_0based[wrc_store1].v[base93];
			  for (int iperm3 = 0; iperm3 < 6; iperm3++) {// 6 perms per bloc
				  int * pw3 = iii[iperm3];
				  int aig3 = 1;
				  for (int j3 = 0, ji3 = 2; j3 < 3; j3++, ji3--) {
					  int iw3 = (bf16_store1.f & (1 << ji3)) ? 1 : 0,
						  ir3 = pw3[j3],
						  cur_pos3 = base93 + ir3;
					  if (iw3 ^ pos_bloc3[ir3]) { // stop if it does not match 
						  aig3 = 0;
						  break;
					  }
					  // if match, put it in output  
					  fibw.perm.v[6 + j3] = cur_pos3;
					  fibw.puz.v[6 + j3] = source[cur_pos3];
				  }
				  if (!aig3) // not a valid perm
					  continue;
				  t_fib[n_fib++] = fibw;  // we store it 
			  }
		  }
	  }
  }

  /* first row of the bloc defined in fibw
	 add r2 r3
	 relabel and compare to maxstore2
  */
  void GO_CAN::Store_Max_Band(int wrc2, int wrc3) {
	  strcpy(w_store2, fibw.puz.v);// first row
	   // next rows reordered
	  char *vi = fibw.perm.v,
		  *r1 = v9puz[wrc2].v,
		  *r2 = v9puz[wrc3].v;
	  for (int i = 0; i < 9; i++) {
		  int col = vi[i];
		  w_store2[9 + i] = r1[col];
		  w_store2[18 + i] = r2[col];
	  }

	  // relabelling before compare
	  relab = relab_init;  // all digits to 0
	  nextdigit = '9';
	  for (int i = 0; i < 27; i++) if (w_store2[i] - '.') {
		  int c = w_store2[i] - '1';
		  if (!relab.v[c])
			  relab.v[c] = nextdigit--;
		  w_store2[i] = relab.v[c];
	  }

	  int ir = strcmp(w_store2, max_store2);
	  if (ir < 0) return;
	  if (ir > 0) { //new max-store2, reset all
		  n_store2 = 0;
		  strcpy(max_store2, w_store2);
	  }
	  // store if >=
	  bsw.fib = fibw;
	  bsw.rows[0] = wrc2;
	  bsw.rows[1] = wrc3;
	  bsw.relab = relab;
	  bsw.nextdigit = nextdigit;
	  if (n_store2 < 2000) // this should never happen, to check
		  bs[n_store2++] = bsw;
  }
  void GO_CAN::Store_Max_Band(int vstore1) {
	  n_fib = 0;
	  wrc_store1 = vstore1 & 31;
	  wp_store1 = vstore1 >> 5;
	  fibw.rc = wrc_store1;
	  USHORT irb1 = wrc_store1 / 3, irb2 = wrc_store1 % 3; // locate the three rows/columns
			 // irb1 is the band indexl irb2 the relative position in the band
	  USHORT w = 3 * irb1, rel2 = w + ((irb2 + 1) % 3), rel3 = w + ((irb2 + 2) % 3);
	  // we start reordering Line 1 and finding all permutations
	  New_Store1();
	  for (int i = 0; i < n_fib; i++) {
		  fibw = t_fib[i];
		  Store_Max_Band(rel2, rel3);
		  Store_Max_Band(rel3, rel2);
	  }
  }
  void GO_CAN::Store_3_Bands(int ira, int irb, int irc) {
	  strcpy(w_store81, max_store54);// first block relabelled
	   // next rows reordered
	  char *vi = bsw.fib.perm.v,
		  *ra = v9puz[ira].v,
		  *rb = v9puz[irb].v,
		  *rc = v9puz[irc].v;
	  for (int i = 0; i < 9; i++) {
		  int col = vi[i];
		  w_store81[54 + i] = ra[col];
		  w_store81[63 + i] = rb[col];
		  w_store81[72 + i] = rc[col];

	  }
	  // relabelling before compare positions 27 to 54
	  relab = w_2b.relab;  // restart from last position
	  nextdigit = w_2b.nextdigit;

	  for (int i = 54; i < 81; i++) if (w_store81[i] - '.') {
		  int c = w_store81[i] - '1';
		  if (!relab.v[c])
			  relab.v[c] = nextdigit--;
		  w_store81[i] = relab.v[c];
	  }

	  int ir = strcmp(w_store81, max_store81);
	  if (ir < 0) return;
	  if (ir > 0)
		  strcpy(max_store81, w_store81);
  }

  /* add block 2 rows ra rb rc see whether it is a new max 54
  */
  int GO_CAN::Store_2_Bands(int ira, int irb, int irc) {
	  strcpy(w_store54, max_store2);// first block relabelled
	   // next rows reordered
	  char *vi = bsw.fib.perm.v,
		  *ra = v9puz[ira].v,
		  *rb = v9puz[irb].v,
		  *rc = v9puz[irc].v;
	  for (int i = 0; i < 9; i++) {
		  int col = vi[i];
		  w_store54[27 + i] = ra[col];
		  w_store54[36 + i] = rb[col];
		  w_store54[45 + i] = rc[col];

	  }

	  // relabelling before compare positions 27 to 54
	  relab = bsw.relab;  // restart from last position
	  nextdigit = bsw.nextdigit;

	  for (int i = 27; i < 54; i++) if (w_store54[i] - '.') {
		  int c = w_store54[i] - '1';
		  if (!relab.v[c])
			  relab.v[c] = nextdigit--;
		  w_store54[i] = relab.v[c];
	  }


	  int ir = strcmp(w_store54, max_store54);
	  if (ir < 0) return 0;
	  if (ir > 0) { //new max-store2, reset all
		  n_store54 = 0;
		  strcpy(max_store54, w_store54);
	  }
	  return 1; // storage done by the caller
  }

  //======================== main steps
 void GO_CAN::InitPuzStd(char *ze) { // ze in the entry 81 or more
	 strncpy(puz_in, ze, 81);
	 for (int i = 0; i < 81; i++)
		 if (puz_in[i] > '0' &&  puz_in[i] <= '9') puz_pat[i] = 1;
		 else {	 puz_in[i] = '.'; puz_pat[i] = 0;	 }
	 max_store1 =n_store1 =0;
	 n_rows = 9;
 }
 void GO_CAN::InitPattern(char *ze) { // ze in the entry 81 or more
	 strncpy(puz_in, ze, 81);
	 for (int i = 0; i < 81; i++)
		 if (puz_in[i] > '0' &&  puz_in[i] <= '9') 
	 {	 puz_pat[i] = 1;  puz_in[i] = '1';	 }
		 else { puz_in[i] = '.'; puz_pat[i] = 0; }
	 max_store1 = 0;
	 n_rows = 9;
 }

 void GO_CAN::Canonical_Band1() {
	 int n_anal;
	 if (n_rows == 9) {
		 // store 18 rows cols in v9pat
		 for (int i = 0; i < 9; i++) {
			 v9puz[i].row(puz_in, i);
			 v9puz[i + 9].col(puz_in, i);
			 v9pat_0based[i].row(puz_pat, i);
			 v9pat_0based[i + 9].col(puz_pat, i);
			 n_anal = 18;
		 }
	 }
	 else {
		 for (int i = 0; i < n_rows; i++) {
			 v9puz[i].row(puz_in, i);
			 v9pat_0based[i].row(puz_pat, i);
			 n_anal = n_rows;
		 }
	 }
	 for (int i = 0; i < n_anal; i++)
		 Store_Max_Pattern(i);
	 if (max_store1 < 20) // should never be with reasonable puzzles
		 return;
	 bf16_store1.f = bf16_directs[max_store1];
	 // process now all entries of the t_store1
	 max_store2[0] = 0; // init comp
	 n_store2 = 0;
	 for (int i = 0; i < n_store1; i++) {
		 Store_Max_Band(t_store1[i]);
	 }
 }
 void GO_CAN::Canonical_Band2() {
	 max_store54[0] = 0;
	 // explore now all active band 1 to sort out bi band actives
	 for (int i = 0; i < n_store2; i++) {
		 bsw = bs[i];
		 int  ib1 = bsw.rows[0] / 3, //first band_stack index
			 g1 = 3 * (ib1 / 3),   // 0 if row, 3 if column
			 g2 = ib1 % 3;   // relative bloc in rows and col
		 for (int ib = 0; ib < 2; ib++) {
			 int ibx = (ib) ? g1 + ((g2 + 2) % 3) : g1 + ((g2 + 1) % 3),
				 rx[3];
			 rx[0] = 3 * ibx;
			 rx[1] = rx[0] + 1;
			 rx[2] = rx[1] + 1;
			 // must test now the six perms on rows or columns
			 for (int ip = 0; ip < 6; ip++) {
				 int ra = rx[iii[ip][0]], rb = rx[iii[ip][1]], rc = rx[iii[ip][2]];

				 if (Store_2_Bands(ra, rb, rc)) {// a new max54 to store
					 w_2b.ibs = i;// pointer to the band_stack source
					 w_2b.ra = ra;
					 w_2b.rb = rb;
					 w_2b.rc = rc;
					 w_2b.relab = relab;
					 w_2b.nextdigit = nextdigit;
					 if (n_store54 < 100)	 t_2b[n_store54++] = w_2b;
				 }
			 }
		 }
	 }
 }
 void GO_CAN::Canonical_Band3() {
	 max_store81[0] = 0;

	 // explore now all active band 1 to sort out bi band actives

	 for (int i = 0; i < n_store54; i++) {
		 w_2b = t_2b[i];
		 bsw = bs[w_2b.ibs];
		 // last band in rx
		 int  ib1 = bsw.rows[0] / 3, //first band_stack index
			 ib2 = w_2b.ra / 3,
			 ib3 = (ib1 > 2) ? (12 - ib1 - ib2) : (3 - ib1 - ib2),
			 rx[3];
		 rx[0] = 3 * ib3;
		 rx[1] = rx[0] + 1;
		 rx[2] = rx[1] + 1;
		 // must test now the six perms on rows or columns
		 for (int ip = 0; ip < 6; ip++) {
			 int ra = rx[iii[ip][0]], rb = rx[iii[ip][1]], rc = rx[iii[ip][2]];
			 Store_3_Bands(ra, rb, rc);
		 }

	 }


 }
 int GO_CAN::Canonical() {// entry is a normalised puzzle
	 Canonical_Band1();
	 if (max_store1 < 20)   // should never be
		 return 0;
	 Canonical_Band2();
	// cout << "n_store1=" << n_store1 << " n_store2="<<n_store2<<" n_store54="<<n_store54<<endl;
	 Canonical_Band3();
	 return 1;
 }

 //===============================================  process
 void Go_c310() {//standard process 
	 cout << "Go_310  canonical standard entry " << sgo.finput_name << endl;
	 finput.open(sgo.finput_name);
	 if (!finput.is_open()) {
		 cerr << "error open file " << sgo.finput_name << endl;
		 return;
	 }
	 char * ze = finput.ze;
	 int nseq = 0;
	 while (finput.GetLigne()) {
		 go_can.InitPuzStd(ze);
		 if (go_can.Canonical()) {
			 if (sgo.vx[0])		 fout1 << ze << ";";
			 fout1 << go_can.max_store81;
			 if (sgo.s_strings[0])
				 fout1 << ";" << sgo.s_strings[0];
			 fout1 << endl;
		 }
	 }

 }
 void Go_c320() {//can on pattern
	 cout << "Go_320  canonical on pattern " << sgo.finput_name << endl;
	 finput.open(sgo.finput_name);
	 if (!finput.is_open()) {
		 cerr << "error open file " << sgo.finput_name << endl;
		 return;
	 }
	 char * ze = finput.ze;
	 int nseq = 0;
	 while (finput.GetLigne()) {
		 go_can.InitPattern(ze);
		 if (go_can.Canonical()) {
			 if (sgo.vx[0])		 fout1 << ze << ";";
			 char *z = go_can.max_store81;
			 for (int i = 0; i < 81; i++)
				 if (z[i] == '9')z[i] = '1';
			 fout1 << z;
			 if (sgo.s_strings[0])
				 fout1 << ";" << sgo.s_strings[0];
			 fout1 << endl;
		 }
	 }

 }
 /*
 
void Go_c20(){ // can (2 longs on clues
	char * ze=myin->ze;
	if(add_fam_id)	if(idpos<82) idpos=0;
	n_rows=9;
	int nseq=0;
	while(myin->GetPuzzle(puz_in )){
		(*myout1)<<puz_in <<";";
		for(int i=0;i<81;i++) 
			if(puz_in[i]-'.')
				puz_in[i]='1';
		if(Canonical( )){
			int p1=0,p2=0,p3=0;  
			for(int i=0;i<27;i++){
				if(max_store81[i]-'.')
					p1 |= 1<<i;
				if(max_store81[i+27]-'.')
					p2 |= 1<<i;   
				if(max_store81[i+54]-'.')
					p3 |= 1<<i;   
			}
			(*myout1)<<p1 <<";"	<<p2<<";"<<p3<<endl;
		}
		else 	(*myout1)<<"0;0;0"<<endl;
	}	
}


void GO_CAN::Do_21(){ // can  in text mode added to output
	// entry must be standard
	char * ze=myin->ze;
	n_rows=9;
	int nseq=0;
	while(myin->GetLigne( )){
		if(strlen(ze)<81) 			continue;
		for(int i=0;i<81;i++) 
			if(ze[i]-'.' && ze[i] - '0')	puz_in[i]='1';
			else puz_in[i]='.';
		if(Canonical( )){
			for(int i=0;i<81;i++) 
				if(max_store81[i] == '9')	max_store81[i]='1';
			if(options.cop4){
				if (!strcmp(ze,max_store81)) // was maxlex
					(*myout1)<<max_store81<<endl;
			}
			else (*myout1)<<ze <<";"	<<max_store81<<endl;
		}
	}	
}
*/