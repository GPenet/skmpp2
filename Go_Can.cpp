



/*
void GO_CAN::Start(){
	
	switch (cop23){
		case 20: Do_20(); 	return;	//can on clues	
		case 21: Do_21(); 	return;	//can on clues	text mode added to entry
		case 30: Do_30();   return; // maxlex or minlex  no change  
		case 31: Do_31();   return; // maxlex or minlex  on pattern
		case 99: Do_99();	return; // test in can mode
	}

}


*/
 
/* recherche des perms possibles pour une même pattern
  On termine par la table des perms autorisées en forme
   .liste de lignes
   .liste de colonnes
   .Indicateur départ lignes ou colonnes

   pate is a 0 based puzzle pattern

   Processing pairs "band stack" same number of clues as first
   Or start with lowest band stack having highest number of clues

*/
/*
void GO_CAN::Study_Perms_Sym(char * pate){

	wsym.nfinal=0; // no perm
// set up datat for later  processing
	patx.SetAll_0();
	for(int i=0;i<81;i++) 
		if(pate[i]-'.') {
			puz_pat[i]=1;
			patx.Set(i);
		}
		else 
			puz_pat[i]=0;

	for(int i=0;i<9;i++){
		v9pat_0based[i].row(puz_pat,i);
		v9pat_0based[i+9].col(puz_pat,i);
	}

// find all bands/stacks with max count.

	int maxband=0,tband[6],nband=0;
	for(int i=0;i<6;i++){
		int cc=(patx & cellsInBandBM[i]).Count();
		if(cc>maxband){
			maxband=cc;
			wsym.bandb=i;
			nband=0;
		}
		if(cc==maxband) 
			tband[nband++]=i;
	}

	wsym.Prepare1(patx);
 
	if(0){
		(*myout1)<<"bandb="<<wsym.bandb
			<<" maxband="<<maxband<<" nband="<<nband<<"->";
		for(int i=0;i<nband;i++)
				(*myout1)<<tband[i]<<"  ";
		(*myout1)<<endl;
		(*myout1)<<"maxrow="<<wsym.maxrow<< " for row" << wsym.ir[0]
		         << " ir[1]" << wsym.ir[1] << " ir[2]"<< wsym.ir[2]<<endl;

		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
			(*myout1)<<wsym.cptrg[i][j]<<" ";
			}
			(*myout1)<< " total="<<wsym.cptr[i] <<endl;
		}
		
				 
		//return;
	}

		// process all pairs with bandb including bandb itself.
	for(int i=0;i<nband;i++){
		wsym.band2=tband[i];
		wsym.Prepare2();
		int ibase2=3*wsym.band2;
		for(int j=0,jx=ibase2;j<3;j++,jx++){
			if((patx & cellsInHouseBM[jx]).Count()-wsym.maxrow)continue;
				//look for perm(s) matching if any
			wsym.ir2[0]=jx;
			wsym.ir2[1]=ibase2+((j+1)%3);
			wsym.ir2[2]=ibase2+((j+2)%3);
			BF81 * tg3r = & zgs.zg3[3*jx];
			int cptw1=wsym.cptrg[0][0],
				cptw2=wsym.cptrg[0][1],
				cptw3=wsym.cptrg[0][2];
			for(wsym.ig1=0;wsym.ig1<3;wsym.ig1++){
				if((patx & tg3r[wsym.ig1]).Count()-cptw1) continue;
					// count matching use it   
				for(wsym.ig2=0;wsym.ig2<3;wsym.ig2++){
					if(wsym.ig2==wsym.ig1) continue;
					if((patx & tg3r[wsym.ig2]).Count()-cptw2) continue;
					wsym.ig3=3-wsym.ig1-wsym.ig2;
					if((patx & tg3r[wsym.ig3]).Count()-cptw3) continue;
					SPS_Row1(); 
				}
			}
		}
	}
}
*/
	/* continued after row1 groups matching
	   ir[0] matches to ir2[0] for groups
	   ig1,ig2,ig3
	   find all valid perms for row 1
	*/
/*
void GO_CAN::SPS_Row1(){
	if(0){
		(*myout1)<<"entry row1 bandb="<<wsym.bandb 
		 << " matches with "<<wsym.band2 <<" row "  << wsym.ir2[0]+1
		 << " for igx " << wsym.ig1 << " " <<  wsym.ig2<< " " <<  wsym.ig3 << " rows order ";
		 for(int i=0;i<9;i++)
			 (*myout1)<< wsym.ir[i]+1<<" ";
		 (*myout1)<<endl;


	}

	char * za=v9pat_0based[wsym.ir[0]].v,
		 * zb=v9pat_0based[wsym.ir2[0]].v;
	int ida=3*wsym.ig1,idb=3*wsym.ig2,idc=3*wsym.ig3;
	USHORT r2a=wsym.ir2[1],r2b=wsym.ir2[2];

	for(int ipa=0;ipa<6;ipa++){
		int aiga=0;
		for(int ia=0;ia<3;ia++){
			USHORT xa=ida+iii[ipa][ia];
			wsym.wp.perm2[ia]=xa;
			if(za[ia]==zb[xa]) continue;
			aiga=1;
			break;
		}
		if(aiga) continue; // does not fit
			// perm a OK look for permb
		for(int ipb=0;ipb<6;ipb++){
			int aigb=0;
			for(int ib=0;ib<3;ib++){
				USHORT xb=idb+iii[ipb][ib];
				wsym.wp.perm2[3+ib]=xb;
				if(za[3+ib]==zb[xb]) continue;
				aigb=1;
				break;
			}
			if(aigb) continue; // does not fit
			// perm a+b  OK look for perm c
			for(int ipc=0;ipc<6;ipc++){
				int aigc=0;
				for(int ic=0;ic<3;ic++){
					USHORT xc=idc+iii[ipc][ic];
					wsym.wp.perm2[6+ic]=xc;
					if(za[6+ic]==zb[idc+iii[ipc][ic]]) continue;
					aigc=1;
					break;
				}
				if(aigc) continue; // does not fit
					// perm   OK try that perm for other (rows)
				wsym.wp.perm1[0]=wsym.ir2[0];
				if(0){
					(*myout1)<<"fit a row1 for unit "<<wsym.ir2[0]<< " ";
					for(int ic=0;ic<9;ic++)
						(*myout1)<<wsym.wp.perm2[ic]+1<<'_';
					
					(*myout1)<<endl;

				}
				SPS_Last(r2a,r2b);
				SPS_Last(r2b,r2a);
			}
		}
	}
//	(*myout1)<<"nfinal="<<wsym.nfinal<<endl;
}
*/
/* after bandb is known
   do all possible assignments
*/
/*
void GO_CAN::WSYM::Prepare1(BF81 & patx){
			// select as start the highest count in the band
	int irel,ibase=3*bandb;
	maxrow=0;
	BF81 *tg3=zgs.zg3;  // 54 triplets row/box col/box

	for(int i=0,ix=ibase;i<3;i++,ix++){
		int cc=(patx & cellsInHouseBM[ix]).Count();
		if(cc>maxrow){
			maxrow=cc;
			irel=i;
		}
	}
	ir[0]=ibase+irel;
	ir[1]=ibase+((irel+1)%3);
	ir[2]=ibase+((irel+2)%3);
	wp.target_is_row=1;
	if(ir[0]>9)
		wp.target_is_row=0;
	for(int i=0;i<3;i++){ // compute all counts row and groups
		cptr[i]=(patx & cellsInHouseBM[ir[i]]).Count();
		BF81 * tg3i= & tg3 [3*ir[i]];
		for(int j=0;j<3;j++)
			cptrg[i][j]=(patx & tg3i[j]).Count();
		
	}
	  // set up all other (rows or columns)  3 to 9
    int irr=3,ibd=(bandb>2)?3:0,ibf=ibd+3;
	for(int i=ibd;i<ibf;i++) if(i-bandb){
		for(int j=0,ij=3*i;j<3;j++,ij++) // ij is row or col 
			ir[irr++]=ij;
	}

}
*/
/* after band2 is known
   do all possible assignments
*/
/*
void GO_CAN::WSYM::Prepare2(){
	  // set up 2 band perms in otherband2
    int ipp=0,ibd=(band2>2)?3:0,ibf=ibd+3;
	for(int i=ibd;i<ibf;i++) if(i-band2){
		otherband2[0][ipp]=i;
		otherband2[1][1-ipp++]=i;
	}


}
*/
/* apply wp.perm2 to the row/col pointed by ir2
   then compare to ir1 stop at first <>
*/
/*
int GO_CAN::WSYM::Check(int a,int b){

	for(int i=a;i<=b;i++){
		char * xa=parent->v9pat_0based[ir[i]].v,
			 * xb=parent->v9pat_0based[ir2[i]].v;
		for(int j=0;j<9;j++)
			if(xa[j]-xb[wp.perm2[j]]) 
				return 1;
	}
	return 0;
}
*/
/* store the current situation as valid perm
   except if no change
   ir1 is the target changed to 0_8
   ir2 will go in wp.perm1
   wp.perm2 is already OK
*/
/*
void GO_CAN::WSYM::Store(){
	int aigperm=0;
	if(wp.perm1[0]>8) 
		wp.target_is_row=0;
	for(int i=0;i<9;i++)  // check perm2 first
		if(wp.perm2[i]-i){
			aigperm=1;
			break;
		}
	for(int i=0;i<9;i++) { // compute and check perm1
		int i1=ir[i],i2=ir2[i];
		if(i1-i2)
			aigperm=1;
		if(i1>8) i1-=9; // adjust index if column
		wp.perm1[i1]=i2;

	}

	if(aigperm && nfinal<100)
		tp[nfinal++]=wp;
}

*/
/* first (row) matches
   build all perms for other (rows)
   and check if all is ok
*/
/*
void  GO_CAN::SPS_Last(USHORT ir2a,USHORT ir2b){
	wsym.ir2[1]=ir2a;
	wsym.ir2[2]=ir2b;
	if(wsym.Check(1,2)) return;
	if(0){
		(*myout1)<<"passe 1-ir2a="<<ir2a <<" ir2b= " <<ir2b << " for ";
		for(int ic=0;ic<9;ic++)
					(*myout1)<<wsym.wp.perm2[ic]+1;
					
		(*myout1)<<endl;
	}
	  // and test all perms 
	for(int ib=0;ib<2;ib++){// 2 perms on band stack
		USHORT * bb=wsym.otherband2[ib];
		USHORT bbx1=3*bb[0];
		for(int ip1=0;ip1<6;ip1++){// 6 perms for first band / stack
			for(int i=0;i<3;i++)
				wsym.ir2[3+i]=bbx1+iii[ip1][i];
			if(wsym.Check(3,5)) continue;

			USHORT bbx2=3*bb[1];

			for(int ip2=0;ip2<6;ip2++){// 6 perms for second band / stack
				for(int i=0;i<3;i++)
					wsym.ir2[6+i]=bbx2+iii[ip2][i];
				if(wsym.Check(6,8)) continue;
					// this is a perm with unchanged pattern store it
				if(0){
					(*myout1)<<"passe 3 for ir1 ";
					for(int ic=0;ic<9;ic++)
						(*myout1)<<wsym.ir[ic]+1;
					(*myout1)<<" ir2 ";
					for(int ic=0;ic<9;ic++)
						(*myout1)<<wsym.ir2[ic]+1;
					(*myout1)<<" order ";
					for(int ic=0;ic<9;ic++)
						(*myout1)<<wsym.wp.perm2[ic]+1;
	
					(*myout1)<<endl;
				}	

				// new perm store it except if no change 
				wsym.Store();

			}			
		}
	}
}



void GO_CAN::Do_99(){
	char * ze=myin->ze;
	int aig1=1;
	while(myin->GetPuzzle(puz_in )){
		  // look for perms for the first one
		(*myout1)<<puz_in <<"study perms"	<<endl;
		if(aig1){
			Study_Perms_Sym(puz_in );
			aig1=0;
		}
			// now find maxtext or mintext

	}	
}

void GO_CAN::Do_31(){
	char * ze=myin->ze;
	PUZC puzf,puzw;
	int aig1=1;
	while(myin->GetPuzzle(puz_in )){
		  // look for perms for the first one
		if(aig1){
			puzf=pze;
			Study_Perms_Sym(puzf.puz );
			aig1=0;
			(*myout2)<<puz_in <<"study pattern"	<<endl;
			(*myout2)<<"nfinal="<<wsym.nfinal<<endl;

		}
			// now find maxtext or mintext
		puzw=pze;
		MaxMinLex(puzw.puz);
		puzf=puzw;
		(*myout2)<< " first "<<puzw.puz <<endl;

		for(int i=0;i<wsym.nfinal;i++){
			puzw=pze;
			 WSYM_PERM wp=wsym.tp[i];
			for(int j=0;j<9;j++){
				VV9 ws0,ws; // ws source reordered
				//create source vector
				int wrc=wp.perm1[j];
				if(wrc<9)
					ws0.row(puz_in,wrc);
				else
					ws0.col(puz_in,wrc-9);
				for(int k=0;k<9;k++) 
					ws.v[k]=ws0.v[wp.perm2[k]];
				for(int k=0;k<9;k++){
					int id;
					if(wp.target_is_row)
						id=9*j+k;
					else
						id=j+9*k;
					puzw.puz[id]=ws.v[k];
				}
			}
			MaxMinLex(puzw.puz);
			(*myout2)<< " perm "<<puzw.puz <<endl;
			int rcomp=strncmp(puzf.puz,puzw.puz,81);
			if((add_fam_id && rcomp>0) ||
				((!add_fam_id) && rcomp<0))
				puzf=puzw;

		}
		(*myout1)<<puzf.puz <<endl;

	}	
}

void GO_CAN::Do_30(){
	char * ze=myin->ze;
	while(myin->GetPuzzle(puz_in )){
			// now find maxtext or mintext
		MaxMinLex(puz_in);
		(*myout1)<<puz_in <<endl;

	}	
}
*/
/* change solution to maxlex 
     or minlex if add_fam_id on
	 puz is text based puzzle nomalized
*/
/*
void GO_CAN::MaxMinLex(char * puz){
	BF16 vus; USHORT nouveau=(add_fam_id)?0:8,tvus[9];
	for(int i=0;i<81;i++){
		if(puz[i]-'.')  {	  
			USHORT ic=puz[i]-'1';
			if(vus.Off(ic)) {
				 vus.Set(ic);
				tvus[ic]=nouveau;	
				 if(add_fam_id)
					nouveau++;
				 else
					 nouveau--;
			}
			puz[i]='1'+tvus[ic]; 
		}
	} 
}
*/