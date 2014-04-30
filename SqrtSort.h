/********* Sqrt sorting *********************************/
/*                                                       */
/* (c) 2014 by Andrey Astrelin                           */
/*                                                       */
/*                                                       */
/* Stable sorting that works in O(N*log(N)) worst time   */
/* and uses O(sqrt(N)) extra memory                      */
/*                                                       */
/* Define SORT_TYPE and SORT_CMP                         */
/* and then call SqrtSort() function                     */
/*                                                       */
/*********************************************************/

#include<memory.h>
#include<malloc.h>

inline void sqrtsort_swap1(SORT_TYPE *a,SORT_TYPE *b){
	SORT_TYPE c=*a;
	*a++=*b;
	*b++=c;
}
inline void sqrtsort_swapN(SORT_TYPE *a,SORT_TYPE *b,int n){
	while(n--) sqrtsort_swap1(a++,b++);
}


static void sqrtsort_MergeRight(SORT_TYPE *arr,int L1,int L2,int M){
	int p0=L1+L2+M-1,p2=L1+L2-1,p1=L1-1;

	while(p1>=0){
		if(p2<L1 || SORT_CMP(arr+p1,arr+p2)>0){
			arr[p0--]=arr[p1--];
		} else{
			arr[p0--]=arr[p2--];
		}
	}
	if(p2!=p0) while(p2>=L1) arr[p0--]=arr[p2--];
}

// arr[M..-1] - free, arr[0,L1-1]++arr[L1,L1+L2-1] -> arr[M,M+L1+L2-1]
static void sqrtsort_MergeLeftWithXBuf(SORT_TYPE *arr,int L1,int L2,int M){
	int p0=0,p1=L1; 
	L2+=L1;
	while(p1<L2){
		if(p0==L1 || SORT_CMP(arr+p0,arr+p1)>0) arr[M++]=arr[p1++];
		else arr[M++]=arr[p0++];
	}
	if(M!=p0) while(p0<L1) arr[M++]=arr[p0++];
}

// arr[0,L1-1] ++ arr2[0,L2-1] -> arr[-L1,L2-1],  arr2 is "before" arr1
static void sqrtsort_MergeDown(SORT_TYPE *arr,SORT_TYPE *arr2,int L1,int L2){
	int p0=0,p1=0,M=-L2;

	while(p1<L2){
		if(p0==L1 || SORT_CMP(arr+p0,arr2+p1)>=0) arr[M++]=arr2[p1++];
		else arr[M++]=arr[p0++];
	}
	if(M!=p0) while(p0<L1) arr[M++]=arr[p0++];
}

static void sqrtsort_SmartMergeWithXBuf(SORT_TYPE *arr,int *alen1,int *atype,int len2,int lkeys){
	int p0=-lkeys,p1=0,p2=*alen1,q1=p2,q2=p2+len2;
	int ftype=1-*atype;  // 1 if inverted
	while(p1<q1 && p2<q2){
		if(SORT_CMP(arr+p1,arr+p2)-ftype<0) arr[p0++]=arr[p1++];
		else arr[p0++]=arr[p2++];
	}
	if(p1<q1){
		*alen1=q1-p1;
		while(p1<q1) arr[--q2]=arr[--q1];
	} else{
		*alen1=q2-p2;
		*atype=ftype;
	}
}

// arr - starting array. arr[-lblock..-1] - buffer (if havebuf).
// lblock - length of regular blocks. First nblocks are stable sorted by 1st elements and key-coded
// keys - arrays of keys, in same order as blocks. key<midkey means stream A
// nblock2 are regular blocks from stream A. llast is length of last (irregular) block from stream B, that should go before nblock2 blocks.
// llast=0 requires nblock2=0 (no irregular blocks). llast>0, nblock2=0 is possible.
static void sqrtsort_MergeBuffersLeftWithXBuf(int *keys,int midkey,SORT_TYPE *arr,int nblock,int lblock,int nblock2,int llast){
	int l,prest,lrest,frest,pidx,cidx,fnext,plast;

	if(nblock==0){
		l=nblock2*lblock;
		sqrtsort_MergeLeftWithXBuf(arr,l,llast,-lblock);
		return;
	}

	lrest=lblock;
	frest=keys[0]<midkey ? 0 : 1;
	pidx=lblock;
	for(cidx=1;cidx<nblock;cidx++,pidx+=lblock){
		prest=pidx-lrest;
		fnext=keys[cidx]<midkey ? 0 : 1;
		if(fnext==frest){
			memcpy(arr+prest-lblock,arr+prest,lrest*sizeof(SORT_TYPE));
			prest=pidx;
			lrest=lblock;
		} else{
			sqrtsort_SmartMergeWithXBuf(arr+prest,&lrest,&frest,lblock,lblock);
		}
	}
	prest=pidx-lrest;
	if(llast){
		plast=pidx+lblock*nblock2;
		if(frest){
			memcpy(arr+prest-lblock,arr+prest,lrest*sizeof(SORT_TYPE));
			prest=pidx;
			lrest=lblock*nblock2;
			frest=0;
		} else{
			lrest+=lblock*nblock2;
		}
		sqrtsort_MergeLeftWithXBuf(arr+prest,lrest,llast,-lblock);
	} else{
		memcpy(arr+prest-lblock,arr+prest,lrest*sizeof(SORT_TYPE));
	}
}

// build blocks of length K
// input: [-K,-1] elements are buffer
// output: first K elements are buffer, blocks 2*K and last subblock sorted
static void sqrtsort_BuildBlocks(SORT_TYPE *arr,int L,int K){
	int m,u,h,p0,p1,rest,restk,p;
		for(m=1;m<L;m+=2){
			u=0;
			if(SORT_CMP(arr+(m-1),arr+m)>0) u=1;
			arr[m-3]=arr[m-1+u];
			arr[m-2]=arr[m-u];
		}
		if(L%2) arr[L-3]=arr[L-1];
		arr-=2;
		for(h=2;h<K;h*=2){
			p0=0;
			p1=L-2*h;
			while(p0<=p1){
				sqrtsort_MergeLeftWithXBuf(arr+p0,h,h,-h);
				p0+=2*h;
			}
			rest=L-p0;
			if(rest>h){
				sqrtsort_MergeLeftWithXBuf(arr+p0,h,rest-h,-h);
			} else {
				for(;p0<L;p0++)	arr[p0-h]=arr[p0];
			}
			arr-=h;
		}
	restk=L%(2*K);
	p=L-restk;
	if(restk<=K) memcpy(arr+p+K,arr+p,restk*sizeof(SORT_TYPE));
	else sqrtsort_MergeRight(arr+p,K,restk-K,K);
	while(p>0){
		p-=2*K;
		sqrtsort_MergeRight(arr+p,K,K,K);
	}
}


static void sqrtsort_SortIns(SORT_TYPE *arr,int len){
	int i,j;
	for(i=1;i<len;i++){
		for(j=i-1;j>=0 && SORT_CMP(arr+(j+1),arr+j)<0;j--) sqrtsort_swap1(arr+j,arr+(j+1));
	}
}

// keys are on the left of arr. Blocks of length LL combined. We'll combine them in pairs
// LL and nkeys are powers of 2. (2*LL/lblock) keys are guarantied
static void sqrtsort_CombineBlocks(SORT_TYPE *arr,int len,int LL,int lblock,int *tags){
	int M,nkeys,b,NBlk,midkey,lrest,u,i,p,v,kc,nbl2,llast;
	SORT_TYPE *arr1;

	M=len/(2*LL);
	lrest=len%(2*LL);
	nkeys=(2*LL)/lblock;
	if(M==0) nkeys=(M-1)/lblock+1;
	if(lrest<=LL){
		len-=lrest;
		lrest=0;
	}
	for(b=0;b<=M;b++){
		if(b==M && lrest==0) break;
		arr1=arr+b*2*LL;
		NBlk=(b==M ? lrest : 2*LL)/lblock;
		u=NBlk+(b==M ? 1 : 0);
		for(i=0;i<=u;i++) tags[i]=i;
		midkey=LL/lblock;
		for(u=1;u<NBlk;u++){
			p=u-1;
			for(v=u;v<NBlk;v++){
				kc=SORT_CMP(arr1+p*lblock,arr1+v*lblock);
				if(kc>0 || (kc==0 && tags[p]>tags[v])) p=v;
			}
			if(p!=u-1){
				sqrtsort_swapN(arr1+(u-1)*lblock,arr1+p*lblock,lblock);
				i=tags[u-1]; tags[u-1]=tags[p]; tags[p]=i;
			}
		}
		nbl2=llast=0;
		if(b==M) llast=lrest%lblock;
		if(llast!=0){
			while(nbl2<NBlk && SORT_CMP(arr1+NBlk*lblock,arr1+(NBlk-nbl2-1)*lblock)<0) nbl2++;
		}
		sqrtsort_MergeBuffersLeftWithXBuf(tags,midkey,arr1,NBlk-nbl2,lblock,nbl2,llast);
	}
	for(p=len;--p>=0;) arr[p]=arr[p-lblock];
}


void sqrtsort_commonSort(SORT_TYPE *arr,int Len,SORT_TYPE *extbuf,int *Tags){
	int lblock,cbuf;

	if(Len<16){
		sqrtsort_SortIns(arr,Len);
		return;
	}
	
	lblock=1;
	while(lblock*lblock<Len) lblock*=2;
	memcpy(extbuf,arr,lblock*sizeof(SORT_TYPE));

	sqrtsort_commonSort(extbuf,lblock,arr,Tags);

	sqrtsort_BuildBlocks(arr+lblock,Len-lblock,lblock);
	cbuf=lblock;
	while(Len>(cbuf*=2)){
		sqrtsort_CombineBlocks(arr+lblock,Len-lblock,cbuf,lblock,Tags);
	}
	sqrtsort_MergeDown(arr+lblock,extbuf,Len-lblock,lblock);
}

void SqrtSort(SORT_TYPE *arr,int Len){
	int L=1;
	SORT_TYPE *ExtBuf;
	int *Tags;

	while(L*L<Len) L*=2;
	int NK=(Len-1)/L+2;
	ExtBuf=(SORT_TYPE*)malloc(L*sizeof(SORT_TYPE));
	if(ExtBuf==NULL) return; // fail
 	Tags=(int*)malloc(NK*sizeof(int));
	if(Tags==NULL) return;

	sqrtsort_commonSort(arr,Len,ExtBuf,Tags);
	free(Tags);
	free(ExtBuf);
}

