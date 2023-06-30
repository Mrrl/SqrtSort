// sqrt-sort: Stable sorting with O(sqrt(N)) external memory.
// 
// MIT License:
// Copyright (c) 2023 The pysoft group.
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this softwareand associated documentation files(the
//    "Software"), to deal in the Software without restriction, including
//    without limitation the rights to use, copy, modify, merge, publish,
//    distribute, sublicense, and /or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so, subject to
//    the following conditions :
//
// The above copyright noticeand this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


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

#include<iterator>
#include<type_traits>


namespace sqrtsort {
	namespace _internal {
		template<typename It>
		using iter_value = typename std::iterator_traits<It>::value_type;

		template<typename It>
		inline void sqrtsort_swap1(It a, It b) {
			iter_value<It> c = *a;
			*a++ = *b;
			*b++ = c;
		}

		template<typename It>
		inline void sqrtsort_swapN(It a, It b, int n) {
			while (n--) sqrtsort_swap1(a++, b++);
		}

		template<typename It, typename Comp>
		static void sqrtsort_MergeRight(It arr, int L1, int L2, int M, Comp comp) {
			int p0 = L1 + L2 + M - 1, p2 = L1 + L2 - 1, p1 = L1 - 1;

			while (p1 >= 0) {
				if (p2 < L1 || comp(arr + p1, arr + p2)>0) {
					arr[p0--] = arr[p1--];
				}
				else {
					arr[p0--] = arr[p2--];
				}
			}
			if (p2 != p0) while (p2 >= L1) arr[p0--] = arr[p2--];
		}

		// arr[M..-1] - free, arr[0,L1-1]++arr[L1,L1+L2-1] -> arr[M,M+L1+L2-1]
		template<typename It, typename Comp>
		static void sqrtsort_MergeLeftWithXBuf(It arr, int L1, int L2, int M, Comp comp) {
			int p0 = 0, p1 = L1;
			L2 += L1;
			while (p1 < L2) {
				if (p0 == L1 || comp(arr + p0, arr + p1) > 0) arr[M++] = arr[p1++];
				else arr[M++] = arr[p0++];
			}
			if (M != p0) while (p0 < L1) arr[M++] = arr[p0++];
		}

		// arr[0,L1-1] ++ arr2[0,L2-1] -> arr[-L1,L2-1],  arr2 is "before" arr1
		template<typename It, typename Comp>
		static void sqrtsort_MergeDown(It arr, It arr2, int L1, int L2, Comp comp) {
			int p0 = 0, p1 = 0, M = -L2;

			while (p1 < L2) {
				if (p0 == L1 || comp(arr + p0, arr2 + p1) >= 0) arr[M++] = arr2[p1++];
				else arr[M++] = arr[p0++];
			}
			if (M != p0) while (p0 < L1) arr[M++] = arr[p0++];
		}

		template<typename It, typename Comp>
		static void sqrtsort_SmartMergeWithXBuf(It arr, int* alen1, int* atype, int len2, int lkeys, Comp comp) {
			int p0 = -lkeys, p1 = 0, p2 = *alen1, q1 = p2, q2 = p2 + len2;
			int ftype = 1 - *atype;  // 1 if inverted
			while (p1 < q1 && p2 < q2) {
				if (comp(arr + p1, arr + p2) - ftype < 0) arr[p0++] = arr[p1++];
				else arr[p0++] = arr[p2++];
			}
			if (p1 < q1) {
				*alen1 = q1 - p1;
				while (p1 < q1) arr[--q2] = arr[--q1];
			}
			else {
				*alen1 = q2 - p2;
				*atype = ftype;
			}
		}

		// arr - starting array. arr[-lblock..-1] - buffer (if havebuf).
		// lblock - length of regular blocks. First nblocks are stable sorted by 1st elements and key-coded
		// keys - arrays of keys, in same order as blocks. key<midkey means stream A
		// nblock2 are regular blocks from stream A. llast is length of last (irregular) block from stream B, that should go before nblock2 blocks.
		// llast=0 requires nblock2=0 (no irregular blocks). llast>0, nblock2=0 is possible.
		template<typename It, typename Comp>
		static void sqrtsort_MergeBuffersLeftWithXBuf(int* keys, int midkey, It arr, int nblock, int lblock, int nblock2, int llast, Comp comp) {
			int l, prest, lrest, frest, pidx, cidx, fnext, plast;

			if (nblock == 0) {
				l = nblock2 * lblock;
				sqrtsort_MergeLeftWithXBuf(arr, l, llast, -lblock, comp);
				return;
			}

			lrest = lblock;
			frest = keys[0] < midkey ? 0 : 1;
			pidx = lblock;
			for (cidx = 1; cidx < nblock; cidx++, pidx += lblock) {
				prest = pidx - lrest;
				fnext = keys[cidx] < midkey ? 0 : 1;
				if (fnext == frest) {
					std::copy_n(arr + prest, lrest, arr + prest - lblock);
					prest = pidx;
					lrest = lblock;
				}
				else {
					sqrtsort_SmartMergeWithXBuf(arr + prest, &lrest, &frest, lblock, lblock, comp);
				}
			}
			prest = pidx - lrest;
			if (llast) {
				plast = pidx + lblock * nblock2;
				if (frest) {
					std::copy_n(arr + prest, lrest, arr + prest - lblock);
					prest = pidx;
					lrest = lblock * nblock2;
					frest = 0;
				}
				else {
					lrest += lblock * nblock2;
				}
				sqrtsort_MergeLeftWithXBuf(arr + prest, lrest, llast, -lblock, comp);
			}
			else {
				std::copy_n(arr + prest, lrest, arr + prest - lblock);
			}
		}

		// build blocks of length K
		// input: [-K,-1] elements are buffer
		// output: first K elements are buffer, blocks 2*K and last subblock sorted
		template<typename It, typename Comp>
		static void sqrtsort_BuildBlocks(It arr, int L, int K, Comp comp) {
			int m, u, h, p0, p1, rest, restk, p;
			for (m = 1; m < L; m += 2) {
				u = 0;
				if (comp(arr + (m - 1), arr + m) > 0) u = 1;
				arr[m - 3] = arr[m - 1 + u];
				arr[m - 2] = arr[m - u];
			}
			if (L % 2) arr[L - 3] = arr[L - 1];
			arr -= 2;
			for (h = 2; h < K; h *= 2) {
				p0 = 0;
				p1 = L - 2 * h;
				while (p0 <= p1) {
					sqrtsort_MergeLeftWithXBuf(arr + p0, h, h, -h, comp);
					p0 += 2 * h;
				}
				rest = L - p0;
				if (rest > h) {
					sqrtsort_MergeLeftWithXBuf(arr + p0, h, rest - h, -h, comp);
				}
				else {
					for (; p0 < L; p0++)	arr[p0 - h] = arr[p0];
				}
				arr -= h;
			}
			restk = L % (2 * K);
			p = L - restk;
			if (restk <= K) std::copy_n(arr + p, restk, arr + p + K);
			else sqrtsort_MergeRight(arr + p, K, restk - K, K, comp);
			while (p > 0) {
				p -= 2 * K;
				sqrtsort_MergeRight(arr + p, K, K, K, comp);
			}
		}

		template<typename It, typename Comp>
		static void sqrtsort_SortIns(It arr, int len, Comp comp) {
			int i, j;
			for (i = 1; i < len; i++) {
				for (j = i - 1; j >= 0 && comp(arr + (j + 1), arr + j) < 0; j--) sqrtsort_swap1(arr + j, arr + (j + 1));
			}
		}

		// keys are on the left of arr. Blocks of length LL combined. We'll combine them in pairs
		// LL and nkeys are powers of 2. (2*LL/lblock) keys are guarantied
		template<typename It, typename Comp>
		static void sqrtsort_CombineBlocks(It arr, int len, int LL, int lblock, int* tags, Comp comp) {
			int M, nkeys, b, NBlk, midkey, lrest, u, i, p, v, kc, nbl2, llast;
			It arr1;

			M = len / (2 * LL);
			lrest = len % (2 * LL);
			nkeys = (2 * LL) / lblock;
			if (M == 0) nkeys = (M - 1) / lblock + 1;
			if (lrest <= LL) {
				len -= lrest;
				lrest = 0;
			}
			for (b = 0; b <= M; b++) {
				if (b == M && lrest == 0) break;
				arr1 = arr + b * 2 * LL;
				NBlk = (b == M ? lrest : 2 * LL) / lblock;
				u = NBlk + (b == M ? 1 : 0);
				for (i = 0; i <= u; i++) tags[i] = i;
				midkey = LL / lblock;
				for (u = 1; u < NBlk; u++) {
					p = u - 1;
					for (v = u; v < NBlk; v++) {
						kc = comp(arr1 + p * lblock, arr1 + v * lblock);
						if (kc > 0 || (kc == 0 && tags[p] > tags[v])) p = v;
					}
					if (p != u - 1) {
						sqrtsort_swapN(arr1 + (u - 1) * lblock, arr1 + p * lblock, lblock);
						i = tags[u - 1]; tags[u - 1] = tags[p]; tags[p] = i;
					}
				}
				nbl2 = llast = 0;
				if (b == M) llast = lrest % lblock;
				if (llast != 0) {
					while (nbl2 < NBlk && comp(arr1 + NBlk * lblock, arr1 + (NBlk - nbl2 - 1) * lblock) < 0) nbl2++;
				}
				sqrtsort_MergeBuffersLeftWithXBuf(tags, midkey, arr1, NBlk - nbl2, lblock, nbl2, llast, comp);
			}
			for (p = len; --p >= 0;) arr[p] = arr[p - lblock];
		}

		template<typename It, typename Buffer, typename Comp>
		void sqrtsort_commonSort(It arr, int Len, Buffer extbuf, int* Tags, Comp comp) {
			int lblock, cbuf;

			if (Len < 16) {
				sqrtsort_SortIns(arr, Len, comp);
				return;
			}

			lblock = 1;
			while (lblock * lblock < Len) lblock *= 2;
			std::copy_n(arr, lblock, extbuf);

			sqrtsort_commonSort(extbuf, lblock, arr, Tags, comp);

			sqrtsort_BuildBlocks(arr + lblock, Len - lblock, lblock, comp);
			cbuf = lblock;
			while (Len > (cbuf *= 2)) {
				sqrtsort_CombineBlocks(arr + lblock, Len - lblock, cbuf, lblock, Tags, comp);
			}
			sqrtsort_MergeDown(arr + lblock, extbuf, Len - lblock, lblock, comp);
		}

		template<typename It, typename Comp>
		void SqrtSort(It arr, int Len, Comp comp) {
			int L = 1;
			iter_value<It>* ExtBuf;
			int* Tags;

			while (L * L < Len) L *= 2;
			int NK = (Len - 1) / L + 2;
			ExtBuf = new iter_value<It>[L];
			if (ExtBuf == NULL) return; // fail
			Tags = new int[NK];
			if (Tags == NULL) return;

			sqrtsort_commonSort(arr, Len, ExtBuf, Tags, comp);
			delete[]ExtBuf;
			delete[]Tags;
		}
	}

	template<typename Compare>
	struct comperator {
		Compare compare;

		explicit comperator(Compare&& comp) :
			compare(std::forward<Compare>(comp))
		{}

		template<typename T, typename U>
		int operator()(T* lhs, U* rhs)
		{
			if (compare(*lhs, *rhs)) {
				return -1;
			}
			if (compare(*rhs, *lhs)) {
				return 1;
			}
			return 0;
		}
	};


	template<typename RandomAccessIterator, typename Compare>
	void sqrtsort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
		_internal::SqrtSort(first, (int)std::distance(first, last), comperator<Compare>(std::move(comp)));
	}

	template<typename RandomAccessIterator>
	void sqrtsort(RandomAccessIterator first, RandomAccessIterator last) {
		sqrtsort(first, last, std::less<_internal::iter_value<RandomAccessIterator>>());
	}
}

