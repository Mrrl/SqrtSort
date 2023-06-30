SqrtSort
========

Stable sorting with O(sqrt(N)) external memory

Idea is the same as in GrailSort ( https://github.com/Mrrl/GrailSort ), but exchange buffer and buffer for block tags are allocated in the heap. Resulting algorithm is 50% faster than std::sort (when number of keys in array more than sqrt(N)) or std::stable_sort.  

Usage
========

C Language:
```c
#define SORT_TYPE int
int SORT_CMP(int *a, int *b){
    return *a - *b;
}
#include "SqrtSort.h"
#include <stdio.h>
int main(){
    int a[]={4,6,43,64,6546,22,7,46,8};
    int sz=sizeof(a)/sizeof(int);
    SqrtSort(a, sz);
    for(int i=0;i<sz;i++) printf("%d ", a[i]);
    return 0;
}
```

C++ Language:
```cpp
#include <iostream>
#include <vector>
#include "sqrtsort.hpp"
int main(){
    std::vector<int> a{4,6,43,64,6546,22,7,46,8};
    sqrtsort::sqrtsort(a.begin(), a.end());
    for(auto p:a) std::cout<<p<<" ";
```
