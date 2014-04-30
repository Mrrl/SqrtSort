SqrtSort
========

Stable sorting with O(sqrt(N)) external memory

Idea is the same as in GrailSort ( https://github.com/Mrrl/GrailSort ), but exchange buffer and buffer for block tags are allocated in the heap. Resulting algorithm is 50% faster than std::sort (when number of keys in array more than sqrt(N)) or std::stable_sort.  

Results for 100M array with mostly different keys:

    SqrtSort              20.106 sec
    qsort (from stdlib)   25.547 sec
    std::stable_sort      28.682 sec
    std::sort             37.264 sec
