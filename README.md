# Corder-TPDS-21
TPDS'21: Workload Balancing via Graph Reordering on Multicore Systems

Two variants of Corder and contemporary degree-based reordering techniques can be found in this folder.

One variant is implemented for [GPOP](https://github.com/souravpati/GPOP)-like frameworks, which adopts CSR partitioning strategy under the partition-centric paradigm. 

The other variant is implemented for [Ligra](https://github.com/faldupriyank/dbg)-like frameworks that follow the classic vertex-centric paradigm.

The processing speed of GPOP is significantly faster than that of Ligra. However, its data structure (CSR) consumes far more time to construct a new graph after the vertex is reordered.  

Detailed descriptions and usages regarding the two frameworks are offered in respective Github repositories. 

## Dataset
The data should be formatted in CSR, following GPOP's instruction. A toy dataset is provided for play. However, its size is too small to exploit the L2-cache-sized partitioning strategy. Please use large graphs as our paper mentioned. 

## Environment && Dependencies
* Intel Skylake Processor 
    * Intel(R) Xeon(R) Silver 4210 CPU @ 2.20GHz 
    * 1MB L2 Cache
* g++ 9.3.0
* Boost 1.75


## Usage
### GPOP
```make clean && make```

```./corder -d /path/to/your/data -a 7 -s 1024 -o /path/to/store/data```

    ./corder
    -d [ --data ] arg           input data path
    -o [ --output ] arg         output file path
    -s [ --size ] arg (=1024)   partition size
    -a [ --algorithm ] arg (=0) reordering algorithm
    -t [ --thread ] arg (=20)   threads

    [ --algorithm ]
    original = 0,
    random = 1,
    sort = 2,
    fbc = 3,
    hc = 4,
    dbg = 5, 
    corder = 6, (a slow version of Corder, which offers better readability)
    fastCorder = 7 (a highly parallized Corder)


### DBG
Replace the ```dbg.h``` file in the DBG's workspace with the patch code from our ```./ligra-patch``` folder. 

Then, follow DBG's instruction to execute the applications.


## Nofication
Notice that the optimal partition size for GPOP varies with the different processor micro-architectures as well as the graph datasets. For convenience, we directly borrow the heuristics of contemporary papers and fix its partition size to be L2 cache size. The optimal combination of GPOP's partition size and Corder's partition size can be explored via exhaustive experiments. 

Another notification is that, we can disrupt the relative order of hot vertices and cold vertices, and obtain suprisingly high performance gain, e.g., speedup > 4x on highly skewed graph. So far, I have not found any disadvantage in disrupting the order, except for losing the theoretical "goodness" of structural information.