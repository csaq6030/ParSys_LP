### Part 2

Complete the following steps:
- Assume that the right-hand-side matrix in our matrix multiplication example is always an **upper triangular matrix**. How would you use this knowledge to sequentially optimize the code? Apply your optimizations and document their impact using your benchmarking regimen.
- What impact does this optimization have on the effectiveness of a simple OpenMP parallelization? How would you optimize your parallelization strategy for this case? Implement your ideas and document their performance impact.


For the sequential optimization, we can simply reduce the workload -> Since it is know that the right hand side matrix is an upper triangular matrix, we can stop the calculation of the vector dot product after the calculation is done for the upper half vector (in respect to the i-th iteration).

|             | wall time [s]       |  user mode cpu [s]  |  kernel mode [s] |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|---------------------|---------------------|------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulSeq10   | 0.0                 | 0.0                 | 0.0              | 0.0               | 120.7             | 3125.2                          | 0.0                            | 
| mmulSeq100  | 0.0                 | 0.0                 | 0.0              | 0.0               | 181.6             | 3361.6                          | 0.0                            | 
| mmulSeq500  | 0.06200000000000001 | 0.06100000000000001 | 0.0              | 0.0               | 1599.3            | 9032.0                          | 0.0                            | 
| mmulSeq1000 | 0.5330000000000001  | 0.527               | 0.0              | 0.0               | 6009.5            | 26670.4                         | 0.0                            | 
| mmulSeq1500 | 1.7570000000000001  | 1.741               | 0.01             | 0.0               | 13347.7           | 56022.8                         | 0.0                            | 



For the parallel version, one has to keep in mind that each iteration has different workload. The work rises lineraly with ith columb of the matrix, therefore the scheduling of the workload to the threads has a huge imapct on the performance of the program. With static scheduling the loop is divided into equal sized chunks and assigned into the available number of threads. Thus, the workload is distribuited unevenly -> some threads might be idle while others still have a lot of computing to do. In this case, a better scheduling approach is dynamic scheduling -> a chunck sized block(default 1) of loop iteration is assigned to the threads, when a thread is finished it recieves a next block of loop iterations.

## Static Scheduling

## Dynamic Scheduling
|             | wall time [s]       |  user mode cpu [s]   |  kernel mode [s]     |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|---------------------|----------------------|----------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulPar10   | 0.0                 | 0.0                  | 0.0                  | 0.0               | 147.7             | 3421.2                          | 0.0                            | 
| mmulPar100  | 0.003               | 0.007000000000000001 | 0.0                  | 0.0               | 208.5             | 3668.4                          | 0.0                            | 
| mmulPar500  | 0.05700000000000001 | 0.213                | 0.0                  | 0.0               | 1626.0            | 9351.6                          | 0.0                            | 
| mmulPar1000 | 0.466               | 1.7800000000000005   | 0.001                | 0.0               | 6036.5            | 27008.0                         | 0.0                            | 
| mmulPar1500 | 1.5610000000000002  | 6.015000000000001    | 0.012999999999999998 | 0.0               | 13374.7           | 56355.6                         | 0.0                            | 

