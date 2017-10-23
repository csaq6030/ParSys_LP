### Part 2

Complete the following steps:
- Assume that the right-hand-side matrix in our matrix multiplication example is always an **upper triangular matrix**. How would you use this knowledge to sequentially optimize the code? Apply your optimizations and document their impact using your benchmarking regimen.
- What impact does this optimization have on the effectiveness of a simple OpenMP parallelization? How would you optimize your parallelization strategy for this case? Implement your ideas and document their performance impact.


For the sequential optimization, we can simply reduce the workload -> Since it is know that the right hand side matrix is an upper triangular matrix, we can stop the calculation of the vector dot product after the calculation is done for the upper half vector (in respect to the i-th iteration).

#### Mean

|             | wall time [s]       |  user mode cpu [s]  |  kernel mode [s] |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|---------------------|---------------------|------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulSeq10   | 0.0                 | 0.0                 | 0.0              | 0.0               | 120.7             | 3125.2                          | 0.0                            | 
| mmulSeq100  | 0.0                 | 0.0                 | 0.0              | 0.0               | 181.6             | 3361.6                          | 0.0                            | 
| mmulSeq500  | 0.06200000000000001 | 0.06100000000000001 | 0.0              | 0.0               | 1599.3            | 9032.0                          | 0.0                            | 
| mmulSeq1000 | 0.5330000000000001  | 0.527               | 0.0              | 0.0               | 6009.5            | 26670.4                         | 0.0                            | 
| mmulSeq1500 | 1.7570000000000001  | 1.741               | 0.01             | 0.0               | 13347.7           | 56022.8                         | 0.0                            | 

#### SD

|             | wall time [s]        |  user mode cpu [s]    |  kernel mode [s]     |  major page fault |  minor page fault  |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|----------------------|-----------------------|----------------------|-------------------|--------------------|---------------------------------|--------------------------------| 
| mmulSeq10   | 0.0                  | 0.0                   | 0.0                  | 0.0               | 1.1874342087037917 | 43.4161260363013                | 0.0                            | 
| mmulSeq100  | 0.0                  | 0.0                   | 0.0                  | 0.0               | 1.2                | 53.52233178776875               | 0.0                            | 
| mmulSeq500  | 0.005999999999999999 | 0.0030000000000000022 | 0.0                  | 0.0               | 1.268857754044952  | 32.93630216038224               | 0.0                            | 
| mmulSeq1000 | 0.006403124237432854 | 0.009000000000000008  | 0.0                  | 0.0               | 1.857417562100671  | 55.48909802835148               | 0.0                            | 
| mmulSeq1500 | 0.01100000000000001  | 0.009433981132056611  | 0.006324555320336759 | 0.0               | 1.4866068747318506 | 48.79508171937003               | 0.0                            | 

### Comparison to -O3 Seq version

#### Mean

|               | wall time [s]       |  user mode cpu [s]  |  kernel mode [s]     |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|---------------|---------------------|---------------------|----------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulSeqO310   | 0.0                 | 0.0                 | 0.0                  | 0.0               | 121.3             | 3138.0                          | 0.0                            | 
| mmulSeqO3100  | 0.0                 | 0.0                 | 0.0                  | 0.0               | 182.1             | 3376.0                          | 0.0                            | 
| mmulSeqO3500  | 0.06100000000000001 | 0.06000000000000001 | 0.001                | 0.0               | 1600.3            | 9052.8                          | 0.0                            | 
| mmulSeqO31000 | 0.5740000000000001  | 0.567               | 0.002                | 0.0               | 6008.6            | 26664.4                         | 0.0                            | 
| mmulSeqO31500 | 1.9509999999999998  | 1.934               | 0.009999999999999998 | 0.0               | 13347.2           | 56015.6                         | 0.0                            | 

#### SD

|               | wall time [s]         |  user mode cpu [s]     |  kernel mode [s]      |  major page fault |  minor page fault  |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|---------------|-----------------------|------------------------|-----------------------|-------------------|--------------------|---------------------------------|--------------------------------| 
| mmulSeqO310   | 0.0                   | 0.0                    | 0.0                   | 0.0               | 1.6155494421403513 | 44.15427499121687               | 0.0                            | 
| mmulSeqO3100  | 0.0                   | 0.0                    | 0.0                   | 0.0               | 1.22065556157337   | 49.63869458396343               | 0.0                            | 
| mmulSeqO3500  | 0.0030000000000000022 | 1.3877787807814457E-17 | 0.0029999999999999996 | 0.0               | 1.7916472867168916 | 44.928387462716714              | 0.0                            | 
| mmulSeqO31000 | 0.004898979485566361  | 0.006403124237432806   | 0.004                 | 0.0               | 1.7435595774162693 | 58.88157606586291               | 0.0                            | 
| mmulSeqO31500 | 0.005385164807134508  | 0.006633249580710805   | 0.00447213595499958   | 0.0               | 1.16619037896906   | 50.03838526571376               | 0.0                            | 
  

&nbsp;

For the parallel version, one has to keep in mind that each iteration has different amount workload. The work rises lineraly with ith columb of the matrix, therefore the scheduling of the workload to the threads has a huge imapct on the performance of the program. With static scheduling the loop is divided into equal sized chunks and assigned into the available number of threads. Thus, the workload is distribuited unevenly -> some threads might be idle while others still have a lot of computing to do. In this case, a better scheduling approach is dynamic scheduling -> a chunck sized block(default 1) of loop iteration is assigned to the threads, when a thread is finished it recieves a next block of loop iterations.

## Static Scheduling

#### Mean


|                   | wall time [s]       |  user mode cpu [s]   |  kernel mode [s]     |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------------|---------------------|----------------------|----------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulParStatic10   | 0.001               | 0.007000000000000001 | 0.0                  | 0.0 				| 147.4   			| 3437.2  						  | 0.0 						   | 
| mmulParStatic100  | 0.001               | 0.003                | 0.0                  | 0.0				| 209.1   			| 3693.6						  | 0.0 						   | 
| mmulParStatic500  | 0.09599999999999999 | 0.266                | 0.0                  | 0.0 				| 1626.7  			| 9332.4						  | 0.0 						   | 
| mmulParStatic1000 | 0.7729999999999999  | 1.816                | 0.003                | 0.0 				| 6035.2  			| 26975.2 						  | 0.0 						   |
| mmulParStatic1500 | 2.593               | 6.008                | 0.009999999999999998 | 0.0 				| 13374.8 			| 56355.6 						  | 0.0 						   |


#### SD
|             		| wall time [s]         |  user mode cpu [s]    |  kernel mode [s]      |  major page fault |  minor page fault  |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------------|-----------------------|-----------------------|-----------------------|-------------------|--------------------|---------------------------------|--------------------------------| 
| mmulParStatic10   | 0.0029999999999999996 | 0.009                 | 0.0                   | 0.0 				| 0.9165151389911681 | 42.4471436023674  			   | 0.0 							| 
| mmulParStatic100  | 0.0029999999999999996 | 0.0045825756949558405 | 0.0                   | 0.0 				| 0.9433981132056604 | 55.51792503327191 			   | 0.0 							| 
| mmulParStatic500  | 0.006633249580710802  | 0.010198039027185565  | 0.0                   | 0.0 				| 1.0049875621120892 | 48.181324182716274 			   | 0.0 							| 
| mmulParStatic1000 | 0.0100498756211209    | 0.015620499351813323  | 0.0045825756949558405 | 0.0 				| 0.9797958971132713 | 45.2831094338717   		       | 0.0 							| 
| mmulParStatic1500 | 0.019519221295943117  | 0.026758176320519078  | 0.006324555320336759  | 0.0 				| 1.16619037896906   | 57.89507751095943 			   | 0.0 							| 


### Dynamic Scheduling

#### Mean

|             | wall time [s]       |  user mode cpu [s]   |  kernel mode [s]     |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|---------------------|----------------------|----------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulPar10   | 0.0                 | 0.0                  | 0.0                  | 0.0               | 147.7             | 3421.2                          | 0.0                            | 
| mmulPar100  | 0.003               | 0.007000000000000001 | 0.0                  | 0.0               | 208.5             | 3668.4                          | 0.0                            | 
| mmulPar500  | 0.05700000000000001 | 0.213                | 0.0                  | 0.0               | 1626.0            | 9351.6                          | 0.0                            | 
| mmulPar1000 | 0.466               | 1.7800000000000005   | 0.001                | 0.0               | 6036.5            | 27008.0                         | 0.0                            | 
| mmulPar1500 | 1.5610000000000002  | 6.015000000000001    | 0.012999999999999998 | 0.0               | 13374.7           | 56355.6                         | 0.0                            | 

#### SD

|             | wall time [s]        |  user mode cpu [s]   |  kernel mode [s]      |  major page fault |  minor page fault  |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|-------------|----------------------|----------------------|-----------------------|-------------------|--------------------|---------------------------------|--------------------------------| 
| mmulPar10   | 0.0                  | 0.0                  | 0.0                   | 0.0               | 1.2688577540449522 | 42.59765251748034               | 0.0                            | 
| mmulPar100  | 0.009                | 0.02100000000000001  | 0.0                   | 0.0               | 1.02469507659596   | 47.343848597257065              | 0.0                            | 
| mmulPar500  | 0.004582575694955838 | 0.006403124237432849 | 0.0                   | 0.0               | 1.3416407864998738 | 61.40879415849166               | 0.0                            | 
| mmulPar1000 | 0.010198039027185555 | 0.013416407864998751 | 0.0029999999999999996 | 0.0               | 0.9219544457292888 | 50.84879546262625               | 0.0                            | 
| mmulPar1500 | 0.019723082923316038 | 0.01910497317454267  | 0.006403124237432849  | 0.0               | 0.7810249675906655 | 48.74258918030514               | 0.0                            | 

