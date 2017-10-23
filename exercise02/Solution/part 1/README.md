
### Part 1

- The three versions use a different implementation of the Matrix `struc`. The `NESTED_VECTOR` is a simple implementation using a `vector<vector<double>>`. `CONTIGUOUS_WITH_MULTIPLICATION` allocates `n*n` as a block and uses address calculation as a way to access it. `CONTIGUOUS_WITH_INDIRECTION` also allocates a block, but uses a `vector<double*>` to access it with already calculated pointers. This acquires additional memory. The `NESTED_VECTOR` should produce the most overhead and therefore be the slowest implementation. `CONTIGUOUS_WITH_MULTIPLICATION` should be faster. `CONTIGUOUS_WITH_INDIRECTION` uses extra memory for the `vector<double*>`, but should have the best performance because of the faster access. 
- We choose the following counters using the gnu time command:
wall time [s], user mode cpu [s], kernel mode [s], major page fault, minor page fault, max resident set size [kbytes] and avg resident set size [kbytes]

- Our suspicion from analyzing the code got confirmed in our measurements. In addition to this the overhead from the `NESTED_VECTOR` also resulted in a bigger memory usage then the `CONTIGUOUS_WITH_INDIRECTION` implementation.

##### Average

|                | wall time [s]        |  user mode cpu [s]   |  kernel mode [s]     |  major page fault |  minor page fault |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|----------------|----------------------|----------------------|----------------------|-------------------|-------------------|---------------------------------|--------------------------------| 
| mmulInd1000    | 15.377               | 15.303999999999998   | 0.032                | 0.0               | 6214.7            | 9.9704832E7                     | 0.0                            | 
| mmulInd100     | 0.009999999999999998 | 0.009999999999999998 | 0.0                  | 0.0               | 403.8             | 4502323.2                       | 0.0                            | 
| mmulInd10      | 0.0                  | 0.0                  | 0.0                  | 0.0               | 344.0             | 3522560.0                       | 0.0                            | 
| mmulInd1500    | 50.718               | 50.577000000000005   | 0.074                | 0.0               | 13540.9           | 2.197356544E8                   | 0.0                            | 
| mmulInd500     | 1.891                | 1.879                | 0.0                  | 0.0               | 1816.9            | 2.76496384E7                    | 0.0                            | 
| mmulMulti1000  | 15.920000000000002   | 15.853000000000003   | 0.032                | 0.0               | 6205.2            | 9.95540992E7                    | 0.0                            | 
| mmulMulti100   | 0.012                | 0.009999999999999998 | 0.0                  | 0.0               | 401.5             | 4463001.6                       | 0.0                            | 
| mmulMulti10    | 0.0                  | 0.0                  | 0.0                  | 0.0               | 342.0             | 3489792.0                       | 0.0                            | 
| mmulMulti1500  | 53.835               | 53.614               | 0.10800000000000001  | 0.0               | 13527.8           | 2.195259392E8                   | 0.0                            | 
| mmulMulti500   | 1.9459999999999997   | 1.934                | 0.001                | 0.0               | 1810.2            | 2.75447808E7                    | 0.0                            | 
| mmulNested1000 | 17.559               | 17.5                 | 0.028000000000000004 | 0.0               | 6391.4            | 1.026031616E8                   | 0.0                            | 
| mmulNested100  | 0.015                | 0.013000000000000001 | 0.0                  | 0.0               | 404.7             | 4517068.8                       | 0.0                            | 
| mmulNested10   | 0.0                  | 0.0                  | 0.0                  | 0.0               | 343.9             | 3522560.0                       | 0.0                            | 
| mmulNested1500 | 59.775               | 59.596000000000004   | 0.091                | 0.0               | 13938.6           | 2.262564864E8                   | 0.0                            | 
| mmulNested500  | 2.168                | 2.1590000000000003   | 0.0                  | 0.0               | 1864.1            | 2.84164096E7                    | 0.0                            | 


##### Standart deviation

|                | wall time [s]         |  user mode cpu [s]    |  kernel mode [s]      |  major page fault |  minor page fault  |  max resident set size [kbytes] | avg resident set size [kbytes] | 
|----------------|-----------------------|-----------------------|-----------------------|-------------------|--------------------|---------------------------------|--------------------------------| 
| mmulInd1000    | 1.2457612130741584    | 1.2286024580799113    | 0.009797958971132713  | 0.0               | 0.6403124237432849 | 8192.0                          | 0.0                            | 
| mmulInd100     | 1.734723475976807E-18 | 1.734723475976807E-18 | 0.0                   | 0.0               | 0.9797958971132712 | 16052.975978303835              | 0.0                            | 
| mmulInd10      | 0.0                   | 0.0                   | 0.0                   | 0.0               | 0.0                | 0.0                             | 0.0                            | 
| mmulInd1500    | 4.039237056672956     | 4.031664792613592     | 0.026907248094147424  | 0.0               | 0.8306623862918073 | 8026.487989151918               | 0.0                            | 
| mmulInd500     | 0.1470680114776833    | 0.1472718574609555    | 0.0                   | 0.0               | 1.7578395831246942 | 15016.184037231298              | 0.0                            | 
| mmulMulti1000  | 0.8062877898120496    | 0.8072180622359734    | 0.006000000000000001  | 0.0               | 1.661324772583615  | 26469.162060027516              | 0.0                            | 
| mmulMulti100   | 0.004                 | 1.734723475976807E-18 | 0.0                   | 0.0               | 1.02469507659596   | 15016.184037231298              | 0.0                            | 
| mmulMulti10    | 0.0                   | 0.0                   | 0.0                   | 0.0               | 0.0                | 0.0                             | 0.0                            | 
| mmulMulti1500  | 2.5399931102268773    | 2.4958733942249554    | 0.05509990925582364   | 0.0               | 0.4                | 6553.6                          | 0.0                            | 
| mmulMulti500   | 0.14793241700181872   | 0.14298251641372103   | 0.0029999999999999996 | 0.0               | 1.9899748742132404 | 32603.74833910973               | 0.0                            | 
| mmulNested1000 | 1.7627447347815275    | 1.7303121105742743    | 0.02039607805437114   | 0.0               | 5.0039984012787215 | 81985.50980655057               | 0.0                            | 
| mmulNested100  | 0.005                 | 0.0045825756949558405 | 0.0                   | 0.0               | 1.4177446878757827 | 23228.328966156823              | 0.0                            | 
| mmulNested10   | 0.0                   | 0.0                   | 0.0                   | 0.0               | 0.3                | 0.0                             | 0.0                            | 
| mmulNested1500 | 5.97190296304285      | 5.826162030015985     | 0.07175653280364096   | 0.0               | 5.0039984012787215 | 81985.50980655057               | 0.0                            | 
| mmulNested500  | 0.22067170185594714   | 0.21741435095227735   | 0.0                   | 0.0               | 4.437341546466759  | 78301.12269335608               | 0.0                            | 

