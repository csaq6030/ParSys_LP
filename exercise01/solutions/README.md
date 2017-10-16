*Part 1*

Study the usage of the gprof tool. Apply it to the matrix multiplication example and study the report. (Also commit it to your git repo) Try this with both a version compiled using -O0 and one with -O3 and observe the differences.


method		prob size N		cpu time/s		binary size

-O0			10				0				38,4kB
			100				0,1
			500				1,78
			1000			19,04
			2000			181,9

-O3			10				0				14kB	
			100				0
			500				0,23
			1000			2,98
			2000			46,28


Think about issues such as statistical variation and external load, and how to minimize their impact on your results.
- statistical variation: can be encountered by having multiple samples or the sample frequncy can be increased (or longer execution of the problem)
- external load: is not an issue with the profiler gprof since it only measures the cpu time

*Part 2*

Try out at least two methods of performing sequential optimization on this matrix multiplication program. For each of them, apply your benchmarking regimen and document your findings. Also try to explain the results you observe.

Optimization 1:

In the constructor of the identity matrix the inner loop can be removed. The vector can be inicialized instead of being resized. Subsequently the 1s can be written. The value n is dacalred as a const. 

Optimization 2:

In the operator*, the most inner loop an (several n³) unnecessary access to the mat is done (exhaustive indexing). The result can be calculated in a temp var and stored into the right place after its calculation.

|           | avg 1   | avg 12 | avg 2  | 
|-----------|---------|--------|--------| 
| 10        |         |        |        | 
| wall time | 0,00    | 0,00   | 0,00   | 
| user time | 0,00    | 0,00   | 0,00   | 
| sys time  | 0,00    | 0,00   | 0,00   | 
| 100       |         |        |        | 
| wall time | 0,01    | 0,01   | 0,01   | 
| user time | 0,00    | 0,00   | 0,00   | 
| sys time  | 0,01    | 0,01   | 0,01   | 
| 500       |         |        |        | 
| wall time | 2.091   | 1,439  | 1,465  | 
| user time | 2.082   | 1.432  | 1,456  | 
| sys time  | 0,00    | 0,00   | 0,00   | 
| 1000      |         |        |        | 
| wall time | 26,118  | 19,448 | 18,351 | 
| user time | 25,995  | 19,303 | 18,304 | 
| sys time  | 0,053   | 0,03   | 0,025  | 
| 1500      |         |        |        | 
| wall time | 105,462 | 65,522 | 69,611 | 
| user time | 105,252 | 65,442 | 69,52  | 
| sys time  | 0,112   | 0,049  | 0,05   | 

*Part 3*

Study OpenMP parallel and for and attempt to parallelize the matrix multiplication program using OpenMP. Once again, apply your benchmarking regimen and document your findings.

omp parallel: 		creates a team of N threads
omp for:			splits the work among available threads
omp parallel for: 	combines both


|           | avg OMP | 
|-----------|---------| 
| 10        |         | 
| wall time | 0,00    | 
| user time | 0,00    | 
| sys time  | 0,00    | 
| 100       |         | 
| wall time | 0,001   | 
| user time | 0,023   | 
| sys time  | 0,00    | 
| 500       |         | 
| wall time | 0,552   | 
| user time | 4,211   | 
| sys time  | 0,001   | 
| 1000      |         | 
| wall time | 5,476   | 
| user time | 42,313  | 
| sys time  | 0,048   | 
| 1500      |         | 
| wall time | 21,662  | 
| user time | 166,982 | 
| sys time  | 0,219   | 




method		prob size N		cpu time/s		binary size

omp			10				0				44,7kB
			100				0
			500				1,32
			1000			12,24
			2000			85,59
