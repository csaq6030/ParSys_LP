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

method		prob size N		cpu time/s		binary size

opt1		10				0				44,5kB
			100				0,02
			500				1,85
			1000			19,54
			2000			174,66

opt2		10				0				38,4kB	
			100				0
			500				1,33
			1000			16,96
			2000			138

opt1+2		10				0				44,5kB	
			100				0
			500				1,6
			1000			15,06
			2000			135,29

*Part 3*

Study OpenMP parallel and for and attempt to parallelize the matrix multiplication program using OpenMP. Once again, apply your benchmarking regimen and document your findings.

omp parallel: 		creates a team of N threads
omp for:			splits the work among available threads
omp parallel for: 	combines both


method		prob size N		cpu time/s		binary size

omp			10				0				44,7kB
			100				0
			500				1,32
			1000			12,24
			2000			85,59
