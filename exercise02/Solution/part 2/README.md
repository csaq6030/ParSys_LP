### Part 2

Complete the following steps:
- Assume that the right-hand-side matrix in our matrix multiplication example is always an **upper triangular matrix**. How would you use this knowledge to sequentially optimize the code? Apply your optimizations and document their impact using your benchmarking regimen.
- What impact does this optimization have on the effectiveness of a simple OpenMP parallelization? How would you optimize your parallelization strategy for this case? Implement your ideas and document their performance impact.


For the sequential optimization, we can simply reduce the workload -> Since it is know that the right hand side matrix is an upper triangular matrix, we can stop the calculation of the vector dot product after the calculation is done for the upper half vector (in respect to the i-th iteration).
