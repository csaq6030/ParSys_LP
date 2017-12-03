## Results
for a 2D 5 point stencil with the dimension of 512 * 512 with North/East/South/West boundaries fixed at 1.0, 0.5, 0.0 and -0.5 respectively.

### Mean

|     | 1     | 2     | 4     | 8     | 
|-----|-------|-------|-------|-------| 
| icc | 547.0 | 279.0 | 145.0 | 117.4 | 
| gcc | 942.0 | 500.1 | 279.0 | 199.2 | 

### Standrat deviation

|     | 1   | 2                   | 4   | 8                   | 
|-----|-----|---------------------|-----|---------------------| 
| icc | 0.0 | 0.0                 | 0.0 | 47.379742506687386  | 
| gcc | 0.0 | 0.30000000000000004 | 0.0 | 0.39999999999999997 | 
