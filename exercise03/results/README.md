### Auto vector

|                    | gcc 7.2 | clang 5.0 | icc 17 | MSVC 19 | 
|--------------------|---------|-----------|--------|---------| 
| c_static double    | 0       | 4         | 4      | 4       | 
| c_static int       | 0       | 8         | 8      | 8       | 
| c_static float     | 0       | 8         | 8      | 8       | 
| c_dynamic double   | 0       | 4         | 4      | 0       | 
| c_dynamic int      | 0       | 8         | 8      | 0       | 
| c_dynamic float    | 0       | 8         | 8      | 0       | 
| cpp_dynamic double | 0       | 0         | 4      | 0       | 
| cpp_dynamic int    | 0       | 0         | 0      | 0       | 
| cpp_dynamic float  | 0       | 0         | 8      | 0       | 


### Loop unrolling

|                    | gcc 7.2 | clang 5.0 | icc 17 | MSVC 19 | 
|--------------------|---------|-----------|--------|---------| 
| c_static double    | no      | yes       | yes    | yes     | 
| c_static int       | no      | yes       | yes    | yes     | 
| c_static float     | no      | yes       | yes    | yes     | 
| c_dynamic double   | no      | yes       | yes    | yes     | 
| c_dynamic int      | no      | yes       | yes    | no      | 
| c_dynamic float    | no      | yes       | yes    | yes     | 
| cpp_dynamic double | no      | yes       | yes    | yes     | 
| cpp_dynamic int    | no      | yes       | yes    | no      | 
| cpp_dynamic float  | no      | yes       | yes    | yes     | 

