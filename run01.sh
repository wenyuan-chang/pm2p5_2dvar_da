gfortran step01.2dvar.f90 -o s01.out -I/Users/wenyuan/bin/netcdf/include -L/Users/wenyuan/bin/netcdf/lib -lnetcdff -L$HOME/lib/lapack-3.8.0 -llapack -lrefblas
./s01.out
rm s01.out
