all: aout

aout: mod_afssh.o AFSSH.o
	gfortran -o aout mod_afssh.o AFSSH.o -O2 -g ~/lapack-3.8.0/liblapack.a ~/lapack-3.8.0/librefblas.a

%.o: %.f90
	gfortran -c $<

clean:
	rm *.o aout

