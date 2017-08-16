#CC=/uufs/landscapearch.arches/sys/pkg/mpich-mx/std/bin/mpicc
#CC=/usr/local/bin/mpicc
#CC=/opt/aci/sw/openmpi/1.10.1_gcc-5.3.1/bin/mpicc
CC=mpicc
trt.x:
	$(CC) tomo2d_mpi_new.c -o $@ -lm
	mv $@ ../executables

clean:
	rm -f *.o core *.x

