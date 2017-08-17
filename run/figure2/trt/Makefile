run:
	cd ../../../src/trt; make clean; make
	#mpirun -np 8 -machinefile ${PBS_NODEFILE} ../../../src/executables/trt.x tomo.deck
	mpirun -np 2 ../../../src/executables/trt.x tomo.deck
	#mpicc ../../../src/executables/trt.x tomo.deck
figure:
	matlab -nodisplay -nosplash < figure2c.m

clean:
	rm -f vel* log* ray* rms* *.png


