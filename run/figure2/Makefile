trt:
	cd trt; make

#######################################################
# run:
#        cd ../../../src/trt; make clean; make
#        mpirun -np 8 -machinefile ${PBS_NODEFILE} ../../../src/executables/trt.x tomo.deck
#        #mpicc ../../../src/executables/trt.x tomo.deck
# figure:
#        matlab -nodisplay -nosplash < figure2c.m

# clean:
#        rm -f vel* log* ray* rms* *.png
#######################################################


modeling:
	mpirun -np 40 -machinefile ${PBS_NODEFILE} ../../../src/executables/a2d.x

ewt:
	cd ewt; make

#######################################################
# run:
#        cp wtw_model.para ../../../src/inversion
#        cd ../../../src/inversion; make clean; make
#        mpirun -np 40 -machinefile ${PBS_NODEFILE} ../../../src/executables/ewt.x
########################################################


mwt:
	# 5-Hz inversion
	mpirun -np 40 -machinefile ${PBS_NODEFILE} ../../../src/executables/ewt.x
	# 10-Hz inversion
	mpirun -np 40 -machinefile ${PBS_NODEFILE} ../../../src/executables/ewt.x

figure:
	matlab -nodisplay -nosplash < figure2.m

.PHONY: trt modeling ewt mwt

