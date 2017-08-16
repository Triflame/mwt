#FTN=/uufs/landscapearch.arches/sys/pkg/mpich-mx/std/bin/mpif90
FTN=mpif90
FLAGS=-c -w
LDFLAGS=-w

OBJS = rdwt.o fft.o model.o mute.o forward.o residual.o string.o \
       get.o wt.o source.o segy.o

all: ewt.x

ewt.x: $(OBJS)
	$(FTN) $(OBJS) $(LDFLAGS) -o $@
	mv $@ ../executables

%.o: %.f
	$(FTN) $< $(FLAGS) -o $@

clean:
	rm -f *.o *.x

