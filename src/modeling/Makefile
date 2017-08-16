#FTN=/usr/local/bin/mpif90
FTN=mpif90
FLAGS=-c -w
LDFLAGS=-w

OBJS = rdwt.o model.o mute.o forward.o string.o \
       get.o source.o a2d.o segy.o

all: a2d.x 

a2d.x: $(OBJS)
	$(FTN) $(OBJS) $(LDFLAGS) -o $@
	mv $@ ../executables

%.o: %.f
	$(FTN) $< $(FLAGS) -o $@

clean:
	rm -f *.o *.x

