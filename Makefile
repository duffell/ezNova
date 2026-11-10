
OBJ = readfile.o report.o timestep.o eos.o fld.o decay.o

default: eznova

%.o: %.c structs.h
	gcc -O3 -c $<

eznova: main.c $(OBJ) structs.h
	gcc -O3 -o eznova main.c $(OBJ)

clean:
	rm -f *.o eznova
