
default: eznova

decay.o: decay.c
	gcc -O3 -c decay.c

fld.o: fld.c
	gcc -O3 -c fld.c

eos.o: eos.c
	gcc -O3 -c eos.c

timestep.o: timestep.c structs.h
	gcc -O3 -c timestep.c

report.o: report.c structs.h
	gcc -O3 -c report.c

readfile.o: readfile.c structs.h
	gcc -O3 -c readfile.c

eznova: main.c readfile.o report.o timestep.o eos.o fld.o decay.o structs.h
	gcc -O3 -o eznova main.c readfile.o report.o timestep.o eos.o fld.o decay.o

clean:
	rm -f *.o eznova
