
default: eznova

timestep.o: timestep.c structs.h
	gcc -c timestep.c

report.o: report.c structs.h
	gcc -c report.c

readfile.o: readfile.c structs.h
	gcc -c readfile.c

eznova: main.c readfile.o report.o timestep.o structs.h
	gcc -o eznova main.c readfile.o report.o timestep.o

clean:
	rm -f *.o eznova
