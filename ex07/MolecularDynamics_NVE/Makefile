CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o
	$(CC) MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm -rf *.o MolDyn_NVE.exe;\
	./clean.sh
run : MolDyn_NVE.exe
	./MolDyn_NVE.exe

restart : MolDyn_NVE.exe
	./restart.sh;\
	./MolDyn_NVE.exe

reset :
	make clean
	cp config.0_bak config.0;
	rm -rf old.0;