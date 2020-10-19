CC = g++
CFLAGS = -Wall -O3 --std=c++11
FOLDER = MolecularDynamics_NVE/

MolDyn_NVE.exe : $(FOLDER)MolDyn_NVE.o
	$(CC) $(FOLDER)MolDyn_NVE.o -o $(FOLDER)MolDyn_NVE.exe
$(FOLDER)MolDyn_NVE.o : $(FOLDER)MolDyn_NVE.cpp $(FOLDER)MolDyn_NVE.h
	$(CC) -c $(FOLDER)MolDyn_NVE.cpp -o $(FOLDER)MolDyn_NVE.o $(CFLAGS)
clean :
	rm -rf $(FOLDER)*.o $(FOLDER)MolDyn_NVE.exe;\
	cd $(FOLDER) && ./clean.sh

run : MolDyn_NVE.exe
	cd $(FOLDER);\
	./MolDyn_NVE.exe

restart : MolDyn_NVE.exe
	cd $(FOLDER);\
	./restart.sh;\
	./MolDyn_NVE.exe

reset :
	make clean
	rm -rf $(FOLDER)old.0;
resset :
	make clean
	cp $(FOLDER)config.0_bak $(FOLDER)config.0;
	rm -rf $(FOLDER)old.0;