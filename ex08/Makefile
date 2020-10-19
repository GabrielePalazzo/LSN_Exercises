CC = g++
CFLAGS = -Wall -O3 --std=c++11
FOLDER = Hydrogen/

main.exe : $(FOLDER)main.o $(FOLDER)random.o $(FOLDER)metropolis.o $(FOLDER)vmc.o
	$(CC) $(FOLDER)vmc.o $(FOLDER)random.o $(FOLDER)metropolis.o $(FOLDER)main.o -o $(FOLDER)main.exe
$(FOLDER)main.o : $(FOLDER)main.cpp $(FOLDER)metropolis.h $(FOLDER)random.h $(FOLDER)vmc.h
	$(CC) -c $(FOLDER)main.cpp -o $(FOLDER)main.o $(CFLAGS)
$(FOLDER)random.o : $(FOLDER)random.cpp $(FOLDER)random.h
	$(CC) -c $(FOLDER)random.cpp -o $(FOLDER)random.o $(CFLAGS)
$(FOLDER)metropolis.o : $(FOLDER)metropolis.cpp $(FOLDER)metropolis.h $(FOLDER)vmc.h
	$(CC) -c $(FOLDER)metropolis.cpp -o $(FOLDER)metropolis.o $(CFLAGS)
$(FOLDER)vmc.o : $(FOLDER)vmc.cpp $(FOLDER)vmc.h
	$(CC) -c $(FOLDER)vmc.cpp -o $(FOLDER)vmc.o $(CFLAGS)
clean :
	rm -rf $(FOLDER)*.o $(FOLDER)main.exe $(FOLDER)seed.out
run: main.exe
	cd $(FOLDER); \
	./main.exe