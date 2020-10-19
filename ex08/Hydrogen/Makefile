CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o metropolis.o vmc.o
	$(CC) vmc.o random.o metropolis.o main.o -o main.exe
main.o : main.cpp metropolis.h random.h vmc.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h vmc.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
vmc.o : vmc.cpp vmc.h
	$(CC) -c vmc.cpp -o vmc.o $(CFLAGS)
clean :
	rm -rf *.o main.exe seed.out
