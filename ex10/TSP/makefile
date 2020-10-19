CC = mpic++
CFLAGS = -Wall -O3 --std=c++11
LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

T_S_P.exe : main.o random.o city.o genetic.o
	$(CC) random.o main.o city.o genetic.o -o T_S_P.exe ${LIBS}
main.o : main.cpp random.h city.h genetic.h
	$(CC) -c main.cpp -o main.o $(CFLAGS) ${INCS}
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
city.o : city.cpp city.h
	$(CC) -c city.cpp -o city.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h city.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
clean :
	rm -rf *.o T_S_P.exe