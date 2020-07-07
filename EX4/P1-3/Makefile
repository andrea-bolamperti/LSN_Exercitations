CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o block.o clean.sh
	./clean.sh
	$(CC) MolDyn_NVE.o block.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
func.o : block.cpp block.h
	$(CC) -c block.cpp -o block.o $(CFLAGS)
clean :
	rm *.o MolDyn_NVE.exe seed.out
