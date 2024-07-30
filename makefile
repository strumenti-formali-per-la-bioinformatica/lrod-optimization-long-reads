CC=g++

CPPFLAGS = -g -Wall -O3 -std=c++11

LROD: main.o read.o kmer.o aligning.o bitarray.o
	$(CC) -o $@ $^ -lpthread -L/usr/local/lib -lnthash
	
main.o: main.cpp read.h kmer.h aligning.h bitarray.h
	$(CC) $(CPPFLAGS) -I/usr/local/include -c main.cpp

bitarray.o: bitarray.cpp 
	$(CC) $(CPPFLAGS) -c bitarray.cpp 

read.o: read.cpp read.h
	$(CC) $(CPPFLAGS) -c read.cpp
	
sortContigSet.o: kmer.cpp read.h bitarray.h
	$(CC) $(CPPFLAGS) -c kmer.cpp
	
aligning.o: aligning.cpp kmer.h read.h bitarray.h
	$(CC) $(CPPFLAGS) -c aligning.cpp

all: LROD
	rm -f *.o

clean:
	rm -f *.o
	rm LROD
