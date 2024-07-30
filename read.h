#ifndef Read_H_INCLUDED 
#define Read_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#ifdef _WIN32
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include <string.h>
#include <math.h>

using namespace std;
//nodi della lista
typedef struct ReadSet{
	char * read;
	long int readLength;
//	long int allReadLen;
	
}ReadSet;
//lista di readSet
typedef struct ReadSetHead{
	ReadSet * readSet; //puntatore alla prossima read
	long int readCount; //numero di read
}ReadSetHead;


typedef struct ReadToKmerIndex{
	long int * index;
	long int indexCount;
}ReadToKmerIndex;

typedef struct ReadToKmerIndexSet{
	ReadToKmerIndex * readToKmerIndex;
	long int readCount;
}ReadToKmerIndexSet;


ReadSetHead * GetReadSetHead(char *filename,char *StrLine, long int maxSize);

long int max (long int x, long int y);

long int min (long int x, long int y);

#endif