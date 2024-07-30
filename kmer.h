#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED
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
#include <cmath>
#include "read.h"
#include "bitarray.h"
using namespace std;

typedef struct KmerHashNode{
    unsigned int kmer;
    int startPositionInArray;
} KmerHashNode;

typedef struct KmerHashTableHead{
    KmerHashNode* kmerHashNode;
    unsigned long int allocationCount;
    long int min;
    long int max;
} KmerHashTableHead;

typedef struct KmerReadNode{
    unsigned int kmer;
    unsigned int readIndex;
    unsigned int position;
    bool orientation;
} KmerReadNode;

typedef struct KmerReadNodeHead{
    KmerReadNode* kmerReadNode;
    long int realCount;
    long int allocationCount;
    long int kmerLength;
    long int startReadIndex;
    long int endReadIndex;
} KmerReadNodeHead;

std::string strup(const std::string& str);
unsigned int hash32shift(unsigned int key);
unsigned int Hash(unsigned int kmer, unsigned int max);
void ReverseComplementKmer(char* kmer, long int kmerLength);
long int SearchKmerHashTable(KmerHashTableHead* kmerHashTableHead, unsigned int kmer);
void sort(KmerReadNode* a, long int left, long int right);
bool DetectSameKmer(const char* kmerf, long int kmerLength);
std::vector<std::string> ReadKmerFile(const char* address, long int maxSize);
KmerReadNodeHead* InitKmerReadNodeHead(const char* address, ReadSetHead* readSetHead, long int kmerLength, long int step, KmerHashTableHead* kmerHashTableHead, int frequencyCutOff);
KmerHashTableHead* GetKmerHashTableHead(const char* address, ReadSetHead* readSetHead, long int kmerLength, long int step, long int min, float maxRatio);
int GetKmerHashTableHead_UnitTest(KmerHashTableHead* kmerHashTableHead);
KmerReadNodeHead* GetKmerReadNodeHeadSub(ReadSetHead* readSetHead, long int kmerLength, long int step, long int intervalCount);
int GetKmerReadNodeHeadSub_UnitTest(KmerReadNodeHead* kmerReadNodeHead);
void InitKmerReadNodeHeadSub(ReadSetHead* readSetHead, KmerReadNodeHead* kmerReadNodeHead, KmerHashTableHead* kmerHashTableHead, long int kmerLength, long int step, long int startReadIndex, long int endReadIndex);

#endif
