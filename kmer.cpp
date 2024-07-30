#ifndef KMER_CPP_INCLUDED
#define KMER_CPP_INCLUDED

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include "kmer.h"
#include "bitarray.h"

std::string strup(const std::string& str) {
    std::string ret = str;
    std::transform(ret.begin(), ret.end(), ret.begin(), ::toupper);
    return ret;
}

unsigned int hash32shift(unsigned int key) {
    key = ~key + (key << 15);
    key = key ^ (key >> 12);
    key = key + (key << 2);
    key = key ^ (key >> 4);
    key = key * 2057;
    key = key ^ (key >> 16);
    return key;
}

unsigned int Hash(unsigned int kmer, unsigned int max) {
    return (hash32shift(kmer) * kmer) % max;
}

void ReverseComplementKmer(char* kmer, long int kmerLength) {
    std::reverse(kmer, kmer + kmerLength);
    for (int i = 0; i < kmerLength; i++) {
        switch (kmer[i]) {
            case 'A': case 'a': kmer[i] = 'T'; break;
            case 'T': case 't': kmer[i] = 'A'; break;
            case 'G': case 'g': kmer[i] = 'C'; break;
            case 'C': case 'c': kmer[i] = 'G'; break;
            case 'N': case 'n': kmer[i] = 'N'; break;
        }
    }
}

long int SearchKmerHashTable(KmerHashTableHead* kmerHashTableHead, unsigned int kmer) {
    long int hashIndex = Hash(kmer, kmerHashTableHead->allocationCount);
    while (true) {
        if (kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0) {
            return -1;
        }
        if (kmerHashTableHead->kmerHashNode[hashIndex].kmer == kmer + 1) {
            return hashIndex;
        } else {
            hashIndex = (hashIndex + 5) % kmerHashTableHead->allocationCount;
        }
    }
    return -1;
}

void sort(KmerReadNode* a, long int left, long int right) {
    if (left >= right) {
        return;
    }
    long int i = left;
    long int j = right;
    KmerReadNode pivot = a[left];
    while (i < j) {
        while (i < j && pivot.kmer <= a[j].kmer) {
            j--;
        }
        if (i < j) {
            a[i] = a[j];
            i++;
        }
        while (i < j && pivot.kmer >= a[i].kmer) {
            i++;
        }
        if (i < j) {
            a[j] = a[i];
            j--;
        }
    }
    a[i] = pivot;
    sort(a, left, i - 1);
    sort(a, i + 1, right);
}

bool DetectSameKmer(const char* kmerf, long int kmerLength) {
    for (long int i = 1; i < kmerLength; i++) {
        if (kmerf[0] != kmerf[i]) {
            return true;
        }
    }
    return false;
}

KmerReadNodeHead* InitKmerReadNodeHead(const char* address, ReadSetHead* readSetHead, long int kmerLength, long int step, KmerHashTableHead* kmerHashTableHead, int frequencyCutOff) {
    std::ifstream file(address);
    if (!file.is_open()) {
        std::cerr << address << ", does not exist!" << std::endl;
        exit(EXIT_FAILURE);
    }

    long int arrayCount = 1000000;
    std::vector<int> freArray(arrayCount, 0);

    long int kmerCount = 0;
    long int allKmerFrequency = 0;
    std::string line;

    // Prima lettura del file per calcolare frequenze e range
    while (std::getline(file, line)) {
        std::string kmer = line.substr(0, kmerLength);
        if (!DetectSameKmer(kmer.c_str(), kmerLength)) {
            continue;
        }
        int frequency = std::stoi(line.substr(kmerLength + 1));
        if (frequency >= arrayCount - 10) {
            continue;
        }
        freArray[frequency]++;
        kmerCount++;
        allKmerFrequency += frequency;
    }

    std::cout << "The number of kmer types: " << kmerCount << ";" << std::endl;
    std::cout << "Sum of frequencies for different kmer: " << allKmerFrequency << ";" << std::endl;

    float acc = 0;
    long int max = 0;
    for (long int i = 0; i < arrayCount; i++) {
        if (freArray[i] != 0) {
            acc += static_cast<float>(freArray[i] * i) / allKmerFrequency;
            if (acc > 0.9 && max == 0) {
                max = i;
                break;
            }
        }
    }

    long int min = std::max(2L, static_cast<long int>(frequencyCutOff));
    if (min > max) {
        std::cerr << "min is larger than max!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "The range of kmer frequency is: [" << min << "," << max << "];" << std::endl;

    // Reinizializza variabili per la seconda lettura
    file.clear();
    file.seekg(0, std::ios::beg);

    kmerCount = 0;
    allKmerFrequency = 0;

    kmerHashTableHead->allocationCount = kmerCount * 1.2;
    kmerHashTableHead->kmerHashNode = new KmerHashNode[kmerHashTableHead->allocationCount]();

    unsigned long int kmerInteger = 0;
    while (std::getline(file, line)) {
        std::string kmer = line.substr(0, kmerLength);
        if (!DetectSameKmer(kmer.c_str(), kmerLength)) {
            continue;
        }
        int frequency = std::stoi(line.substr(kmerLength + 1));
        if (frequency >= arrayCount - 10 || frequency < min || frequency > max) {
            continue;
        }

        SetBitKmer(&kmerInteger, kmerLength, kmer.c_str());
        long int hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount);
        while (true) {
            if (kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0) {
                kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;
                break;
            } else {
                hashIndex = (hashIndex + 1) % kmerHashTableHead->allocationCount;
            }
        }
    }

    auto kmerReadNodeHead = new KmerReadNodeHead();
    kmerReadNodeHead->realCount = 0;
    kmerReadNodeHead->allocationCount = allKmerFrequency * 1.1;
    kmerReadNodeHead->kmerReadNode = new KmerReadNode[kmerReadNodeHead->allocationCount]();

    long int readLength = 0;
    char* kmer1 = new char[kmerLength + 1];

    for (long int i = 0; i < readSetHead->readCount; i++) {
        readLength = readSetHead->readSet[i].readLength;
        for (int j = 0; j < readLength - kmerLength + 1 - step; j += step) {
            strncpy(kmer1, readSetHead->readSet[i].read + j, kmerLength);
            kmer1[kmerLength] = '\0';
            if (!DetectSameKmer(kmer1, kmerLength)) {
                continue;
            }
            SetBitKmer(&kmerInteger, kmerLength, kmer1);
            long int hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
            if (hashIndex != -1) {
                auto& node = kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount++];
                node.kmer = kmerInteger;
                node.readIndex = i + 1;
                node.position = j;
                node.orientation = true;
            } else {
                ReverseComplementKmer(kmer1, kmerLength);
                SetBitKmer(&kmerInteger, kmerLength, kmer1);
                hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
                if (hashIndex != -1) {
                    auto& node = kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount++];
                    node.kmer = kmerInteger;
                    node.readIndex = i + 1;
                    node.position = j;
                    node.orientation = false;
                }
            }
        }
    }

    sort(kmerReadNodeHead->kmerReadNode, 0, kmerReadNodeHead->realCount - 1);

    kmerInteger = kmerReadNodeHead->kmerReadNode[0].kmer + 1;
    for (long int i = 0; i < kmerReadNodeHead->realCount; i++) {
        if (kmerReadNodeHead->kmerReadNode[i].kmer != kmerInteger) {
            kmerInteger = kmerReadNodeHead->kmerReadNode[i].kmer;
            long int hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
            kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray = i;
        }
    }

    delete[] kmer1;
    kmerReadNodeHead->kmerLength = kmerLength;

    return kmerReadNodeHead;
}

KmerHashTableHead* GetKmerHashTableHead(const char* address, ReadSetHead* readSetHead, long int kmerLength, long int step, long int min, float maxRatio) {
    std::cout << "\nFunzione GetKmerHashTableHead\n";

    std::ifstream file(address);
    if (!file.is_open()) {
        std::cerr << address << ", does not exist!" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto kmerHashTableHead = new KmerHashTableHead();

    long int arrayCount = 1000000;
    std::vector<int> freArray(arrayCount, 0);

    long int kmerCount = 0;
    long int allKmerFrequency = 0;
    std::string line;

    while (std::getline(file, line)) {
        std::string kmer = line.substr(0, kmerLength);
        if (!DetectSameKmer(kmer.c_str(), kmerLength)) {
            continue;
        }
        int frequency = std::stoi(line.substr(kmerLength + 1));
        if (frequency >= arrayCount - 10) {
            continue;
        }
        freArray[frequency]++;
        kmerCount++;
        allKmerFrequency += frequency;
    }

    std::cout << "all kmer frequency: " << allKmerFrequency << "\n";
    std::cout << "The number of kmer types: " << kmerCount << ";\n";

    float acc = 0;
    long int max = 0;

    if (maxRatio > 1) {
        maxRatio = 1;
    }

    for (long int i = 0; i < arrayCount; i++) {
        if (freArray[i] != 0) {
            acc += static_cast<float>(freArray[i] * i) / allKmerFrequency;
            if (acc > maxRatio && max == 0) {
                max = i;
                break;
            }
        }
    }

    std::cout << "il max di acc è: " << max << "\n";

    if (min < 2) {
        min = 2;
    }

    if (min > max) {
        std::cerr << "The parameter minimumKmerFrequency is larger than maxKmerFrequencyRatio. Please increase the value of maxKmerFrequencyRatio or decrease the value of minimumKmerFrequency!" << std::endl;
        exit(EXIT_FAILURE);
    }

    kmerCount = 0;
    allKmerFrequency = 0;

    file.clear();
    file.seekg(0, std::ios::beg);

    while (std::getline(file, line)) {
        std::string kmer = line.substr(0, kmerLength);
        if (!DetectSameKmer(kmer.c_str(), kmerLength)) {
            continue;
        }
        int frequency = std::stoi(line.substr(kmerLength + 1));
        if (frequency >= arrayCount - 10 || frequency < min || frequency > max) {
            continue;
        }
        kmerCount++;
        allKmerFrequency += frequency;
    }

    std::cout << "There are " << kmerCount << " kmer for overlap detection;\n";
    if (kmerCount <= 0) {
        exit(EXIT_FAILURE);
    }

    kmerHashTableHead->allocationCount = kmerCount * 2;
    kmerHashTableHead->kmerHashNode = new KmerHashNode[kmerHashTableHead->allocationCount]();

    unsigned long int kmerInteger = 0;

    file.clear();
    file.seekg(0, std::ios::beg);

    while (std::getline(file, line)) {
        std::string kmer = line.substr(0, kmerLength);
        if (!DetectSameKmer(kmer.c_str(), kmerLength)) {
            continue;
        }
        int frequency = std::stoi(line.substr(kmerLength + 1));
        if (frequency >= arrayCount - 10 || frequency < min || frequency > max) {
            continue;
        }

        SetBitKmer(&kmerInteger, kmerLength, kmer.c_str());

        long int hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount);

        while (true) {
            if (kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0) {
                kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;
                break;
            } else {
                hashIndex = (hashIndex + 5) % kmerHashTableHead->allocationCount;
            }
        }
    }

    std::cout << "\nstampa della kmerHashtable: \n";
    for (unsigned int i = 0; i < kmerHashTableHead->allocationCount; i++) {
        std::cout << "elemento " << i << " = " << kmerHashTableHead->kmerHashNode[i].kmer << "\n";
    }

    kmerHashTableHead->min = min;
    kmerHashTableHead->max = max;

    std::cout << "kmer = " << kmerHashTableHead->kmerHashNode->kmer << " e startPositionInArray = " << kmerHashTableHead->kmerHashNode->startPositionInArray << "\n";

    return kmerHashTableHead;
}

int GetKmerHashTableHead_UnitTest(KmerHashTableHead* kmerHashTableHead) {
    long int kmerCount = 0;
    for (unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++) {
        if (kmerHashTableHead->kmerHashNode[i].kmer < 0) {
            return 0;
        } else if (kmerHashTableHead->kmerHashNode[i].kmer > 0) {
            kmerCount++;
        }
    }

    std::cout << "KMERCOUNT = " << kmerCount << "\n";

    if (kmerCount <= 0) {
        return 0;
    }
    return 1;
}

KmerReadNodeHead* GetKmerReadNodeHeadSub(ReadSetHead* readSetHead, long int kmerLength, long int step, long int intervalCount) {
    std::cout << "\n Funzione GetKmerReadNodeHeadSub\n";

    long int max = 0;
    long int allKmerFrequency = 0;
    long int startReadIndex = 0;
    long int endReadIndex = startReadIndex + intervalCount;

    std::cout << "readSetHead->readCount: " << readSetHead->readCount << "\n";

    while (true) {
        if (startReadIndex >= readSetHead->readCount) {
            break;
        }
        if (endReadIndex >= readSetHead->readCount) {
            endReadIndex = readSetHead->readCount - 1;
        }
        allKmerFrequency = 0;
        for (long int i = startReadIndex; i <= endReadIndex; i++) {
            allKmerFrequency += readSetHead->readSet[i].readLength - kmerLength;
        }
        if (allKmerFrequency > max) {
            max = allKmerFrequency;
        }
        startReadIndex = endReadIndex + 1;
        endReadIndex += intervalCount;
    }

    std::cout << "End allkmerfrequency: " << allKmerFrequency << "\n";
    max = max / step;

    auto kmerReadNodeHead = new KmerReadNodeHead();
    kmerReadNodeHead->realCount = 0;
    kmerReadNodeHead->allocationCount = max * 1.1;
    kmerReadNodeHead->kmerReadNode = new KmerReadNode[kmerReadNodeHead->allocationCount]();

    return kmerReadNodeHead;
}

void InitKmerReadNodeHeadSub(ReadSetHead* readSetHead, KmerReadNodeHead* kmerReadNodeHead, KmerHashTableHead* kmerHashTableHead, long int kmerLength, long int step, long int startReadIndex, long int endReadIndex) {
    std::cout << "funzione: InitKmerReadNodeHeadSub\n";
    std::cout << "\n startIndex = " << startReadIndex << " endReadIndex = " << endReadIndex << "\n ";

    kmerReadNodeHead->realCount = 0;

    for (unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++) {
        kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1;
    }

    for (long int i = 0; i < kmerReadNodeHead->allocationCount; i++) {
        kmerReadNodeHead->kmerReadNode[i].kmer = 0;
    }

    long int readLength = 0;

    char* kmer1 = new char[kmerLength + 1];

    long int hashIndex = 0;
    unsigned long int kmerInteger = 0;

    std::cout << "\n stampa della struttura kmerReadNode\n";

    for (long int i = startReadIndex; i <= endReadIndex; i++) {
        readLength = readSetHead->readSet[i].readLength;

        for (int j = 0; j < readLength - kmerLength + 1 - step; j += step) {
            strncpy(kmer1, readSetHead->readSet[i].read + j, kmerLength);
            kmer1[kmerLength] = '\0';

            SetBitKmer(&kmerInteger, kmerLength, kmer1);

            hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);

            if (hashIndex != -1) {
                std::cout << "\n kmer: " << kmer1 << " kmerInteger: " << kmerInteger << " readIndex: " << i << "\n";

                auto& node = kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount++];
                node.kmer = kmerInteger;
                node.readIndex = i + 1;
                node.position = j;
                node.orientation = true;
            } else {
                ReverseComplementKmer(kmer1, kmerLength);
                SetBitKmer(&kmerInteger, kmerLength, kmer1);
                hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
                if (hashIndex != -1) {
                    std::cout << "\n kmer: " << kmer1 << " kmerInteger: " << kmerInteger << " readIndex: " << i << "\n";

                    auto& node = kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount++];
                    node.kmer = kmerInteger;
                    node.readIndex = i + 1;
                    node.position = j;
                    node.orientation = false;
                }
            }
        }
    }
    sort(kmerReadNodeHead->kmerReadNode, 0, kmerReadNodeHead->realCount - 1);

    kmerInteger = kmerReadNodeHead->kmerReadNode[0].kmer + 1;
    for (long int i = 0; i < kmerReadNodeHead->realCount; i++) {
        if (kmerReadNodeHead->kmerReadNode[i].kmer != kmerInteger) {
            kmerInteger = kmerReadNodeHead->kmerReadNode[i].kmer;
            hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
            kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray = i;
        }
    }
    delete[] kmer1;
    kmerReadNodeHead->kmerLength = kmerLength;
    kmerReadNodeHead->startReadIndex = startReadIndex;
    kmerReadNodeHead->endReadIndex = endReadIndex;
}

int GetKmerReadNodeHeadSub_UnitTest(KmerReadNodeHead* kmerReadNodeHead) {
    if (kmerReadNodeHead->realCount <= 0) {
        std::cout << kmerReadNodeHead->realCount << std::endl;
        return 0;
    }
    for (long int i = 0; i < kmerReadNodeHead->realCount - 1; i++) {
        std::cout << kmerReadNodeHead->kmerReadNode[i].kmer << std::endl;
        std::cout << kmerReadNodeHead->kmerReadNode[i + 1].kmer << std::endl;
        if (kmerReadNodeHead->kmerReadNode[i].kmer < kmerReadNodeHead->kmerReadNode[i + 1].kmer) {
            return 0;
        }
    }
    return 1;
}

#endif
