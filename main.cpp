#include <unistd.h>
#include <stdio.h>
#ifdef _WIN32
#include <malloc.h>
#else
#include <stdlib.h>
#endif
#include <string.h>
#include <iostream>
#include <ctype.h>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <chrono>
#include <sstream>
#include <sys/resource.h>
#include "kmer.h"
#include "read.h"
#include "aligning.h"
#include "bitarray.h"
#include <nthash/nthash.hpp>

using namespace std;
using namespace nthash;
using namespace chrono;

// Funzione per ottenere l'uso della memoria
long getMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // La memoria massima residente in KB
}

// Funzione per scrivere i benchmark su un file CSV
void writeBenchmarkToCSV(double elapsedSeconds, long memoryUsage) {
    ofstream file("benchmark_results_w303.csv");
    if (file.is_open()) {
        file << "Running time,Memory (Mb)\n";
        file << elapsedSeconds << " s," << (memoryUsage / 1024) << " Mb\n";
        file.close();
    }
}

void print_usage() {
    printf("\nLROD può rilevare regioni di sovrapposizione tra letture lunghe.\n");
    printf("\nIl formato del comando è il seguente: ");
    printf("\nLROD -r <long-read-file> -c <kmer-frequency-file> -o result-file [opzioni]\n");
    printf("\nUso:\n");
    printf("\t-r long-read-file: file di input in formato FASTA\n");
    printf("\t-c kmer-frequency-file: ogni riga nel file delle frequenze dei kmer deve essere \"kmer kmer-frequency\"\n");
    printf("\t-o result-file: file di risultati\n");
    printf("\t-t count: numero di thread (predefinito 1)\n");
    printf("\t-k kmerLength: lunghezza del kmer (predefinito 15)\n");
    printf("\t-q smallKmerLength: lunghezza del piccolo kmer (predefinito 9)\n");
    printf("\t-f minimumKmerFrequency: frequenza minima del kmer (predefinito 2)\n");
    printf("\t-m maxKmerFrequencyRatio: rapporto massimo di frequenza del kmer (deve essere inferiore a 1, predefinito 0.9)\n");
    printf("\t-s kmerStep: passo del kmer (predefinito 1)\n");
    printf("\t-d distance: piccola distanza usata per determinare se due kmer comuni sono coerenti (predefinito 400)\n");
    printf("\t-e distance: grande distanza usata per determinare se due kmer comuni sono coerenti (predefinito 1500)\n");
    printf("\t-a min-overlap-length: lunghezza minima di sovrapposizione tra due letture lunghe (predefinito 500)\n");
    printf("\t-b length-ratio: rapporto massimo di lunghezza tra due regioni allineate (predefinito 0.3)\n");
    printf("\t-h, -help\n");
    printf("\t--generate-kmer-file-only: genera solo il file delle frequenze dei kmer e termina\n");
}

std::unordered_map<std::string, int> computeKmerFrequencies(const std::string& sequence, int kmerLength, int step) {
    std::unordered_map<std::string, int> kmerFrequencies;
    NtHash nth(sequence, 1, kmerLength);

    while (nth.roll()) {
        std::string kmer = sequence.substr(nth.get_pos(), kmerLength);
        kmerFrequencies[kmer]++;
    }

    return kmerFrequencies;
}

std::string readLongReadFile(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string sequence;
    std::string line;

    while (std::getline(file, line)) {
        if (line[0] != '>') {
            sequence += line;
        }
    }

    return sequence;
}

void generateKmerFrequencyFile(const std::string& inputFilePath, const std::string& outputFilePath, int kmerLength, int step) {
    std::ifstream inputFile(inputFilePath);
    std::ofstream outputFile(outputFilePath);
    std::string sequence = readLongReadFile(inputFilePath);
    std::unordered_map<std::string, int> kmerFrequencies = computeKmerFrequencies(sequence, kmerLength, step);

    for (const auto& entry : kmerFrequencies) {
        outputFile << entry.first << " " << entry.second << "\n";
    }
}

void convertToFASTA(const std::string& inputFilePath, const std::string& outputFilePath) {
    std::ifstream inputFile(inputFilePath);
    std::ofstream outputFile(outputFilePath);
    
    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    std::string line;
    int kmerCount = 1;

    while (std::getline(inputFile, line)) {
        size_t spaceIndex = line.find(' ');
        if (spaceIndex != std::string::npos) {
            std::string kmer = line.substr(0, spaceIndex);
            outputFile << ">kmer" << kmerCount << std::endl;
            outputFile << kmer << std::endl;
            kmerCount++;
        }
    }

    inputFile.close();
    outputFile.close();
    std::cout << "Conversion completed. Output file: " << outputFilePath << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage();
        return 1;
    }

    long int maxSize = 1000000;
    char* StrLine = (char*)malloc(sizeof(char) * maxSize);
    char* readFile = NULL;
    char* kmerFrequencyFile = NULL;
    char* outputKmerFile = (char*)malloc(sizeof(char) * 500);
    FILE* fp4;

    strcpy(outputKmerFile, "output");

    long int step = 1;
    long int kmerLength = 15;
    long int smallKmerLength = 9;
    long int threadCount = 1;

    long int smallIntervalDistance = 400;
    long int largeIntervalDistance = 1500;
    long int overlapLengthCutOff = 500;
    float lengthRatio = 0.3;
    int frequencyCutOff = 3;
    long int minimumKmerFrequency = 2;
    float maxKmerFrequencyRatio = 0.9;
    bool generateKmerFileOnly = false;

    struct option long_options[] = {
        {"readFile", required_argument, NULL, 'r'},
        {"kmerFrequencyFile", required_argument, NULL, 'c'},
        {"outputKmerFile", optional_argument, NULL, 'o'},
        {"step", optional_argument, NULL, 's'},
        {"kmerLength", optional_argument, NULL, 'k'},
        {"smallKmerLength", optional_argument, NULL, 'q'},
        {"minimumKmerFrequency", optional_argument, NULL, 'f'},
        {"maxKmerFrequencyRatio", optional_argument, NULL, 'm'},
        {"threadCount", optional_argument, NULL, 't'},
        {"smallIntervalDistance", optional_argument, NULL, 'd'},
        {"largeIntervalDistance", optional_argument, NULL, 'e'},
        {"overlapLengthCutOff", optional_argument, NULL, 'a'},
        {"lengthRatio", optional_argument, NULL, 'b'},
        {"generate-kmer-file-only", no_argument, NULL, 'g'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}
    };

    int ch = 0;

    while ((ch = getopt_long(argc, argv, "c:r:o:m:n:d:k:e:a:s:t:f:q:b:gh", long_options, NULL)) != -1) {
        switch (ch) {
            case 'r': readFile = (char*)(optarg); break;
            case 'c': kmerFrequencyFile = (char*)(optarg); break;
            case 'o': outputKmerFile = (char*)optarg; break;
            case 'k': kmerLength = atoi(optarg); break;
            case 'q': smallKmerLength = atoi(optarg); break;
            case 'f': minimumKmerFrequency = atoi(optarg); break;
            case 'm': maxKmerFrequencyRatio = atof(optarg); break;
            case 's': step = atoi(optarg); break;
            case 'd': smallIntervalDistance = atoi(optarg); break;
            case 'e': largeIntervalDistance = atoi(optarg); break;
            case 'a': overlapLengthCutOff = atoi(optarg); break;
            case 't': threadCount = atoi(optarg); break;
            case 'b': lengthRatio = atof(optarg); break;
            case 'g': generateKmerFileOnly = true; break;
            case 'h':
                print_usage();
                return 0;
            default:
                return -1;
        }
    }

    if (minimumKmerFrequency >= maxKmerFrequencyRatio) {
        printf("Regolazione dei parametri: impostazione di maxKmerFrequencyRatio su minimumKmerFrequency + 0.1\n");
        maxKmerFrequencyRatio = minimumKmerFrequency + 0.1;
    }

    if ((fp4 = fopen(readFile, "r")) == NULL) {
        printf("File delle letture lunghe mancante!\n");
        printf("Per un uso dettagliato di LROD, utilizzare il comando: -h o -help!\n");
        return 2;
    }
    fclose(fp4);

    if ((fp4 = fopen(kmerFrequencyFile, "r")) == NULL) {
        printf("Generazione del file delle frequenze dei kmer...\n");
        generateKmerFrequencyFile(readFile, kmerFrequencyFile, kmerLength, step);
    } else {
        fclose(fp4);
    }

    // Se l'opzione --generate-kmer-file-only è stata specificata, termina l'esecuzione qui
    if (generateKmerFileOnly) {
        printf("File delle frequenze dei kmer generato: %s\n", kmerFrequencyFile);
        return 0;
    }

    strcat(outputKmerFile, ".csv");

    if ((fp4 = fopen(outputKmerFile, "w")) == NULL) {
        printf("%s, non esiste!\n", outputKmerFile);
        return 5;
    }
    fclose(fp4);

    printf("\nInizio caricamento delle letture lunghe!\n");
    printf("funzione GetReadSetHead\n");

    auto start = high_resolution_clock::now();
    long initialMemoryUsage = getMemoryUsage();

    ReadSetHead* readSetHead = GetReadSetHead(readFile, StrLine, maxSize);

    if (readSetHead->readCount <= 1) {
        printf("Il numero di letture è inferiore a uno!\n");
        return 6;
    }

    printf("Inizio costruzione della tabella hash dei kmer!\n");

    KmerHashTableHead* kmerHashTableHead = GetKmerHashTableHead(kmerFrequencyFile, readSetHead, kmerLength, step, minimumKmerFrequency, maxKmerFrequencyRatio);

    if (GetKmerHashTableHead_UnitTest(kmerHashTableHead) == 0) {
        printf("Errore nella costruzione della tabella hash dei kmer!\n");
        return 7;
    }

    long int subReadCount = 50000;

    KmerReadNodeHead* kmerReadNodeHead = GetKmerReadNodeHeadSub(readSetHead, kmerLength, step, subReadCount);

    printf("Preparazione per rilevare sovrapposizioni tra letture lunghe!\n");
    GetCommonKmerHeadAllThreadNew(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputKmerFile, step, threadCount, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio, subReadCount);

    if (GetCommonKmerHeadAllThreadNew_UnitTest(outputKmerFile, readSetHead->readCount) == 0) {
        return 8;
    }

    auto end = high_resolution_clock::now();
    long finalMemoryUsage = getMemoryUsage();
    double elapsedSeconds = duration_cast<seconds>(end - start).count();
    long memoryUsage = finalMemoryUsage - initialMemoryUsage;

    writeBenchmarkToCSV(elapsedSeconds, memoryUsage);

    printf("Fatto!\n");
    printf("Nome del file di risultati è: %s!\n", outputKmerFile);

    return 0;
}
