#ifndef Read_CPP_INCLUDED
#define Read_CPP_INCLUDED
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef _WIN32
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include "read.h"

using namespace std;

// Funzione per restituire il valore massimo tra due numeri
long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

// Funzione per restituire il valore minimo tra due numeri
long int min(long int x, long int y) {
    return (x < y) ? x : y;
}

// @param StrLine -> array di size maxsize
// @param maxSize -> taglia dell'array (1000000)
// @param fp -> puntatore al file che contiene la read (long_read.fa)
ReadSetHead * GetReadSetHead(char *filename, char *StrLine, long int maxSize)
{
    printf("\nfunzione GetReadSetHead\n");
    ReadSetHead * readSetHead = (ReadSetHead *)malloc(sizeof(ReadSetHead));
    if (!readSetHead) {
        printf("Memory allocation error!\n");
        return NULL;
    }

    readSetHead->readSet = NULL;
    readSetHead->readCount = 0;

    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("error opening file!\n");
        free(readSetHead);
        return NULL;
    }

    long int readIndex = -1;
    long int allocatedReads = 0;

    // Conta le reads e popola la struttura readSetHead in un unico ciclo
    while ((fgets(StrLine, maxSize, fp)) != NULL) { // legge maxSize caratteri in ogni ciclo
        if (StrLine[0] == '>') { // se il primo carattere è > -> è una read -> incrementa il contatore
            readSetHead->readCount++;
            readIndex++;
            if (readIndex >= allocatedReads) {
                // Rialloca memoria per readSet
                allocatedReads += 100; // incrementa di 100000 alla volta per ridurre le riallocazioni
                readSetHead->readSet = (ReadSet *)realloc(readSetHead->readSet, sizeof(ReadSet) * allocatedReads);
                if (!readSetHead->readSet) {
                    printf("Memory allocation error!\n");
                    fclose(fp);
                    free(readSetHead);
                    return NULL;
                }
                // Inizializza le nuove reads allocate
                for (long int i = readIndex; i < allocatedReads; i++) {
                    readSetHead->readSet[i].readLength = 0;
                    readSetHead->readSet[i].read = NULL; // Inizializza i puntatori a NULL
                }
            }
            continue;
        }

        // Elimina il carattere di terminazione \n o \r
        long int readLength = strlen(StrLine);
        while (readLength > 0 && (StrLine[readLength - 1] == '\n' || StrLine[readLength - 1] == '\r')) {
            readLength--;
        }
        StrLine[readLength] = '\0'; // Aggiunge il terminatore nullo

        // Popolare la struttura readSetHead
        if (readIndex >= 0 && readIndex < readSetHead->readCount) {
            readSetHead->readSet[readIndex].readLength = readLength;
            readSetHead->readSet[readIndex].read = (char *)malloc(sizeof(char) * (readLength + 1));
            if (!readSetHead->readSet[readIndex].read) {
                printf("Memory allocation error!\n");
                fclose(fp);
                for (long int i = 0; i <= readIndex; i++) {
                    free(readSetHead->readSet[i].read);
                }
                free(readSetHead->readSet);
                free(readSetHead);
                return NULL;
            }
            strncpy(readSetHead->readSet[readIndex].read, StrLine, readLength);
            readSetHead->readSet[readIndex].read[readLength] = '\0';
        }
    }

    fclose(fp);

    printf("Number of long reads: %ld;\n", readSetHead->readCount);

    return readSetHead;
}


#endif
