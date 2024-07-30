#ifndef Aligning_CPP_INCLUDED 
#define Aligning_CPP_INCLUDED 
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
#include <ctype.h>
#include <time.h>

#include "aligning.h"
#include "read.h"

using namespace std;
//Reallocate space for CommonKmerHead
//Rialloca dinamicamente lo spazio di memoria per la struttura dati CommonKmerHead
//In sostanza, questa funzione è utilizzata per gestire la crescita dinamica della struttura dati CommonKmerHead, consentendo di adattarsi a una quantità crescente di dati senza sprechi di memoria o perdite di informazioni.
void ReAllocateCommonKmer(CommonKmerHead * commonKmerHead){
	unsigned long int maxCount = commonKmerHead->allocationCount*1.5; //Calcola la nuova dimensione massima per l'allocazione come 1,5 volte la dimensione attuale 
	CommonKmer * commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*maxCount); //Alloca memoria per un nuovo array di CommonKmer con la nuova dimensione (maxCount).
	for(long int i = 0; i < commonKmerHead->allocationCount; i++){
		//Copia i dati dall'array esistente commonKmerHead->commonKmer al nuovo array commonKmer
		commonKmer[i].readIndex = commonKmerHead->commonKmer[i].readIndex;
		commonKmer[i].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
		commonKmer[i].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
		commonKmer[i].orientation = commonKmerHead->commonKmer[i].orientation;
	}
	free(commonKmerHead->commonKmer); //Libera la memoria occupata dall'array precedente commonKmerHead->commonKmer
	commonKmerHead->commonKmer = commonKmer; //Assegna il puntatore al nuovo array commonKmer all'interno della struttura commonKmerHead.
	commonKmerHead->allocationCount = maxCount; //Aggiorna la variabile allocationCount nella struttura commonKmerHead con la nuova dimensione massima.
	
}

//When two reads have the same kmer, store the position and direction of the kmer on the two reads to commonKmerHead
//Questa funzione inserisce i k-mer comuni tra due read nell'oggetto CommonKmerHead. Le read considerate sono quella su cui stiamo iterando (R1) e quella nella kmerReadNode.
void InsertCommonToTwoReadAligningHead(CommonKmerHead * commonKmerHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int hashIndex, unsigned long int readIndex, unsigned long  long int position, bool orien){ //readIndex è i+1
	//printf("\n funzione InsertCommonToReadAligningHead\n");
	
	/*if (hashIndex != -1)
		printf("\n hashIndex in InsetCommonTo: %ld\n",hashIndex);*/

	//Controlla se l'indice di inizio della posizione nell'array dell'hash table è diverso da -1. Se lo è, significa 
	//che ci sono k-mer comuni da gestire; altrimenti, non fa nulla e termina la funzione.
	if(kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray == -1){ 
		return; //controlla se nella read c'è il kmer
	}

	//Altro errore (?), anche se l'indice della tabella hash è -1, lui esegue il codice. TEST:
	/*if (hashIndex == -1){
		printf("kmer in hash table in poszione %ld :%ld\n",hashIndex, kmerHashTableHead->kmerHashNode[hashIndex].kmer);
	}*/

	unsigned long int i = kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray; //la posizione del kmer nella struttura kmerRreadNode,
	unsigned long int kmer = kmerHashTableHead->kmerHashNode[hashIndex].kmer - 1; //kmerInteger in posizione hashIndex, nella read i-esima
	//printf("\n kmer = %ld, start position del kmer in kmerHashTable: %ld\n",kmer,i);

	//il kmer è presente nella hash table, ora scorro kmerReadNode fino a quando non trovo il kmer
	//nella kmerReadNode ci sono i kmer, la posizione all'interno della read e la readIndex di quelli presenti nella tabella hash
	for(; i < kmerReadNodeHead->realCount; i++){ 
		//se il k-mer nella kmerReadNode è diverso dal k-mer che stiamo valutando, si passa al k-mer successivo (vedere funzione chiamante)
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmer){ 
			//printf("\n kmerReadNode i = %ld\n", kmerReadNodeHead->kmerReadNode[i].kmer);
			break;
		}

		//controlla se readIndex si trova nell'intervallo delle reads contenute nella kmerReadNode
		if(kmerReadNodeHead->kmerReadNode[i].readIndex <= readIndex && readIndex <= kmerReadNodeHead->endReadIndex){             
			continue; //se readIndex non si trova nell intervallo compreso tra kmerReadNode[i].readIndex e kmerReadNodeHead -> endReadIndex, passo al kmerReadNode[i] successivo
		}
		//Altrimenti 
		//Assegna l'indice della read corrente all'array di k-mer comuni.
		//gli assegnamenti vengono fatti solo se il kmer non rientra nell'intervallo specificato in precedenza
		//l'indice della read in kmerReadNode lo salva in commonKmerHead (ovvimanete stiamo analizzando kmer uguali, presenti nell'itervallo che ci interessa)
		commonKmerHead->commonKmer[commonKmerHead->realCount].readIndex = kmerReadNodeHead->kmerReadNode[i].readIndex;
		
		// printf("\n kmer in kmer readNode %ld readIndex in kmerReadNodeHead: %ld\n",kmerReadNodeHead->kmerReadNode[i].kmer, kmerReadNodeHead->kmerReadNode[i].readIndex);
		// printf("\n readIndex in commonkmer: %ld\n",commonKmerHead->commonKmer[commonKmerHead->realCount].readIndex);

		//Assegna la posizione del k-mer nella read destra all'array di k-mer comuni.
		//rigthPosition è la posizione del kmer comune della read considerata da kmerReadNode (chiamiamola R2)
		commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = kmerReadNodeHead->kmerReadNode[i].position; //position in kmerReadNode è la posizione del kmer all'interno della readConsiderata in kmerReadNode

		//Assegna la posizione del k-mer nella read sinistra all'array di k-mer comuni.
		//leftPosition è la posizione j della read corrente (chiamiamola R1)
		commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = position; //position è posizione del kmer nella read corrente

		//Verifica l'orientazione del k-mer e imposta il valore di orientazione dell'array di k-mer comuni di conseguenza.
		if(orien == kmerReadNodeHead->kmerReadNode[i].orientation){ //Controlla se l'orientazione del k-mer attuale è uguale all'orientazione del k-mer rappresentato dall'elemento corrente nell'array kmerReadNode.
			commonKmerHead->commonKmer[commonKmerHead->realCount].orientation = 0; //Se l'orientazione è la stessa, viene impostato il valore di orientazione dell'array di k-mer comuni a 0. Questo potrebbe indicare che il k-mer è in una posizione normale nella read.
		}else{
			commonKmerHead->commonKmer[commonKmerHead->realCount].orientation = 1; //Se l'orientazione è diversa, viene impostato il valore di orientazione dell'array di k-mer comuni a 1. Questo potrebbe indicare che il k-mer è in una posizione invertita o complementare nella read.
		}
		commonKmerHead->realCount++; //Incrementa il contatore realCount dell'oggetto CommonKmerHead per tener traccia del numero totale di k-mer comuni inseriti.
		
		//commonKmerHead -> size++; //variabile aggiunta da noi per controllare l'effettiva size della commonKmerHead 
		if(commonKmerHead->realCount >= commonKmerHead->allocationCount){ //Verifica se il numero di k-mer comuni supera la capacità attuale dell'oggetto CommonKmerHead.
			ReAllocateCommonKmer(commonKmerHead); //Se sì, chiama la funzione ReAllocateCommonKmer per riallocare dinamicamente la memoria e aumentare la capacità dell'array di k-mer comuni.
		}
	}

	
}

// funzione utilizzata per effettuare un test unitario sui risultati generati dalla funzione GetCommonKmerHeadAllThreadNew
int GetCommonKmerHeadAllThreadNew_UnitTest(char * result, long int readCount){
	FILE * fp = NULL;
	if((fp = fopen(result,"r")) == NULL){ //Apre il file risultante in modalità lettura.
        return 0;
   	}
	long int maxSize = 20000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * p =NULL;
	const char * split = ","; //Per ogni riga, suddivide la riga utilizzando la virgola come delimitatore.
	
	//Legge ogni riga del file e la analizza.
	while((fgets(line, maxSize, fp)) != NULL){ 

		p = strtok(line,split); //La riga appena letta viene suddivisa in token utilizzando la virgola come delimitatore. Il puntatore p viene assegnato al primo token estratto dalla riga.
		//Il token estratto viene convertito in intero con atoi()
		if(atoi(p) < 0 || atoi(p) > readCount){ //Verifica che gli indici delle read siano compresi tra 0 e il numero totale di read (readCount). Questo controllo garantisce che i riferimenti alle read siano validi e non superino l'indice massimo.
			return 0;
		}
		
		//Dopo il primo controllo dell'indice della read, la funzione strtok() viene chiamata nuovamente per estrarre il token successivo dalla stessa riga. Questo passaggio avanza nel processo di tokenizzazione della riga, consentendo di accedere agli altri valori presenti nella riga e di effettuare ulteriori controlli su di essi.
		//Questo processo di estrazione dei token continua fino a quando non si esauriscono i token nella riga. Ogni token viene convertito in un intero con atoi() e controllato per assicurarsi che cada nell'intervallo valido degli indici delle letture.
		p = strtok(NULL,split);
		if(atoi(p) < 0 || atoi(p) > readCount){
			return 0;
		}
		//Dopo aver controllato gli indici delle letture, il codice procede ad estrarre il token successivo con strtok(NULL, split). Questo token viene quindi convertito in un intero con atoi() e controllato per garantire che rappresenti un valore di orientazione valido, cioè 0 o 1. Se il valore non corrisponde a nessuna di queste opzioni, la funzione restituisce 0, indicando un fallimento nel test dell'unità.
		p = strtok(NULL,split);
		if(atoi(p) != 0 && atoi(p) != 1){ //Assicura che l'orientazione sia 0 o 1.
			return 0;
		}
		
		p = strtok(NULL,split);
		if(atoi(p) < 0){ //Verifica che le posizioni dei k-mer siano non negative. 
			return 0;
		}
		
		p = strtok(NULL,split);
		if(atoi(p) < 0){
			return 0;
		}
		
		p = strtok(NULL,split);
		if(atoi(p) < 0){
			return 0;
		}
		
		p = strtok(NULL,split);
		if(atoi(p) < 0){
			return 0;
		}
		
	
	}
	
	// Se tutti i controlli su una riga hanno successo, il ciclo continua con la prossima riga. Se uno qualsiasi dei controlli fallisce per una qualsiasi riga, il test ritorna 0, indicando un problema nei risultati ottenuti. Se tutti i controlli hanno successo per tutte le righe, il test ritorna 1, indicando che i risultati sono validi.
	fclose(fp);
	return 1;
	
}

//@param KmerHashTableHead -> tabella hash dei k-mers
//@param KmerReadNodeHead -> lista dei k-mers nelle read (Ancora da riempire)
//@param ReadSetHead -> lista di reads 
//@param char* outputFile -> file dove viene scritto l'overlap
//@param int largeIntervalDistance = 1500
//@param int kmerLenght -> lunghezza dei k-mers
//@param int step -> valore per leggere i k-mers, settato a 1
CommonKmerHead * GetCommonKmerHeadAllThreadNew(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int subReadCount){
	printf("\n funzione GetCommonKmerHeadAllThreadNew\n");
	long int startReadIndex = 0;
	//subReadCount = 50.000
	long int endReadIndex = subReadCount;
	FILE * fp = NULL;
	if((fp = fopen(outputFile,"w")) == NULL){ //Viene aperto un file di output specificato dal parametro outputFile in modalità scrittura ("w"). Se il file non può essere aperto correttamente, viene stampato un messaggio di errore e il programma termina.
        printf("%s, does not exist!", outputFile);
        exit(0);
   	}
	fclose(fp);
	printf("\n ciclo in GetCommonKmerHeadAllThreadNew\n");
	while(true){ //Viene eseguito un ciclo while(true) per processare le read in batch, cioè un sottoinsieme alla volta. Questo ciclo continua finché non vengono processate tutte le read.
		//All'interno del ciclo, vengono definiti gli indici di inizio e fine del batch corrente in base al parametro subReadCount, che specifica il numero di read da elaborare in ogni batch.
		//startReadIndex = 0
		if(startReadIndex >= readSetHead->readCount){
			break;
		}
		//endReadIndex = 50.000
		if(endReadIndex >= readSetHead->readCount){
			//endReadIndex = 3
			endReadIndex = readSetHead->readCount - 1;
		}
		
		//Viene inizializzata la struttura dati KmerReadNodeHead per contenere i k-mer relativi al batch corrente. Questa inizializzazione coinvolge la chiamata a InitKmerReadNodeHeadSub
		InitKmerReadNodeHeadSub(readSetHead, kmerReadNodeHead, kmerHashTableHead, kmerLength, step, startReadIndex, endReadIndex);
		//Viene chiamata la funzione GetCommonKmerHeadAllThread per trovare i k-mer comuni tra le read nel batch corrente. Questa funzione sembra essere responsabile dell'esecuzione del processo di rilevamento dei k-mer comuni utilizzando più thread.
		GetCommonKmerHeadAllThread(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputFile, step, totalThreadNumber, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio,startReadIndex);
		
		//Gli indici di inizio e fine vengono aggiornati per puntare al prossimo batch di read.
		startReadIndex = endReadIndex + 1;
		endReadIndex = endReadIndex + subReadCount;

	} // Il ciclo continua finché non vengono processate tutte le read.
}

//Alignment function between sequences
//Questa funzione sembra essere progettata per coordinare il processo di rilevamento dei k-mer comuni tra le read utilizzando più thread. 
//La funzione prende diversi parametri, tra cui puntatori a strutture dati come KmerHashTableHead, KmerReadNodeHead e ReadSetHead, insieme ad altre informazioni come la lunghezza dei k-mer, i file di input e output, i parametri per il controllo degli intervalli, ecc.
CommonKmerHead * GetCommonKmerHeadAllThread(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int startReadIndex){
	printf("\n Funzione GetCommonKmerHeadAllThread\n");
	printf ("totalThreadNumber = %ld\n", totalThreadNumber);
	//TotalThreadNumber = 1
	pthread_t tid[totalThreadNumber]; //Viene creato un array di thread tid con una dimensione pari al numero totale di thread specificato dal parametro totalThreadNumber.
    
    long int i = 0;

    GetCommonKmerHeadP * getCommonKmerHeadP = new GetCommonKmerHeadP[totalThreadNumber]; //è una struttura che contiene i parametri necessari per eseguire il processo di rilevamento dei k-mer comuni in un singolo thread. 
    long int * threadSignal = new long int[totalThreadNumber]; //non usato
    
	//Viene stampato un messaggio che indica il numero di thread abilitati, insieme ad altre informazioni sui parametri passati alla funzione.
	printf("The number of enabled threads is: %ld;\n",totalThreadNumber); //1
	printf("The length of kmer is: %ld;\n",kmerLength); //15
	printf("The length of smallkmer is: %ld;\n",smallKmerLength); //9
	printf("The minimum overlap length of each pair reads is: %ld;\n",overlapLengthCutOff); //500
		
	//Viene inizializzato un ciclo for per creare e avviare i thread. Per ogni thread, vengono passati i parametri appropriati attraverso la struttura GetCommonKmerHeadP.
    for(i = 0; i< totalThreadNumber; i++){
        getCommonKmerHeadP[i].kmerHashTableHead = kmerHashTableHead;
        getCommonKmerHeadP[i].kmerReadNodeHead = kmerReadNodeHead;
		getCommonKmerHeadP[i].readSetHead = readSetHead;
		getCommonKmerHeadP[i].kmerLength = kmerLength;
		getCommonKmerHeadP[i].readFile = readFile;
		getCommonKmerHeadP[i].smallIntervalDistance = smallIntervalDistance;
		getCommonKmerHeadP[i].largeIntervalDistance = largeIntervalDistance;
		getCommonKmerHeadP[i].overlapLengthCutOff = overlapLengthCutOff;
		getCommonKmerHeadP[i].smallKmerLength = smallKmerLength;
		getCommonKmerHeadP[i].lengthRatio = lengthRatio;
		getCommonKmerHeadP[i].startReadIndex = startReadIndex;
		
		//Questo frammento di codice è responsabile della creazione di nomi univoci per i file di output di ciascun thread.
		//viene creato un nome di file temporaneo per il thread corrente, in modo che ciascun thread scriva i suoi risultati in un file separato. Questi file temporanei vengono successivamente combinati in un unico file di output.
		char * outputFileTemp = (char *)malloc(sizeof(char)*(strlen(outputFile)+10));

		//La funzione sprintf viene utilizzata per scrivere il nome del file di output univoco nel formato "outputFile-threadIndex", 
		//dove outputFile è il nome del file di output originale e threadIndex è l'indice del thread corrente. Questo assicura che ogni thread scriva i suoi risultati in un file separato.
		//stampa in outputFileTemp
		sprintf(outputFileTemp, "%s-%ld", outputFile, i);
		getCommonKmerHeadP[i].outputFile = outputFileTemp;
		
		//Successivamente, vengono impostati altri parametri nella struttura getCommonKmerHeadP[i], come step, threadIndex e totalThreadNumber, che vengono utilizzati per il corretto funzionamento della logica del thread.
		getCommonKmerHeadP[i].step = step;
        getCommonKmerHeadP[i].threadIndex = i;
        getCommonKmerHeadP[i].totalThreadNumber = totalThreadNumber;

        //Call multiple threads for sequence overlap detection
		//Viene chiamata la funzione pthread_create per creare un thread. Se questa operazione fallisce, il programma esce con un messaggio di errore.
		//&tid[i]: un puntatore alla variabile tid[i], che conterrà l'ID del nuovo thread creato.
		//NULL: un puntatore a un oggetto pthread_attr_t che specifica gli attributi del nuovo thread. In questo caso, viene utilizzato NULL per indicare l'uso degli attributi predefiniti.
		//GetCommonKmerHeadThread: la funzione che il thread eseguirà. Questa funzione è responsabile della logica di elaborazione dei k-mer comuni per il thread corrente.
		//(void *)&getCommonKmerHeadP[i]: un puntatore generico al parametro che verrà passato alla funzione del thread. In questo caso, getCommonKmerHeadP[i] è un puntatore alla struttura GetCommonKmerHeadP specifica per il thread corrente.
        if(pthread_create(&tid[i], NULL, GetCommonKmerHeadThread, (void *)&getCommonKmerHeadP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
       }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL); //Dopo aver creato tutti i thread, la funzione attende che ciascun thread termini utilizzando pthread_join.
    }
	
	
	FILE * fp = NULL; //fp è un puntatore al file di output principale

	long int Size = 10000;
	//Viene allocata memoria per una stringa str utilizzata per leggere i contenuti dai file temporanei dei singoli thread.
	char * str = (char *)malloc(sizeof(char)*Size);
	for(i = 0; i < totalThreadNumber; i++){
        FILE * fpTemp = NULL; //fpTemp viene utilizzato per aprire e leggere i file temporanei generati dai singoli thread.
		
		//Per ogni thread, il relativo file temporaneo viene aperto in modalità di lettura ("r") utilizzando il percorso specificato dalla variabile outputFile all'interno della struttura getCommonKmerHeadP[i]. Se il file temporaneo non può essere aperto correttamente, il programma stampa un messaggio di errore indicando il nome del file e termina l'esecuzione.
		if((fpTemp = fopen(getCommonKmerHeadP[i].outputFile,"r")) == NULL){
        	printf("%s, does not exist!", getCommonKmerHeadP[i].outputFile);
        	exit(0);
    	}
		
		if((fp = fopen(outputFile,"a")) == NULL){
        	printf("%s, does not exist!", outputFile);
        	exit(0);
   		}
		setbuf(fp,NULL); //Una volta aperto il file principale, il buffer del file viene disabilitato con setbuf(fp, NULL) per garantire che i dati vengano scritti direttamente nel file senza essere bufferizzati.
		
		//Questo ciclo legge righe dal file temporaneo del thread (fpTemp) e le scrive nel file di output principale (fp)
		//Viene iterato su ciascun thread, aprendo il relativo file temporaneo e scrivendo i suoi contenuti nel file di output principale. Successivamente, i file temporanei vengono chiusi e rimossi.
		while((fgets(str, Size, fpTemp)) != NULL){ //Utilizza fgets per leggere una riga dal file temporaneo (fpTemp) e la memorizza nella stringa str.
			setbuf(fp,NULL);
			long int len = strlen(str);
			if(str[len - 1]=='\n' || str[len - 1]=='\r'){ //Verifica se l'ultima carattere della riga (tranne il carattere di terminazione) è un newline (\n) o un carriage return (\r).
				str[len - 1] = '\0'; //Se sì, rimuove tale carattere sostituendolo con il terminatore di stringa \0. Questo serve per garantire che non ci siano newline aggiuntivi tra le righe del file di output.
			}
			fprintf(fp, "%s\n", str);
			fflush(fp);

		}
		fflush(fp);
		fclose(fp);
		fclose(fpTemp);
		remove(getCommonKmerHeadP[i].outputFile);
	}
	
}

void * GetCommonKmerHeadThread(void * arg){ 
    printf("Funzione GetCommonKmerHeadThread\n");
    GetCommonKmerHeadP * getCommonKmerHeadP = (GetCommonKmerHeadP *)arg; //il puntatore arg viene reinterpretato come un puntatore a GetCommonKmerHeadP per accedere ai parametri passati alla funzione.
     
    //Estrazione dei parametri dalla struttura getCommonKmerHeadP
    KmerHashTableHead * kmerHashTableHead = getCommonKmerHeadP->kmerHashTableHead;
    KmerReadNodeHead * kmerReadNodeHead = getCommonKmerHeadP->kmerReadNodeHead;
    ReadSetHead * readSetHead = getCommonKmerHeadP->readSetHead;
    long int kmerLength = getCommonKmerHeadP->kmerLength;
    char * readFile = getCommonKmerHeadP->readFile;
    char * outputFile = getCommonKmerHeadP->outputFile;
    long int step = getCommonKmerHeadP->step;
    long int threadIndex = getCommonKmerHeadP->threadIndex;
    long int totalThreadNumber = getCommonKmerHeadP->totalThreadNumber;
    long int smallIntervalDistance = getCommonKmerHeadP->smallIntervalDistance;
    long int largeIntervalDistance = getCommonKmerHeadP->largeIntervalDistance;
    long int overlapLengthCutOff = getCommonKmerHeadP->overlapLengthCutOff;
    long int startReadIndex = getCommonKmerHeadP->startReadIndex;

    //Allocazione della memoria per la struttura AdjGraphHead e inizializzazione dei suoi membri
    AdjGraphHead * G = (AdjGraphHead *)malloc(sizeof(AdjGraphHead)); //Viene allocata la memoria per la struttura AdjGraphHead.
    G->allocationCountGraph = 20000; //numero massimo di grafi adiacenti che possono essere memorizzati.
    G->graph = (AdjGraph*)malloc(sizeof(AdjGraph)* G->allocationCountGraph); //Viene allocata la memoria per l'array graph, che conterrà i grafi adiacenti.
    G->realCountGraph = 0; //realCountGraph viene inizializzato a 0, indicando che non ci sono ancora grafi adiacenti presenti.
    G->reverseAllocationCountGraph = 20000; //Viene inizializzato reverseAllocationCountGraph con il valore 20000, che rappresenta il numero massimo di grafi adiacenti invertiti che possono essere memorizzati.
    G->reverseGraph = (AdjGraph*)malloc(sizeof(AdjGraph)* G->reverseAllocationCountGraph); //Viene allocata la memoria per l'array reverseGraph, che conterrà i grafi adiacenti invertiti.
    G->reverseRealCountGraph = 0; //reverseRealCountGraph viene inizializzato a 0, indicando che non ci sono ancora grafi adiacenti invertiti presenti.
    G->allocationCountArc = 20000; //Viene inizializzato allocationCountArc con il valore 20000, che rappresenta il numero massimo di archi che possono essere memorizzati.
    G->arcIndex = (ArcIndex*)malloc(sizeof(ArcIndex)* G->allocationCountArc); //Viene allocata la memoria per l'array arcIndex, che conterrà gli indici degli archi.
    G->realCountArc = 0; //realCountArc viene inizializzato a 0, indicando che non ci sono ancora archi presenti.
    G->nodeCount = 0; //nodeCount viene inizializzato a 0, indicando che non ci sono ancora nodi presenti.
    G->kmerLength = getCommonKmerHeadP->kmerLength; //Il valore di kmerLength viene inizializzato utilizzando il valore presente in getCommonKmerHeadP.
    G->smallKmerLength = getCommonKmerHeadP->smallKmerLength; //Il valore di smallKmerLength viene inizializzato utilizzando il valore presente in getCommonKmerHeadP.
    G->lengthRatio = getCommonKmerHeadP->lengthRatio; //Il valore di lengthRatio viene inizializzato utilizzando il valore presente in getCommonKmerHeadP.
    G->overlapLengthCutOff = getCommonKmerHeadP->overlapLengthCutOff; //Il valore di overlapLengthCutOff viene inizializzato utilizzando il valore presente in getCommonKmerHeadP.
    
    G->largestIntervalDistance = smallIntervalDistance; //Il valore di largestIntervalDistance viene inizializzato utilizzando il valore di smallIntervalDistance (400)
    
    //Inizializzazione degli elementi degli array graph, reverseGraph e arcIndex:
    for(long int n=0;n<G->allocationCountGraph;n++){ //Itera su tutti gli elementi degli array graph, reverseGraph e arcIndex.
        G->graph[n].dataLeft = -1; 
        G->graph[n].dataRight = -1;
        
        G->graph[n].visit = 0; 
        G->reverseGraph[n].dataLeft = -1; 
        G->reverseGraph[n].dataRight = -1;

        G->reverseGraph[n].visit = 0; 
        G->arcIndex[n].startIndex = 0; 
        G->arcIndex[n].endIndex = 0;
    }
    
    //Allocazione e inizializzazione della struttura dati CommonKmerHead:
    
    long int maxCount = 100000; //Definisce il massimo numero di elementi che la struttura CommonKmerHead può contenere
    
    CommonKmerHead * commonKmerHead = (CommonKmerHead *)malloc(sizeof(CommonKmerHead)); //Alloca memoria per la struttura CommonKmerHead
    commonKmerHead->commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*maxCount); //Alloca memoria per l'array commonKmer all'interno della struttura CommonKmerHead.
    commonKmerHead->realCount = 0; //tiene traccia del numero reale di elementi presenti nell'array commonKmer.
    commonKmerHead->allocationCount = maxCount; //Imposta il membro allocationCount con il valore maxCount, che rappresenta la capacità massima dell'array commonKmer.
    commonKmerHead ->size = 0;

    //Allocazione e inizializzazione della struttura dati localG:
    AdjGraphHead * localG = (AdjGraphHead *)malloc(sizeof(AdjGraphHead)); //Alloca memoria per la struttura localG di tipo AdjGraphHead.
    localG->largestIntervalDistance = largeIntervalDistance; //Imposta il membro largestIntervalDistance di localG con il valore largeIntervalDistance.
    localG->allocationCountGraph = 2*largeIntervalDistance; //Imposta il membro allocationCountGraph di localG con il valore 2 * largeIntervalDistance.
    localG->graph = (AdjGraph*)malloc(sizeof(AdjGraph)* localG->allocationCountGraph); //Alloca memoria per l'array graph all'interno della struttura localG. L'array avrà una lunghezza pari a localG->allocationCountGraph e conterrà istanze di AdjGraph.
    localG->realCountGraph = 0; //Inizializza il contatore realCountGraph a 0, che tiene traccia del numero reale di elementi presenti nell'array graph.
    localG->allocationCountArc = 20000; //Imposta il membro allocationCountArc di localG con il valore 20000.
    localG->arcIndex = (ArcIndex*)malloc(sizeof(ArcIndex)* localG->allocationCountArc); //Alloca memoria per l'array arcIndex all'interno della struttura localG. L'array avrà una lunghezza pari a localG->allocationCountArc e conterrà istanze di ArcIndex.
    localG->realCountArc = 0; //Inizializza il contatore realCountArc a 0, che tiene traccia del numero reale di elementi presenti nell'array arcIndex.
    localG->nodeCount = 0; //Inizializza il contatore nodeCount a 0, che tiene traccia del numero di nodi nell'insieme dei nodi dell'grafo.
    //Alloca memoria per localLeftRead di localG, una stringa che sarà utilizzata per memorizzare le informazioni sulle sequenze di lettura.
    localG->localLeftRead = (char *)malloc(sizeof(char)*(localG->allocationCountGraph + 2*kmerLength));
    //Alloca memoria per localRightRead di localG, una stringa che sarà utilizzata per memorizzare le informazioni sulle sequenze di lettura.
    localG->localRightRead = (char *)malloc(sizeof(char)*(localG->allocationCountGraph + 2*kmerLength));
    //Imposta i membri lengthRatio, kmerLength e smallKmerLength di localG con i valori corrispondenti ottenuti da getCommonKmerHeadP.
    localG->lengthRatio = getCommonKmerHeadP->lengthRatio;
    localG->kmerLength = getCommonKmerHeadP->kmerLength;
    localG->smallKmerLength = getCommonKmerHeadP->smallKmerLength;

    //Allocazione e inizializzazione della struttura dati localCommonKmerHead
    CommonKmerHead * localCommonKmerHead = (CommonKmerHead *)malloc(sizeof(CommonKmerHead)); //Alloca memoria per la struttura localCommonKmerHead di tipo CommonKmerHead.
    localCommonKmerHead->allocationCount = 2*largeIntervalDistance; 

    //Alloca memoria per l'array commonKmer all'interno della struttura localCommonKmerHead.
    localCommonKmerHead->commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*localCommonKmerHead->allocationCount*2);
    localCommonKmerHead->realCount = 0; //Inizializza il contatore realCount a 0, che tiene traccia del numero reale di elementi presenti nell'array commonKmer.
    
    //Dichiarazione di variabili e apertura dei file di input e output:
    unsigned long int kmerInteger = 0; //Sarà utilizzata per memorizzare il valore intero corrispondente al k-mer durante la conversione.
    long int hashIndex; //utilizzata per memorizzare l'indice calcolato dal funzionamento dell'hashing.
    unsigned long int readIndex = 0; //Sarà utilizzata per tenere traccia dell'indice delle read.
    long int Size = 200000; //Sarà utilizzata per specificare la dimensione del buffer per la lettura dei file.
    long int readLength = 0; //Sarà utilizzata per memorizzare la lunghezza della read durante la lettura del file.
    
    char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1)); //memorizzare il k-mer
    char * kmer2 = (char *)malloc(sizeof(char)*(kmerLength + 1)); 

    FILE * fp = fopen(readFile,"r"); 
    FILE * fp1 = fopen(outputFile,"w");
    
    //Alloca memoria per un array di lunghezza readSetHead->readCount di interi long int e assegna il puntatore a forwardKmerCount. Questo array verrà utilizzato per memorizzare i conteggi dei k-mer in avanti.
    long int * forwardKmerCount = (long int *)malloc(sizeof(long int)*readSetHead->readCount);
    //Alloca memoria per un array di lunghezza readSetHead->readCount di interi long int e assegna il puntatore a reverseKmerCount. Questo array verrà utilizzato per memorizzare i conteggi dei k-mer rovesciati.
    long int * reverseKmerCount = (long int *)malloc(sizeof(long int)*readSetHead->readCount);
    
    char buf[1024]; //Dichiara una stringa di buffer buf di dimensione 1024. Questa stringa sarà utilizzata per leggere temporaneamente i dati dai file.
    time_t timep; // Sarà utilizzata per memorizzare l'istante temporale durante l'esecuzione.
    double sencond; //Sarà utilizzata per memorizzare il tempo trascorso in secondi durante l'esecuzione.
    int pi = 0; //Questa variabile verrà utilizzata per calcolare l'indice di partizione delle read per i thread.
    float per = 0.0; //Questa variabile verrà utilizzata per memorizzare la percentuale di avanzamento durante l'esecuzione.
    
    //Questa riga calcola il valore di pi, che rappresenta il numero di read per ogni thread. Diviosne del carico di lavoro
    pi = readSetHead->readCount/totalThreadNumber;
    
    //itera sul batch di reads
    printf("\n startReadIndex: %ld fino a readCount %ld\n",startReadIndex, readSetHead->readCount);
    long int i = 0;
    //il for più esterno itera sulle read presenti nel file
    for(i = startReadIndex; i < readSetHead->readCount; i++){
    
        printf("Agisco sulla read: %ld, prendo i kmer e li inserisco in CommonKmer\n",i);
        //Assegna l'indice della read corrente incrementato di 1 a commonKmerHead->readIndex e readIndex, utilizzati per tenere traccia dell'indice della read in commonKmerHead e in una variabile separata.
        //commonKmerHead contiene tutti i kmer comuni. readIndex serve per far vedere di essere arrivati alla kmer
        //i+1 nel calcolo dei common kmer
        commonKmerHead->readIndex = i + 1;
        readIndex = i + 1; //punta alla read successiva rispetto a quella che stiamo analizzando
        printf("\n readIndex di commonKmerHead: %ld\n", commonKmerHead->readIndex);

        //Controlla se l'indice della read è divisibile per il numero totale di thread (totalThreadNumber) utilizzando l'operatore modulo (%). Se non lo è, salta il resto del ciclo per quella specifica iterazione.
        if(readIndex % totalThreadNumber != threadIndex){ //serve per vedere se un thread puo analizzare quella read
            continue;
        }
        
        //Se l'indice della read è un multiplo di pi e ci sono più di un thread attivo (indicato da totalThreadNumber > 1)
        if(readIndex % pi == 0 && totalThreadNumber > 1){
            per = (float(readIndex)/readSetHead->readCount); //calcola la percentuale di completamento (per) rispetto al numero totale di read
            printf("Thread %ld starts to detect overlap!\n", threadIndex);
        }
        
        readLength = readSetHead->readSet[i].readLength; //lunghezza della read corrente 

        //Itera attraverso la read corrente con una finestra scorrevole di dimensione step, estraendo così i k-mer dalla read.
        //itera sui kmer della read i
        for(long int j = 0; j < readLength - kmerLength + 1;j = j + step){
            strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength); //Utilizza strncpy per copiare il k-mer dalla read corrente (readSetHead->readSet[i]) in kmer1.
            
            kmer1[kmerLength] = '\0'; //Aggiunge il carattere terminatore

            //se il k-mer ha caratteri tutti uguali, ritorna false e passa al k-mer successivo
            if(DetectSameKmer(kmer1, kmerLength) != true){ // Delete kmer that contains only one kind of base
                continue;
            }
            
            SetBitKmer(&kmerInteger, kmerLength, kmer1); //Converte il k-mer in un intero utilizzando SetBitKmer.
            hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger); //Cerca il k-mer nell'hash table dei k-mer per ottenere l'indice hash.
            
            if(hashIndex!=-1){ //Se l'indice hash è valido (hashIndex!=-1)

                // printf("\nkmer di riferimento: %s\n",kmer1);

                //inserisce il k-mer nella struttura dati commonKmerHead, specificando che il k-mer è stato trovato nella posizione j della read corrente.
                printf(" kmer (integer) da inserire in common kmer= %ld se la posizione in KHash è != -1\n", kmerHashTableHead->kmerHashNode[hashIndex].kmer);
        
                InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex, j, 1);//readIndex,posizione k-mer,orientamento
            
            }else{ //Se l'indice hash non è valido
                ReverseComplementKmer(kmer1, kmerLength); //cerca il complemento inverso del k-mer
                SetBitKmer(&kmerInteger, kmerLength, kmer1); //lo converte in un intero
                hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger); //lo cerca nell'hash table dei k-mer.

                // MODIFICA: Aggiunto controllo per verificare se il complemento inverso è presente nella hash table
                if(hashIndex != -1){ //Se viene trovato, inserisce il k-mer nella struttura dati commonKmerHead.
                    InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex, j, 0); //lo zero indica il reverse
                }
            }   
        }

        //qui ottiene la lista di commonKmer dell'i-esima read. Su questa lista chiama le altre funzioni
        
        /* STAMPA STRUTTURA COMMON KMER
        printf("stampa della common kmer della read %ld con realcount: %ld\n",i,commonKmerHead->realCount);
        for (int k = 0; k < commonKmerHead -> realCount; k++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[k].readIndex);
        }

        NOTA: IL REAL COUNT VIENE POSTO DI VOLTA IN VOLTA A 0, ANCHE SE NELLA STRUTTURA SONO PRESENTI DEI COMMON K-MER. QUINDI HO INSERITO UNA VARIABILE
        SIZE PER STAMPARE TUTTA LA STRUTTURA E NON FINO A REAL COUNT.

        printf("\n stampa struttura totale:\n");
        for (int i = 0; i < commonKmerHead -> size; i++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[i].readIndex);
        }
         */

        // printf("\n RemoveLowNumberKmer\n"); 

        //Dopo aver estratto tutti i k-mer dalla read corrente e inserito quelli validi nella struttura dati commonKmerHead, rimuove i k-mer con conteggio basso dalla struttura dati utilizzando RemoveLowNumberKmer.
        RemoveLowNumberKmer(commonKmerHead, forwardKmerCount, reverseKmerCount, readSetHead->readCount);
        
        /*printf("\n stampa commonKmerHead dopo RemoveLowNumber (realcount) \n");
        for (int i = 0; i < commonKmerHead -> realCount; i++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[i].readIndex);
        }*/

        //Riferimento algoritmo 1 nel paper: Sorting common k-mers in CKS;
        if(commonKmerHead->realCount > 100000){ //Se il numero di k-mer nella struttura dati commonKmerHead supera una soglia di 100.000
            //printf("\n heapSort \n");
            heapSort(commonKmerHead->commonKmer, 0, commonKmerHead->realCount - 1);
        }else{
            //printf("\n sort \n");
            sort(commonKmerHead->commonKmer, 0, commonKmerHead->realCount - 1);
        }
        
        //printf("\n RemoveMultipleSameKmer\n");
        RemoveMultipleSameKmer(commonKmerHead); //Rimuove i k-mer duplicati dalla struttura dati commonKmerHead.
        
        /*
        printf("\n stampa struttura dopo rimozione multimple:\n");
        for (int i = 0; i < commonKmerHead -> realCount; i++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[i].readIndex);
        }
        */

        /*printf("\n stampa struttura totale:\n");
        for (int i = 0; i < commonKmerHead -> size; i++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[i].readIndex);
        }*/

        //ci permette di ottenere i risultati delle sovrapposizioni tra le read utilizzando i k-mer comuni estratti dalla read corrente (riferendosi alle altre read, tramite la kmerReadNode)
        GetOverlapResult(G, commonKmerHead, readSetHead, localG, localCommonKmerHead, fp1);
        
        
        /*printf("Stampa di G\n");
    
        for (int i = 0; i < G ->allocationCountGraph; i++){
            if(G->graph[i].dataLeft != -1 || G->graph[i].dataRight != -1){
                printf("\n dataleft: %ld, dataright: %ld\n",G->graph[i].dataLeft,G->graph[i].dataRight);
            }
        }*/

        commonKmerHead->realCount = 0; //Resetta il conteggio dei k-mer nella struttura dati commonKmerHead per la successiva iterazione del ciclo.

        /*printf("\n riporto realcount a 0\n");
        printf("\n stampo la struttura per vedere che non è vuota\n");

        for (int i=0; i<commonKmerHead->size; i++){
            printf("\n commonKmer readIndex %ld\n",commonKmerHead->commonKmer[i].readIndex);
        }*/
    }

    /*
    printf("Stampa di G totale\n");
    
    for (int i = 0; i < G ->allocationCountGraph; i++){
        printf("\n dataleft: %ld, dataright: %ld\n",G->graph[i].dataLeft,G->graph[i].dataRight);
    }
    */
    //chiusura files
    fclose(fp);
    fflush(fp1);
    fclose(fp1);
}


//La funzione SubRemoveMultipleSameKmer si occupa di rimuovere i k-mer duplicati in un intervallo specificato
//da startIndex e endIndex all'interno della struttura CommonKmerHead. Questa funzione utilizza due variabili
//booleane (token e token1) per tenere traccia dei duplicati trovati e per marcare i k-mer duplicati con una
//posizione di inizio di -1, che indica che devono essere rimossi.
//La funzione, come descritta, è efficace per rimuovere k-mer duplicati e garantire che ogni k-mer sia unico
//entro l'intervallo specificato.
//If some of the same kmer positions of the two reads are close, remove part of the kmers to reduce memory
void SubRemoveMultipleSameKmer(CommonKmerHead * commonKmerHead, long int startIndex, long int endIndex){
	//La variabile booleana token viene utilizzata per tenere traccia se è stato trovato almeno un duplicato 
	//per il k-mer corrente (commonKmerHead->commonKmer[j]). Se token diventa vero, il k-mer corrente viene 
	//contrassegnato per la rimozione.
	bool token = false; //Indica se è stato trovato almeno un duplicato per il k-mer corrente.
	bool token1 = false; //Indica se è stato trovato un k-mer con differenze nelle posizioni di inizio e fine che soddisfano una determinata condizione.
	// Loop attraverso tutti i k-mer nell'intervallo specificato.
	for(long int j = startIndex; j < endIndex; j++){
		token = false;
		// Confronta il k-mer corrente con tutti i k-mer successivi nell'intervallo.
		for(long int m = j + 1; m <= endIndex; m++){
			// Salta i k-mer già contrassegnati per la rimozione (posizione di inizio impostata a -1).
			if(commonKmerHead->commonKmer[m].leftPosition == -1){ 
				continue;
			}

			//confronta le posizioni di inizio dei due k-mer. Se le posizioni di inizio coincidono		
			if(commonKmerHead->commonKmer[j].leftPosition == commonKmerHead->commonKmer[m].leftPosition){
				token = true; //viene impostata la variabile token su true
				//e il secondo k-mer (commonKmerHead->commonKmer[m]) viene contrassegnato per la rimozione 
				//impostando la sua posizione di inizio a -1.
				commonKmerHead->commonKmer[m].leftPosition = -1;
				continue;
			}
			
			//Controlla se le differenze nelle posizioni di inizio e fine dei due k-mer sono uguali e inferiori a 10.
			if(abs(commonKmerHead->commonKmer[j].leftPosition - commonKmerHead->commonKmer[m].leftPosition) == abs(commonKmerHead->commonKmer[j].rightPosition - commonKmerHead->commonKmer[m].rightPosition) 
			  && abs(commonKmerHead->commonKmer[j].leftPosition - commonKmerHead->commonKmer[m].leftPosition) < 10){
				token1 = true; //Imposta la variabile token1 su true.
				// Contrassegna il secondo k-mer per la rimozione.
				commonKmerHead->commonKmer[m].leftPosition = -1;
			}
		}
		// Se è stato trovato almeno un duplicato per il k-mer corrente, 
		//contrassegna anche il k-mer corrente per la rimozione.
		if(token != false){
			commonKmerHead->commonKmer[j].leftPosition = -1;
		}
	}
}

//La funzione RemoveMultipleSameKmer ha il compito di rimuovere k-mer duplicati dalla struttura CommonKmerHead.
//La funzione utilizza una funzione di supporto SubRemoveMultipleSameKmer per effettuare la rimozione all'interno
//di intervalli specifici, quindi compatta la struttura rimuovendo gli elementi marcati per la rimozione.
void RemoveMultipleSameKmer(CommonKmerHead * commonKmerHead){
	long int startIndex = 0; //l'indice di inizio dell'intervallo corrente di k-mer con lo stesso readIndex.
	//readIndex: l'indice di lettura corrente, inizialmente impostato sull'indice di lettura del primo k-mer.
	long int readIndex = commonKmerHead->commonKmer[0].readIndex;
	//endIndex: l'indice di fine dell'intervallo corrente di k-mer con lo stesso readIndex.
	long int endIndex = -1;
	 // Scansiona tutti i k-mer in commonKmerHead per determinare gli intervalli di duplicati
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		
		//Quando viene trovato un k-mer con un readIndex diverso da quello corrente, viene considerato come fine di
		//un intervallo di k-mer duplicati.
		if(commonKmerHead->commonKmer[i].readIndex != readIndex){//startIndex=0 i = 2 endIndex=2-1
			endIndex = i - 1;
			 // Rimuove i duplicati nell'ultimo intervallo
			SubRemoveMultipleSameKmer(commonKmerHead, startIndex, endIndex); //considera l'intervallo 0
			
			//L'intervallo successivo inizia dall'indice corrente, e readIndex viene aggiornato.
			startIndex = i;
			readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}
	endIndex = commonKmerHead->realCount - 1;
	
	//Dopo il ciclo, viene chiamata SubRemoveMultipleSameKmer sull'ultimo intervallo di k-mer.
	SubRemoveMultipleSameKmer(commonKmerHead, startIndex, endIndex);
	
	long int shiftCount = 0;
	
	// Compatta la struttura eliminando gli elementi marcati per la rimozione
	//Un altro ciclo for scorre nuovamente tutti i k-mer.
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].leftPosition == -1){
			//Se un k-mer è stato marcato per la rimozione 
			//(indicato da leftPosition impostato a -1), shiftCount viene incrementato.
			shiftCount++;
		}else if(shiftCount != 0){
			//Se shiftCount non è zero, i k-mer vengono spostati in avanti per riempire gli spazi vuoti lasciati 
			//dai k-mer rimossi.
			commonKmerHead->commonKmer[i - shiftCount].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
			commonKmerHead->commonKmer[i - shiftCount].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
			commonKmerHead->commonKmer[i - shiftCount].orientation = commonKmerHead->commonKmer[i].orientation;
			commonKmerHead->commonKmer[i - shiftCount].readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}

	// Aggiorna il conteggio reale dei k-mer. Viene ridotto di shiftCount, che rappresenta il numero di k-mer rimossi.
	commonKmerHead->realCount = commonKmerHead->realCount - shiftCount;
}

//Questa funzione, RemoveLowNumberKmer, esclude la possibilità di sovrapposizione se il numero di k-mer uguali
//per una coppia di letture è relativamente basso.
//If the number of the same kmer for a pair of reads is relatively small, the possibility of overlap is excluded
void RemoveLowNumberKmer(CommonKmerHead * commonKmerHead, long int * forwardKmerCount, long int * reverseKmerCount, long int readCount){
	
	// Inizializza i contatori di k-mer per le letture.
	for(long int i = 0; i < readCount; i++){
		forwardKmerCount[i] = 0;
		reverseKmerCount[i] = 0;
	}
	// Conta i k-mer per ogni lettura e orientazione.
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].orientation == 0){
			forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]++;
		}else{
			reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]++;
		}
	}

	// Esclude i k-mer con conteggio basso.
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] < 15 && commonKmerHead->commonKmer[i].orientation == 0){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		if(forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= 15 && commonKmerHead->commonKmer[i].orientation == 0 && forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] <= reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
		if(reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] < 15 && commonKmerHead->commonKmer[i].orientation == 1){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
		if(reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= 15 && commonKmerHead->commonKmer[i].orientation == 1 && forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
	}
	
	long int shiftCount = 0;
	
	// Rimuove i k-mer esclusi e aggiorna il conteggio.
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].leftPosition == -1){
			shiftCount++;
		}else if(shiftCount != 0){
			commonKmerHead->commonKmer[i - shiftCount].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
			commonKmerHead->commonKmer[i - shiftCount].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
			commonKmerHead->commonKmer[i - shiftCount].orientation = commonKmerHead->commonKmer[i].orientation;
			commonKmerHead->commonKmer[i - shiftCount].readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}
	// Aggiorna il numero reale di k-mer dopo la rimozione.	
	commonKmerHead->realCount = commonKmerHead->realCount - shiftCount;
}


void sortGraph(AdjGraph * graph, long int left, long int right){
	if(left >= right){
        return ;
    }
    long int i = left;
    long int j = right;
    long int key = graph[left].dataLeft;
	long int dataRight = graph[left].dataRight;
	
    while(i < j){
        while(i < j && key <= graph[j].dataLeft){
            j--;
        }
		
		if(i < j){ // Scambia gli elementi a[i] e a[j]
			graph[i].dataLeft = graph[j].dataLeft;
			graph[i].dataRight = graph[j].dataRight;
			i++;
		}

         
        while(i < j && key > graph[i].dataLeft){
            i++;
        }
		
		if(i < j){ // Scambia gli elementi a[i] e a[j]
			graph[j].dataLeft = graph[i].dataLeft;
			graph[j].dataRight = graph[i].dataRight;
			j--;
		}
        
    }
    
    graph[i].dataLeft = key;
	graph[i].dataRight = dataRight;
	
    sortGraph(graph, left, i - 1);
    sortGraph(graph, i + 1, right);
}

void swapCommonKmer(CommonKmer *a, long int left, long int right){
	// Scambia i valori dell'elemento readIndex
	long int temp = a[left].readIndex;
 	a[left].readIndex = a[right].readIndex;
	a[right].readIndex = temp;
	
	// Scambia i valori dell'elemento leftPosition
	temp = a[left].leftPosition;
 	a[left].leftPosition = a[right].leftPosition;
	a[right].leftPosition = temp;
	
	// Scambia i valori dell'elemento rightPosition
	temp = a[left].rightPosition;
 	a[left].rightPosition = a[right].rightPosition;
	a[right].rightPosition = temp;
	
	// Scambia i valori dell'elemento orientation
	bool temp1 = a[left].orientation;
 	a[left].orientation = a[right].orientation;
	a[right].orientation = temp1;
}

//La funzione downToMaxHeap implementa l'operazione "down-heap" su un heap massimo, partendo da un nodo padre
//specificato e continuando fino alla fine dell'heap. 
void downToMaxHeap(CommonKmer *a, long int bgn, long int end){
    long int child;
    long int parent = bgn;
	// Continua finché il nodo padre ha almeno un figlio
    while ((child = parent * 2 + 1) < end)
    {
		// Trova il figlio più grande
        if ((child < end - 1) && ((a[child].readIndex < a[child + 1].readIndex) || (a[child].readIndex == a[child + 1].readIndex && a[child].leftPosition < a[child + 1].leftPosition)))
            ++child;   
		// Se il figlio più grande è maggiore del padre, scambia il padre con il figlio
        if ((a[child].readIndex > a[parent].readIndex) || (a[child].readIndex == a[parent].readIndex && a[child].leftPosition > a[parent].leftPosition))
            swapCommonKmer(a, child, parent);
        else
            break;
        parent = child;  // Passa al figlio appena scambiato e continua il processo
    }
}

void buildMaxHeap(CommonKmer * a, long int bgn, long int end){
    if (bgn >= end - 1)
        return;

    int parent = end / 2 - 1;  // Calcola l'indice del genitore dell'ultimo elemento nell'array
    while (parent >= 0) // Scorri dall'indice del genitore fino all'inizio dell'array
    {
		 // Ripristina la proprietà di heap massimo partendo dal genitore corrente
        downToMaxHeap(a, parent, end);
        --parent; // Passa al genitore precedente
    }
}

void heapSort(CommonKmer *a, long int bgn, long int end){

    buildMaxHeap(a, bgn, end); // Costruzione dell'heap massimo

    while (end > 1)
    {
        swapCommonKmer(a, 0, --end); //Scambia il primo elemento (massimo) con l'ultimo
        downToMaxHeap(a, 0, end); //Ripristina la proprietà di heap massimo nell'array ridotto
    }
}

//La funzione sort ordina gli elementi di un array CommonKmer in base al campo readIndex, e in caso di parità,
//ordina per leftPosition in ordine ascendente.
//Sort the kmer positions on the same read in ascending order
void sort(CommonKmer *a, long int left, long int right)
{

	if(left >= right){
        return ;
    }

    long int i = left;
    long int j = right;
    long int key = a[left].readIndex;
	long int leftPosition = a[left].leftPosition;
	long int rightPosition = a[left].rightPosition;
    bool orientation = a[left].orientation;
	
    while(i < j){
        while(i < j && (key < a[j].readIndex || (key == a[j].readIndex && leftPosition <= a[j].leftPosition))){
            j--;
        }
		
		if(i < j){ // Scambia gli elementi a[i] e a[j]
			a[i].readIndex = a[j].readIndex;
			a[i].leftPosition = a[j].leftPosition;
			a[i].rightPosition = a[j].rightPosition;
			a[i].orientation = a[j].orientation;
			i++;
		}
      
        while(i < j && (key > a[i].readIndex || (key == a[i].readIndex && leftPosition > a[i].leftPosition))){
            i++;
        }
		
		if(i < j){ // Scambia gli elementi a[i] e a[j]
			a[j].readIndex = a[i].readIndex;
			a[j].leftPosition = a[i].leftPosition;
			a[j].rightPosition = a[i].rightPosition;
			a[j].orientation = a[i].orientation;
			j--;
		}
    }

	// Posiziona l'elemento chiave nella sua posizione corretta nell'array ordinato
    a[i].readIndex = key;
	a[i].leftPosition = leftPosition;
	a[i].rightPosition = rightPosition;
	a[i].orientation = orientation;

	// Ordina ricorsivamente le due parti dell'array
    sort(a, left, i - 1);
    sort(a, i + 1, right); 
}


void DestroyGraph(AdjGraphHead * G){
	G->realCountGraph = 0;
	G->realCountArc = 0;
}

//determina se l'intervallo tra i piccoli k-mer soddisfa la condizione di consistenza. 
//Se soddisfa, restituisce 1; altrimenti, restituisce 0.
//In sintesi, questa funzione verifica se le posizioni di inizio e di fine delle letture di sinistra e destra 
//sono entro una certa distanza fissa (300) dalle rispettive estremità delle letture. 
//Se una di queste posizioni supera la distanza consentita, la funzione restituisce 0, altrimenti restituisce 1.
//According to the graph of the local area, determine whether the interval between the small kmers meets the consistency condition. 
//If it meets, it returns 1; if it does not, it returns 0.
long int Overlap_DisplayLocalRegion(AdjGraphHead * G,long int leftLen,long int rightLen){
	if(G->realCountArc <= 0){
		// Se non ci sono archi reali, restituisce 0.
		printf("assenza di archi reali\n");
		return 0;
	}
    // Inizializza le posizioni di inizio e fine per le letture di sinistra e destra
	long int leftStartpos=-1,leftEndpos=-1,rightStartpos=-1,rightEndpos=-1;
	
	// Ottiene le posizioni di inizio e fine dal grafico degli archi.
	leftStartpos = G->graph[G->arcIndex[0].startIndex].dataLeft;
	leftEndpos = G->graph[G->arcIndex[G->realCountArc - 1].endIndex].dataLeft;
	rightStartpos = G->graph[G->arcIndex[0].startIndex].dataRight;
	rightEndpos = G->graph[G->arcIndex[G->realCountArc - 1].endIndex].dataRight;

	// Calcola la distanza tra le lunghezze delle letture di sinistra e destra.
	long int distance = abs(leftLen - rightLen);
	
	distance = 300;  // Imposta una distanza fissa di 300.

	// Controlla se le posizioni di inizio o le posizioni di fine sono oltre la distanza consentita.
	if((leftStartpos > distance || rightStartpos > distance) || (leftLen - leftEndpos > distance || rightLen - rightEndpos > distance)){
		// Se le posizioni di inizio o di fine superano la distanza consentita, restituisce 0.
		return 0;
	}else{
		// Se le posizioni di inizio e di fine sono entro la distanza consentita, restituisce 1.
		return 1;
	}
}

//La funzione Overlap_Display_Graph determina il risultato finale dell'overlap tra due sequenze di lettura 
//(leftRead e rightRead) basandosi sul grafico del k-mer locale. Restituisce 1 se l'overlap soddisfa le 
//condizioni di consistenza, altrimenti restituisce 0. 
//Secondo il percorso del grafo (l'array di grafi), ottiene il risultato finale del rilevamento dell'overlap.
//According to the graph path, get the final overlap detection result
//Riferimento Algoritmo 6 del paper 
long int Overlap_Display_Graph(AdjGraphHead * G, long int leftIndex,long int rightIndex,bool orien,long int leftLen,long int rightLen,FILE * fp, long int leftStartpos, long int leftEndpos, long int rightStartpos, long int rightEndpos, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead){
	// Se gli indici sono uguali, non c'è overlap da verificare.
	//leftIndex = commonKmerHead->readIndex; (indice della read in cui si trova il common kmer()
	//rightIndex = readIndex; (readIndex = commonKmerHead->commonKmer[0].readIndex;)
	if(leftIndex == rightIndex){
		return 0;
	}
	long int kmerLength = G->kmerLength;
	
	long int i;
	long int t = 0;

	long int overlenleft,overlenright; // Variabili per calcolare le lunghezze dell'overlap.
	long int MinOverLen,MaxOverLen;

	// Inizializza le variabili temporanee per le posizioni di start e end.
	//leftStartpos = graph[iniStartIndex].dataLeft (la read che stiamo considerando R1)
	long int tempLeftStart = leftStartpos;
	//leftEndpos = graph[lastEndIndex].dataLeft
	long int tempLeftEnd = leftEndpos;
	//rightStartpos = graph[iniStartIndex].dataRight (l'altra read R2)
	long int tempRightStart = rightStartpos;
	//rightEndpos = graph[lastEndIndex].dataRight
	long int tempRightEnd = rightEndpos;

	printf("TLS = %ld, TLE = %ld, TRS = %ld, TRE = %ld\n", tempLeftStart, tempLeftEnd, tempRightStart, tempRightEnd);
	
	// Calcola il minimo e il massimo overlap tra le letture di sinistra e destra.
	//Riferimento algoritmo 6
	//MinOverlapLen: Minimo di OverlapLenR1 e OverlapLenR2
	MinOverLen = min(leftEndpos -leftStartpos,rightEndpos - rightStartpos);
	//MaxOverlapLen: Massimo di OverlapLenR1 e OverlapLenR2
	MaxOverLen = max(leftEndpos -leftStartpos,rightEndpos - rightStartpos);
	// La lunghezza minima dell'allineamento deve essere maggiore di 300.
	//questo controllo non è scritto nel paper
	//The minimum alignment length is greater than 300
	if(MinOverLen < 300){return 0;} 
	// Controlla il rapporto di lunghezza per l'overlap.
	//Verifica che il rapporto tra la differenza delle lunghezze massime e minime e la lunghezza massima non superi un certo valore (G->lengthRatio).
	//lengthRatio = 0.3
	//riferimento all'algoritmo 6 nel paper (condizione 3 da rispettare)
	if((float)(MaxOverLen - MinOverLen)/MaxOverLen > G->lengthRatio ){return 0;} 

	// Variabili per le distanze massime e gli intervalli.
	t=G->largestIntervalDistance; //t=400 (alfa)
	long int maxIntervalDistance = localG->largestIntervalDistance; //1500 (beta)
	
	// Variabili per le posizioni locali di start e end.
	long int localLeftStart = 0;
	long int localLeftEnd = 0;
	long int localRightStart = 0;
	long int localRightEnd = 0;

	long int dd = 0;
	
	long int alignLen = 30;
	float alignRtio = 0.6;
	/* Regardless of whether the two readings R1 and R2 are forward aligned and reverse aligned, there are four alignment situations.
	The head of R1 is aligned to the tail of R2, the tail of R1 is aligned to the head of R2, R1 contains R2 and R2 contains R1.
	*/ //-->spiegazione sotto
	//Indipendentemente dal fatto che le due letture R1 e R2 siano allineate in avanti o all'indietro, esistono quattro situazioni di allineamento.
	//La testa di R1 è allineata alla coda di R2, la coda di R1 è allineata alla testa di R2, R1 contiene R2 e R2 contiene R1.
	//Case of positive alignment;
	// Caso di allineamento positivo.
	if(orien==0){
		// Controlla se le posizioni di start sono entro i limiti accettabili.
		//The length of the aligned position area on the two kmers should not be too large

	//leftStartpos = graph[iniStartIndex].dataLeft (la read che stiamo considerando R1)
	//leftEndpos = graph[lastEndIndex].dataLeft
	//rightStartpos = graph[iniStartIndex].dataRight (l'altra read R2)
	//rightEndpos = graph[lastEndIndex].dataRight, 

		//nelle seguenti righe di codice considera i 4 possibili casi di allineamento spiegati nel paper
		//riferimento algoritmo 6 pag. 5 condizioni (1) (2) (3) (4)
		if((rightStartpos > leftStartpos && leftStartpos > maxIntervalDistance) || (rightStartpos <= leftStartpos && rightStartpos > maxIntervalDistance)){ 
			return 0;
		}
		if((leftLen-leftEndpos > rightLen-rightEndpos && rightLen-rightEndpos > maxIntervalDistance) || (leftLen-leftEndpos <= rightLen-rightEndpos && leftLen-leftEndpos > maxIntervalDistance)){ 
			return 0;
		}
		if(rightStartpos > leftStartpos){
			if(leftStartpos > maxIntervalDistance){
				return 0;
			
			//riferimento algoritmo 6 condizione (1) da rispettare
			//se è maggiore di 400 (t) la condizione non è rispettata e analizza gli shorter kmer
			}else if(leftStartpos > t){ //t=400
				localLeftStart = 0; // Calcola le posizioni locali di start e end.
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightStartpos - leftStartpos;
				localRightEnd = rightStartpos - 1;
				// Verifica la consistenza del k-mer più corto.
				//If the distance between two kmers is less than maxIntervalDistance, but greater than t, then use a shorter kmer to further determine whether the two kmers are consistent kmer
				//t = 400 e maxIntervalDistance = 1500
				//smallKmerLength = 9
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				//alignRtio = 0.6
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightStartpos = rightStartpos - leftStartpos;
			leftStartpos = 1;
		}else{
			
			if(rightStartpos > maxIntervalDistance){
				return 0;
			//riferimento algoritmo 6 condizione (1) da rispettare
			//se è maggiore di 400 (t) la condizione non è rispettata e analizza gli shorter kmer
			}else if(rightStartpos > t){
				localLeftStart = leftStartpos - rightStartpos;
				localLeftEnd = leftStartpos - 1;
				localRightStart = 0;
				localRightEnd = rightStartpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftStartpos = leftStartpos-rightStartpos;
			rightStartpos = 1;
		}
		
		// Calcola le nuove posizioni di end.
		if(leftLen-leftEndpos > rightLen-rightEndpos){
			if(rightLen-rightEndpos > maxIntervalDistance){
				return 0;
			//riferimento algoritmo 6 condizione (1) da rispettare
			//se è maggiore di 400 (t) la condizione non è rispettata e analizza gli shorter kmer
			}else if(rightLen-rightEndpos > t){
				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftEndpos + rightLen - rightEndpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightLen - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftEndpos = rightLen - rightEndpos + leftEndpos;
			rightEndpos = rightLen;
			
		}else{
			if(leftLen-leftEndpos > maxIntervalDistance){
				return 0;
			//riferimento algoritmo 6 condizione (1) da rispettare
			//se è maggiore di 400 (t) la condizione non è rispettata e analizza gli shorter kmer
			}else if(leftLen-leftEndpos > t){
				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftLen - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightEndpos + leftLen - leftEndpos - 1;
				

				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightEndpos = leftLen - leftEndpos + rightEndpos;
			leftEndpos = leftLen;
		}
		
	}
//Reverse alignment
// Caso di allineamento inverso.
	if(orien==1){
		
		if((rightLen-rightEndpos > leftStartpos && leftStartpos > maxIntervalDistance) || (rightLen-rightEndpos <= leftStartpos && rightLen-rightEndpos > maxIntervalDistance)){
			return 0;
		}
		if((leftLen-leftEndpos > rightStartpos && rightStartpos > maxIntervalDistance) || (leftLen-leftEndpos <= rightStartpos && leftLen-leftEndpos > maxIntervalDistance)){
			return 0;
		}
		
		if(rightLen-rightEndpos > leftStartpos){
			if(leftStartpos > maxIntervalDistance){
				return 0;
			}else if(leftStartpos > t){

				localLeftStart = 0;
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightEndpos + kmerLength + leftStartpos;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightEndpos = rightEndpos + leftStartpos;
			leftStartpos = 1;
		}else{
			if(rightLen-rightEndpos > maxIntervalDistance){
				return 0;
			}else if(rightLen-rightEndpos > t){

				localLeftStart = leftStartpos - rightLen + rightEndpos + kmerLength;
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightLen - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftStartpos = leftStartpos - rightLen + rightEndpos;
			rightEndpos = rightLen;
		}
		
		if(leftLen-leftEndpos > rightStartpos){
			if(rightStartpos > maxIntervalDistance){ 
				return 0;
			}else if(rightStartpos > t){

				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftEndpos + kmerLength + rightStartpos - 1;
				localRightStart = 0;
				localRightEnd = rightStartpos - 1;

				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftEndpos = leftEndpos + rightStartpos;
			rightStartpos = 1;
		}else{
			if(leftLen-leftEndpos > maxIntervalDistance){
				return 0;
			}else if(leftLen-leftEndpos > t){

				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftLen - 1;
				localRightStart = rightStartpos - leftLen + leftEndpos + kmerLength;
				localRightEnd = rightStartpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);

				if(ss < alignRtio){
					return 0;
				}
				
			}
			
			rightStartpos = rightStartpos - leftLen + leftEndpos;
			leftEndpos = leftLen;
			
		}

	}

	// Calcola le nuove lunghezze dell'overlap.
	overlenleft = leftEndpos - leftStartpos;
	overlenright = rightEndpos - rightStartpos;
	
	long int tempMaxLen = MaxOverLen;
	long int tempMinLen = MinOverLen;
	
	MinOverLen = min(overlenleft,overlenright);
	MaxOverLen = max(overlenleft,overlenright);

	// Verifica se l'overlap soddisfa i criteri di lunghezza e rapporto.
	//riferimento algoritmo 6
	//overlapLengthCutOff is the minimum overlap length between two long reads (default 500)
	//sarebbe la epsilon nel paper
	//condizione 2 dell'algoritmo 6: MinOverlapLen > ε;
	if(MaxOverLen<G->overlapLengthCutOff){ return 0;}
	//G->lengthRatio = gamma (γ) = 0.3
	//riferimento algoritmo 6 condizione 3
	//(3) (MaxOverlapLen − MinOverlapLen) / MaxOverlapLen < γ .
	if((float)(MaxOverLen - MinOverLen)/MaxOverLen > G->lengthRatio ){ return 0;} 
	// Scrive i risultati dell'allineamento sul file.
	//Output alignment results to file
	if(leftIndex > rightIndex){	
		fprintf(fp,"%ld,%ld,%d,%ld,%ld,%ld,%ld,%ld,%ld\n",leftIndex,rightIndex,orien,
												leftStartpos,leftEndpos,leftLen,rightStartpos,rightEndpos,rightLen);
	}else{
		fprintf(fp,"%ld,%ld,%d,%ld,%ld,%ld,%ld,%ld,%ld\n",rightIndex,leftIndex,orien,
												rightStartpos,rightEndpos,rightLen,leftStartpos,leftEndpos,leftLen);
	}

	return 1;

}

/*Questa funzione aumenta dinamicamente la capacità di memorizzazione del grafo (diretto o inverso) 
moltiplicandola per 1,5 e copia i dati esistenti nella nuova area di memoria allocata.
*/
void ReAllocateAdjGraph(AdjGraphHead * G, long int a){
	
	if(a == 0){
		G->allocationCountGraph = G->allocationCountGraph*1.5;
		AdjGraph * graph = (AdjGraph *)malloc(sizeof(AdjGraph)*G->allocationCountGraph);
		for(long int i = 0; i < G->realCountGraph; i++){
			graph[i].dataLeft = G->graph[i].dataLeft;
			graph[i].dataRight = G->graph[i].dataRight;

			graph[i].visit = G->graph[i].visit;
		}
		free(G->graph);
		G->graph = graph;
	}else{
		G->reverseAllocationCountGraph = G->reverseAllocationCountGraph*1.5;
		AdjGraph * graph = (AdjGraph *)malloc(sizeof(AdjGraph)*G->reverseAllocationCountGraph);
		for(long int i = 0; i < G->reverseRealCountGraph; i++){
			graph[i].dataLeft = G->reverseGraph[i].dataLeft;
			graph[i].dataRight = G->reverseGraph[i].dataRight;

			graph[i].visit = G->reverseGraph[i].visit;
		}
		free(G->reverseGraph);
		G->reverseGraph = graph;
	}
		
}


void ReAllocateArcIndex(AdjGraphHead * G){
	
	G->allocationCountArc = G->allocationCountArc*2;
	
	ArcIndex * arcIndex = (ArcIndex *)malloc(sizeof(ArcIndex)*G->allocationCountArc);
	
	for(long int i = 0; i < G->realCountArc; i++){
		arcIndex[i].startIndex = G->arcIndex[i].startIndex;
		arcIndex[i].endIndex = G->arcIndex[i].endIndex;
	}
	
	free(G->arcIndex);
	
	G->arcIndex = arcIndex;
}

/*
La funzione AddEdgeInGraph aggiunge archi al grafo basandosi sui k-mer comuni. L'aggiunta degli archi dipende
da diversi fattori, inclusi gli indici dei nodi nel grafo, le lunghezze dei k-mer, le distanze tra i k-mer, 
e altri parametri specifici del grafo. Aggiunge archi ai kmer comuni consistenti
Parametri della funzione
*AdjGraphHead G: Puntatore alla struttura del grafo di adiacenza.
bool orientation: Indica se si sta lavorando sul grafo diretto (0) o inverso (1).
long int leftIndex, rightIndex: Indici del nodo sinistro e destro del grafo.
long int largestIntervalDistance, maxIntervalDistance: Distanze massime consentite tra i nodi (caso peggiore e vaso normale).
*AdjGraphHead localG: Puntatore alla struttura del grafo locale.
*CommonKmerHead localCommonKmerHead: Puntatore alla struttura dei k-mer comuni locali.
**char leftRead, rightRead: Puntatori ai read sinistro e destro.
*/
//With the common kmer as the vertex, add edges to the common consistent kmer.
long int AddEdgeInGraph(AdjGraphHead * G, bool orientation, long int leftIndex , long int rightIndex, long int largestIntervalDistance, long int maxIntervalDistance, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead){
	//rightIndex = endIndex
	//leftIndex = startIndex
	AdjGraph * graph = NULL; //Viene inizializzato un puntatore al grafo su cui lavorare (graph).
	
	//Se l'orientamento è 0, si lavora sul grafo diretto; altrimenti, si lavora sul grafo inverso.
	if(orientation == 0){
		//se il commonKmer in un'altra read (non quella corrente) è stato già visitato
		//rightIndex = endIndex
		if(G->graph[rightIndex].visit == true){
			//Se il nodo destro nel grafo è già stato visitato, la funzione restituisce 0 indicando che l'aggiunta dell'arco non è possibile.
			return 0;
		}
		graph = G->graph;
	}else{
		//Se l'orientamento è 1, si lavora sul grafo inverso;
		if(G->reverseGraph[rightIndex].visit == true){
			return 0;
		}
		graph = G->reverseGraph;
	}
	
	long int kmerLength = 13; //non viene usato

	//Riferimento algoritmo 4 nel paper: Determine_consistent_1 (CKS[start],CKS[end], k, ks , α,γ )	
	//Getting the positions of the start and end common k-mers (commento mio)
	//m1>n1
	//leftIndex = startIndex
	//rightIndex = endIndex
	//n1 è la posizione di partenza del kmer della long read1
	long int n1 = graph[leftIndex].dataLeft;
	//n2 è la posizione finale del kmer della long read1
	long int n2 = graph[leftIndex].dataRight;
	//m1 è la posizione di partenza del kmer della long read2
	long int m1 = graph[rightIndex].dataLeft;
	//m2 è la posizione finale del kmer della long read2
	long int m2 = graph[rightIndex].dataRight;
					
	//Viene calcolata la distanza tra i k-mer sinistro e destro (distanza tra due common kmers)
	//maxIntervalDistance = 1500
	//For two common k-mers (n1, On1, n2, On2) and (m1, Om1, m2, Om2 ), two distances D1 = |m1 - n1| 
	//and D2 = |m2 − n2| can be calculated. (commento mio)
	//qui viene usato gia beta (1500) e non alfa (400) come specificato nel paper algoritmo 4
	//condizione c3 nel paperà
	//probabilmente questa funzione viene chiamata anche dai kmer più piccoli
	if(abs(m1-n1) > maxIntervalDistance || abs(m2-n2) > maxIntervalDistance){
		//Se la distanza tra i k-mer supera la distanza massima consentita, la funzione restituisce 0.
		return 0; //ritorna false
	}
	
	//Vengono calcolati i valori minimi e massimi tra le distanze dei k-mer.
	//servono per la condizione c4 algoritmo 4
	long int minvalue = min(abs(m1-n1), abs(m2-n2));
	long int maxvalue = max(abs(m1-n1), abs(m2-n2));
	
	long int temp = 0; //Viene impostato il valore di temp a 0. Questo valore verrà utilizzato per determinare se l'aggiunta dell'arco è avvenuta con successo.
	float ss1 = 0; //ss1 potrebbe rappresentare un punteggio di similarità tra due segmenti di sequenze. Se il punteggio è 
		//superiore a una certa soglia (ad esempio, 0.5), allora le sequenze sono considerate sufficientemente
		//simili per aggiungere un arco tra i nodi.
		//ss1 non viene mai modificato quindi è sempre 0
	long int alignLength = 30;
	
	//In base all'orientamento dei k-mer e alla loro lunghezza, viene verificato se è possibile aggiungere l'arco al grafo.
	//G->lengthRatio = 0.3 (sarebbe gamma (γ) nel paper)
	//riferimento algoritmo 4: se l'orientamento è 0
	//&& float(maxvalue - minvalue)/maxvalue < G->lengthRatio  (condizione C4 nel paper)
	//&& ((n1 < m1 && n2 < m2) || (n1>m1 && n2>m2)) (condizione C1 nel paper)
	//orientation = 0 --> forward
	//se tutte queste condizioni sono soddisfatte il kmer è consistente con quello preso in considerazione
	//riferimento algoritmo 4 (Determine_Consistent_1)
	if(orientation == 0 && float(maxvalue - minvalue)/maxvalue < G->lengthRatio && ((n1 < m1 && n2 < m2) || (n1>m1 && n2>m2))){
		
		//largestIntervalDistance = 400 (come nel paper) alfa (α)
		//riferimento all'algoritmo 4 nel paper (condizione C3 nel paper)
		if(abs(m1 - n1) <= largestIntervalDistance && abs(m2 - n2) <= largestIntervalDistance){
			//Se le condizioni sono soddisfatte, viene aggiunto un arco tra i nodi sinistro e destro con un 
			//peso calcolato. Viene aggiornata la struttura del grafo con il nuovo arco.
			temp = 1;
		
		//altrimenti se anche una delle due condizioni non è soddisfatta significa che non è possibile determinare
		//la consistenza nella prima fase, viene quindi utilizzata la seconda fase, che analizza ulteriormente la 
		//consistenza basandosi su k-mer più piccoli (ks < k).
		//verifica se abs (m1 - n1) <= maxIntervalDistance (1500) 
		//maxIntervalDistance sarebbe beta (β) che fornisce più k-mer comuni candidati
		//Riferimento algoritmo 5: Determine_consistent _2 (CKS[start],CKS[end], k, ks , β)
		}else if(abs(m1 - n1) <= maxIntervalDistance && abs(m2 - n2) <= maxIntervalDistance){
			//a riga 1303 fa < 0.3; a riga 1313 fa < 0.1
			//The values of these parameters will need to be changed according to the sequencing error rate of
			//TGS. If the error rate is low, reducing the values of these parameters could improve the
			//precision of the result. (commento mio)
			//rifà questo calcolo con un valore più piccolo (0.1 invece di 0.3)
			if(float(maxvalue - minvalue)/maxvalue < 0.1){ 
				temp = 1;
			//se la posizione di partenza del kmer della read1 è < della posizione finale del kmer della read1 e la 
			//posizione di partenza del kmer della read2 è < della posizione finale del kmer della read2
			//significa che i kmer sono entrambi in forward
			}else if(n1 < m1 && n2 < m2){ //condizione C1 nel paper
				//Se le condizioni sono soddisfatte, viene chiamata la funzione GetCommonShorterKmer per ottenere i k-mer comuni più corti.
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, n1, m1, n2, m2, G->smallKmerLength, 0);
			//significa che i kmer sono entrambi in forward
			}else if(n1>m1 && n2>m2){ //condizione C1 nel paper)
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, m1, n1, m2, n2, G->smallKmerLength, 0);
			}
		}
		//L'orientamento del grafo e il confronto delle distanze tra i k-mer sono fondamentali per determinare se aggiungere o meno un arco al grafo.
		//ss1 potrebbe rappresentare un punteggio di similarità tra due segmenti di sequenze. Se il punteggio è 
		//superiore a una certa soglia (ad esempio, 0.5), allora le sequenze sono considerate sufficientemente simili per aggiungere un arco tra i nodi.

		//if two common k-mers are consistent, the two distances between them should be similar (dal paper)
		//Nota: ss1 è sempre zero, non lo modifica
		if(temp == 1 || ss1 > 0.5){
			//Se l'aggiunta dell'arco è avvenuta con successo, viene memorizzato l'arco nel grafo (G->arcIndex) insieme al peso dell'arco calcolato.
			G->arcIndex[G->realCountArc].startIndex = leftIndex;
			G->arcIndex[G->realCountArc].endIndex = rightIndex;
			G->arcIndex[G->realCountArc].weight = ((((float)minvalue/maxvalue) + (1 - (float)maxvalue/maxIntervalDistance))/2) * (maxvalue); //no riferimenti nel paper
			G->realCountArc++; //Viene incrementato il contatore degli archi (G->realCountArc).
			if(G->realCountArc >= G->allocationCountArc){
				//Se il numero di archi nel grafo supera la capacità allocata, viene eseguita la riallocazione della memoria per gli archi.
				ReAllocateArcIndex(G);
			}
			//Se l'aggiunta dell'arco è avvenuta con successo, la funzione restituisce 1.
			return 1;
		}else{
			//Se l'aggiunta dell'arco non è avvenuta con successo, la funzione restituisce 2.
			return 2;
		}
		
		
	}
	
	//orientamento = 1 --> reverse
	//riferimento algoritmo 4, condizione c2 nell'algoritmo 4
	if(orientation == 1 && float(maxvalue - minvalue)/maxvalue < G->lengthRatio && ((n1<m1 && n2>m2) || (n1>m1 && n2<m2))){
		
		if(abs(m1 - n1) <= largestIntervalDistance && abs(m2 - n2) <= largestIntervalDistance){
			temp = 1;
		}else if(abs(m1 - n1) <= maxIntervalDistance && abs(m2 - n2) <= maxIntervalDistance){
			if(float(maxvalue - minvalue)/maxvalue < 0.1){
				temp = 1;
			}else if(n1<m1 && n2>m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, n1, m1, m2, n2, G->smallKmerLength, 1);
			}else if(n1>m1 && n2<m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, m1, n1, n2, m2, G->smallKmerLength, 1);
			}
		}
		//ss1 è sempre zero, non lo modifica
		if(temp == 1 || ss1 > 0.5){
			G->arcIndex[G->realCountArc].startIndex = leftIndex;
			G->arcIndex[G->realCountArc].endIndex = rightIndex;
			G->arcIndex[G->realCountArc].weight = ((((float)minvalue/maxvalue) + (1 - (float)maxvalue/maxIntervalDistance))/2) * (maxvalue);
			G->realCountArc++;
			if(G->realCountArc >= G->allocationCountArc){
				ReAllocateArcIndex(G);
			}
			return 1;
		}else{
			return 2;
		}
		
	}
	
	return 0; //Se non è possibile aggiungere l'arco, la funzione restituisce 0.
}

/*La funzione CreatGraphSinglePath è responsabile di creare un singolo percorso nel grafo di adiacenza, partendo da un certo indice di inizio 
e arrivando a un certo indice di fine. Questo percorso viene creato esplorando i nodi del grafo in base alla direzione specificata (grafo diretto o inverso) 
e tenendo conto delle distanze massime consentite tra i nodi.
parametri:
*AdjGraphHead G: Puntatore alla struttura che contiene il grafo di adiacenza.
*CommonKmer commonKmer: Array di k-mer comuni.
unsigned long int startIndex, endIndex: Indici di inizio e fine del percorso nel grafo.
long int a: Indica se il percorso viene creato nel grafo diretto (0) o inverso (1).
long int leftIndex, rightIndex: Indici dei read sinistro e destro.
long int leftLen, rightLen: Lunghezza dei read sinistro e destro.
*FILE fp: Puntatore al file in cui scrivere i risultati.
*AdjGraphHead localG: Puntatore alla struttura che contiene il grafo locale.
*CommonKmerHead localCommonKmerHead: Puntatore alla struttura che contiene i k-mer comuni locali.
**char leftRead, rightRead: Puntatori ai read sinistro e destro.
*/
//riferimento algoritmo 2: Determine_Consistent_1
long int CreatGraphSinglePath(AdjGraphHead * G, CommonKmer * commonKmer, unsigned long int startIndex, unsigned long int endIndex, long int a, long int leftIndex, long int rightIndex, long int leftLen, long int rightLen, FILE * fp, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead) {

	//Viene inizializzato un puntatore al grafo (graph) e una variabile per contare il numero di nodi nel grafo (realCountGraph). La scelta del grafo (diretto o inverso) avviene in base al valore di a.
	AdjGraph * graph = NULL;
	long int realCountGraph = 0;
	//a indica se il percorso viene creato nel grafo diretto o inverso
	if(a == 0){ 
		//percorso creato nel grafo diretto
		graph = G->graph;
		realCountGraph = G->realCountGraph;
	}else{
		//percorso creato nel grafo inverso
		graph = G->reverseGraph;
		realCountGraph = G->reverseRealCountGraph;
	}
	
	printf("\n stampa del grafo graph\n");

	for(int i = 0; i < realCountGraph; i++){
		printf("\n dataLeft = %ld, dataRight = %ld\n", graph[i].dataLeft, graph[i].dataRight);
	}

	//Vengono inizializzati i conteggi delle distanze massime tra i nodi, sia per il grafo corrente (largestIntervalDistance) che per il grafo locale (maxIntervalDistance).
	long int largestIntervalDistance = G->largestIntervalDistance; //400
	long int maxIntervalDistance = localG->largestIntervalDistance; //1500
	
	for(long int n = 0; n < realCountGraph; n++){
		graph[n].visit = 0; //Vengono impostati tutti i nodi del grafo come non visitati.
	}

	startIndex = 0;
	endIndex = 0;

	//leftRead è la read corrente che stiamo analizzando 
	//rigthRead è quella che abbiamo ottenuto grazie alla commonKmer
	long int rightReadLength = strlen(rightRead); //si calcola la lunghezza della rightRead
	
	//Viene eseguito un loop sui nodi del grafo fino al penultimo nodo.
	for(long int j = 0; j < realCountGraph - 1; j++){
		if(graph[j].visit == 1){ //se il nodo è già stato visitato, si passa al prossimo nodo.
			continue;
		}
		
		if(graph[j].visit == 0 && a == 0){ //forward perchè a==0
			if(graph[j].dataLeft > graph[j].dataRight){
				//Viene verificato se il nodo corrente soddisfa i requisiti di distanza massima tra i nodi in base alla direzione del percorso (a).
				if(graph[j].dataRight > maxIntervalDistance){ //vede se è maggiore di 1500
					continue;
				}
			}else{
				if(graph[j].dataLeft > maxIntervalDistance){
					continue;
				}
			}
		}
		
		if(graph[j].visit == 0 && a == 1){ //reverse
			if(graph[j].dataLeft > rightReadLength - graph[j].dataRight){
				if(rightReadLength - graph[j].dataRight > maxIntervalDistance){
					continue;
				}
			}else{
				if(graph[j].dataLeft > maxIntervalDistance){
					continue;
				}
			}
		}

		startIndex = j;

		//se realCountGraph è 10, il ciclo for va da 0 a 9 quindi quando startIndex = 8
		//10-8 = 2 e ritorna 0. Questo significa che è rimasto un solo kmer quindi non si può effettuare 
		//il confronto. startIndex va da 0 a 9 ma quando è 8 ritorna zero
		if(realCountGraph - startIndex <= 2){
			return 0;
		}
		endIndex = startIndex + 1; // = j+1
		long int lastEndIndex = endIndex;
		long int iniStartIndex = startIndex;
		bool token = false;
		
		int edgeCount = 0; //contatore archi
		
		while(endIndex < realCountGraph && startIndex < realCountGraph){
			
			//Viene chiamata la funzione AddEdgeInGraph per aggiungere un arco al grafo.
			//riferimento algoritmo 4: verifica della consistenza
			long int edge = AddEdgeInGraph(G, a, startIndex, endIndex, largestIntervalDistance, maxIntervalDistance, localG, localCommonKmerHead, leftRead, rightRead);
			
			//Se viene aggiunto un arco, i nodi corrispondenti vengono contrassegnati come visitati e il 
			//conteggio degli archi (edgeCount) viene incrementato.
			//se il kmer è consistente, quindi edge == 1, esegue il codice dell'if
			if(edge == 1){
				graph[startIndex].visit = 1;
				graph[endIndex].visit = 1;
				lastEndIndex = endIndex;
				startIndex = endIndex;
				endIndex = startIndex + 1;
				token = true;
				edgeCount++;
			}else if(edge == 0){
				//sarebbe la funzione chaining_from_start: Aggiunge kmer[0] alla catena, poi confronta il kmer[0] con 
				//il kmer[1]; se non sono consistenti, confronta il kmer[0] con il kmer[2]. Quindi incrementa solo 
				//endIndex da 1 a 2 mentre startIndex rimane 0. Se sono consistenti,
				//aggiunge kmer[2] alla catena. Poi confronta kmer[2] con kmer[3] e cosi via..
				endIndex++;
				continue;
			}else{
				//i kmer non sono consistenti e l'arco non è stato aggiunto
				break;
			}
		}

		//Creata la catena di k-mer consistenti, ovvero aggiunti gli archi, si arriva all'ultima fase: determinare l'overlap sfruttando la catena.

		//Viene verificato se il percorso creato ha almeno 3 archi (indicativo di un percorso significativo).
		//riferimento algoritmo 2 (nella catena ci devono essere almeno 2 elementi, cioè 2 archi)
		//riferimento algoritmo 2: if |chain| > 2 then return chain else return null;
		if(token == true && edgeCount > 2){
			long int result = 0;
			
			//Se il percorso soddisfa i criteri, viene chiamata la funzione Overlap_Display_Graph per 
			//visualizzare il risultato e il valore di ritorno viene controllato.
			//a indica se il percorso viene creato nel grafo diretto (0) o inverso (1)
			if(a == 0){                                                                                                                                                          
				result = Overlap_Display_Graph(G, leftIndex,rightIndex,a,leftLen,rightLen,fp,graph[iniStartIndex].dataLeft, graph[lastEndIndex].dataLeft, graph[iniStartIndex].dataRight, graph[lastEndIndex].dataRight, localG,localCommonKmerHead, leftRead, rightRead);
			}else{
				result = Overlap_Display_Graph(G, leftIndex,rightIndex,a,leftLen,rightLen,fp,graph[iniStartIndex].dataLeft, graph[lastEndIndex].dataLeft, graph[lastEndIndex].dataRight, graph[iniStartIndex].dataRight, localG,localCommonKmerHead, leftRead, rightRead);
			}
			
			if(result == 1){
				return 1;
			}
		}
	}
	
	return 0;	
	
}

//La funzione CreatGraphLocalRegion costruisce grafi locali utilizzando piccoli k-mer per sequenze locali. 
//Questa funzione considera la distanza tra k-mer per determinare se esistono collegamenti validi tra nodi 
//adiacenti nel grafo. La funzione controlla la distanza e le condizioni di relazione tra i k-mer per creare
//archi tra i nodi. Funzione per costruire grafi locali di piccoli k-mer per sequenze locali
//Constructing graphs of small kmers for local sequences
long int CreatGraphLocalRegion(AdjGraphHead * G, long int distance) {
	long int count,throld;
	unsigned long  long int i,j,k,m,n;
	long int n1,n2,m1,m2;
	unsigned long int f1,f2;
	count = 0;
	long int minvalue;
	long int maxvalue;
	//Viene inizializzata una variabile graph per riferirsi al grafo contenuto nella struttura AdjGraphHead.
	AdjGraph * graph = G->graph;
	long int realCountGraph = G->realCountGraph; //contiene il numero effettivo di nodi nel grafo.

	long int firstIntervalDistance = distance;
	long int largestIntervalDistance = 2*distance;
	long int edgeCount = 0;

	//Un doppio ciclo for viene utilizzato per iterare attraverso i nodi del grafo.
	for(j=0; j<realCountGraph; ++j){ //Il primo ciclo itera attraverso ogni nodo nel grafo
		edgeCount = 0;
		 //Il secondo ciclo controlla i nodi successivi fino a un massimo di 5 nodi
		for(k=j+1; k<realCountGraph && k < j + 5; ++k){
			
			n1 = graph[j].dataLeft;
			n2 = graph[j].dataRight;
			m1 = graph[k].dataLeft;
			m2 = graph[k].dataRight;

			 // Se la distanza tra i k-mer supera la distanza massima consentita, interrompe il ciclo
			if(abs(m1-n1) > largestIntervalDistance && abs(m2-n2) > largestIntervalDistance){
				break;
			}
			// Calcola i valori minimi e massimi tra le distanze dei k-mer
			minvalue = min(abs(m1-n1),abs(m2-n2));
			maxvalue = max(abs(m1-n1),abs(m2-n2));
			// Se il rapporto tra la differenza dei valori massimo e minimo supera il rapporto di lunghezza, salta l'iterazione
			if(float(maxvalue - minvalue)/maxvalue >= G->lengthRatio){
				continue;
			}
			 // Aggiunge un arco tra i nodi se le condizioni sono soddisfatte
			if(n1<m1 && n2<m2){
				if(abs(m1-n1) < firstIntervalDistance && abs(m2-n2) < firstIntervalDistance || float(maxvalue - minvalue)/maxvalue < G->lengthRatio){	
					//Se le condizioni sono soddisfatte, un arco viene aggiunto tra i nodi j e k.
					//L'arco viene registrato nell'array arcIndex della struttura AdjGraphHead.
					G->arcIndex[G->realCountArc].startIndex = j;
					G->arcIndex[G->realCountArc].endIndex = k;
					G->arcIndex[G->realCountArc].weight =1;
					G->realCountArc++;
					//Se il numero di archi supera la capacità attuale (allocationCountArc), viene effettuata una riallocazione
					if(G->realCountArc >= G->allocationCountArc){
						ReAllocateArcIndex(G);
					}
				}
			}else if(n1>m1 && n2>m2){
				if(abs(m1-n1) < firstIntervalDistance && abs(m2-n2) < firstIntervalDistance || float(maxvalue - minvalue)/maxvalue < G->lengthRatio){
					G->arcIndex[G->realCountArc].startIndex = j;
					G->arcIndex[G->realCountArc].endIndex = k;
					G->arcIndex[G->realCountArc].weight =1;
					G->realCountArc++;  
					if(G->realCountArc >= G->allocationCountArc){
						ReAllocateArcIndex(G);
					}
				}
			}else{continue;}
				
		}
	}
	
	
}


/* Parametri
@param AdjGraphHead *G: Puntatore alla struttura che rappresenta il grafo di adiacenza globale.
@param CommonKmerHead *commonKmerHead: Puntatore alla struttura che contiene i k-mer comuni tra le letture.
@param ReadSetHead *readSetHead: Puntatore alla struttura che contiene le letture delle sequenze.
@param AdjGraphHead *localG: Puntatore alla struttura che rappresenta il grafo di adiacenza locale.
@param CommonKmerHead *localCommonKmerHead: Puntatore alla struttura che contiene i k-mer comuni locali.
@param FILE *fp: Puntatore al file dove verranno scritti i risultati.
@param realCommonKmerCount: Numero reale di k-mer comuni
@param readIndex: Indice della lettura corrente nel k-mer: viene inizializzato con il valore dell'indice di lettura del primo k-mer.
*/
// La funzione GetOverlapResult elabora un insieme di k-mer comuni per identificare le sovrapposizioni
// Utilizza grafi di adiacenza per rappresentare queste sovrapposizioni e seleziona il grafo ottimale per trovare un percorso singolo che rappresenti la sovrapposizione.
// getOverlapResult viene chiamato sulla read i-esima usando la commonKmerHead i-esima (che ottiene con kmerReadNode)
// viene chiamata tante volte quante sono le read
void GetOverlapResult(AdjGraphHead * G, CommonKmerHead * commonKmerHead, ReadSetHead * readSetHead, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, FILE * fp){
	printf("\n Funzione GetOverlapResult1\n");
	long int leftIndex, rightIndex, orien;
	long int leftLen,rightLen;
	long int realCommonKmerCount;
	long int s = 0 ;
	//readIndex è uguale alla readIndex di kmerReadNode ovvero la read dove è stato trovato il primo kmer che era presente anche nella hashTable
	unsigned long  long int readIndex = commonKmerHead->commonKmer[0].readIndex;
	unsigned long  long int startIndex = 0; //assume il valore di m 
	unsigned long  long int endIndex = 0; //
	bool orientation;

	
	printf("\n readIndex di dove siamo arrivati (CommonKmerHead->readIndex): %ld realCount: %ld allocationCount %ld\n",commonKmerHead->readIndex,commonKmerHead->realCount,commonKmerHead->allocationCount);
	
	printf("read index della prima cella di commonKmer: %ld\n",readIndex);

	printf("CommonKmerHead realCount: %ld",commonKmerHead->realCount);

	printf("\n Stampa della struttura CommmonKmerHead in getOverlapResult\n");

	for (int i = 0; i < commonKmerHead->realCount; i++){
		printf("\n common kmer read index: %ld\n",commonKmerHead->commonKmer[i].readIndex);
	}
		
	
	//Itera attraverso tutti i k-mer comuni (commonKmerHead->realCount).
	//commonKmerHead->realCount è il contatore dei kmer comuni nella read corrente perchè ogni volta viene
	//posta a zero (ma gli elementi nella commonKmer ci sono sempre)

	//star for
	for(long int m = 0; m < commonKmerHead->realCount; m++){
		//readIndex = commonKmerHead->commonKmer[0].readIndex;
		//inizialmente sono uguali quindi passa ad m = 1
		//se m=1 e entra nell'if, endIndex=0
		//printf("m=%ld\t", m);
		if(commonKmerHead->commonKmer[m].readIndex != readIndex){ //se l'indice della read del kmer comune[m] è											  
																// diverso dall indice della read del kmer comune [0]
			//entra nell'if quando m=5, quindi endIndex=4																
			endIndex = m - 1;

			//commonKmerHead è la struttura padre di commonKmer
			//leftIndex è l'indice della read corrente che sta valutando
			//readIndex è l'indice di dove siamo arrivati
			leftIndex = commonKmerHead->readIndex; //indice della read in cui si trova il common kmer

			//rightIndex punta alla prima cella di commonKmer
			//inizialmente readIndex = commonKmerHead->commonKmer[0].readIndex;
			//poi viene modificato in readIndex = commonKmerHead->commonKmer[m].readIndex;
			rightIndex = readIndex; 

			//printf("endIndex = %lu, leftIndex di commonKmerHead = %ld e rightIndex = %ld\n", endIndex, leftIndex, rightIndex);

			//si prende la lunghezza della read in cui si trova il commonKmer in posizione leftIndex-1
			leftLen = readSetHead->readSet[leftIndex - 1].readLength;

			//si prende la lunghezza della read in posizione rightIndex-1
			rightLen = readSetHead->readSet[rightIndex - 1].readLength;
			
			//quindi il contatore dei grafi dovrebbe essere uguale al contatore dei CKS
			//basandoci sull'algoritmo 1 presente nel paper
			//qui iniziano i grafi
			//riferimento all'algoritmo 1 nel paper: se count < 5 ritorna null perchè considera solo due reads
			//nel codice invece abbiamo più reads quindi non ritorna null
			if(G->realCountGraph < 5 && G->reverseRealCountGraph < 5){ //Algoritmo 1 c'è l'if di m < 5
				
				//reinizializza i contatori del grafo
				G->realCountGraph = 0;
				G->realCountArc = 0;
				G->reverseRealCountGraph = 0;
				
				if(commonKmerHead->commonKmer[m].orientation == 0){ //oriientation == 0 -> kmer in forward
					//aggiunge il kmer corrente al grafo o al grafo inverso a seconda dell'orientamento

					//dataLeft è la posizione j del kmer relativo alla read i-esima che sta analizzando la posizione del kmer salvato in kmerReadNode
					//dataRight è la posizione del kmer nella read puntata da kmerReadNode[m]: questo contiene la posizione del kmer della read[i]

					G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition; //posizione del kmer salvato in kmerReadNode
					G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition; //posizione del kmer della read[i]
					G->graph[G->realCountGraph].visit = 1;
					G->realCountGraph++;

					if(G->realCountGraph >= G->allocationCountGraph){
						/*Questa funzione aumenta dinamicamente la capacità di memorizzazione del grafo (diretto o inverso) 
							moltiplicandola per 1,5 e copia i dati esistenti nella nuova area di memoria allocata.
						*/
						ReAllocateAdjGraph(G, 0);
					}
				}else{ //se orientation == 1 ->reverse
					G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
					G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
					G->reverseGraph[G->reverseRealCountGraph].visit = 1;
					G->reverseRealCountGraph++;
					if(G->reverseRealCountGraph >= G->reverseAllocationCountGraph){
						ReAllocateAdjGraph(G, 1);
					}
				}
				
				//aggiorna startIndex e readIndex
				startIndex = m; // read valutata nel ciclo for
				readIndex = commonKmerHead->commonKmer[m].readIndex; //indice della read del commonKmer
				continue; //continua con il prossimo kmer
			}
			//Riferimento algoritmo 1 e pag. 3-4 nel paper: If M > N and M > count (count = 5), 
			
			//determina  l'orientamento del grafo
			if(G->realCountGraph > G->reverseRealCountGraph){
				orientation = false; //forward
				//dal paper: If N > M and N > count, LROD keeps the opposite k-mers and ignores the positive k-mers.
			}else{
				orientation = true;// reverse
			}
			
			//recupera le letture leftRead e rightRead
			//leftIndex è l'indice della read in cui si trova il common kmer
			//leftRead è la read corrente che stiamo analizzando 
			char * leftRead = readSetHead->readSet[leftIndex - 1].read;

			//rigthRead è quella che abbiamo ottenuto grazie alla commonKmer
			char * rightRead = readSetHead->readSet[rightIndex - 1].read;
			
			printf("left read %s\n",leftRead);
			printf("right read %s\n",rightRead);

			//Questa funzione viene utilizzata per aggiungere un arco al grafo di adiacenza, basandosi su criteri 
			//specifici relativi ai k-mer comuni e alle distanze tra di essi. Il risultato indica se l'aggiunta dell'arco è avvenuta con successo o meno.
			//riferimento algoritmo 1: creazione della catena nel paper
			//CreatGraphSinglePath è il Determine_consistent_1 nel paper
			long int result = CreatGraphSinglePath(G,commonKmerHead->commonKmer, startIndex, endIndex, orientation,leftIndex, rightIndex, leftLen, rightLen, fp, localG,localCommonKmerHead, leftRead, rightRead);
			
			//Reinizializza i contatori del grafo.
			G->realCountGraph = 0;
			G->realCountArc = 0;
			G->reverseRealCountGraph = 0;
			
			//Aggiorna startIndex e readIndex.
			startIndex = m;
			readIndex = commonKmerHead->commonKmer[m].readIndex;
		}
		// se l'orientamento del commonKmer[m] == 0
		if(commonKmerHead->commonKmer[m].orientation == 0){ //orientation == 0 -> kmer in forward
			//dataLeft è la posizione (position) del kmer nella kmerReadNode
			G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition; //posizione del kmer nella kmerRead
			//dataRight è la posizione del kmer nella read i-esima che stiamo analizzando
			G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition; //posizione del kmer nnella read

			G->graph[G->realCountGraph].visit = 1;
			G->realCountGraph++;
			if(G->realCountGraph >= G->allocationCountGraph){
				ReAllocateAdjGraph(G, 0);
			}
		}else{ //kmer in reverse
			G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition; 
			G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;

			G->reverseGraph[G->reverseRealCountGraph].visit = 1;
			G->reverseRealCountGraph++;
			if(G->reverseRealCountGraph >= G->reverseAllocationCountGraph){
				ReAllocateAdjGraph(G, 1);
			}
		}
		
	}
	//end for
	
	//Dopo il loop, aggiorna endIndex per includere l'ultimo k-mer.
	//nel nostro caso commonKmerHead->realCount = 8, quindi endIndex = 7
	endIndex = commonKmerHead->realCount - 1;
	printf("endIndex successivo = %ld\n", endIndex);
	//Se il numero di nodi nel grafo (reverse o forward) è maggiore o uguale a 5
	if(G->realCountGraph >= 5 || G->reverseRealCountGraph >= 5){
		if(G->realCountGraph > G->reverseRealCountGraph){
			orientation = false; //considera il grafo in forward
		}else{
			orientation = true; //reverse
		}
		leftIndex = commonKmerHead->readIndex; //indice della read che stiamo analizzando nel for della funzione chiamante
		rightIndex = readIndex; //indice della read index in common kmer

		leftLen = readSetHead->readSet[leftIndex - 1].readLength;
		rightLen = readSetHead->readSet[rightIndex - 1].readLength;
		char * leftRead = readSetHead->readSet[leftIndex - 1].read;
		char * rightRead = readSetHead->readSet[rightIndex - 1].read;

		// Questa funzione viene utilizzata per aggiungere un arco al grafo di adiacenza, basandosi su criteri
		// specifici relativi ai k-mer comuni e alle distanze tra di essi. Il risultato indica se l'aggiunta dell'arco è avvenuta con successo o meno.
		// riferimento algoritmo 1: creazione della catena nel paper
		// CreatGraphSinglePath è il Determine_consistent_1 nel paper
		long int result = CreatGraphSinglePath(G,commonKmerHead->commonKmer, startIndex, endIndex, orientation, leftIndex, rightIndex, leftLen, rightLen, fp, localG,localCommonKmerHead, leftRead, rightRead);
	
		//pone i contatori a 0 per il prossimo ciclo
		G->realCountGraph = 0;
		G->realCountArc = 0;
		G->reverseRealCountGraph = 0;	
	}

	//creazione della catena consistente, che corrisponde all'array di grafi contenuti in G.
}

//La funzione DetectCommon si occupa di trovare k-mer comuni tra due sequenze, considerando una certa distanza 
//massima dal k-mer attuale nella sequenza di lettura. Se trova un k-mer comune, lo aggiunge alla struttura CommonKmerHead.
void DetectCommon(CommonKmerHead * commonKmerHead, long int position, char * kmer, char * read, long int readLength, long int kmerLength, long int distance){
	bool t = true;
	//Calcola la distanza massima da considerare come il massimo tra la distanza fornita e 100
	//più un ulteriore margine di 100 per estendere la ricerca su un intervallo maggiore
	distance = max(distance, 100) + 100; 
	// Calcola la posizione di partenza per la ricerca nella lettura
	long int i = max(0, position - distance); //Questo garantisce che la ricerca non inizi da una posizione negativa.
	//Scansiona la sequenza di lettura dalla posizione di partenza fino alla posizione position+distance o fino
	//alla fine della sequenza cercando k-mer che corrispondano al k-mer fornito.
	for(; i < readLength - kmerLength + 1 && i < position + distance; i++){
		t = true;
		// Confronta il k-mer corrente con la sequenza di lettura
		//Per ogni posizione i, confronta i caratteri del k-mer con i caratteri della sequenza di lettura. 
		//Se tutti i caratteri corrispondono, il k-mer è considerato comune.
		for(long int j = 0; j < kmerLength; j++){
			if(kmer[j] != read[i + j]){
				t = false;
				break;
			}
		}
		// Se il k-mer è comune, le posizioni corrispondenti nella sequenza di lettura (rightPosition) e 
		//nella posizione originale (leftPosition) vengono salvate nella struttura commonKmerHead.
		if(t == true){
			commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = i;
			commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = position;		
			commonKmerHead->realCount++; //il contatore viene incrementato
			if(commonKmerHead->realCount >= commonKmerHead->allocationCount){
				ReAllocateCommonKmer(commonKmerHead);  // Se necessario, rialloca la memoria per commonKmerHead
			}
		}
	}
}

//La funzione GetCommonShorterKmer si occupa di trovare k-mer comuni tra due sequenze di lettura, costruire un 
//grafo locale basato su questi k-mer e determinare se esiste un sufficiente numero di k-mer comuni per giustificare un allineamento.
//When the distance between the two common kmers k1 and k2 is too large, the sequence between the two kmers is further compared with a smaller kmer. 
//If there are multiple small kmers with the same direction and adjacent positions, the k1 and k1 k2 belongs to the same consistent kmer.
long int GetCommonShorterKmer(AdjGraphHead * G, CommonKmerHead * commonKmerHead, char * leftRead, char * rightRead, long int leftStartPosition, long int leftEndPosition, long int rightStartPosition, long int rightEndPosition, long int kmerLength, bool orientation)
{
	//kmerLength = 9
	//Aggiornamento delle posizioni di inizio per includere la lunghezza del k-mer
	//leftStartPosition = n1 (posizione di partenza del kmer nella long read1)
	leftStartPosition = leftStartPosition + kmerLength;
	//rightStartPosition = n2 (posizione di partenza del kmer nella long read2)
	rightStartPosition = rightStartPosition + kmerLength;
	//riferimento algoritmo 2: il primo viene aggiunto (Adding CKS[start] to the chain;)
	commonKmerHead->realCount = 1; //Inizializzazione delle strutture dati per contenere i kmer trovati
	G->realCountGraph = 0;
	G->realCountArc = 0;
	commonKmerHead->commonKmer[0].rightPosition = 0;
	commonKmerHead->commonKmer[0].leftPosition = 0;	
	
	//Calcolo delle lunghezze temporanee delle letture
	//leftEndPosition = m1
	long int tempLeftReadLength = leftEndPosition - leftStartPosition + 1;
	//leftEndPosition = m2
	long int tempRightReadLength = rightEndPosition - rightStartPosition + 1;
	
	//Copia delle sottosequenze di lettura nelle strutture locali del grafico (in buffer locali)
	strncpy(G->localLeftRead, leftRead + leftStartPosition, tempLeftReadLength);
	strncpy(G->localRightRead, rightRead + rightStartPosition, tempRightReadLength);
	
	G->localLeftRead[tempLeftReadLength] = '\0';
	G->localRightRead[tempRightReadLength] = '\0';

	char * tempLeftRead = G->localLeftRead;
	char * tempRightRead = G->localRightRead;
	
	//Se l'orientamento è 1, calcola il complemento inverso della lettura destra	
	if(orientation == 1){ //orientation == 1 --> reverse	
		ReverseComplementKmer(tempRightRead, tempRightReadLength);
	}
	
	//Inizializzazione del k-mer
	char kmer[kmerLength + 1];
	long int kmerCount = tempLeftReadLength - kmerLength + 1;
	
	//Scansione della sequenza di sinistra per trovare k-mer di lunghezza kmerLength e confrontarli 
	//con la sequenza di destra
	for(long int i = 0; i < kmerCount; i++){
		strncpy(kmer, tempLeftRead + i, kmerLength);
		kmer[kmerLength] = '\0';
		//Se il k-mer non è un k-mer valido, continua con il successivo
		if(DetectSameKmer(kmer,kmerLength) != true){
			continue;
		}

		//Ogni k-mer viene confrontato con la sequenza di destra per trovare k-mer comuni.
		//I k-mer comuni vengono aggiunti alla struttura commonKmerHead.
		//Cerca k-mer comuni tra le due sequenze
		DetectCommon(commonKmerHead, i, kmer, tempRightRead, tempRightReadLength, kmerLength, abs(tempLeftReadLength - tempRightReadLength));
	}
	
	//Aggiungi un k-mer comune fittizio per garantire che ci sia almeno un k-mer alla fine di entrambe le sequenze
	commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = tempRightReadLength - kmerLength;
	commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = tempLeftReadLength - kmerLength;	
	commonKmerHead->realCount++;
	
	RemoveMultipleSameKmer(commonKmerHead); //Rimuove eventuali k-mer duplicati
	
	//Calcola il numero minimo di k-mer richiesti per considerare un allineamento
	long int num = 2 + tempLeftReadLength/300;
	if(commonKmerHead->realCount < num){
		//Se il numero di k-mer comuni trovati è inferiore a una soglia calcolata, la funzione restituisce 0.
		return 0;
	}
	
	//Aggiunge i k-mer trovati al grafico
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[i].leftPosition;
		G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[i].rightPosition;
		G->realCountGraph++;
		if(G->realCountGraph >= G->allocationCountGraph){
			ReAllocateAdjGraph(G, 0); //Rialloca la memoria se necessario
		}
	}
	
	CreatGraphLocalRegion(G,300); //Crea un grafo locale basato sui k-mer trovati
	
	//Verifica la sovrapposizione tra le due sequenze basata sul grafico locale
	long int result = Overlap_DisplayLocalRegion(G,tempLeftReadLength,tempRightReadLength);
	
	return result; //Il risultato della verifica viene restituito
}



#endif
