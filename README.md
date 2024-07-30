# Introduzione

Questo repository contiene il codice sorgente e la documentazione per l'algoritmo LROD (Long Read Overlap Detection), migliorato per il rilevamento delle sovrapposizioni nelle letture lunghe di DNA. Le modifiche apportate all'algoritmo includono l'integrazione del calcolo delle frequenze dei k-mer direttamente all'interno di LROD, eliminando così la dipendenza da strumenti esterni come DSK.

# Per installare LROD:

Crea una directory principale, importando tutti i codici sorgenti, o clona questa repository (es: LROD). <br>

Esegui i seguenti comandi:

```
cd LROD
make all
```

Per eseguire LROD, utilizza il seguente comando:

```
./LROD -r <long-read-file> -c <kmer-frequency-file> -o <result-file> [opzioni]
```

Se non viene specificato il `<kmer-frequency-file>`, esso verrà creato automaticamente.

Bisogna installare la libreria [ntHash](https://github.com/bcgsc/ntHash?tab=readme-ov-file) per il corretto funzionamento. 

# Parametri


```
-r <long-read-file>: file di input in formato fasta.
-c <kmer-frequency-file>: ogni riga nel file delle frequenze dei kmer deve essere "kmer kmer-frequency".
-o <result-file>: file di output.
-t <count>: numero di thread (predefinito 1).
-k <kmerLength>: lunghezza del kmer (predefinito 15).
-q <smallKmerLength>: lunghezza del piccolo kmer (predefinito 9).
-f <minimumKmerFrequency>: frequenza minima del kmer (predefinito 2).
-m <maxKmerFrequencyRatio>: rapporto massimo di frequenza del kmer (deve essere inferiore a 1, predefinito 0.9).
-s <kmerStep>: passo del kmer (predefinito 1).
-d <distance>: piccola distanza usata per determinare se due kmer comuni sono coerenti (predefinito 400).
-e <distance>: grande distanza usata per determinare se due kmer comuni sono coerenti (predefinito 1500).
-a <min-overlap-length>: lunghezza minima di sovrapposizione tra due letture lunghe (predefinito 500).
-b <length-ratio>: rapporto massimo di lunghezza tra due regioni allineate (predefinito 0.3).
-h: mostra le regole di utilizzo di LROD.
--generate-kmer-file-only: genera solo il file delle frequenze dei kmer e termina.

```
Esempio di Comando

```
./LROD -r w303_illumina.fa -c kmer_file_test_w303.txt -o result_w303
```

Esempio di Comando per la sola generazione del file delle frequenze

```
./LROD -r ecoli_illumina.fa -c kmer_file.txt -o result --generate-kmer-file-only
```

# Output


Il file di output `<result-file>` contiene i risultati delle sovrapposizioni. La struttura del file è la seguente:

```
Prima colonna: numero della prima lettura.
Seconda colonna: numero della seconda lettura.
Terza colonna: orientamento dell'allineamento (0 per allineamento in avanti, 1 per allineamento inverso).
Quarta colonna: posizione iniziale nella prima lettura.
Quinta colonna: posizione finale nella prima lettura.
Sesta colonna: lunghezza della prima lettura.
Settima colonna: posizione iniziale nella seconda lettura.
Ottava colonna: posizione finale nella seconda lettura.
Nona colonna: lunghezza della seconda lettura.

```
Esempio di una riga nel file di output:

```

36423,1,0,5326,9923,9923,1,4364,10479

```

Significa che la regione [5326,9923] nella prima lettura è sovrapposta in avanti con la regione [1,4364] nella seconda lettura.

Nella directory Testing, sono presenti i risultati dei benchmark.<br>
Nella directory Datasets, i dataset utilizzati per il testing.

I file delle frequenze vengono automaticamente generati dal codice se non presenti.

