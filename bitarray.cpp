#include <stdio.h>
#ifndef BITS_H_INCLUDED 
#define BITS_H_INCLUDED 

#define BASE 2 
#define SHIFT 2 
#define MASK 3 


void print_array(unsigned long int *vect, int length){
    for(int i = 0; i < length; i++){
        printf("[%ld] ",vect[i]);
    }
    printf("\n");
}

//converte il valore i-esimo del kmer (value A,T,C,G) in un valore binario
// &= ~ + un AND bitwise -> ~ fa il complemento del bit (0 diventa 1)
// & -> operatore AND 
// << -> shift a sinistra ; x<<y sposta i bit di x a sinistra di y posizioni es 1011 << 1 = 0110
//quindi setto specifici bit dell'array, in base a value -> 
//A' potrebbe essere rappresentato come 00,
//'T' potrebbe essere rappresentato come 11,
//'G' potrebbe essere rappresentato come 10,
//'C' potrebbe essere rappresentato come 01,
//'N' potrebbe essere rappresentato come 00,

//Questa funzione imposta i bit all'interno di un array di bit in base al valore del carattere value. L'argomento bit_number indica l'indice del bit all'interno dell'array di bit che deve essere impostato.
//In breve, questa funzione traduce i nucleotidi in valori binari all'interno dell'array di bit.
//PRENDE l'array, l'indice i e il carattere i-esimo del singolo kmer
//bit_array sarebbe l'indirizzo 140733961459888
//*bit_array = 0
void SetBitTemp(unsigned long int * bit_array, unsigned long int bit_number, char value){

    //4printf("VALORE DI BIT_ARRAY IN SetBitTemp: %lu\n", bit_array);
    //5printf("VALORE DI *BIT_ARRAY IN SetBitTemp: %lu\n", *bit_array);
    //NOTA: gli interi, quando si utilizzano queste operazioni, vengono visti come cifre binarie
    //end (&): se entrambi i bit sono 1, allora il risultato è 1; altrimenti (in tutti gli altri casi), il risultato è 0

    if(value == 'A')
    {
        unsigned long int temp = 3; // 00000000011
        //bit_array è l'indirizzo 140733961459888
        //effettua l'end (&), poi il complemento (~), poi lo shift a sinistra (<<)
        //inzialmente temp = 3
        //bit_number è la i: alla prima iterazione è 1
        //quindi 2*bit_number = 2
        //quindi temp<<(2*bit_number) effettua lo shift a sinistra di due posizioni di temp
        //inizialmente temp=011=3, poi temp=1100=12, poi temp=110000=48, poi temp=11000000=192 
        //quindi lo shift a sinistra di due posizioni è una moltiplicazione di un intero * 4

        //1printf("temp<<(2*bit_number): %lu\n", temp<<(2*bit_number));
        //2printf("~(temp<<(2*bit_number)): %lu\n", ~(temp<<(2*bit_number)));

        (*bit_array) &= ~(temp<<(2*bit_number)); 
        //printf("valore di A = %lu\n", (*bit_array) &= ~(temp<<(2*bit_number)));
        //printf("%c = %ld\n",value,bit_array[bit_number] );
        //printf("%c = %ld\n",value,*bit_array);

    }
    if(value == 'T')
    {
        unsigned long int temp = 3;
        (*bit_array) |= (temp<<(2*bit_number));
        //printf("LETTERA=: %c, *bit_array=%d\n", value, (*bit_array));

        //printf("VALORE DI T: %lu\n", *bit_array);
        //printf("%c = %ld\n",value,bit_array[bit_number] );
        //printf("%c = %ld\n",value,*bit_array);

    }
    if(value == 'G')
    {
        unsigned long int temp = 1;
        (*bit_array) |= (temp<<(2*bit_number));
        temp = 1;
        (*bit_array) &= ~(temp<<(2*bit_number+1));
        //printf("LETTERA=: %c, *bit_array=%d\n", value, (*bit_array));

       // printf("%c = %ld\n",value,bit_array[bit_number] );
        //printf("%c = %ld\n",value,*bit_array);
    }
    if(value == 'C')
    {
        unsigned long int temp = 1;
        (*bit_array) |= (temp<<(2*bit_number+1));
        temp = 1;
        (*bit_array) &= ~(temp<<(2*bit_number)); 
        //printf("LETTERA=: %c, *bit_array=%d\n", value, (*bit_array));

        //printf("%c = %ld\n",value,bit_array[bit_number] );
        //printf("%c = %ld\n",value,*bit_array);
    }
    if(value == 'N')
    {
        unsigned long int temp = 3;
        (*bit_array) &= ~(temp<<(2*bit_number)); 
        //printf("LETTERA=: %c, *bit_array=%d\n", value, (*bit_array));

        //printf("%c = %ld\n",value,bit_array[bit_number] );
        //printf("%c = %ld\n",value,*bit_array);

    }

}
/*@param unsigned long int * bit_array -> array a 32 bit senza segno positivi
  @param long int len -> lunghezza del kmer
  @param char *value -> il kmer
*/

//Questa funzione imposta i bit di un array di bit in base ai valori di un k-mer rappresentato come una stringa di caratteri.
//La funzione prende tre parametri: un puntatore a un array di bit, la lunghezza del k-mer e il valore del k-mer rappresentato come una stringa di caratteri.
//bit_array sarebbe &kmerInteger, cioè l'indirizzo: 140733961459888
//kmerLength = 15
//kmer è uguale proprio alla stringa del kmer. La prima è AAAAAAAAAAAAAAC
void SetBitKmer(unsigned long int * bit_array, long int len, const char * value){
    //printf(" In set bit kmer = %ld\n",*bit_array);
    int i = 0;
    //per i chee va da 0 a 15-1
    for(i = 0; i<len; i++){ //Utilizzando un ciclo for, la funzione scorre attraverso ogni carattere della stringa value.
        //bit_array passato a SetBitTemp = 140733961459888
        //i sarebbe l'indice corrente e value[i] è l'i-esimo carattere del kmer
        SetBitTemp(bit_array, i, value[i]); // imposta il bit corrispondente nell'array di bit in base al valore del carattere ('0' o '1').
        //printf ("bit array passato a setBitTemp: %lu\n", bit_array);
    } 
    //printf("VALORE ARRAY DI BIT: %ld\n",*bit_array); //stampa il valore dell'array di bit dopo aver impostato i bit del k-mer.
}

char GetBit(unsigned long long int * bit_array, unsigned int bit_number){
    unsigned long long int temp = 3;
    int a = ((*bit_array)>>(2*bit_number)) & temp;
    printf("tyuioidh3gyiufohinbyugroifhb3vycg879ifhbivgyreohijcj VALORE DI A IN GetBit: %d\n", a);
    if(a == 0){
        return 'A';
    }
    if(a == 3){
        return 'T';
    }
    if(a == 1){
        return 'G';
    }
    if(a == 2){
        return 'C';
    }
}

void GetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value){
    int i = 0;
    for(i = 0; i<len; i++){
        value[i] = GetBit(bit_array, i);
    }
    value[i] = '\0';
}

void SetBit(char bit_array[], unsigned int bit_number, char value) 
{ 
    unsigned int shift = 6 - BASE*(bit_number & MASK);
    
    if(value == 'A')
    {
        bit_array[bit_number >> SHIFT] &= ~ (3 << shift); 
    }
    if(value == 'T')
    {
        bit_array[bit_number >> SHIFT] |= 3 << shift; 
    }
    if(value == 'G')
    {
        bit_array[bit_number >> SHIFT] &= ~ (2 << shift);
        bit_array[bit_number >> SHIFT] |= 1 << shift;
    }
    if(value == 'C')
    {
        bit_array[bit_number >> SHIFT] &= ~ (1 << shift);
        bit_array[bit_number >> SHIFT] |= 2 << shift; 
    } 
    if(value == 'N')
    {
        bit_array[bit_number >> SHIFT] &= ~ (3 << shift); 
    }
} 

char GetBit(char bit_array[], unsigned int bit_number)
{
    unsigned int shift = 6 - BASE*(bit_number & MASK);
    char base = (bit_array[bit_number >> SHIFT] & (3 << shift)) >> shift; 
    if(base == 0)
    {
        return 'A';
    }
    if(base == 3)
    {
        return 'T';
    }
    if(base == 1)
    {
        return 'G';
    }
    if(base == 2)
    {
        return 'C';
    }
}

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value) 
{ 
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        SetBit(bit_array, bit_start_number + i, value[i]);
    }
} 

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value, int index) 
{ 
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        SetBit(bit_array, bit_start_number + i, value[i]);
    }
    int j = ((bit_length/4) + 1)*4;
    for(i = i; i<j; i++){
        SetBit(bit_array, bit_start_number + i, 'A');
    }
} 

char * GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length)
{
    unsigned int i = 0;
    char * value = new char[bit_length + 1];
    for(i = 0; i < bit_length; i++)
    {
        value[i] = GetBit(bit_array, bit_start_number + i);
    }
    value[i] = '\0';
    return value;
}  

void GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value)
{
    unsigned int i = 0;
    for(i = 0; i < bit_length; i++)
    {
        value[i] = GetBit(bit_array, bit_start_number + i);
    }
    value[i] = '\0';
} 


#endif 
