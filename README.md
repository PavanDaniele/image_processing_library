# Image processing library

- [Scppo](#Scopo)
- [Introduzione al progetto](#Introduzione-al-progetto)
- [Il Progetto](#Il-Progetto)
   - [PARTE 1: strutture dati, operazioni matematiche e gestione memoria](#PARTE-1:-strutture-dati,-operazioni-matematiche-e-gestione-memoria)
   - [PARTE 2: Operazioni semplici con immagini](#PARTE-2:-Operazioni-semplici-con-immagini)
   - [PARTE 3: Convoluzione e filtraggio](#PARTE-3:-Convoluzione-e-filtraggio)
- [Green-Screen (o chroma-key)](#Green-Screen-(o-chroma-key))
- [Equalizzazione Immagine](#Equalizzazione-Immagine)
- [Gestione Errori](#Gestione-Errori)
- [Green-Screen (o chroma-key)](#Green-Screen-(o-chroma-key))

- Equalizzazione Immagine

# Scopo 
Realizzare una libreria in linguaggio C++ per effettuare operazioni di elaborazione di immagini (image processing).
Il progetto è diviso in tre parti: 
-	Per la prima parte si dovrà implementare un insieme di metodi per la gestione di matrici a 3 dimensioni utile per comprendere alcune operazioni matematiche e di gestione della memoria.
-	nella seconda parte dovranno essere sviluppati alcuni semplici metodi per l’elaborazione di immagini (conversione da colore a scala di grigi, blending).
-	nell’ultima parte verrà richiesto di implementare l’algoritmo di convoluzione a cui si applicheranno filtri per ottenere delle trasformazioni dell’immagine.

# Introduzione al progetto
Le immagini sono codificate nel computer tramite griglie ordinate di pixel. Un pixel è l’entità più piccola per descrivere il colore di un’immagine in un determinato punto.

I pixel di un’immagine tipicamente hanno valori compresi tra 0 e 255, pertanto saranno necessari 8 bit per rappresentare un pixel.

Le immagini sono rappresentate in memoria tramite matrici (l’ordine dei pixel è fondamentale). Nel caso di immagini in scala di grigi abbiamo solo un canale (il canale di luminosità). 

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/0c56fdca-07b8-4505-bf5b-50087187ba07)

Nel caso di immagini a colori (il nostro caso), per rappresentare i colori visibili vengono utilizzati 3 canali (rosso, verde e blu). Ogni pixel quindi è composto da 3 valori interi (rosso, verde e blu). Pertanto la matrice avrà tre dimensioni: larghezza x altezza x 3. Tipicamente si parla di immagini a 24 bit (8 bit * 3 canali = 24 bit).

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1a96ee56-b6bb-48a3-964b-38dd8edb70f2)

La libreria in C++ per leggere e scrivere immagini in formato BITMAP vi viene fornita nel repository del progetto.

Eseguendo il comando “make testbmp”, verrà generato un file eseguibile “test_bmplib” che, una volta eseguito, creerà una scacchiera colorata.

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/919d5e37-1a20-45f8-be44-d74948216267)

# Il Progetto
L’image processing, o elaborazione delle immagini, è una branca della computer science nata negli anni ‘60 con il fine di applicare delle trasformazioni alle immagini più o meno in modo automatico. Applicazioni tipiche sono la rimozione del rumore, l’equalizzazione, la conversione di formati, scalatura, ritaglio o il rilevamento di bordi. In sostanza le operazioni basilari che trovate all’interno di Photoshop o GIMP. 

## PARTE 1: strutture dati, operazioni matematiche e gestione memoria
Per le operazioni successive dobbiamo realizzare una nuovo tipo di dato chiamato tensor . Il tipo tensor altri non è che una matrice a 3 dimensioni (altezza x larghezza x canali, vedi figura sotto) in cui saranno memorizzati i pixel dell’immagine e sui quali applicheremo varie trasformazioni. Nella conversione da Bitmap a Tensor, ai pixel viene fatto un casting a float in quanto le operazioni di elaborazione avranno bisogno di lavorare con dati in virgola, positivi e negativi. Dovremo quindi abbandonare i pixel [0,255] per un po’.
La classe a cui ci appoggeremo (da implementare) sarà la seguente:

class Tensor{

private:
    float * data;
    int c;
    int r; 
    int d;

public:
    Tensor();
    ~Tensor(); 
. . . }

La matrice a tre dimensioni sarà memorizzata nella variabile data e il suo accesso avverrà tramite overloading dell’operatore ( ) passando tre indici: riga, colonna e profondità. 

Viene lasciata libertà sulla struttura di memoria (la variabile data della struttura di cui sopra) purchè non vengano utilizzate oggetti Vector e Array

Vanno implementati gli operatori matematici (+ - * / ), di selezione (i ,j , k) e quello di assegnamento ( = )

Vedere il file Tensor.h per maggiori dettagli.


## PARTE 2: Operazioni semplici con immagini
Le operazioni su immagini vengono eseguite per mezzo della classe DAISGram che contiene all’interno un attributo data di tipo Tensor.

Si dovranno implementare i seguenti metodi (maggiori dettagli in DAISGram.h):

### Image Brightening:
Aumenta il valore di tutti i pixel di una costante bright. L’effetto è alzare la luminosità complessiva dell’immagine.

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/22d52b6d-e5ad-43fc-8af7-2468463b8ed3)

### Conversione a scala di grigi:
Converte un’immagine a colori in una a scala di grigi

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/fde45e80-dff6-4e5c-9fae-90871543f6da)

### Image Blending:
Fonde due immagini insieme tramite combinazione convessa. 
Blend = alpha * A + (1-alpha)* B;

alpha è in [0,1]

Le immagini “a” e “b” devono avere le stesse dimensioni, altrimenti l’operazione non è possibile.

- Immagine A: ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/f577c68b-ed64-4e2b-944a-b34e7a37ffc3)
- Immagine B: ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/9d62a2b6-6334-4655-aa2d-f457ea2556fa)
- Blend (alpha = 0): ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/7c562b3a-a691-440d-badc-51ee76e4c32a)
- Blend (alpha = 0.25): ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1a2c2873-284b-44f7-b7b6-82787b26db7c)
- Blend (alpha = 0.5): ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/bb49c9e5-03cf-4659-a871-c8bd314d69fa)
- Blend (alpha = 0.75): ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/547c915e-d2f9-4725-af61-3b33d993614c)
- Blend (alpha = 1): ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1c8a1528-d7fb-4c2a-8d40-5a1196a932a4)

### Andy Warhol:
Crea una nuova DAISGram contenente l’immagine in stile “Andy Warhol”.

Input: ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/e56233f8-10fa-4c33-b604-3eded46df859)     
Output: ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/07942773-637d-4452-924b-859b2c9ef2a5)

L’immagine risultante è una composizione di 4 immagini così realizzate:
-	In alto a sinistra viene replicata l’immagine originale
-	In alto a destra, a partire dall’immagine originale, viene invertito il canale Rosso con il canale Verde
- In basso a sinistra, a partire dall’immagine originale, viene invertito il canale Verde con il canale Blu
-	In basso a destra, a partire dall’immagine originale, viene invertito il canale Rosso con il canale Blu

L’immagine finale avrà quindi dimensione doppia dell’originale (a parte nel numero di canali) in quanto le 4 immagini saranno concatenate per formarne una unica come da esempio sopra.

## PARTE 3: Convoluzione e filtraggio
La convoluzione è una delle tecniche di analisi di segnali maggiormente applicate, ed è inoltre largamente applicata nel contesto di immagini. E’ uno dei metodi di processamento di immagini più semplici ma potenti e rappresenta uno dei fondamenti della visione artificiale, nell’analisi di segnali, nelle moderne reti di intelligenza artificiale e nelle reti profonde (deep learning). 

La convoluzione, non è altro che una somma pesata dei valori dell’immagine rispetto a quelli di un filtro, chiamato kernel. Il kernel è una matrice di valori reali (tipicamente sono matrici quadrate a dimensioni dispari non molto grandi, come una 3x3 o 5x5) che definiscono “l’importanza” dei pixel sottostanti. In funzione della configurazione dei valori del filtro abbiamo diversi effetti sull’immagine di input.

Il kernel abbiamo detto essere una matrice di float, possiamo quindi utilizzare la stessa struttura Tensor anche per memorizzare il filtro.

Come funziona al lato pratico? prendiamo ad esempio la figura qui sotto. Abbiamo un’immagine 5x5x1 e un filtro 3x3x1: 
Fig 1: esempio di calcolo della convoluzione tra immagine (sinistra), filtro (centro). Il risultato è a destra

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/c22b9704-5779-4163-9ceb-b73b60a5231a)
Sovrapponiamo il kernel all’immagine in alto a sinistra (prima posizione) e calcoliamo la convoluzione (somma dei prodotti):

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/e006e952-fe11-4858-a86b-df64ecd44240)

I valori dei pixel dell’immagine sotto il filtro sono: IMG:| 7 | 2 | 3 | 4 | 5 | 3 | 3 | 3 | 2 |

Quelli del filtro: K:| 1 | 0 | -1 | 1 | 0 | -1 | 1 | 0 | -1 |

Il risultato della convoluzione è la somma dei prodotti, quindi:   7 * 1 + 2 * 0 + 3 * -1 + 4 * 1 + 5 * 0 + 3 * -1 + 3 * 1 + 3 * 0 + 2 * -1 = 6 

Questo sarà il valore del pixel in alto a sinistra nell’immagine risultante, ovvero quello al centro del kernel. La stessa operazione viene effettuata su tutta l’immagine (vedi animazione) facendo scorrere il filtro di una posizione alla volta su tutte le righe e colonne.

L’operazione di convoluzione è implementata in Tensor.h 
Tensor convolve(const Tensor &f); 

e riceve in input un filtro (Tensor f), applica il padding all’oggetto (vedi sotto) e ne calcola la convoluzione. 

Si dovranno implementare 4 filtri da applicare all’immagine per mezzo della convoluzione:

#### Sharpen
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/d352a345-91cf-45d9-b023-4913efb308ea)

Matrice 3x3, permette di enfatizzare i dettagli
Valori:
- 0  -1  0
- -1  5 -1
- 0  -1  0

#### Edge
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/acf79105-c4a2-4d78-80a5-6c82d65a8f26)

Matrice 3x3, calcola i contorni dell’immagine
Valori:
- -1  -1  -1
- -1   8  -1
- -1  -1  -1

#### Emboss
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/e7f805cb-87dd-4024-a7a0-d6ec185c5546)

Matrice 3x3, aggiunge un senso di profondità o rilievo all’immagine
Valori:
- -2 -1  0
- -1  1  1
-  0  1  2


#### Smoothing
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/b449ca65-2ae8-4a52-b18d-3f829092d4ea)

Matrice hxh, permette di rimuove il rumore.  
c = 1/(h*h)
Ad esempio w = 3 e h = 3 => c = 1/9 
Valori:
- c c c
- c c c
- c c c

## Problema 1: le immagini sono rappresentate tramite una matrice a 3 dimensioni ma il filtro è a 2 dimensioni.
possiamo quindi:
1) avere filtri diversi per ogni canale (uno per il canale rosso, uno per il verde ed uno per il blu)
2) utilizziamo lo stesso filtro per tutti i canali. 

Nel nostro caso vogliamo un filtro diverso per ogni canale, nel caso in cui lo stesso filtro sia applicato a più canali, questo si dovrà copiare sui relativi canali del filtro. In sostanza, filtro e immagine devono avere la stessa depth.

## Problema 2: Cosa fare sui bordi? 
Come potete notare, l’immagine risultante dalla convoluzione è più piccola dell’immagine originale perdendo l’informazione sui bordi. Ad esempio nella figura sopra otteniamo un’immagine con dimensione 3x3x1 (vedi parte destra della Figura 1).

Ci sono varie tecniche per risolvere questo problema e ottenere un’immagine delle dimensioni originali dopo la convoluzione. Nel nostro caso useremo una tecnica chiamata padding. Nel padding si estende l’immagine con un bordo di larghezza opportuna di pixel tutti a zero. In questo modo, a convoluzione applicata, otterremo un’immagine della stessa dimensione dell’input. 

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1541e65e-6566-4802-a0da-6344cd884a6f)

Nella figura a sinistra abbiamo: in blu l’immagine (a singolo canale), in bianco il padding, in grigio il filtro e in verde l’immagine risultante. Notate come la dimensione dell’immagine verde sia la medesima dell’immagine blu.

Ma quanto bordo aggiungere? Dipende dalle dimensioni del filtro! 

In generale il padding si calcola come:
P = (F-1)/2 mantenendo la parte intera. Dove P è il valore del padding risultante ed F è la dimensione del filtro.

Quindi un filtro 3x3 avrà come padding orizzontale (3-1)/2 = 1 e verticale (3-1)/2 = 1. Un filtro 5x5 avrà come padding orizzontale (5-1)/2 = 2 e verticale (5-1)/2 =2.

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/78ca7ee0-760a-4cfb-95ea-9e2ce04e172c)

Se il filtro è di dimensione 3x3 estenderemo quindi l’immagine di 1 pixel per lato, ad esempio se un’immagine è a 283x362x3, l’immagine risultante con il padding sarà di 285*364x3 (1 pixel sul lato sinistro, 1 sul lato destro, 1 pixel in alto e 1 sotto). 
#### (NB: L’operazione di padding viene effettuata, ovviamente, prima della convoluzione.)

## Problema 3: Una volta calcolata la convoluzione i valori hanno un range diverso da [0,255]!
Come potete vedere nella figura 1, i valori risultanti dalla convoluzione non sono in [0,255]. Prima di salvare un’immagine bisogna convertire i valori del campo data nel range di cui sopra. Pertanto dobbiamo “scalare” il nostro Tensore. 

Si dovranno implementare due metodi (vedete il dettaglio nel file Tensor.h ):
void clamp(float low, float high);
void rescale(float new_max=1.0);

Il metodo di rescale porta i valori del Tensore in [0,new_max] tramite questa operazione:
T = (T - min(T))/(max(T)-min(T));

dovranno quindi essere moltiplicati per new_max in modo da essere nel range [0,new_max].

Il metodo clamp vincola l’intervallo di valori tra un valore minimo e uno massimo (maggiori dettagli in tensor.h).

# Green-Screen (o chroma-key):
Il green screen è una tecnica cinematografica che permette di sostituire uno sfondo a tinta uniforme con uno desiderato. E’ tutt’ora largamente utilizzato in tantissime produzioni (Matrix, Avengers etc…).

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1bceb0a7-0815-499d-a9c8-7303361b92ba)  ![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/3e3ae752-379f-47d5-a63f-7cd8a07bc6aa)

L’idea alla base è quella che, date due immagini (una di primo piano e una di sfondo), vengano sostituiti i pixel di un certo colore nell’immagine in primo piano con i corrispondenti pixel dell’immagine di sfondo.

Si dovrà implementare il metodo:  DAISGram greenscreen(DAISGram & bkg, int rgb[], float threshold[]);
che riceve in input un background (delle stesse dimensioni dell’oggetto attuale), il colore di sfondo (come vettore rgb[3]={0,128,00}) e un vettore di soglia (threshold[3]={10,20,10}) che decide per ogni canale l’intervallo intorno al colore selezionato.

Nell’esempio di cui sopra andremo a scambiare i pixel tra primo piano e sfondo se questa condizione è soddisfatta per tutti e 3 i canali:
pixel(i,j,0)>=(rgb[0]-threshold[0]) && pixel(i,j,0)<=(rgb[0]+threshold[0]) …

![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/2f683ff5-27f3-4728-a78c-08cc6fd048ef)
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/35681aa5-ca18-48fc-8256-2bd46e18cead)
![image](https://github.com/PavanDaniele/image_processing_library/assets/127297363/1a798a7a-bd1a-4f7f-acbb-6e5ee0f626e0)

# Equalizzazione Immagine:
L’equalizzazione dell’immagine è una tecnica che permette di utilizzare tutto il range di valori che i pixel possono acquisire. In questo modo il range di luminosità viene “stirato” rendendo le immagini piu vive. DAISGram equalize();

Il metodo è ben descritto in questo link . 

Consideriamo l’equalizzazione indipendente per canale, quindi dato un canale:
1.	Creare un istogramma della distribuzione delle intensità (vettore di 256 elementi)
2.	Calcolare la funzione cumulativa (cdf) dell’istogramma.
3.	Calcolare il minimo della funzione cumulativa (cdfmin)

Dato un valore v di luminosità, la sua versione equalizzata h(v) è:

# Gestione Errori
Nel caso in cui vengano passati dei parametri sbagliati, ad esempio le dimensioni di due Tensor non sono concordi, oppure Tensor è a NULL, si dovrà generare un errore e bloccare l’esecuzione. Gli errori sono gestiti tramite eccezioni che trovate dentro il file dais_exc.h.



