# R-Soft-Inverted-Pendulum
Modellazione e controllo del sistema Soft Inverted Pendulum con un giunto rotoidale alla base.

## Contenuti:
* [1. Requisiti](#1-requisiti)
* [2. Guida al Codice](#2-guida-al-codice)
* [3. Modellazione](#3-modellazione)
* [4. Proprietà Strutturali](#4-proprietà-strutturali)
* [5. Controllo](#5-controllo)
* [6. Controllo Adattivo](#6-controllo-adattivo)
* [7. Simulazioni](#7-simulazioni)

## 1) Requisiti
Per eseguire gli scripts, è necessario scaricare il Robotics Toolbox di Peter Corke, reperibile
al seguente [link](https://petercorke.com/toolboxes/robotics-toolbox/). Inoltre, per eseguire i file
simulink, è necessario avere una versione superiore al R2021b.

## 2) Guida al Codice
Il pkg è organizzato con 3 directory, in cui sono depositate le funzioni utili per l'analisi e il controllo
del modello. Di seguito una descrizione di queste 3:

- `Della Santina`: Directory contenente le funzioni fornite alla consegna del progetto. Queste funzioni servono per 
confrontare il modello originale dell'articolo con l'implementazione in MATLAB.

- `my_functions`: Directory contenente funzioni di math-utility.

- `origin_soft_pendulum`: Directory contenente funzioni autogenerate da MATLAB che implementano la dinamica
del Soft Inverted Pendulum.

Il resto degli scripts e delle funzioni sono descritte nelle sezioni successive.

## 3) Modellazione
- La modellazione del sistema Soft-Inverted Pendulum è svolta dallo script `Soft_origin.m`. Esso, dopo aver calcolato
gli elementi dinamici di interesse, genera delle funzioni utili per le simulazioni. 

- La modellazione del sistema R-Soft Inverted Pendulum è invece implementata dallo script `Rsoft_model.m`. Analogamente allo script
precedente, esso implementa il calcolo degli elementi dinamici di interesse e genera delle funzioni utili per l'implementazione in simulazione.
Inoltre, lo script effettua una linearizzazione del sistema.

- L'analisi degli equilibri è svolta dalo script `equilibria_analysis.m`, il quale analizza la stabilità, la raggiungibilità e l'osservabilità del sistema linearizzato.
Gli equilibri sono ricavati dai file `.mat` generati dallo script `Roft_model.m`.

- Gli scripts `validazione_coriolis.m` e `validazione_dinamica.m` confrontano i valori delle matrici dinamiche
calcolate dagli scripts precedenti, particolarizzandole per un set di valori delle variabili di stato soddisfacente.

- Una simulazione puramente cinematica è implementata dallo script `test_kin.m`.

## 4) Proprietà Strutturali
Per un controllo più approfondito delle proprietà strutturali del sistema, si effettua un calcolo iterativo
della distribuzione di accessibilità e della codistribuzione di osservabilità. Le distribuzioni sono calcolate nello script `acc_obs_dist.m`. 
Inoltre, sono state implementate diverse funzioni utili:
- `lieBracket.m`: Essa implementa il calcolo della Lie Bracket tra vettori.
- `rowLieBracket.m`: Essa implementa il calcolo della Lie Bracket tra un vettore e un covettore.
- `filtration.m`: Essa implementa il calcolo iterativo della filtrazione, utile per la distribuzione di accessibilità.
- `rowFiltration.m`: Essa implementa il calcolo iterativo della filtrazione per la codistribuzione di osservabilità.

## 5) Controllo
Per il controllo, è stato implementato lo script `feedback_lin.m` che calcola gli ingressi e il nuovo stato, impostata l'output da controllare.
Inoltre, sono state scritte le funzioni `collocatedFL.m` e `collocatedFL2.m` che implementano rispettivamente la lin. in feedback su \theta_r e \alpha_{s}.
L'implementazione del controllore sarà poi presente nel file simulink `R_soft_sim.slx`, il quale ha una sezione apposita.

## 6) Controllo Adattivo
Il controllo adattivo è implementato dalla funzione `regressorSoftInverted.m`. Essa si occupa di restituire il regressore del Soft Inverted Pendulum indicato in questo [articolo](https://ieeexplore.ieee.org/abstract/document/9482817).

## 7) Simulazioni
Il file simulink `R_soft_sim.slx` contiene diverse simulazioni.
- **Collocated Feedback Linearization Spong Like**: Controllo sulla variabile di giunto alla base. 

- **Feedback Linearization Alpha**: Controllo sull'inclinazione della tip rispetto al sistema di riferimento fisso.

- **Autonomous Soft Inverted Pendulum**: Sistema Soft Inverted Pendulum Autonomo o con ingressi lentamente variabili.

- **Adaptive Controller Soft Inverted Pendulum**: Controllo adattivo del sistema Soft Inverted Pendulum. I risultati sono del tutto paragonabili con quelli
dell'[articolo](https://ieeexplore.ieee.org/abstract/document/9482817). 

Lo script `test_din.m` non fa altro che dichiarare i parametri della simulazione e registrarne i valori di simulazione. E' stata poi implementata la funzione `plot_Rosft.m` che visualizza il sistema e la traiettoria. 
