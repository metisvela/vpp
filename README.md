# vpp
Velocity Prediction Program altamente modulare per prevedere il comportamento in acqua degli skiff a vela e confrontare gli effetti di setup diversi (set di vele, piani di deriva, disposizione foil) sulle prestazioni.
# Funzionamento
## Moduli da importare
Per partire, è necessario prima installare i pacchetti:
- numpy
- matplotlib
- time
- pandas
- scipy.interpolate
- scipy.optimize
- scipy.misc
- sys
- [xfoil](https://pypi.org/project/xfoil/)

Purtroppo l'ultimo modulo è necessario anche per una barca senza foil (viene usato anche per simulare deriva e timone) e sembra dare problemi di installazione su MacOS e Linux.

## Come funziona
Per analogia con altri programmi per l'ingegneria, questo vpp è suddiviso in pre-processing, processing e post-processing.
### Pre-processing
**Nota**: Se si vuole solo ottenere dei risultati senza cambiare le funzioni sottostanti il calcolo, è sufficiente modificare i moduli **main** e **input_data**.

Nel modulo **input_data** vanno inseriti tutti i dati riguardanti barca, foil, equipaggio, vele, campo di regata, deriva, profili idrodinamici. Seguire le istruzioni per essere sicuri di usare il formato corretto, soprattutto per i profili.

Il modulo **main** può essere considerato il *front-end*. Nella prima sezione vengono inizializzate le classi corrispondenti a barca, vele, equipaggio, foil, deriva, timone (non ancora implementato), campo di regata. Ogni classe prende come input il dizionario uscente da *input_data_dictionary* che contiene i suoi dati. 

Vi è la possibilità di inizializzare più di una classe dello stesso tipo, utile nel caso si vogliano fare confronti fra setup diversi.

La variabile **AWAvector** definisce gli angoli di vento apparente iniziali e finali, così come lo step. Si consiglia di non modificare.

**TWS** è la velocità del vento reale sul campo di regata, in nodi.

Vanno infine definite tutte le polari che si vogliono calcolare tramite l'inizializzazione delle classi vpp a cui vanno date in input gli oggetti inizializzati in precedenza, nell'ordine:
- Barca
- Campo di regata
- Equipaggio
- Vele
- Deriva
- Foil
- Vettore vento apparente **AWAvector**
- Vento reale **TWS**

### Processing
Nella fase di processing è sufficiente inserire dentro la variabile **vpplist** le classi Vpp() inizializzate (che può anche essere una sola). 

### Post-processing
Per vizualizzare correttamente i grafici finali, inserire nella variabile **legend** una stringa per ogni Vpp() inizializzata (ad esempio "No Foil", "Foil wortmann", "Deriva grande", ecc.)
