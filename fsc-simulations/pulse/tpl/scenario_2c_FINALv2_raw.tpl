//Number of population samples (demes)
5
//Population effective sizes (number of genes) (haploid!)
NNK //NK 0
NSK //SK 1
NWP //WP 2
NEP //EP 3
NNP //NP 4
//Sample sizes, sampling time, average inbreeding level of the samples
10 0 0
10 0 0
10 0 0
10 0 0
10 0 0
//Growth rates	: negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time (in generations), source, sink, migrants, new size, new growth rate, migr. matrix 
24  historical event
TADNKSK 0 1 RNKSK 1 0 0 
TADNKSK 1 0 RSKNK 1 0 0
TADNKNP 0 4 RNKNP 1 0 0 
TADNKNP 4 0 RNPNK 1 0 0
TADSKNP 1 4 RSKNP 1 0 0 
TADSKNP 4 1 RNPSK 1 0 0
TADNPWP 4 2 RNPWP 1 0 0 
TADNPWP 2 4 RWPNP 1 0 0
TADNPEP 4 3 RNPEP 1 0 0 
TADNPEP 3 4 REPNP 1 0 0
TADWPEP 2 3 RWPEP 1 0 0 
TADWPEP 3 2 REPWP 1 0 0
TK1 1 0 1 FANK 0 0
TADAKNP 0 4 RAK1NP 1 0 0 
TADAKNP 4 0 RNP1AK 1 0 0
TKNP 4 0 1 FAN1KNP 0 0
TADWPAKNP 2 0 RWPAKNP 1 0 0
TADWPAKNP 0 2 RAKNPWP 1 0 0
TADEPAKNP 3 0 REPAKNP 1 0 0
TADEPAKNP 0 3 RAKNPEP 1 0 0
TP1 3 2 1 FANP 0 0
TADAPAKNP 0 2 R1AKNPAP 1 0 0
TADAPAKNP 2 0 RAP1AKNP 1 0 0
TKPNP 2 0 1 FAN1KPNP 0 0
//Number of independent loci [chromosome], flag to say whether they have same structure (0) or not (1) 
100 0
//Per chromosome: Number of linkage blocks (blocks might differ in terms of type of markers, recombination rate, or mutation rate)
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 1000000 1.0E-08 1.25E-08 0.33
