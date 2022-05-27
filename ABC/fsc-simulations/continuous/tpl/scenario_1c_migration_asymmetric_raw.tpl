//Number of population samples (demes)
5
//Population effective sizes (number of genes) (haploid!)
NNK //northern Khoe-San 0
NSK //southern Khoe-San 1
NWP //western RHG 2
NEP //eastern RHG 3
NNP //western-eastern RHG neigbours 4
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
5
// mutation matrix 0 -> most recent
0 MNKSK 0 0 MNKNP
MSKNK 0 0 0 MSKNP
0 0 0 MWPEP MWPNP
0 0 MEPWP 0 MEPNP
MNPNK MNPSK MNPWP MNPEP 0
// mutation matrix 1 -> before WP-EP split
0 MNKSK 0 0 MNKNP
MSKNK 0 0 0 MSKNP
0 0 0 0 MAP1NP
0 0 0 0 0
MNPNK MSKNP MNP1AP 0 0
// mutation matrix 2 -> before AP-NP split
0 MNKSK MNKAPNP 0 0
MSKNK 0 MSKAPNP 0 0
MAPNPNK MAPNPSK 0 0 0
0 0 0 0 0
0 0 0 0 0
// mutation matrix 3 -> before NK-SK split
0 0 MAKAPNP 0 0
0 0 0 0 0
MAPNPAK 0 0 0 0
0 0 0 0 0
0 0 0 0 0
// mutation matrix 4 -> before first split -> no migration
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
//historical event: time (in generations), source, sink, migrants, new size, new growth rate, migr. matrix 
4  historical event
TP1 3 2 1 FANP 0 1
TPNP 4 2 1 FAN1PNP 0 2
TK1 1 0 1 FANK 0 3
TKPNP 2 0 1 FAN1KPNP 0 4
//Number of independent loci [chromosome], flag to say whether they have same structure (0) or not (1) 
100 0
//Per chromosome: Number of linkage blocks (blocks might differ in terms of type of markers, recombination rate, or mutation rate)
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 2000000 1.0E-08 1.25E-08 0.33
