//Parameters for the coalescence simulation program :
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP0$
NPOP1$
//Samples sizes and samples age
34
76
//Growth rates  : negative growth implies population expansion
R0$
R1$
//Number of migration matrices
4
//Migration matrix 0
0.0000 Mix01$
Mix10$ 0.0000
//Migration matrix 1
0.0000 0.0000
0.0000 0.0000
//Migration matrix 2
0.0000 Mig01$
Mig10$ 0.0000
//Migration matrix 3
0.0000 0.0000
0.0000 0.0000
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
5 historical event
TDIV0$ 0 0 0 1 0 1
TDIV0$ 1 1 0 1 0 1
TDIV1$ 0 0 0 RES0$ 0 2
TDIV1$ 1 1 0 RES1$ 0 2
TDIV2$ 0 1 1 RES2$ 0 3
//Number of independent loci [chromosome]
1 0
//Per chromosome:
1
//per Block: data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1.918e-8 OUTEXP