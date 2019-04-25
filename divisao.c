#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double h;
int n;

struct inf{
	int ini;
	int fim;
};



int main (int argc, char *argv[]) {



	int numprocs, proc, ini, fim, res;
	double pi, *sum, sum_tot;


	n = 100;
	numprocs = 7;


	struct inf xh[numprocs];



	printf("\n");
	printf("Numero de processos: %d", numprocs);
	printf("\n");

	ini = 0;
	fim = 0;
	res = n;

	for(proc=1; proc<=numprocs; proc++){

		ini  = fim + 1;
		fim += (int)n / numprocs;

		res -= (int)n / numprocs;

		if(proc == numprocs){
			fim += res;
		}

		xh[proc].ini = ini;
		xh[proc].fim = fim;

		printf("\n");
		printf("ini: %d   fim: %d",xh[proc].ini,xh[proc].fim);

	}



	return(0);
 }
