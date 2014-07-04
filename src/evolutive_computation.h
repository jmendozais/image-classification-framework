/*
 * evolutive_computation.h
 *
 *  Created on: Jan 19, 2013
 *      Author: jmendoza
 */

#ifndef EVOLUTIVE_COMPUTATION_H_
#define EVOLUTIVE_COMPUTATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*******************************************************/
/* ARTIFICIAL BEE COLONY OPTIMIZATION : Binary Case :D */
/*******************************************************/

#define NP 40  /* Numero de abejas ( poblacion = exploradoras + supervisoras ) */
#define FoodNumber NP/2 /*Numero de fuentes de alimento, la mitad de la poblacion */
#define limit 25  /* una fuente de alimento, que no es mejorada en "limite" iteraciones es abandonada por su abeja empleada*/
#define maxCycle 3000 /* numero de ciclos de recoleccion ( criterio de parada */
#define MAXD 500
/* Variables del problema especifico*/
#define D 20 /*numero de parametros del problema a optimizar*/
#define lb -5.12 /* limite interior de los parametros */
#define ub 5.12 /* limite superior de los parametros */


#define runtime 30  /* numero de veces que puede ser ejecutado el algoritmo*/

class ObjectiveFunction {
public:
	// dirty trick : ... { return -1; };
	virtual double perform ( double *v ) = 0;
};
typedef ObjectiveFunction* FunctionCallbackPtr;


class BinaryABCOptimizacion {
private:
	double Foods[FoodNumber][MAXD]; /*Fuentes de alimento, cada fuente de alimento contiene una posible solucion que esta siendo optimizada*/
	double f[FoodNumber];  /*f es un arreglo que contiene los valores de la funcion objetivo a ser optimizada*/
	double fitness[FoodNumber]; /*Este arreglo fitnes contiene el valor de aptitud de cada fuente de alimento*/
	double trial[FoodNumber]; /*Este arreglo contiene el numero de pruebas realizadas en cada fuente de alimento, si pasa el limite la abeja empleada lo abandona*/
	double prob[FoodNumber]; /*Este arreglo contiene las probabilidades de cada solucion de ser escogida*/
	double solution [MAXD]; /*Nueva solucion*/
	double ObjValSol; /*Valor objetivo */
	double FitnessSol; /*Valor de aptitud de la nueva solucion*/
	int neighbour, param2change; /*param2change corrresponde a j, neighbour corresponde a k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
	double GlobalMin; /*Solucion optima obtenida por el algoritmo ABC*/
	double GlobalParams[MAXD]; /*Contiene los parametros de la solucion optima*/
	double GlobalMins[runtime]; /*Contiene el minimo global de cada iteracion*/
	FunctionCallbackPtr function; /* referencia a la funcion objetivo */
	double r; /* numero aleatorio entre [0,1)*/

	int dim_;
public:
	BinaryABCOptimizacion( int dim, FunctionCallbackPtr _function ) {
		function = _function;
		dim_ = dim;
	}
	double* getOptimalParams() { return GlobalParams; }
	void setObjectiveFunction( FunctionCallbackPtr _function ) {
		function = _function;
	}
	double calcularFitness(double fun) {
		double result=0;
		if(fun>=0)
			result=1/(fun+1);
		else
			result=1+fabs(fun);
		return result;
	}
	void memorizarMejoresFuentes() {
	  int i,j;
		for(i=0;i<FoodNumber;i++)	{
			//print_vector( Foods[i] );
			if (f[i]<GlobalMin)	{
				GlobalMin=f[i];
				for(j=0;j<dim_;j++)
					GlobalParams[j]=Foods[i][j];
			}
		}
	}
	void inicializarFuente(int index)	{
	   int j;
	   for (j=0;j<dim_;j++) {
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			Foods[index][j]= round(r);//r*(ub-lb)+lb;
			solution[j]=Foods[index][j];
		}
		f[index]=function->perform(solution);
		fitness[index]=calcularFitness(f[index]);
		trial[index]=0;
	}
	/* Inicializar todas las fuentes de alimento */
	void inicializar() {
		int i;
		for(i=0;i<FoodNumber;i++) {
			inicializarFuente(i);
		}
		GlobalMin=f[0];
		for(i=0;i<dim_;i++)
			GlobalParams[i]=Foods[0][i];
	}
	void enviarAbejasEmpleadas() {
	  int i,j;
	  /*fase de las abejas empleadas*/
		for (i=0;i<FoodNumber;i++) {
			/* El parametro a ser cambiado es determinado aleatoriamente */
			r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			param2change=(int)(r*dim_);

			/* Se selecciona aleatoriamente un vecino para obtener una solucion mutante de i */
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			neighbour=(int)(r*FoodNumber);

			/* La solucion aleatoriamente escojida debe ser diferente de i*/
			while(neighbour==i)	{
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				neighbour=(int)(r*FoodNumber);
			}
			for(j=0;j<dim_;j++)
				solution[j]=Foods[i][j];

			/* Se genera una nueva solucion: v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			solution[param2change]= 1 - solution[param2change]; //Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

			/* Si el valor de parametro genrado esta fuera de los limites se corrige*/
			ObjValSol=function->perform(solution);
			FitnessSol=calcularFitness(ObjValSol);

			/* Se selecciona la solucion con mejor valor de aptitud ( criterio voraz ) */
			if (FitnessSol>fitness[i]) {
				/*Si la solucion mutante es mejor que la solucion actual, se reemplaza y se resetea el contador de prueba a 0*/
				trial[i]=0;
				for(j=0;j<dim_;j++)
					Foods[i][j]=solution[j];
				f[i]=ObjValSol;
				fitness[i]=FitnessSol;
			}	else {   /*Si la solucion no puede ser mejorada se incrementa el contador de prueba*/
				trial[i]=trial[i]+1;
			}
		}
	}

	/* Una fuente de alimente es escogida con una probabilidad proporcional a su valor de aptitud */
	void calcularProbabilidades()	{
		int i;
		double maxfit;
		maxfit=fitness[0];
		for (i=1;i<FoodNumber;i++) {
			if (fitness[i]>maxfit)
			maxfit=fitness[i];
		}
	 	for (i=0;i<FoodNumber;i++) {
			prob[i]=(0.9*(fitness[i]/maxfit))+0.1;
		}
	}

	void enviarAbejasSupervisoras() {
	  int i,j,t;
	  i=0;
	  t=0;
	  /* Fase de abejas supervisoras */
	  while(t<FoodNumber)	{
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			if(r<prob[i]) {/*Escoger una fuende de alimento de acuerdo a su probabilidad de ser escogida*/
				t++;

				/* Se escoge el parametro a se modificado aleatoriamente */
				r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				param2change=(int)(r*dim_);

				/* Se selecciona aleatoriamente una solucion a ser usada para producir una solucion mutante de la solucion i*/
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				neighbour=(int)(r*FoodNumber);

				/* La solucion aleatoriamente seleccionada debe ser diferente de la solucion i*/
				while(neighbour==i)	{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					neighbour=(int)(r*FoodNumber);
				}
				for(j=0;j<dim_;j++)
					solution[j]=Foods[i][j];

				/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				solution[param2change]= 1 - Foods[i][param2change];//+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

				/* Si el parametro generado esta fuera de los limites se corrige */
				ObjValSol=function->perform(solution);
				FitnessSol=calcularFitness(ObjValSol);

				/* Se selecciona la solucion con mejor aptitud ( criterio voraz ) */
				if (FitnessSol>fitness[i]) {
				/* Si la solucion mutante es mejor que la solucion actual, se resetea el contador de pruebas*/
					trial[i]=0;
					for(j=0;j<dim_;j++)
						Foods[i][j]=solution[j];
					f[i]=ObjValSol;
					fitness[i]=FitnessSol;
				}	else {   /* Si la solucion i no es mejorada se incrementa */
					trial[i]=trial[i]+1;
				}
			} /* Si se ha iterado sobre todas las fuentes sin completar el alimento */
			i++;
			if (i==FoodNumber)i=0;
		}
	}
	/* Determinar que fuentes de alimento exeden el limite de pruebas*/
	void enviarAbejasExploradoras()	{
		int maxtrialindex,i;
		maxtrialindex=0;
		for (i=1;i<FoodNumber;i++) {
			if (trial[i]>trial[maxtrialindex])
				maxtrialindex=i;
		}
		if(trial[maxtrialindex]>=limit)	{
			inicializarFuente(maxtrialindex);
		}
	}
	double ejecutar() {
		printf ( "ABC begining to execute\n");
		int iter,run,j;
		double mean, temp;
		mean=0;
		srand(time(NULL));
		for(run=0;run<runtime;run++)		{
			inicializar();
			memorizarMejoresFuentes();
			for (iter=0;iter<maxCycle;iter++)	{
				enviarAbejasEmpleadas();
				calcularProbabilidades();
				enviarAbejasSupervisoras();
				memorizarMejoresFuentes();
				enviarAbejasExploradoras();
				printf ( "Best solution in cycle %d: %f\n", iter, GlobalMin );
				for ( int i = 0; i < dim_; ++ i ) {
					float some = (float)GlobalParams[i];
					printf ( "%f ", some );
				}
				printf ( "\n" );
			}
			printf ( "Global sol\n" );
			for ( int i = 0; i < dim_; ++ i ) {
				float some = (float)GlobalParams[i];
				printf ( "%f ", some );
			}
			printf ( "\n" );
			printf("%d. Iteracion: %f \n",run+1,GlobalMin );
			GlobalMins[run]=GlobalMin;
			mean +=  GlobalMin;
		}
		mean = mean/runtime;
		printf("Umbral promedio de %d iteraciones: %e\n", runtime, temp);
		return temp;
	}
};

#endif /* EVOLUTIVE_COMPUTATION_H_ */
