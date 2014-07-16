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
#include <set>
#include <assert.h>

/*******************************************************/
/* ARTIFICIAL BEE COLONY OPTIMIZATION : Binary Case :D */
/*******************************************************/

#define NP 60  /* Numero de abejas ( poblacion = exploradoras + supervisoras ) */
#define FOOD_NUMBER NP/2 /*Numero de fuentes de alimento, la mitad de la poblacion */
#define MAXD 630
#define LIMIT MAXD*0.50*FOOD_NUMBER  /* (dim * food_number) una fuente de alimento, que no es mejorada en "limite" iteraciones es abandonada por su abeja empleada*/
#define maxCycle 2000 /* numero de ciclos de recoleccion ( criterio de parada */
/* Variables del problema especifico*/
#define MU_MAX 0.9
#define MU_MIN 0.005
//#define DP 20 /*numero de parametros del problema a optimizar*/
#define LB 0.0 /* limite interior de los parametros */
#define UB 1.0 /* limite superior de los parametros */
#define P_LOCAL 0.1 // default: 0.01
#define N_LOCAL 40

#define runtime 30  /* numero de veces que puede ser ejecutado el algoritmo*/

inline bool one(double value) {
	return value > 0.5;
}
class ObjectiveFunction {
public:
	// dirty trick : ... { return -1; };
	virtual double perform ( double *v ) = 0;
};
typedef ObjectiveFunction* FunctionCallbackPtr;

/* DistABC */
class BinaryABCOptimizacion {
private:
	double Foods[FOOD_NUMBER][MAXD]; /*Fuentes de alimento, cada fuente de alimento contiene una posible solucion que esta siendo optimizada*/
	double f[FOOD_NUMBER];  /*f es un arreglo que contiene los valores de la funcion objetivo a ser optimizada*/
	double fitness[FOOD_NUMBER]; /*Este arreglo fitnes contiene el valor de aptitud de cada fuente de alimento*/
	double trial[FOOD_NUMBER]; /*Este arreglo contiene el numero de pruebas realizadas en cada fuente de alimento, si pasa el limite la abeja empleada lo abandona*/
	double prob[FOOD_NUMBER]; /*Este arreglo contiene las probabilidades de cada solucion de ser escogida*/
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
	int cycle_iter_;
public:
	BinaryABCOptimizacion( int dim, FunctionCallbackPtr _function ) {
		function = _function;
		//assert ( dim < MAXD );
		printf("ABC dim: %d\n", dim);
		dim_ = dim;
	}
	double* getOptimalParams() { return GlobalParams; }
	void setObjectiveFunction( FunctionCallbackPtr _function ) {
		function = _function;
	}
	double jackard(int i, int j) {
		int m01, m10, m11;
		m01 = m10 = m11 = 0;
		for(int idx = 0; idx < dim_; ++ idx) {
			if(!one(Foods[i][idx]) && one(Foods[j][idx]) )
				++ m01;
			if(one(Foods[i][idx]) && !one(Foods[j][idx]))
				++ m10;
			if(one(Foods[i][idx]) && one(Foods[j][idx]))
				++ m11;
		}
		return 1.0 - m11 * 1.0 / ( m01 + m10 + m11 );
	}
	void newSolution(double r, int i_, int nb, double* new_solution) {
		double A = r*jackard(i_, nb);
		int m01, m10, m11;
		int n1, n0; n1 = n0 = 0;
		int idx1[MAXD];
		int idx0[MAXD];
		for(int j = 0; j < dim_; ++ j)
			if(one(Foods[i_][j]))
				idx1[n1++] = j;
			else
				idx0[n0++] = j;
		//printf("n0 %d, n1 %d\n", n0, n1);
		double mind = 1e10, dist;
		int bm01, bm10, bm11;
		for(m01 = 0; m01 <= n1; ++ m01) {
			//printf("it\n");
			m11 = n1 - m01;
			for(m10 = 0; m10 <= n0; ++ m10) {
				dist = fabs(1 - m11 * 1.0 / ( m01 + m10 + m11 ) - A );
				//printf("%d %d %f\n",m01, m10, dist);
				if ( dist < mind ) {
					mind = dist;
					bm01 = m01;
					bm10 = m10;
					bm11 = m11;
				}
			}
		}
		for(int j = 0; j < dim_; ++ j)
			new_solution[j] = 0;
		//printf("m11 m10 %d %d\n", bm11, bm10);

		int pos = -1, t;
		for(int j = 0; j < bm11; ++ j) {
			pos = rand()%(n1-j);
			if(pos != n1-j-1) {
				t = idx1[n1-j-1];
				idx1[n1-j-1] = idx1[pos];
				idx1[pos] = t;
			}
		}
		for(int j = n1-bm11; j < n1; ++ j)
			new_solution[idx1[j]] = 1;

		for(int j = 0; j < bm10; ++ j) {
			pos = rand()%(n0-j);
			if(pos != n0-j-1) {
				t = idx0[n0-j-1];
				idx0[n0-j-1] = idx0[pos];
				idx0[pos] = t;
			}
		}
		for(int j = n0-bm10; j < n0; ++ j)
			new_solution[idx0[j]] = 1;

		/*TEST
		printf("Old: ");
		for(int j = 0; j < dim_; ++ j)
			printf("%d",(int)Foods[i_][j]);
		printf("\nNew: ");
		for(int j = 0; j < dim_; ++ j)
			printf("%d",(int)Foods[i_][j]);
		printf("\n");*/
	}
	void localSearch(double* solution, double *newSolution) {
		int idx1[MAXD];
		int idx0[MAXD];
		int n0, n1; n0 = n1 = 0;
		for(int i = 0; i < dim_; ++i) {
			if(one(solution[i]))
				idx1[n1++] = i;
			else
				idx0[n0++] = i;
			newSolution[i] = solution[i];
		}
		int p1 = idx1[rand()%n1];
		int p0 = idx0[rand()%n0];
		double t = newSolution[p1];
		newSolution[p1] = newSolution[p0];
		newSolution[p0] = t;
	}
	void localSearchModule() {
		int i,j, t;
		t = 0;
		i = 0;
		r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
		if(r > P_LOCAL)
			return;
	  /*fase de las abejas empleadas*/
		while(t < N_LOCAL) {
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			if(r<prob[i]) {/*Escoger una fuende de alimento de acuerdo a su probabilidad de ser escogida*/
				t++;
				localSearch(Foods[i], solution);
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
			if (i==FOOD_NUMBER)
				i=0;
		}
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
		for(i=0;i<FOOD_NUMBER;i++)	{
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
		for(i=0;i<FOOD_NUMBER;i++) {
			inicializarFuente(i);
		}
		GlobalMin=f[0];
		for(i=0;i<dim_;i++)
			GlobalParams[i]=Foods[0][i];
	}
	void enviarAbejasEmpleadas() {
	  int i,j;
	  /*fase de las abejas empleadas*/
		for (i=0;i<FOOD_NUMBER;i++) {
			/* El parametro a ser cambiado es determinado aleatoriamente */
			r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			param2change=(int)(r*dim_);

			/* Se selecciona aleatoriamente un vecino para obtener una solucion mutante de i */
			r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			neighbour=(int)(r*FOOD_NUMBER);

			/* La solucion aleatoriamente escojida debe ser diferente de i*/
			while(neighbour==i)	{
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				neighbour=(int)(r*FOOD_NUMBER);
			}
			for(j=0;j<dim_;j++)
				solution[j]=Foods[i][j];

			double base = 1.005; // 0
			/* Se genera una nueva solucion: v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			double mu = MU_MAX -  (MU_MAX - MU_MIN) * (1 - pow(base, -1 * cycle_iter_));
			//solution[param2change]= 1 - Foods[i][param2change];//Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2; //1 - solution[param2change];
			double var = (1 - pow(base, -1 * (cycle_iter_+1))) - (1 - pow(base, -1 * cycle_iter_));
			var = var*(r-0.5);
			if(i==0) printf("mu + var: %lf %lf\n", mu, var);
			newSolution(mu + var, i, neighbour, solution);
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
	void calcularProbabilidadesMaximo()	{
		int i;
		double maxfit;
		maxfit=fitness[0];
		for (i=1;i<FOOD_NUMBER;i++) {
			if (fitness[i]>maxfit)
			maxfit=fitness[i];
		}
	 	for (i=0;i<FOOD_NUMBER;i++) {
			prob[i]=(0.9*(fitness[i]/maxfit))+0.1;
			printf("%.2f ", prob[i]);
		}
	 	printf("\n");
	}
	void calcularProbabilidadesSuma()	{
		int i;
		double sumfit;
		sumfit=fitness[0];
		for (i=1;i<FOOD_NUMBER;i++) {
			sumfit += fitness[i];
		}
		for (i=0;i<FOOD_NUMBER;i++) {
			prob[i]=(0.9*(fitness[i]/sumfit))+0.1;
			printf("%.2f ", prob[i]);
		}
		printf("\n");
	}

	void enviarAbejasSupervisoras() {
	  int i,j,t;
	  i=0;
	  t=0;
	  /* Fase de abejas supervisoras */
	  while(t<FOOD_NUMBER)	{
			r = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			if(r<prob[i]) {/*Escoger una fuende de alimento de acuerdo a su probabilidad de ser escogida*/
				t++;
				/* Se escoge el parametro a se modificado aleatoriamente */
				r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				param2change=(int)(r*dim_);

				/* Se selecciona aleatoriamente una solucion a ser usada para producir una solucion mutante de la solucion i*/
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				neighbour=(int)(r*FOOD_NUMBER);

				/* La solucion aleatoriamente seleccionada debe ser diferente de la solucion i*/
				while(neighbour==i)	{
					r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
					neighbour=(int)(r*FOOD_NUMBER);
				}
				for(j=0;j<dim_;j++)
					solution[j]=Foods[i][j];

				double base = 1.005;
				/*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
				r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
				//double mu = MU_MAX - (MU_MAX - MU_MIN) * (log2(1 + (maxCycle - cycle_iter_) * 1.0 /maxCycle));
				double mu = MU_MAX -  (MU_MAX - MU_MIN) * (1 - pow(base, -1 * cycle_iter_));
				//solution[param2change]= 1 - Foods[i][param2change];//Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2; //1 - solution[param2change];
				double var = (1 - pow(base, -1 * (cycle_iter_+1))) - (1 - pow(base, -1 * cycle_iter_));
				var = var*(r-0.5);
				if(t==0)
					printf(" 2. mu + var: %lf %lf\n", mu, var);
				newSolution(mu + var, i, neighbour, solution);

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
			if (i==FOOD_NUMBER)i=0;
		}
	}
	/* Determinar que fuentes de alimento exeden el limite de pruebas*/
	/* A comparative study ... Karagoba suggest 5% - 10% of scouts */
	void enviarAbejasExploradoras()	{
		int num_scouts = 0.1 * NP, maxtrialindex,i,s;
		if (num_scouts == 0)
			num_scouts = 1;
		for(s=0;s<num_scouts;++s) {
			maxtrialindex=0;
			for (i=1;i<FOOD_NUMBER;i++) {
				if (trial[i]>trial[maxtrialindex])
					maxtrialindex=i;
			}
			if(trial[maxtrialindex]>=LIMIT)	{
				inicializarFuente(maxtrialindex);
			}
		}
	}
	double ejecutar() {
		printf ( "ABC begining to execute\n");
		int iter, run,j;
		double mean, temp;
		mean=0;
		srand(time(NULL));
		for(run=0;run<runtime;run++)		{
			inicializar();
			memorizarMejoresFuentes();
			for (iter=0;iter<maxCycle;iter++)	{
				cycle_iter_ = iter;
				enviarAbejasEmpleadas();
				//calcularProbabilidadesSuma();
				calcularProbabilidadesMaximo();
				enviarAbejasSupervisoras();
				localSearchModule();
				memorizarMejoresFuentes();
				enviarAbejasExploradoras();
				printf ( "Best solution in cycle %d: %f\n", iter, -1*(log2(GlobalMin)/100.0 - 1));
				for ( int i = 0; i < dim_; ++ i ) {
					float some = (float)GlobalParams[i];
					printf ( "%d", (int)some );
				}
				printf ( "\n" );
			}
			printf ( "Global sol\n" );
			for ( int i = 0; i < dim_; ++ i ) {
				float some = (float)GlobalParams[i];
				printf ( "%d ",  (int)some );
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
