#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX 2
#define ERRO 1e-6

void pontilhado() {
	printf("\n---------------------------------------\n");
}

//----------------------Funções que compõem o sistema----------------
double f1(double x[MAX]) {
	return(log(pow(x[0],2)+pow(x[1],2))-sin(x[0]*x[1])-log(2)-log(M_PI));
}


double f2(double x[MAX]) {
	return(exp(x[0]-x[1])+cos(x[0]*x[1]));
}
//----------------------Fim das funções que compõem o sistema -------------

// --------------------Funções que compõem o jacobiano ---------------------
double f1_x1(double x[MAX]) {
	return(2*x[0]/(pow(x[0],2)+pow(x[1],2))-x[1]*cos(x[0]*x[1]));
}

double f1_x2(double x[MAX]) {
	return(2*x[1]/(pow(x[0],2)+pow(x[1],2))-x[0]*cos(x[0]*x[1]));
}

double f2_x1(double x[MAX]) {
	return(exp(x[0]-x[1])-x[1]*sin(x[0]*x[1]));
}

double f2_x2(double x[MAX]) {
	return(-exp(x[0]-x[1])-x[0]*sin(x[0]*x[1]));
}
//--------------------------Fim das funções que compõem o jacobiano-------------------

void imprimeResultado(double *raizes, int dim) { //Apenas imprime resultado do sistema
	int i;
	
	pontilhado();
	printf("Resultado:");
	pontilhado();
	for(i = 0; i < dim; i++)
        	printf("x_%d = %.3lf\n",i + 1,raizes[i]);
}

double maxElementoVetor(double * v, int n) { // Retorna o elemento máximo de um vetor
	int i;
	double maior = v[0];

	for(i = 0; i < n;i++) {
		if(maior < v[i])
			maior = v[i];
	}

    return(maior);
}

double **multPorEscalar(double num, double **matriz, int linha,int coluna) { // Multiplica os elementos de uma matriz por um escalar
	int i,j;								// e retorna a nova matriz
	double **A = (double**)malloc(linha*sizeof(double*));
	
	for(i = 0; i < linha; i++)
		A[i] = (double*)malloc(coluna*sizeof(double));

	for(i = 0; i < linha; i++)
		for( j = 0; j < coluna; j++)
			A[i][j] = num * matriz[i][j];
	return(A);
}

void pivoteamento(double **M, int k, int dim) { // Uma função auxiliar da subsreversa, ajuda a resolver um sistema linear
	int i;
	double maior = M[k][k], *aux_linha; // a variável maior guarda o elemento da diagonal principal que a principio utilizaríamos
	// para zerar os elementos da mesma coluna nas linhas inferiores, a linha auxiliar será utilizada na hora de trocar as linhas da
	// matriz, caso o maior seja superado por outro elemento abaixo de sua coluna
	
	for(i = k; i < dim; i++) {
		if(fabs(maior) < fabs(M[i][k])) {
			maior = M[i][k];
			
			aux_linha = M[i];
			M[i] = M[k];
			M[k] = aux_linha;
		}
	}
}
  			
void triangsup(double **M, int dim) { //Triangulariza uma matriz
        int i, j, k=0;
        double lambda;
        
        for(k = 0; k < (dim-1); k++) // k é a linha que nos estamos trabalhando
        	for(i = (k+1); i < dim; i++) { // i são as linhas abaixo de k onde iremos operar
        		pivoteamento(M,k, dim);
        		
        		(M[k][k] != 0)? (lambda = (-M[i][k]/M[k][k])) : (lambda = 0); // lambda é a constante que iremos
        		// multiplicar os elementos da linha k e somar com os respectivos elementos da linha i para o escalonamento
        	
        		for(j = 0; j <= dim; j++) // aqui faz a operação na linha i
        			M[i][j] += lambda*M[k][j];  
              	} 
}

double *subsreversa(double **M, int dim) { // Resolve um sistema linear, substituição reversa
	int i, j;
	double somatorio=0, *raizes;
	
	raizes = (double*)malloc(2*sizeof(double));
	
	for(i = (dim-1); i >= 0; i--) {
		for(j = (i+1); j < dim; j++)
			somatorio += M[i][j]*raizes[j];
		raizes[i] = (M[i][dim] - somatorio)/M[i][i];

		somatorio = 0;
	}
	
	return(raizes);
}

void *Broyden(double (*equacao[MAX])(), double x[MAX]) {
	int i,j,c = 0;
	double F_x[MAX], *y, **jacob_aum; // matriz jacobina aumentada
	
    	jacob_aum = (double**)malloc(MAX * sizeof(double*));

	for(i = 0; i < MAX; i++)
		jacob_aum[i] = (double*)malloc((MAX + 1) * sizeof(double));

	do {
		for(i = 0; i < MAX; i++)
			F_x[i] = -equacao[i](x);
		
		//------------ Calcula o Jacobiano e insere elementos de F(x)
		jacob_aum[0][0] = f1_x1(x);
		jacob_aum[0][1] = f1_x2(x);
		jacob_aum[0][2] = F_x[0];
		jacob_aum[1][0] = f2_x1(x);
		jacob_aum[1][1] = f2_x2(x);
		jacob_aum[1][2] = F_x[1];
		//-------------Resolve o sistema composto pelo jacobiano aumentado	
		
		

		triangsup(jacob_aum,MAX);
		y = subsreversa(jacob_aum,MAX);
		
		for(i = 0; i < MAX; i++)
			x[i] += y[i];


	} while(maxElementoVetor(y,MAX) > ERRO);
	
	imprimeResultado(x,MAX); //Apresenta o resultado final na tela
}

int main() {
	double (*funcoes[MAX])() = {f1,f2}; 
	double x1[2]; // valor inicial
	
	x1[0] = x1[1] = 2;

	Broyden(funcoes,x1);

	return(0);
}
