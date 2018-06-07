#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//------------------- Funções que compõem o sistema -----------------------
double f1(double x1, double x2, double x3) {
	return(3*x1-cos(x2*x3)-0.5);
}

double f2(double x1, double x2, double x3) {
	return(4*pow(x1,2)-625*pow(x2,2)+2*x2-1);
}

double f3(double x1, double x2, double x3) {
	return(exp(-x1*x2)+20*x3+(10*M_PI-3)/3);
}

//-------------------Fim das funções que compõem o sistema --------------------

//------------------Funções que compõem o Jacobiano ----------------------------
double f1_x1(double x1, double x2, double x3) {
	return(3);
}

double f1_x2(double x1, double x2, double x3) {
	return(sin(x2*x3)*x3);
}

double f1_x3(double x1, double x2, double x3) {
	return(sin(x2*x3)*x2);
}

double f2_x1(double x1, double x2, double x3) {
	return(8*x1);
}

double f2_x2(double x1, double x2, double x3) {
	return(-1250*x2+2);
}

double f2_x3(double x1, double x2, double x3) {
	return(0);
}

double f3_x1(double x1, double x2, double x3) {
	return(-x2*exp(-x1*x2));
}

double f3_x2(double x1, double x2, double x3) {
	return(-x1*exp(-x1*x2));
}

double f3_x3(double x1, double x2, double x3) {
	return(20);
}

// ----------------------------Fim das funções que compõem o Jacobiano ------------

 
double **inicia_J(double x1, double x2, double x3) { //Calcula o jacobiano
	int i;
	double **J;
	
	J = (double**)malloc(3*sizeof(double*));
	for(i = 0; i < 3; i++)
		J[i] = (double*)malloc(3*sizeof(double));
	
	J[0][0] = f1_x1(x1, x2, x3);
	J[0][1] = f1_x2(x1, x2, x3);
	J[0][2] = f1_x3(x1, x2, x3);
	J[1][0] = f2_x1(x1, x2, x3);
	J[1][1] = f2_x2(x1, x2, x3);
	J[1][2] = f2_x3(x1, x2, x3);
	J[2][0] = f3_x1(x1, x2, x3);
	J[2][1] = f3_x2(x1, x2, x3);
	J[2][2] = f3_x3(x1, x2, x3);
	
	return(J);
}


double **inicia_F(double x1, double x2, double x3) { // calcula o F(x)
	int i;
	double **F;
	
	F = (double**)malloc(3*sizeof(double*));
	for(i = 0; i < 3; i++)
		F[i] = (double*)malloc(sizeof(double));
	
	F[0][0] = f1(x1, x2, x3);
	F[1][0] = f2(x1, x2, x3);
	F[2][0] = f3(x1, x2, x3);
	
	return(F);
}

void geraArquivo(double **J, double **F) { // Gera um arquivo com a matriz expandida para que outra função leia e resolva
					// o sistema linear
	int i,j;
	FILE *fp = fopen("matriz.txt", "w");
	
	for(i = 0; i < 3; i++)
		fprintf(fp,"%lf %lf %lf %lf\n", J[i][0], J[i][1], J[i][2], -F[i][0]);
	
	fclose(fp);
}
 	
double **ler(char *nome_arqu_matriz, int dim) { // Ler matriz expandida do arquivo e guarda na matriz M
        int i, j;
        double **M, **aux, num;
        FILE *fp = fopen(nome_arqu_matriz, "r");
        
        M = malloc(dim*sizeof(double*)); // aloca memoria para as linhas da matriz
        aux = M;       
       
        for(i = 0; i < dim; i++) {
        	*aux = malloc((dim+1)*sizeof(double)); // aloca memoria para as colunas da linha i da matriz
        	
        	for(j = 0; j <= dim; j++)
        		fscanf(fp,"%lf", &M[i][j]);
        	        	
        	aux++;
        }	
        
       fclose(fp);
       return(M);
}    

void triangsup(double **M, int dim) { // Triangulariza matrizes
        int i, j, k=0;
        double lambda;
        
        for(k = 0; k < (dim-1); k++) // k é a linha que nos estamos trabalhando
        	for(i = (k+1); i < dim; i++) { // i são as linhas abaixo de k onde iremos operar
        		(M[k][k] != 0)? (lambda = (-M[i][k]/M[k][k])) : (lambda = 0); // lambda é a constante que iremos
        		// multiplicar os elementos da linha k e somar com os respectivos elementos da linha i para o escalonamento
        	
        		for(j = 0; j <= dim; j++) // aqui faz a operação na linha i
        			M[i][j] += lambda*M[k][j];  
              	} 
}


void subsreversa(double **M,double *raizes, int dim) { // substituição reversa
	int i, j;
	double somatorio=0;
	
	for(i = (dim-1); i >= 0; i--) {
		for(j = (i+1); j < dim; j++)
			somatorio += M[i][j]*raizes[j];
		raizes[i] = (M[i][dim] - somatorio)/M[i][i];

		somatorio = 0;
	}
}

void pontilhado() { // Função para visual estético
        printf("\n-----------------------------------\n");
}
void imprime(double **M, int dim) { // Imprime matrizes 
        int i, j;
        
        pontilhado(); // função que torna mais apresentável os resultados na tela
        for(i = 0; i < dim; i++) {
                for(j = 0; j < (dim+1); j++) {
                        printf("%.2lf\t", M[i][j]);
                }
                printf("\n");
        }
        pontilhado();
}

int main() {
	double **J, **F, **M, *s, *x_k;
	int i, k = 0, k_max = 2;
	
	s = (double*)malloc(3*sizeof(double));
	x_k = (double*)malloc(3*sizeof(double));
	
	for(i = 0; i < 3; i++)
		x_k[i] = 0;
		
	while(k < k_max) { // Algoritmo de Newton
		k++;
			
		J = inicia_J(x_k[0],x_k[1],x_k[2]);
		F = inicia_F(x_k[0],x_k[1],x_k[2]);
		
		geraArquivo(J,F);
		M = ler("matriz.txt", 3);
		triangsup(M,3);
		subsreversa(M,s,3);
		
		for(i = 0; i < 3; i++) 
			x_k[i] += s[i];
	}
	
	
	printf("\nSolução:\n");
	for(i = 0; i < 3; i++)
		printf("x_%d: %lf\n", i, x_k[i]); 
	 
	return(0);
}	
		
		
		