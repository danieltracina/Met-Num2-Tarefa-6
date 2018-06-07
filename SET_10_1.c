#include <stdio.h>
#include <math.h>
#define MAX 3
#define ERRO 1e-5

//--------------- Funções que compõem o sistema------------------------
double f1(double v[MAX]) {
	return(pow(v[0],2)+2*pow(v[1],2)-v[1]-2*v[2]);
}

double f2(double v[MAX]){
	return(pow(v[0],2)-8*pow(v[1],2)+10*v[2]);
}

double f3(double v[MAX]) {
	return(pow(v[0],2)/(7.0*v[1]*v[2])-1);
}

//------------------------------Fim--------------------------------

double norma(double *v,int dim) { // calcula norma euclidiana de um vetor
	int i;
	double sum = 0;
	for(i = 0; i < dim; i++)
		sum += pow(v[i],2);
	sum = sqrt(sum);

	return sum;
}

double calcTol(double norm, double norma_vetor) { // calcula o Tol
	return(fabs(norm - norma_vetor)/norm);
}

int main()
{
	double x_0[MAX] = {1,1,1}, tol = 1, norm, norma_vetor; 
	double (*equacao[MAX])() = {f1,f2,f3};
	int i;

	while(tol > ERRO) { 
		
		for(i = 0; i < MAX;i++)
		{
			norma_vetor = norma(x_0,MAX); 

			x_0[i] = equacao[i](x_0);

			norm = norma(x_0,MAX);
		}
		
		tol = calcTol(norm,norma_vetor);

	}

	printf("Resultado:\n");
	for(i = 0; i < MAX;i++)
		printf("x_%d =  %lf\n",i + 1, x_0[i]);
		
	return(0);

}
