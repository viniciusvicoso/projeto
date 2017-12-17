//Os valores de xa foram {0.1, 0.1, -0.1}, {0.5, 0.5, -0.5}, {0.8, 0.8, -0.8}

#include<stdio.h>
#include<math.h>

#define T 20
#define N 3
#define h 1e-06
#define tol 1e-07


typedef double (*sistfunc)();


double f1(double x[N])
{
	return(pow(x[0],2) - 81*pow((x[1]+0.1), 2) + sin(x[2]) + 1.06);
}


double f2(double x[N])
{
	return(exp(x[0]*x[1]) + 20*x[2] + ((10*M_PI) - 3)/3.0);
}


double f3(double x[N])
{
	return((3*x[0]) - cos(x[1] *x[2]) - 1.0/2);
}


double df(double f(), double x[N], int k)
{
	int i;
	double aux[N], v;
	
	
	for(i = 0; i < N; i++)
		aux[i] = x[i];
	
	aux[k] = x[k] + h/2.0;
	
	v = f(aux);
	
	aux[k] = x[k] - h/2.0;	
	
	v = (v - f(aux))/h;
	
	return(v);
}


void imprime(double M[N][N+1])
{
	int i, j;
	
	
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N+1; j++)
			
			printf("%.2lf\t", M[i][j]);
		
		printf("\n");
	}
	
	printf("\n\n");
}


void diagonalizando(double M[20][20])
{
	double c;
	int i, j, k;
	
	
	for(j = 0; j < N; j++)
    {
        for(i = 0; i < N; i++)
        {
            if(i != j)
            {
                c = M[i][j]/M[j][j];
                 
                for(k = 0; k < (N+1); k++)
                
                    M[i][k] = M[i][k] - (c*M[j][k]);
	    }
	}
    }

}


void jacobiano(sistfunc equacao[], double M[N][N+1], double x[N], double xa[N], double mt[N][N+1])
{
	int i, j;
	
	
	for(i = 0; i < N; i++)
	{
		x[i] = -equacao[i](xa);
					
		for(j = 0; j < N; j++)
		
			M[i][j] = df(equacao[i], xa, j);
				
	}
	
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		
			mt[i][j] = M[i][j];
	}
	
	for(i = 0; i < N; i++)
	
		mt[i][N] = x[i];
		
	
	printf("\nMatriz total:\n");
	imprime(M);
		
}

int main(int argc, char **argv)
{
	FILE *in;
	double x[N] = {0}, xa[N] = {0.1, 0.1, -0.1}, Y[N], aux;	//xa é o vetor chute
	double norm, norma, norm1, norma1, M[N][N+1], mt[N][N+1];
	double A[T][T], c;
	int i, j, k, q = 0;
	sistfunc equacao[N]= {f1, f2, f3};

	
	in = fopen(argv[1], "w");
	
	
	do
	{
		
		norm = 0.0;
		
			
		for(i = 0; i < N; i++)
		{		
			x[i] = equacao[i](xa);
 								
			norm += pow(fabs(x[i]),2);
		}
		
				
		norma = sqrt(norm);
		
		
		jacobiano(equacao, M, x, xa, mt);
		
		
		for(i = 0; i < N; i++)
	   	{
			for(j = 0; j < (N+1); j++)
			
				A[i][j] = mt[i][j];        		
				
			printf("\n");   	
		}


		
		diagonalizando(A);
		
			  			
		for(i = 0; i < N; i++)
 
    	   	 	Y[i] = A[i][N]/A[i][i];

	
		for(i = 0; i < N; i++)
		
			xa[i] = xa[i]+Y[i];
		
			
		norm1 = 0.0;
		
			
		for(i = 0; i < N; i++)
		{	
			x[i] = equacao[i](xa);
 							
			norm1 += pow(fabs(x[i]),2);
		}
	
			
		norma1 = sqrt(norm1);
			
		
		fprintf(in, "%d\t%lf\n", itt, norma1);
		
		q++;
		
		   
	}while(fabs(norma1-norma) > tol);
	
	
	printf("\n\nIterações: %d\n", q);

}
