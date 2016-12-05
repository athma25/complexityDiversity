#include"pdesolver.h"

void main()
{
	struct parameters par;
	par.a=1;
	par.k=1;
	par.K=1;
	par.na=2;
	par.nk=2;
	par.r=1;
	par.d=.01;				//bin size
	
	float trunc=.001;			//truncation carrying capacity value
	float t,dt=.1;
	
	int i,j,N,check;			//N: number of bins
	int maxIter=100;
	float x=0,xMax;				//dx: bin size; xMax: phenotype space truncation
	
	float **pop,*binPt;
	FILE *f,*g,*h;				//output; f has time and then pop distribution of one iteration in a line

	f=fopen("popEvol","w");
	g=fopen("evolTime","w");
	h=fopen("binPts","w");

	xMax=powf(-log(trunc/par.K)/par.k,1/par.nk)+trunc;
	N=2*ceil(xMax/par.d);

	pop=(float **)malloc(sizeof(float *)*2);
	pop[0]=(float *)malloc(sizeof(float)*N);
	pop[1]=(float *)malloc(sizeof(float)*N);
	binPt=(float *)malloc(sizeof(float)*N);

	binPt[0]=-xMax+par.d/2;
	fprintf(h,"%f\n",binPt[0]);
	for(i=1;i<N;i++)
	{
		binPt[i]=binPt[i-1]+par.d;
		fprintf(h,"%f\n",binPt[i]);
	}
	
//initial population distribution
	for(i=0;i<N;i++)
		pop[0][i]=1;

//checks
//	printf("xMax=%f\t N=%d\n",xMax,N);

//simulation
	i=1;
	t=0;
	while(i<maxIter)
	{
		check=0;
		j=0;
//		printf("CHECK POINT 1 \t i=%d\n",i);
		while(check==0 && j<N)
		{
			pop[1][j]=pop[0][j]+par.r*dt*pop[0][j]*(1-trapComp(N,j,pop[0],binPt,par)/carry(binPt[j],par));
			if(pop[1][j]<0)
				check=1;
			j=j+1;
//			printf("CHECK POINT 2\t j=%d\n",j);
		}
		printf("i=%d \t check=%d \t pop[1]=%f \t pop[N]=%f\n",i,check,pop[1][0],pop[1][N]);

		if(check==0)
		{
			fprintf(f,"%f\t",t);
			for(j=0;j<N;j++)
			{
				fprintf(f,"%f\t",pop[1][j]);
				pop[0][j]=pop[1][j];
			}
			fprintf(f,"\n");
			i=i+1;
			t=t+dt;
			fprintf(g,"%f\n",t);
			printf("%d\t %f\n",i,t);
		}
		if(check==1)
			dt=dt/2;
	}

	fclose(f);
	fclose(g);
	fclose(h);
}
