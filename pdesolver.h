#include<stdio.h>
#include<stdlib.h>
#include<math.h>

struct parameters
{
	float nk;
	float na;
	float a;
	float k;
	float K;
	float d;
	float r;
};

float alpha(float x, float y, struct parameters par)
{
	return exp(-par.a*pow(fabs(x-y),par.na));
}

float carry(float x, struct parameters par)
{
	return par.K*exp(-par.k*pow(fabs(x),par.nk));
}

float trapComp(float n, int i, float *pop, float *binPt, struct parameters par)	//trapezoid method with n integrals at p(i); p denotes population size vector
{
	int j;
	float x=0;
	for(j=0;j<n-1;j++)
		x=x+(alpha(binPt[i],binPt[j],par)*pop[j]+alpha(binPt[i],binPt[j+1],par)*pop[j+1])*par.d/2;
	return x;
}
