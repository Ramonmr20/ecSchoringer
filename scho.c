#include <stdio.h>
#include <math.h>
#include "complex.h"

//Parámetros
#define N 1000
#define ND 100
#define NCICLOS 250
#define LAMBDA 0.3
#define PASOS 10000

//Constantes
#define PI 3.141592
#define K0 2*PI*NCICLOS/(1.*N)
#define STIL 1/(4.*K0*K0)

fcomplex funcionInicial(int j);
void inicializarFuncioninicial(fcomplex phi[], fcomplex alpha[], fcomplex gamma[], double v[], int n);
void calculoBeta(fcomplex phi[], fcomplex beta[],fcomplex gamma[], int n);
void calculoX(fcomplex alpha[], fcomplex beta[], fcomplex x[], int n);
void paso(fcomplex phi[], fcomplex alpha[], fcomplex beta[], fcomplex x[], fcomplex gamma[], int n, int *t);

int main(){
    
    fcomplex phi[N+1];//Función de onda
    fcomplex alpha[N], beta[N], x[N], gamma[N];
    double v[N+1];
    double norma, pd;
    int n, t;
    FILE *f1, *f2, *f3;
    
    f1 = fopen("ecScho.txt","w");
    f2 = fopen("potencial.txt","w");
    f3 = fopen("norma.txt", "w");
    
    //Inicializo variables
    n = N;
    t = 0;

    inicializarFuncioninicial(phi, alpha, gamma, v, n); //Inicializa valores y calcula el valor de alpha y gamma
    
    //Output potencial
    for(int i=0;i<=n;i++)
        fprintf(f2,"%lf\t",v[i]);
    
    //Programa  y output
    for(;t<PASOS;){
        
        inicializarFuncioninicial(phi, alpha, gamma, v, n);
        
        for(int a=0;a<ND;a++) paso(phi, alpha, beta, x, gamma, n, &t);
        
        norma = 0;
        //Output
        for(int i=0;i<=n;i++){
            fprintf(f1,"%lg\t",Cabs(phi[i])*Cabs(phi[i]));
            norma += Cabs(phi[i])*Cabs(phi[i]);
        }
        fprintf(f1,"\n");
        
        fprintf(f3,"%lg\t",norma);
    }
    
    fclose(f1);
    fclose(f2);
    fclose(f3);
}

void inicializarFuncioninicial(fcomplex phi[], fcomplex alpha[], fcomplex gamma[], double v[], int n){
    
    fcomplex aux0, aux1, mUno;
    
    //Cálculo phi0 e inicializacion potencial V
    for(int i=0;i<=n;i++){
        phi[i] = funcionInicial(i);
        if((i<=3*N/5)&&(i>=2*N/5)){
            v[i] = LAMBDA*K0*K0;
        }else v[i] = 0;
    }
    
    //Cálculo alphas y gammas
    aux0 = Complex(-2,2/(1.*STIL));
    mUno = Complex(-1,0);
    alpha[n-1] = Complex(0,0);
    gamma[n-1] = aux0;
    
    for(int i=n-2;i>=0;i--){ 
        //Calculo alpha
        alpha[i] = Cdiv(mUno,gamma[i+1]);
        
        //Calculo gamma
        gamma[i] = Cadd(aux0,alpha[i]);
        aux1 = Complex(v[i],0);
        gamma[i] = Csub(gamma[i],aux1);
    }
}

fcomplex funcionInicial(int j){
    fcomplex phi, expoC;
    double expoR;
    
    expoC = Cgauss(K0*j,1);// e^ik0x
    expoR = exp(-8*(4*j-N)*(4*j-N)/(1.*N*N));
    
    phi = RCmul(expoR,expoC);
    
    return phi;
}

void calculoBeta(fcomplex phi[], fcomplex beta[],fcomplex gamma[], int n){
    
    fcomplex b[n];
    fcomplex aux;
    
    beta[n-1] = Complex(0,0);
    
    for(int i=n-2;i>=0;i--){
        //Calculo b
        aux = Complex(0,4);
        aux = RCmul(1/(1.*STIL),aux);
        b[i+1] = Cmul(aux,phi[i+1]);
        
        //Calculo beta
        beta[i] = Csub(b[i+1],beta[i+1]);
        beta[i] = Cdiv(beta[i],gamma[i+1]);
        
    }
    
}

void calculoX(fcomplex alpha[], fcomplex beta[], fcomplex x[], int n){
    
    x[0] = Complex(0,0);
    
    for(int i=1;i<=n;i++){
        x[i] = Cmul(alpha[i-1],x[i-1]);
        x[i] = Cadd(x[i],beta[i-1]);
    }
}

void paso(fcomplex phi[], fcomplex alpha[], fcomplex beta[], fcomplex x[], fcomplex gamma[], int n, int *t){
    
    calculoBeta(phi, beta, gamma, n);
    calculoX(alpha,beta,x,n);
    
    //Calculo phis
    phi[0] = Complex(0,0);
    phi[n] = Complex(0,0);
    
    for(int i=1;i<n;i++){
        phi[i] = Csub(x[i],phi[i]);
    }
    
    *t += 1;
}

double calculoPd(fcomplex phi[]){
    double pd = 0;
    
    for(int i=4*N/5.; i<N; i++){
        pd += Cabs(phi[i])*Cabs(phi[i]);
    }
    
    return pd;
}
