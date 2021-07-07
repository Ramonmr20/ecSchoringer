#include <stdio.h>
#include <math.h>
#include <time.h>
#include "complex.h"
#include "gsl_rng.h"

//Parámetros
#define N 2000
#define ND 100
#define NCICLOS N/10
#define LAMBDA 0.1
#define PASOS 1000

//Constantes
#define PI 3.141592
#define K0 2*PI*NCICLOS/(1.*N)
#define STIL 1/(4.*K0*K0)

gsl_rng *tau;

fcomplex funcionInicial(int j);
void inicializarFuncioninicial(fcomplex phi[], fcomplex alpha[], fcomplex gamma[], double v[], int n, double lam);
void calculoBeta(fcomplex phi[], fcomplex beta[],fcomplex gamma[], int n);
void calculoX(fcomplex alpha[], fcomplex beta[], fcomplex x[], int n);
void paso(fcomplex phi[], fcomplex alpha[], fcomplex beta[], fcomplex x[], fcomplex gamma[], int n, int *t);
double calculoPd(fcomplex phi[], int n);
double calculoPi(fcomplex phi[], int n);

int main(){
    
    //Inicializar num aleatorios
    extern gsl_rng *tau;
    int semilla = time(NULL);
    
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);
    
    fcomplex phi[N+1];//Función de onda
    fcomplex alpha[N], beta[N], x[N], gamma[N];
    double v[N+1];
    double norma, pd, pi, al, lam;
    int n, t, mt;
    FILE *f1;
    
    f1 = fopen("prueba.txt","w");
    
    //Inicializo variables
    n = N;
    t = 0;
    lam = LAMBDA;

    inicializarFuncioninicial(phi, alpha, gamma, v, n, lam); //Inicializa valores y calcula el valor de alpha y gamma
    
    
    //Programa  y output
    
    for(int p=0;p<PASOS;p++){
        printf("%i\n",p);
        
        for(int a=0;a<ND;a++) paso(phi, alpha, beta, x, gamma, n, &t);
        
        //Calculo PD
        pd = calculoPd(phi, n);
        al = gsl_rng_uniform(tau);
        
        if(al<pd){
            mt++;
            
            inicializarFuncioninicial(phi, alpha, gamma, v, n, lam);
            continue;
        }
        
        for(int i=4*n/5.;i<=n; i++){
            phi[i] = Complex(0,0);
        }
        
        //Normalización phi
        norma = 0;
        //Output
        for(int i=0;i<=n;i++){
            norma += Cabs(phi[i])*Cabs(phi[i]);
        }
        
        norma = sqrt(norma);
        
        for(int i=0; i<=n; i++){
            phi[i] = RCmul(1/(1.*norma),phi[i]);
        }
        
        //Calculo PI
        pi = calculoPi(phi, n);
        al = gsl_rng_uniform(tau);
        
        if(al<pi){
            inicializarFuncioninicial(phi, alpha, gamma, v, n, lam);
            continue;
        }
        
        for(int i=0.;i<=n/5.; i++){
            phi[i] = Complex(0,0);
        }
        
        //Normalización phi
        norma = 0;
        //Output
        for(int i=0;i<=n;i++){
            norma += Cabs(phi[i])*Cabs(phi[i]);
        }
        
        norma = sqrt(norma);
        
        for(int i=0; i<=n; i++){
            phi[i] = RCmul(1/(1.*norma),phi[i]);
        }
        
    }
    
    printf("COEF:: %lf\n",mt/(1.*PASOS));
    fclose(f1);
}

void inicializarFuncioninicial(fcomplex phi[], fcomplex alpha[], fcomplex gamma[], double v[], int n, double lam){
    
    fcomplex aux0, aux1, mUno;
    
    //Cálculo phi0 e inicializacion potencial V
    for(int i=0;i<=n;i++){
        phi[i] = funcionInicial(i);
        if((i<=3*n/5)&&(i>=2*n/5)){
            v[i] = lam*K0*K0;
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

double calculoPd(fcomplex phi[], int n){
    double pd = 0;
    
    for(int i=4*n/5.; i<n; i++){
        pd += Cabs(phi[i])*Cabs(phi[i]);
    }
    
    return pd;
}

double calculoPi(fcomplex phi[], int n){
    double pi = 0;
    
    for(int i=0.; i<n/5.; i++){
        pi += Cabs(phi[i])*Cabs(phi[i]);
    }
    
    return pi;
}
