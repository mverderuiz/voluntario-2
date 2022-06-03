//Comprobar que la norma se conserva
//Modificar el valor de lambda

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <complex>
#define PI 3.14159265

using namespace std;

int main (void)
{
    srand (time(NULL));

    //Variables modificables
    int N=2000;
    int nciclos=50;
    int duracion=1000;
    double lambda=0.1;
    int m=1000;
    int nd=500;
    double s=0.01; //Discretización del tiempo
    double h=0.01; //Discretización del espacio

    ofstream datos("datos.txt");

    //Otras variables
    complex<double> phi[N+1], phi0[N+1], xi[N+1], aux3;
    complex<double> A0[N], b[N], beta[N], alpha[N], gamma[N];
    double modulo;
    double Aplus, Aminus;
    double k0virg, V[N+1], svirg;
    double aux1, aux2, aleatorio;
    int j, k, l, n, omega;  
    int  mt=0;  
    double transmision, reflexion, Pd;
    bool encontrado = false;

    k0virg=2.*PI*nciclos/N;
    svirg=1/(4*k0virg*k0virg);

    //Calculo del potencial
    for (j=0; j<N+1; j++)
    {
        aux1=(double) j;
        aux2=(double)N/5.;
        if ((aux1>=2*aux2)&&(aux1<=3*aux2)) V[j]=lambda*k0virg*k0virg;
        else V[j]=0.;
    }

    //Cálculo de la función de onda inical
    phi0[0]=complex<double>(0.,0.);
    phi0[N]=complex<double>(0.,0.);
    for (j=1; j<N; j++)
    {
        aux1=exp(-8.*(4*j-N)*(4*j-N)/(N*N));
        phi0[j]=complex<double>(cos(k0virg*j),sin(k0virg*j))*aux1;
    }
    //Normalizacion de la funcion de onda
    modulo=0.;
    for (j=1; j<N; j++) modulo= modulo + real(phi0[j])*real(phi0[j]) + imag(phi0[j])*imag(phi0[j]);
    modulo=sqrt(modulo);
    for (j=1; j<N; j++) phi0[j]=phi0[j]/modulo;

    //Calculamos los coeficientes gamma y alpha
    Aplus=1.;
    Aminus=1.;
    for (j=1; j<N; j++) A0[j]=complex<double>(-2.-V[j],2*Aminus/svirg);
    beta[N-1]=0.;
    alpha[N-1]=0.;
    gamma[N-1]=Aminus/A0[N-1];
    for (j=N-2; j>=0; j--) 
    {
        alpha[j]=-Aminus*gamma[j+1];
        gamma[j]=Aplus/(A0[j]+alpha[j]);
    }

    xi[0]=0;

    //Iniciamos el bucle de pasos
    for (l=0; l<m; l++)
    {
        for (j=0; j<=N; j++) phi[j]=phi0[j];
        encontrado=false;
        omega=0;

        do{

            for (n=0; n<=nd; n++)
            {
                for (j=1; j<N; j++) b[j]=complex<double>(0, 4)*phi[j]/svirg;
                for (j=N-2; j>=0; j--) beta[j]=gamma[j+1]*(b[j+1]-beta[j+1]);    

                for (j=1; j<N; j++)
                {
                    xi[j]=alpha[j-1]*xi[j-1]+beta[j-1];
                    phi[j]=xi[j]-phi[j];
                } 
            }

            transmision=0.;
            for (j=4*aux2; j<=N; j++) transmision = transmision + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
            aleatorio = (rand() % 1001)/1000.;
            if (aleatorio<transmision) 
            {
                mt++;
                encontrado=true;
            }
            else
            {
                for (j=4*aux2; j<N; j++) phi[j]=0.;
                modulo=0.;
                for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
                modulo=sqrt(modulo);
                for (j=1; j<N; j++) phi[j]=phi[j]/modulo;

                reflexion=0.;
                for (j=1; j<=aux2; j++) reflexion = reflexion + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
                aleatorio = (rand() % 1001)/1000.;
                if (aleatorio<reflexion) encontrado=true;
                else 
                {
                    for (j=0; j<=aux2; j++) phi[j]=0.;
                    modulo=0.;
                    for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
                    modulo=sqrt(modulo);
                    for (j=1; j<N; j++) phi[j]=phi[j]/modulo;
                }
            }
            omega=omega+1;
            
        }while (encontrado==false && omega<10);   
        datos << mt << endl;    
    }
    
    Pd=1.*mt/m;

    cout << Pd << endl;

    datos.close();
    return 0;
}