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
    //Variables modificables
    int N=1000;
    int nciclos=50;
    int duracion=5000;
    double lambda=1.;
    double s=0.01; //DiscretizaciÃ³n del tiempo
    double h=0.01; //DiscretizaciÃ³n del espacio
    double transmision;

    ofstream moment("momentesperado.txt");
    ofstream energ("energiesperado.txt");
    ofstream pos("posesperada.txt");
    ofstream trans("transmision.txt");
    ofstream datos("datos.txt");
    ofstream modulotxt("modulo.txt");
    
    //Otras variables
    complex<double> phi[N+1], xi[N+1], aux3;
    complex<double> A0[N], b[N], beta[N], alpha[N], gamma[N];
    double modulo;
    double Aplus, Aminus;
    double k0virg, V[N+1], svirg;
    double aux1, aux2;
    double posesperada, pos2esperada, momentesperado, moment2esperado, incertpos, incertmoment;
    double potenesperado, poten2esperado, incertpoten, energiaesperada, incertenergia;
    complex<double> d1[N+1], d2[N+1];
    double planck=1.;
    int j, n;    

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

    //CÃ¡lculo de la funciÃ³n de onda inical
    phi[0]=complex<double>(0.,0.);
    phi[N]=complex<double>(0.,0.);
    for (j=1; j<N; j++)
    {
        aux1=exp(-8.*(4*j-N)*(4*j-N)/(N*N));
        phi[j]=complex<double>(cos(k0virg*j),sin(k0virg*j))*aux1;
    }
    //Normalizacion de la funcion de onda
    modulo=0.;
    for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
    modulo=sqrt(modulo);
    for (j=1; j<N; j++) phi[j]=phi[j]/modulo;

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
    for (j=0; j<N+1; j++) datos << j*h << ", " << norm(phi[j]) << "\n"; datos << endl;

    //Iniciamos el bucle de pasos
    for (n=0; n<duracion; n++)
    {
        //Escribimos el modulo
        modulo=0.;
        for (j=1; j<N; j++) modulo= modulo + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
        modulo=sqrt(modulo);
        modulotxt << n << "\t" << modulo << endl;

        for (j=1; j<N; j++) b[j]=complex<double>(0, 4)*phi[j]/svirg;
        for (j=N-2; j>=0; j--) beta[j]=gamma[j+1]*(b[j+1]-beta[j+1]);    

        for (j=1; j<N; j++)
        {
            xi[j]=alpha[j-1]*xi[j-1]+beta[j-1];
            phi[j]=xi[j]-phi[j];
        } 

        for (j=0; j<N+1; j++) datos << j*h << ", " << norm(phi[j]) << "\n";
        datos << endl;

        transmision=0.;
        for (j=4*aux2; j<=N; j++) transmision = transmision + real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]);
        trans << n << "\t" << transmision << endl;

        //Derivadas
        d1[0]=phi[1]/h;
        d1[N]=(phi[N]-phi[N-1])/h;
        d2[0]=0.;
        d2[N]=0.;
        for (j=1; j<N; j++)
        {
            d1[j]=(phi[j+1]-phi[j-1])/(2*h);
            d2[j]=(phi[j+1]-2.*phi[j]+phi[j-1])/(h*h);
        }

        //Valores esperados
        posesperada = 0.;
        pos2esperada =0.;
        momentesperado = 0.;
        moment2esperado = 0.;
        potenesperado = 0.;
        poten2esperado = 0.;
        for (j=0; j<=N; j++)
        {
            posesperada = posesperada + j*h*(real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]));
            pos2esperada = pos2esperada + j*h*j*h*(real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]));
            momentesperado = momentesperado + h*real(conj(phi[j])*(-complex<double>(0,1.)*d1[j]));
            moment2esperado = moment2esperado - h*h*real(conj(phi[j])*d2[j]);
            potenesperado = potenesperado + V[j]*(real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]));
            poten2esperado = poten2esperado + V[j]*V[j]*(real(phi[j])*real(phi[j]) + imag(phi[j])*imag(phi[j]));
        }
        
        incertpoten = sqrt(poten2esperado - potenesperado*potenesperado);
        incertmoment = sqrt(moment2esperado - momentesperado*momentesperado);
        incertpos = sqrt(pos2esperada - posesperada*posesperada);
        energiaesperada = moment2esperado + potenesperado ;
        incertenergia = sqrt(4*incertpoten*incertpoten+incertmoment*incertmoment);

        pos << n*s+0.01 << "\t" << posesperada << "\t" << incertpos  << endl;
        moment << n*s+0.01 << "\t" << momentesperado << "\t" << incertmoment << endl;
        energ << n*s+0.01 << "\t" << energiaesperada << "\t" << incertenergia << "\t" << moment2esperado << "\t" << potenesperado   << endl;

    } 

    energ.close();
    trans.close();
    moment.close();
    pos.close();
    modulotxt.close();
    datos.close();
    return 0;
}