/////////////////////////////////////////////
/////   Programa: transmision.cpp      /////
/////   Lenguaje: C++                   /////
/////   Creador: Adrián Arce Sánchez    /////
/////   Fecha de creación: 25/05/2020   /////
/////   Fecha ult.modific: 16/07/2020   /////
/////////////////////////////////////////////


/////////////////////////////////////////////////////////
/////   Descripción del programa:                   /////
/////   Resolucion de la ecuación de Scrhodinguer   /////
/////   para el potencial cuadrado.                 /////
/////   Calculo de valores esperados                /////
/////   y coef. de transmision.                     /////
/////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "gsl_rng.h"

using namespace std;

const int N = 1000;
const int Nd = 50;
const int n_ciclos = N/10;
const double lambda = 0.3;
const int bucle = 100000;
const int seed = 1003423;
const int tanda = 1000;
const bool detector_activado=true;
const bool parada_tras_deteccion=true;
const bool efecto_detec=true;

///////////////////////////////////////////////////////////////////
//DEfinición de los cálculos para numeros complejos

typedef struct FCOMPLEX {double r;double i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b)
// Output c=a+b
{
	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
// Output c=a-b
{
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}

fcomplex Cmul(fcomplex a, fcomplex b)
// Output c=a*b
{
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(double re, double im)
// Output c=(re, im)
{
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg (fcomplex z)
// Output c=conjudado(z)
{
	fcomplex c;
	c.r=z.r;
	c.i=-z.i;
	return c;
}

fcomplex Cdiv (fcomplex a, fcomplex b)
// Output c=a/b
{
	fcomplex c;
	double r,den;
	if (fabs(b.r)>=fabs(b.i)){
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabs (fcomplex z)
//output c=|z|
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x==0.0)
		ans=y;
	else if (y==0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z)
// output c=sqrt(z)
{
	fcomplex c;
	double x,y,w,r;
	if ((z.r==0.0)&&(z.i==0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x>=y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r>= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i>=0.0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul (double x, fcomplex a)
// output: x=a*x
{
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}


fcomplex Cpow(fcomplex x, int n)
// output x=x^n
{
	fcomplex c;
	int i;
	c=Complex(x.r,x.i);
	for (i=1;i<n;i++)
	{
		c=Cmul(c,x);
	}
	return c;
}


fcomplex Cgauss (double x,double y)
// Devuelve el nœmero complejo correspondiente al modulo y y la fase x 
{
	fcomplex c;
	c=Complex(cos(x),sin(x));
	c=RCmul(y,c);
	return c;
	}

////////////////////////////////////////////////////////////////////

fcomplex inicia_valor_funcion(double K0, fcomplex funcion[]);
void calcula_potencial_Vj (double Vj[], double K0);
void calcula_K0(double & K0);
void calcula_S(double & S,double K0);
void calcula_A(fcomplex & A_menos, fcomplex A_cero[], fcomplex & A_mas, double S, double Vj[]);
void calcula_alpha_y_gamma_inv(fcomplex alpha[],fcomplex gamma_inv[], fcomplex A_menos, fcomplex A_cero[], fcomplex A_mas);
void calcula_b (fcomplex b[], double S, fcomplex funcion[]);
void calcula_beta(fcomplex beta[], fcomplex b[], fcomplex gamma_inv[], fcomplex A_mas);
void calcula_chi (fcomplex chi[], fcomplex alpha[], fcomplex beta[]);
void calcula_funcion (fcomplex funcion[], fcomplex chi[]);
bool detecta_transmision(fcomplex funcion[], int &Mt);
bool detecta_reflexion(fcomplex funcion[], int &Mr);
double calcula_norma(fcomplex funcion[]);
double posicion_esperada(fcomplex funcion[]);
fcomplex momento_esperado(fcomplex funcion[]);
fcomplex energia_esperada(fcomplex funcion[],double V[]);

gsl_rng *tau;

////////////////////////////////////////////////////////////////////////////////

int main (void)
{
    fcomplex funcion[N+1];

    fcomplex chi[N+1];
    fcomplex alpha[N], beta[N];

    fcomplex A_menos, A_cero[N+1], A_mas;

    double Vj[N+1];
    fcomplex b[N+1], gamma_inv[N];

    double K0, S;

    ofstream salida_posicion, salida_momento, salida_energia, salida_coeficientes;
    double norma=0;

    bool deteccion;
    int paso=0;

    int Mt=0, Mr=0; //Contadores para los coeficientes de transmision
    double K=0, R=0;

    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,seed); //Inicializamos la semilla

    //Apertura de ficheros
    salida_posicion.open("resultados/posicion_esperada.dat");
    salida_momento.open("resultados/momento_esperado.dat");
    salida_energia.open("resultados/energia_esperada.dat");
    salida_coeficientes.open("resultados/coeficientes.dat");

    for(int intento=0; intento<tanda; intento++)
    {
        //Inicio de programa

        //Calculo de las constantes necesarias
        calcula_K0(K0);
        calcula_S(S,K0);

        //Calculo del potencial
        calcula_potencial_Vj(Vj,K0);

        calcula_A(A_menos, A_cero, A_mas, S, Vj);

        //Calculo del vector alpha y gamma_inv
        calcula_alpha_y_gamma_inv(alpha, gamma_inv, A_menos, A_cero ,A_mas);

        //Calculo de la condicion inicial de la onda
        inicia_valor_funcion(K0,funcion);


        paso=0;
        deteccion=false;

        //Inicio del bucle
        do
        {
            calcula_b(b,S,funcion);
            calcula_beta(beta, b, gamma_inv, A_mas);
            calcula_chi(chi, alpha, beta);
            calcula_funcion(funcion,chi);
            if(((paso%Nd==0) && (paso!=0)) && detector_activado==1)
            {
                deteccion=detecta_transmision(funcion,Mt);
                if (deteccion==false) // EN caso de no detectar la particula se observa el detector de reflexion
                {
                    deteccion=detecta_reflexion(funcion,Mr);
                }
                if(parada_tras_deteccion==false)
                    deteccion=false; //Para permitir que el programa continue tras la deteccion, necesario para valores esperados
            }
        
            //Escritura valores esperado
            //Solo activa en caso de desactivacion de parada tras deteccion de particula
            if((parada_tras_deteccion==false)&&(paso%10==0))
            {
                //Valor esperado de la posicion
                salida_posicion << paso << " " << posicion_esperada(funcion) << endl;

                //Valor esperado del momento
                salida_momento << paso << " " << momento_esperado(funcion).r << " " << momento_esperado(funcion).i << endl;

                //Valor esperado de la energía
                salida_energia << paso << " " << energia_esperada(funcion,Vj).r << " " << energia_esperada(funcion,Vj).i << endl;
            }

            //Conteo de tiempo
            paso++;

        }while((paso<bucle) && (deteccion==false));
    }

    //Ajuste de la suma de K y R en funcion de las tandas realizadas y escritura en fichero de estos
    salida_coeficientes << "Transmision: " << Mt/(1.*tanda) << " Reflexion: " << Mr/(1.*tanda) << " Valor N: " << N << " Valor lambda: " << lambda << " Tanda: " << tanda << endl;


    //Cierre de ficheros
    salida_posicion.close();
    salida_momento.close();
    salida_energia.close();
    salida_coeficientes.close();
}
///////////////////////////////////////////////////////////////////////////////

//Funcion que recibe como parametros K0, j para obtener la imajen de j,0
fcomplex inicia_valor_funcion(double K0, fcomplex funcion[])
{
    double base;
    double arg;
    double norma;

    funcion[0] = Complex(0,0);
    funcion[N] = Complex(0,0);

    for(int j=1;j<N;j++)
    {
        base = K0*j*1.;
        arg = exp(-8. * pow((4.*j-N), 2) / (N*N*1.));

        funcion[j] = Cgauss(base, arg);
    }

    norma=calcula_norma(funcion);

    for(int i=0; i<=N; i++)
    {
        funcion[i]=RCmul(pow(norma,-0.5),funcion[i]);
    }
}


void calcula_potencial_Vj (double Vj[], double K0)
{
    double min, max;

    min = (double)(2*N/5);
    max = (double)(3*N/5);

    for (int j=0; j<=N;j++)
    {
        if ((j>=min)&&(j<=max))
            Vj[j]=lambda*K0*K0;
        else
            Vj[j]=0.;
    }
}

void calcula_K0(double & K0)
{
    K0 = 2.*M_PI*n_ciclos/(N*1.);
}

void calcula_S(double & S,double K0)
{
    S = 0.25/(K0*K0*1.);
}

void calcula_A(fcomplex & A_menos, fcomplex A_cero[], fcomplex & A_mas, double S, double Vj[])
{
    A_menos = Complex(1,0);
    A_mas = Complex(1,0);

    fcomplex sumando_s=Complex(0,2./S);
    fcomplex aux = Cadd(Complex(-2,0),sumando_s);

    for (int j = 0; j<=N; j++)
        A_cero[j] = Csub(aux, Complex(Vj[j],0));
    
}

//Calcula el vector gamma y el vector alpha
void calcula_alpha_y_gamma_inv(fcomplex alpha[],fcomplex gamma_inv[], fcomplex A_menos, fcomplex A_cero[], fcomplex A_mas)
{
    //Condiciones de contorno
    alpha[N-1] = Complex (0,0);

    for (int j=N-2; j>=0; j--)
    {
        gamma_inv[j+1] = Cadd(A_cero[j+1], Cmul(A_mas,alpha[j+1]));
        alpha[j] = RCmul(-1,Cdiv(A_menos,gamma_inv[j+1]));
    }
    gamma_inv[0] = Cadd(A_cero[0], Cmul(A_mas,alpha[0]));
}

//Calcula el vector b para la iteracion n-esima
void calcula_b (fcomplex b[], double S, fcomplex funcion[])
{
    for(int j=0;j<=N;j++)
        b[j] = Cmul(Complex(0,4./S), funcion[j]);
}

//Calcula el vector beta para la iteracion n-esima
void calcula_beta(fcomplex beta[], fcomplex b[], fcomplex gamma_inv[], fcomplex A_mas)
{
    //Condiciones de contorno
    beta[N-1] = Complex (0,0);
    for (int j=N-2; j>=0; j--)
    {
        beta[j] = Cdiv(Csub(b[j+1],Cmul(A_mas,beta[j+1])),gamma_inv[j+1]);
    }
}

void calcula_chi (fcomplex chi[], fcomplex alpha[], fcomplex beta[])
{
    //Condiciones de contorno
    chi[0] = Complex (0,0);

    for (int j=1; j<=N; j++)
    {
        chi[j]=Cadd( Cmul( chi[j-1], alpha[j-1] ) , beta[j-1]);
    }
}

void calcula_funcion (fcomplex funcion[], fcomplex chi[])
{
    for (int j=1; j<N; j++)
    {
        funcion[j] = Csub(chi[j] , funcion[j]);
    }
}

bool detecta_transmision(fcomplex funcion[], int &Mt)
{
    double probabilidad=0;
    double nueva_norma=0;
    
    for (int k=(4*N)/5; k<=N; k++)
            probabilidad += pow(Cabs(funcion[k]),2);
    
    if(gsl_rng_uniform(tau)<probabilidad)
    {
        Mt++;
        return true;
    }else if(efecto_detec)
    {
        //Colapso de la funcion de onda en el detector
        for (int k=4*N/5; k<=N; k++)
            funcion[k]=Complex(0,0);

        //Calculo de la nueva norma de la onda tras el colapso
        for(int i=0; i<=N; i++)
            nueva_norma += pow(Cabs(funcion[i]),2);

        //Ajuste de los valores normalizados de la onda
        for(int i=0; i<=N; i++)
            funcion[i]=RCmul(pow(nueva_norma,-0.5),funcion[i]);
        return false;
    }
    

}

bool detecta_reflexion(fcomplex funcion[], int &Mr)
{
    double probabilidad=0;
    double nueva_norma=0;

    for (int k=0; k<=N/5; k++)
            probabilidad += pow(Cabs(funcion[k]),2);
    

    if(gsl_rng_uniform(tau)<probabilidad)
    {
        Mr++;
        return true;
    }else if(efecto_detec)
    {
        //Colapso de la funcion de onda en el detector
        for (int k=0; k<=N/5; k++)
            funcion[k]=Complex(0,0);

        //Calculo de la nueva norma de la onda tras el colapso
        for(int i=0; i<=N; i++)
            nueva_norma += pow(Cabs(funcion[i]),2);

        //Ajuste de los valores normalizados de la onda
        for(int i=0; i<=N; i++)
            funcion[i]=RCmul(pow(nueva_norma,-0.5),funcion[i]);
        return false;
    }
    

}

double calcula_norma(fcomplex funcion[])
{
    double norma=0;
    for(int j=0; j<=N; j++)
        norma += pow(Cabs(funcion[j]),2);
    return norma;
}

double posicion_esperada(fcomplex funcion[])
{
    double valor_esperado=0;
    for(int j=0; j<=N; j++)
        valor_esperado += pow(Cabs(funcion[j]),2)*j;
    return valor_esperado;
}

fcomplex momento_esperado(fcomplex funcion[])
{
    fcomplex valor_esperado=Complex(0,0);
    fcomplex aux=Complex(0,0);

    for(int j=1; j<N; j++)
    {
        aux=Csub(funcion[j+1],funcion[j-1]);
        aux=Cmul(Complex(0,-0.5),aux);
        aux=Cmul(Conjg(funcion[j]),aux);
        valor_esperado = Cadd(valor_esperado,aux);
    }

    return valor_esperado;

}

fcomplex energia_esperada(fcomplex funcion[],double V[])
{
    fcomplex valor_esperado=Complex(0,0);
    fcomplex aux=Complex(0,0);

    for(int j=1; j<N; j++)
    {
        aux=RCmul(2,funcion[j]);
        aux=Csub(aux,funcion[j+1]);
        aux=Csub(aux,funcion[j-1]);
        aux=Cadd(aux,RCmul(V[j],funcion[j]));
        aux=Cmul(Conjg(funcion[j]),aux);
        valor_esperado = Cadd(valor_esperado,aux);
    }

    return valor_esperado;

}