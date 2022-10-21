/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
*********************************************/

#include"opt_alg.h"
#include<cmath>
#include<ctime>
#include<cstdlib>


void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

long double f1(long double x);

int f_calls = 0;

int main()
{
	try
	{

	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab1()
{
	srand(time(NULL));
	double x0 = rand()%200 - 100; //staring point in range (-100, 100)
	double d = 0.1, alpha = 1.1;

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}


// long double* expansion(matrix(*ff)(matrix, matrix, matrix), long double x0, long double d, double alpha, int Nmax, int &f_calls, matrix ud2)

long double f1(long double x){
	f_calls++;
	
	long double a = (-1) * std::cos(0.1*x) *
		(
			1 / std::exp(
				pow(0.1*x-2*M_PI,2)
			)
		) +
		0.002* pow(0.1*x,2);

		return a;
}