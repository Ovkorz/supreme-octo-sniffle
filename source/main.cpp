/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
*********************************************/

#include<cmath>
#include<ctime>
#include<cstdlib>

#include"opt_alg.h"

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
	std::cout<< "lab 1" << std::endl;

	srand((time(NULL)));
	long double* temp_tab = new long double[2]{ 0.0, 0.0 };
	double x0; int rand_temp;									//ponieważ chcę losować liczby niecałkowite
	double d = 0.1;												//to jest krok chyba, nie pamiętam jaki powinien być, jak coś to się zmieni
	double alfa; int alfa_temp;									// też chcę niecałkowite więc podobnie - zamiana na double potem
	int N_max = 1000;												//na razie dam 5 jak coś to zmień
	

	for (int j = 0; j < 3; j++) {
		alfa_temp = rand() % 301 + 100;							
		alfa = alfa_temp / 100.0;
		std::cout << "\nalpha = " << alfa << std::endl;
		std::cout << "Lp.\tx0\ta\tb\tf_calls" << std::endl;

		for (int i = 0; i < 100; i++) {
			rand_temp = rand() % 20001 - 10000;					//o tutaj liczby int
			x0 = rand_temp / 100.0;										//dzielenie tak aby były liczby niecałkowite 
			temp_tab = expansion(&f1, x0, d, alfa, N_max, f_calls);
			std::cout 	<< i+1 << "\t"
						<< x0 << "\t"
						<< temp_tab[0] << "\t"
						<< temp_tab[1] << "\t"
						<< f_calls 
			<< std::endl;

			f_calls = 0;
		}
	}

}

void lab2()
{
	//Funkcja testowa
	double s = 0.5, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
	int Nmax = 1000;
	solution opt;
	matrix x0, s0;
	s0 = matrix(2, 1, s);
	x0 = 2 * rand_mat(2, 1) - 1;
	cout << x0 << endl << endl;
	opt = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();
	opt = Rosen(ff2T, x0, s0, alphaR, beta, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Ramie robota
	s = 2;
	x0 = 10 * rand_mat(2, 1);
	cout << x0 << endl << endl;
	opt = HJ(ff2R, x0, s, alphaHJ, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();
	s0 = matrix(2, 1, s);
	opt = Rosen(ff2R, x0, s0, alphaR, beta, epsilon, Nmax);
	cout << opt << endl << endl;
	solution::clear_calls();
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