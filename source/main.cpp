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
#include<fstream>
#include<cstdlib>
#include<ctime>

#include"opt_alg.h"

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

long double f1(long double x);

int f_calls = 0;
const int Nmax = 1e7;

int main()
{
	try
	{
		lab4();
		
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

// void lab1()
// {
// 	std::cout<< "lab 1" << std::endl;

// 	srand((time(NULL)));
// 	long double* temp_tab = new long double[2]{ 0.0, 0.0 };
// 	double x0; int rand_temp;									//ponieważ chcę losować liczby niecałkowite
// 	double d = 0.1;												//to jest krok chyba, nie pamiętam jaki powinien być, jak coś to się zmieni
// 	double alfa; int alfa_temp;									// też chcę niecałkowite więc podobnie - zamiana na double potem
// 	int N_max = 1000;												//na razie dam 5 jak coś to zmień
	

// 	for (int j = 0; j < 3; j++) {
// 		alfa_temp = rand() % 301 + 100;							
// 		alfa = alfa_temp / 100.0;
// 		std::cout << "\nalpha = " << alfa << std::endl;
// 		std::cout << "Lp.\tx0\ta\tb\tf_calls" << std::endl;

// 		for (int i = 0; i < 100; i++) {
// 			rand_temp = rand() % 20001 - 10000;					//o tutaj liczby int
// 			x0 = rand_temp / 100.0;										//dzielenie tak aby były liczby niecałkowite 
// 			temp_tab = expansion(&f1, x0, d, alfa, N_max, f_calls);
// 			std::cout 	<< i+1 << "\t"
// 						<< x0 << "\t"
// 						<< temp_tab[0] << "\t"
// 						<< temp_tab[1] << "\t"
// 						<< f_calls 
// 			<< std::endl;

// 			f_calls = 0;
// 		}
// 	}

// }

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
	//Funkcja testowa
	ofstream S;
	S.open("wyniki.txt");
	double epsilon = 1e-3, h = -1;
	cout << h << endl << endl;
	int Nmax = 5000;
	matrix x0;
	solution opt;
	x0 = 20 * rand_mat(2, 1) - 10;
	//x0(0) = -6.95284;
	//x0(1) = 9.72172;
	//cout << x0(0)<< " " <<x0(1) << endl << endl;
	//opt = SD(ff4T, gf4T, x0, h, epsilon, Nmax);
	//cout << opt << endl << endl;
	//solution::clear_calls();
	//opt = CG(ff4T, gf4T, x0, h, epsilon, Nmax);
	//cout << opt << endl << endl;
	//solution::clear_calls();
	//opt = Newton(ff4T, gf4T, Hf4T, x0, h, epsilon, Nmax);
	//cout << opt << endl << endl;
	//solution::clear_calls();

	for (int i = 0; i < 100; i++) {
		x0 = 20 * rand_mat(2, 1) - 10;
		cout << x0 << endl << endl;
		opt = SD(ff4T, gf4T, x0, h, epsilon, Nmax);
		S << x0(0) << " " << x0(1) << " " << opt.x(0) << " " << opt.x(1) << " " << opt.y << " " << opt.f_calls << " " << opt.g_calls << " ";
		solution::clear_calls();
		opt = CG(ff4T, gf4T, x0, h, epsilon, Nmax);
		S << opt.x(0) << " " << opt.x(1) << " " << opt.y << " " << opt.f_calls << " " << opt.g_calls << " ";
		solution::clear_calls();
		opt = Newton(ff4T, gf4T, Hf4T, x0, h, epsilon, Nmax);
		S<< opt.x(0) << " " << opt.x(1) << " " << opt.y << " " << opt.f_calls << " " << opt.g_calls <<" "<< opt.H_calls << endl;
		solution::clear_calls();
	}

	////Regresja liniowa
	epsilon = 1e-5, h = 0.001;
	Nmax = 20000;
	x0 = matrix(3, 1);
	opt = CG(ff4R, gf4R, x0, h, epsilon, Nmax);
	cout << opt << endl << endl;

	int n = 3, m = 100, P = 0;
	matrix X(n, m), Y(1, m);
	ifstream Sin("XData.txt");
	Sin >> X;
	Sin.close();
	Sin.open("YData.txt");
	Sin >> Y;
	Sin.close();
	for (int i = 0; i < m; ++i)
	{
		h = (trans(opt.x) * X[i])();
		h = 1.0 / (1.0 + exp(-h));
		if (round(h) == Y(0, i))
			++P;
	}
	cout << P << endl << endl;
}

void lab5()
{
	solution temp; temp.clear_calls();
	//starting point
	matrix x0(2,1); x0(0) = 3.4; x0(1) = 2.2;
	matrix ld(2,1); ld(0) = 0.1; ld(1) = 0.2;

	//parameters
	const double epsilon = 1e-8, c = 1e7, dc = 100;

	const string file_name_prefix = "lab5_";
	const int a_values[] = {1, 10, 100};

	ofstream output(file_name_prefix + "theoretical.txt");
	for(int i = 0; i <= 100; i++){

		x0(0) = rand() % 10 - 5; x0(1) = rand() % 10 - 5;
		double w = i * 0.01;
		output 	<< w << "," << x0(0) << "," << x0(1) << ",";

		for(int j = 0; j < 3; j++){	

			// matrix a_w:
			// [[a],
			// 	[w]]
			matrix a_w(2,1); a_w(0) = a_values[j]; a_w(1) = w;
	
			solution Xopt = Powell(ff5T, x0, matrix(), epsilon, Nmax, a_w);
			matrix f_values = ff5_values(Xopt.x, a_w(0));

			output << Xopt.x(0) << "," << Xopt.x(1) << ","
					<< f_values(0) << "," << f_values(1)  << "," << Xopt.f_calls << ",";

			Xopt.clear_calls();

		}
		output << endl;


	} 
	output.close();

	output.open(file_name_prefix + string("real.txt"));
	matrix pen_coef(2,2); pen_coef(0,0) = c; pen_coef(1,0) = dc; 

	// pen_coef
	// [[c, w],
	//  [dc, 0]]

	

	for(int i = 0; i <= 100; i++){
		ld(0) = (rand() % 3) / 4.; ld(1) =( rand() % 2)/4.;

		double w = i * 0.01;
		pen_coef(0,1) = w;

		solution Xopt = pen(ff5rwPT, Powell, ld, pen_coef, epsilon, Nmax);
		matrix f_values = ff5rwP_values(Xopt.x);

		output 	<< w << "\t" << ld(0) << "\t" << ld(1) << "\t" << Xopt.x(0) << "\t" << Xopt.x(1) << "\t"
				<< f_values(0) << "\t" << f_values(1) << "\t" << Xopt.f_calls << endl;

		#ifdef VERBOSE
				cout<< "[lab5] RW problem optimization result: " << print_m_l(ld, "ld0") 
					<< ", w: " << w << "," << print_m_l(Xopt, "Xopt") << endl;
		#endif


		Xopt.clear_calls();

	}

	output.close();

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