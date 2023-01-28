#include"user_funs.h"
#include<cmath>
#include<iostream>
// #include<cstring>

double geometric_mean (const matrix &A){
    double sum_of_squares = 0;
    int *n = get_size(A);

    for(int i = 0; i<n[0]; i++){
		for(int j = 0; j < n[1]; j++){
        	sum_of_squares += pow(A(i,j), 2);
		}
    }

    return sqrt(sum_of_squares);
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5,1,10 });
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	y = abs(max - 40);
	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
	return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(2, 1), Y_ref(2, new double[2]{ 3.14,0 });
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
	int n = get_len(Y[0]);
	y = 0;
	for (int i = 0; i < n; ++i)
		y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) +
		pow(Y_ref(1) - Y[1](i, 1), 2) +
		pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
	y = y * 0.1;
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	double a = 0.98, b = 0.63, g = 9.81, PA = 1, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 10;
	matrix dY(3, 1);
	double FAout = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA) : 0;
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));
	return dY;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	double mr = 1, mc = 5, l = 0.5, b = 0.5;
	double I = mr * l * l / 3 + mc * l * l;
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
	return dY;
}


matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	if (isnan(ud2(0, 0)))
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
		y = ff4T(ud2[0] + x * ud2[1], ud1);
	return y;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix g(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	return g;
}

matrix Hf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix H(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
	return H;
}

matrix ff4R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	static bool read = true;
	if (read)
	{
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
		read = false;
	}
	double h;
	y = 0;
	for (int i = 0; i < m; ++i)
	{
		h = (trans(x) * X[i])();
		h = 1.0 / (1.0 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}
	y = y / m;
	return y;
}

matrix gf4R(matrix x, matrix ud1, matrix ud2)
{
	int m = 100;
	int n = get_len(x);
	matrix g(n, 1);
	static matrix X(n, m), Y(1, m);
	static bool read = true;
	if (read)
	{
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
		read = false;
	}
	double h;
	for (int j = 0; j < n; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			h = (trans(x) * X[i])();
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
	return g;
}



matrix ff5a(matrix x, matrix a, matrix ud2 = NAN){
	try{
		return matrix(
			a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2))
		);
	}
	catch (string ex_info)
	{
		throw ("ff5a(...):\n" + ex_info);
	}
}

matrix ff5b(matrix x, matrix a, matrix ud2 = NAN){
	try{	
		return matrix(
			1/a * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2))
		);
	}
	catch (string ex_info)
	{
		throw ("ff5b(...):\n" + ex_info);
	}
}

matrix ff5Ta(matrix h, matrix x, matrix coef){

	try{
		double a = coef(0,0);
		matrix d = get_col(coef,1);

		return ff5a(x + h * d, matrix(a));

	}
	catch (string ex_info)
	{
		throw ("ff5Ta(...):\n" + ex_info);
	}

}

matrix ff5Tb(matrix h, matrix x, matrix coef){

	try{
		double a = coef(0,0);
		matrix d = get_col(coef,1);
		
		return ff5b(x + h * d, matrix(a));
	}
	catch (string ex_info)
	{
		throw ("ff5Tb(...):\n" + ex_info);
	}
}

matrix ff5_values(matrix x, matrix a, matrix w){
	try{
		int *x_size = get_size(x);
		if(x_size[0] != 2 || x_size [1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
		
		matrix dY(2,1);
		dY(0) = ff5a(x, matrix(a))(0);
		dY(1) = ff5b(x, matrix(a))(0);

		return dY;
	}
	catch (string ex_info)
	{
		throw ("ff5_values(...):\n" + ex_info);
	}
}

matrix ff5(matrix x, matrix a, matrix w){
	try{
	int *x_size = get_size(x);
		if(x_size[0] != 2 || x_size [1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
		
		return w * ff5a(x, matrix(a)) + (1-w) * ff5b(x, matrix(a));
	}
	catch (string ex_info)
	{
		throw ("ff5(...):\n" + ex_info);
	}
}

matrix ff5T(matrix h, matrix x, matrix coef){
	// coef matrix:
	// [[a, d1],
	//  [w, d2]]
	try{
		int *x_size = get_size(x),
			*h_size = get_size(h),
			*a_size = get_size(coef);

		if(	x_size[0] != 2 || x_size[1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
		
		if(	h_size[0] != 1 || h_size[1] != 1)
			throw("wrong \"h\" matrix size: " + to_string(h_size[0]) + ", " + to_string(h_size[1])+"; expected: 1,1");

		if(	a_size[0] != 2 || a_size[1] != 2)
			throw("wrong \"coef\" matrix size: " + to_string(a_size[0]) + ", " + to_string(a_size[1])+"; expected: 2,2");

		double a = coef(0,0), w = coef(1,0);
		w = (w < 0 || w > 1) ? 0.5 : w;
		a = (a == 0) ? 1 : a;

		matrix d = get_col(coef,1);
		
		return ff5(x + h * d, matrix(a), matrix(w));
		// return w * ff5a(x + h * d, matrix(a)) + (1-w) * ff5b(x + h * d, matrix(a));
	}
	catch (string ex_info)
	{
		throw ("ff5T(...):\n" + ex_info);
	}
}

matrix ff5rwP_values(matrix ld){

	try{
		int *ld_size = get_size(ld);

		if(ld_size[0] != 2 || ld_size[1] != 1)
			throw("wrong \"ld\" matrix size: " +  to_string(ld_size[0]) + ", " + to_string( ld_size[1])+"; expected: 2,1");

		const double P_force = 1000, E_young_mod = 2.07e11, material_density = 7800;
		
		double length = ld(0), diam = ld(1);

		double 
		deflection = (
			(64 * P_force * pow(length, 3))/
			(3 * E_young_mod * M_PI * pow(diam, 4))
		),
		stress = (
			(32 * P_force * length)/
			(M_PI * pow(diam, 3))
		),
		mass = (
			M_PI * pow(diam/2, 2) * length * material_density
		);

		matrix dY(3,1);
		dY(0) = mass; 
		dY(1) = deflection; 
		dY(2) = stress;
		return dY;
	}
	catch (string ex_info)
	{
		throw ("ff5rwP_values(...):\n" + ex_info);
	}
}

matrix ff5rwP(matrix ld, matrix w, matrix c){
	//ld matrix
	// [[l],
	// 	[d]]
	try{
		
		// const double P_force = 1000, E_young_mod = 2.07e11, material_density = 7800;
		const double	LEN_MIN = 0.2, LEN_MAX = 1,
						DIAM_MIN = 0.01, DIAM_MAX = 0.05,
						DEFL_MAX = 0.005, STRESS_MAX = 3e8;

		double length = ld(0), diam = ld(1);

		matrix f_values = ff5rwP_values(ld);

		double 
		mass = f_values(0),
		deflection = f_values(1),
		stress = f_values(2);

		double penalty = pow(max(0., LEN_MIN - length),2) + pow( max(0., length - LEN_MAX), 2);

		penalty += pow( max(0., DIAM_MIN - diam), 2) + pow( max(0., diam - DIAM_MAX), 2);
		penalty += pow( max(0., deflection - DEFL_MAX), 2) + pow( max(0., stress - STRESS_MAX), 2);

		#ifdef VERBOSE
			cout<< "[ff5rwP] input " << print_m_l(ld, "ld") << ", " << print_m_l(w, "w") << ", " << print_m_l(c, "c")
				<< ", deflection: "<< deflection << ", stress: " << stress
				<< ", mass: "<< mass << ", penalty: "<< penalty<<endl;
		#endif

		return w * mass + (1-w) * deflection + c * penalty;
	}
	catch (string ex_info)
	{
		throw ("ff5rwP(...):\n" + ex_info);
	}
}

matrix ff5rwPT(matrix h, matrix ld, matrix coef){
	// coef
	// [[c, d1],
	//  [w, d2]]
	try{
		int *h_size = get_size(h),
			*ld_size = get_size(ld),
			*coef_size = get_size(coef);

		if(ld_size[0] != 2 || ld_size[1] != 1)
			throw("wrong \"ld\" matrix size: " +  to_string(ld_size[0]) + ", " + to_string( ld_size[1])+"; expected: 2,1");

		if(h_size[0] != 1 || h_size[1] != 1)
			throw("wrong \"h\" matrix size: " +  to_string(h_size[0]) + ", " + to_string( h_size[1])+"; expected: 1,1");

		if(coef_size[0] != 2 || coef_size[1] != 2)
			throw("wrong \"coef\" matrix size: " +  to_string(coef_size[0]) + ", " + to_string( coef_size[1])+"; expected: 2,2");

		double c = coef(0,0), w = coef(1,0);
		w = (w < 0 || w > 1) ? 0.5 : w;
		c = (c == 0) ? 1 : c;

		matrix d = get_col(coef,1);
		
		#ifdef VERBOSE
			cout<< "[ff5rwPT] c:" << c << ", w: " << w << ", " << print_m_l(d, "d") <<endl;
		#endif


		return ff5rwP(ld + h * d, matrix(w), matrix(c));
	}
	catch (string ex_info)
	{
		throw ("ff5rwPT(...):\n" + ex_info);
	}
}

// matrix ff6T(matrix x, matrix ud1, matrix ud2)
// {
// 	matrix y;
// 	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
// 	return y;
// }


// matrix ff6R(matrix x, matrix ud1, matrix ud2)
// {
// 	matrix y;
// 	int N = 1001;
// 	static matrix X(N, 2);
// 	static bool read = true;
// 	if (read)
// 	{
// 		ifstream S("polozenia.txt");
// 		S >> X;
// 		S.close();
// 		read = false;
// 	}
// 	matrix Y0(4, new double[4]{ 0,0,0,0 });
// 	matrix* Y = solve_ode(df6, 0, 0.1, 100, Y0, ud1, x[0]);
// 	y = 0;
// 	for (int i = 0; i < N; ++i)
// 		y = y + abs(X(i, 0) - Y[1](i, 0)) + abs(X(i, 1) - Y[1](i, 2));
// 	y = y / (2.0 * N);
// 	return y;
// }

// matrix df6(double t, matrix Y, matrix ud1, matrix ud2)
// {
// 	double m1 = 5, m2 = 5, k1 = 1, k2 = 1, F = 1;
// 	double b1 = ud2(0), b2 = ud2(1);
// 	matrix dY(4, 1);
// 	dY(0) = Y(1);
// 	dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
// 	dY(2) = Y(3);
// 	dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;
// 	return dY;
// }