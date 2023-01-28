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