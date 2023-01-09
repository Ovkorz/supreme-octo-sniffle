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

matrix ff5a(matrix x, matrix ud1, matrix ud2){
	try{
		int *x_size = get_size(x);
		if(x_size[0] != 2 || x_size [1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
	
		return matrix(
			a * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2))
		)
	}
	catch (string ex_info)
	{
		throw ("ff5a(...):\n" + ex_info);
	}
}

matrix ff5b(matrix x, matrix ud1, matrix ud2){
	try{
		int *x_size = get_size(x);
		if(x_size[0] != 2 || x_size [1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
	
		return matrix(
			1/a * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2))
		)
	}
	catch (string ex_info)
	{
		throw ("ff5b(...):\n" + ex_info);
	}
}

matrix ff5Ta(matrix h, matrix x, matrix coef){

	try{

		int *x_size = get_size(x),
			*h_size = get_size(h),
			*a_size = get_size(coef);

		if(	x_size[0] != 2 || x_size[1] != 1)
			throw("wrong \"x\" matrix size: " +  to_string(x_size[0]) + ", " + to_string( x_size[1])+"; expected: 2,1");
		
		if(	h_size[0] != 1 || h_size[1] != 1)
			throw("wrong \"h\" matrix size: " + to_string(h_size[0]) + ", " + to_string( h_size[1])+"; expected: 1,1");

		if(	a_size[0] != 2 || a_size[1] != 2)
			throw("wrong \"coef\" matrix size: " + to_string(a_size[0]) + ", " + to_string( a_size[1])+"; expected: 2,2");

		double a = coef(0,0);
		matrix d = get_col(coef,1);
		
		// matrix dY(
		// 	a * (pow(x(0) + h(0) *d(0) - 2, 2) + pow(x(1) + h(0) * d(1) - 2, 2))
		// );

		return ff5a(x + h * d);

	}
	catch (string ex_info)
	{
		throw ("ff5Ta(...):\n" + ex_info);
	}

}

matrix ff5Tb(matrix h, matrix x, matrix coef){

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

		double a = coef(0,0);
		matrix d = get_col(coef,1);
		
		return ff5b(x + h * d);
	}
	catch (string ex_info)
	{
		throw ("ff5Tb(...):\n" + ex_info);
	}
}
