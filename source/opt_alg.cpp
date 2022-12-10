#include"opt_alg.h"
#include"user_funs.h"

long double* expansion(long double(*ff)(long double), long double x0, long double d, double alpha, int Nmax, int &f_calls)
{
	try
	{
		long double* p = new long double[2]{ 0,0 };
		//Tu wpisz kod funkcji
		
		long double x1 = x0 + d;

		long double 
			y0 = (*ff)(x0),
			y1 = (*ff)(x1);
		if(y0 == y1){
			p[0] = x0; p[1] = x1;
			return p;
		}
		else if(y0 > y1){
			d= -d;
			x1 = x0 + d;
			y1 = (*ff)(x1);
			if(y0 <= y1){
				p[0] = x0; p[1] = x1;
				return p;
			}
		}

		int i = 0;
		long double 
			x_prev = x0,
			x_crrnt = x1,
			x_next = 0;

		do {
			if (f_calls > Nmax){
				cerr << "Too many f_calls!"<<endl;
				return NULL;
			}
			i++;
			x_prev = x_crrnt;
			x_crrnt = x_next;
			x_next = x0 + pow(alpha, i) * d;
			
		} while (x_crrnt <= x_next);

		if (d > 0) {
			p[0] = x_prev;
			p[1] = x_next;
		}
		else {
			p[0] = x_next;
			p[1] = x_prev;
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		int n = get_dim(x0);

		matrix *simplex = new matrix[n+1];
		simplex[0] = x0;
		for(int i = 1; i <= n; i++){
			double* temp = new double[n];
			for(int j = 0; j < n; j++) temp[j] = 0;
			
			temp[i-1] = lambda;
			simplex[i] = matrix(x0 + matrix(n, temp));
		}

		int max_index =0, min_index=0;
		while(true){

			// ######## set up of indexes ########
			solution y_min(simplex[0]); y_min.fit_fun(ff, ud1, ud2);
			solution y_max(simplex[0]); y_max.fit_fun(ff, ud1, ud2);
			
			for(int i = 1; i <= n; i++){
				solution X(simplex[i]).fit_fun(ff, ud1, ud2);
				if(X.y < y_min.y){
					y_min.y = X.y;
					min_index = i;
				}
				else if (X.y > y_max.y){
					y_max.y = X.y;
					max_index = i;
				}
			}
			//#############################

			//######## mid-point ########
			matrix sum_filtered(0);
			for(int i = 0; i <= n; i++){
				if(i != max_index) sum_filtered += simplex[i];
			}
			matrix mid_point = sum_filtered / matrix(n);

			//######## reflection ########
			matrix refl_point = simplex[max_index] + alpha * (mid_point - simplex[max_index]);
			solution y_refl(refl_point).fit_fun(ff, ud1, ud2);

			//######## expansion ########
			if(y_refl.y < y_min.y){
				matrix exp_point(mid_point + delta * (refl_point - mid_point));
				solution y_exp(exp_point).fit_fun(ff, ud1, ud2);
				if(y_exp.y < y_refl.y) 
					simplex[max_index] = exp_point;

				else
					simplex[max_index] = refl_point;				
			}
			//######## reflection acceptance ########
			else if(y_refl.y < y_max.y)
				simplex[max_index] = refl_point;

			//######## contraction ########
			else {
				matrix contr_point = mid_point + beta * (simplex[max_index] - mid_point);
				solution y_contr(contr_point).fit_fun(ff, ud1, ud2);
				if(y_contr.y < y_exp.y)
					simplex[max_index] = contr_point;

			//######## reduction ########
				else{
					for(int i = 0; i <=n; i++){
						if(i == min_index) continue;
						else{
							simplex[i] = delta * (simplex[i] + simplex[min_index]);
						}
					}
				}
			}

			if(solution.f_calls > Nmax)
				throw "Too many calls.";

			bool stop_condition = false;
			for(int i = 0; i<=n; i++){
				if(geometric_mean(simplex[i] - simplex[min_index]) < epsilon){
					stop_condition = true;
					break;
				}
			}

			if(stop_condition){
				Xopt.x = simplex[min_index];
				x.fit_fun(ff, ud1, ud2);
				break;
			}

		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
