#include"opt_alg.h"
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
		Xopt.ud = trans(x0);
		solution XB(x0), X;
		int n = get_dim(XB);
		matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
		XB.fit_fun(ff, ud1, ud2);

		while (true)
		{
			// wstepne przeszukiwanie do uzupelnienia
			X.x = XB.x;
			XB.fit_fun(ff, ud1, ud2);

			for(int i = 0; i < n; i++){

				X.x[i] += s[i] * D[i];		
				X.fit_fun(ff, ud1, ud2);
				if(X.y[i] >= XB.y[i]){
					X.x[i] = XB.x[i];
					p[i]++;
					s[i] *= -beta;
				}
				else {
					l[i] += s[i];
					s[i] *= alpha;
				}

			}
			
			Xopt.ud.add_row(trans(XB.x));
			XB.x = X.x;

			bool change = true;
			for (int i = 0; i < n; ++i)
				if (l(i) == 0 || p(i) == 0)
				{
					change = false;
					break;
				}
			if (change)
			{
				matrix Q(n, n), v(n, 1);
				for (int i = 0; i < n; ++i)
					for (int j = 0; j <= i; ++j)
						Q(i, j) = l(i);
				Q = D * Q;
				v = Q[0] / norm(Q[0]);
				D.set_col(v, 0);
				for (int i = 1; i < n; ++i)
				{
					matrix temp(n, 1);
					for (int j = 0; j < i; ++j)
						temp = temp + trans(Q[i]) * D[j] * D[j];
					v = (Q[i] - temp) / norm(Q[i] - temp);
					D.set_col(v, i);
				}
				s = s0;
				l = matrix(n, 1);
				p = matrix(n, 1);
			}

			if(X.f_calls > Nmax) throw "Too many f calls!";

			double max_step = s[0];
			for(int i = 1; i < n; i++){
				if(s[i]> max_step) max_step = s[i];
			}

			if(max_step < epsilon) break;
		// warunki stopu
		}
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
