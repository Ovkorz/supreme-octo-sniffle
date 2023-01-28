#include"opt_alg.h"
#include"user_funs.h"

using namespace std;

// class Fibonacci{
// public:
//     int a = 0, b = 1, c = 1, k = 3;
//     int generate(int);
// };

int Fibonacci::generate(int n){
	if(n <= 0 ) return -1;

	int i = 2;

	if(n == k) return c;
	else if (n == k - 1) return b;
	else if (n == k - 2) return a;
	else if(n > k) i = k;
	else {
		a = 0;
		b = 1;
		i = 2;
	}

	while(i < n){
		c = a + b;
		a = b;
		b = c;

		i++;
	}

    return c;
}

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
		Fibonacci fib;


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

	#ifdef VERBOSE
	cout<< "Rosenbrock optimization..." <<endl;
	#endif

	try
	{
		solution Xopt;
		Xopt.ud = trans(x0);
		solution XB(x0), X;
		int n = get_dim(XB);
		matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);

		#ifdef VERBOSE
			cout << "Initial state: "<< print_m_l(x0, "x0") << ", " << print_m_l(s0, "s0") << ", " << print_m_l(ud1, "ud1") << ", " << print_m_l(ud2, "ud2")
			<< ", epsilon: " << epsilon	<<endl;
		#endif

		
		int counter = 0;
		while (true)
		{
			X.x = XB.x;
			XB.fit_fun(ff, ud1, ud2);

			#ifdef VERBOSE
				cout<<"\nIteration "<<counter << ", ";
				cout << "init: " << print_m_l(XB, "XB") <<", " << print_m_l(s, "step") << "; ";
			#endif
			
			// wstepne przeszukiwanie do uzupelnienia

			for(int i = 0; i < n; i++){

				#ifdef VERBOSE
					cout <<"shift: i = " << i << ", ";
				#endif

				X.x(i) = X.x(i) + s(i);
				X.fit_fun(ff, ud1, ud2);

				if(X.y >= XB.y){
					X.x(i) = XB.x(i);
					p(i) = p(i) + 1;
					s(i) = s(i)*(-beta);

					#ifdef VERBOSE
						cout<<"bad step, "
							<< print_m_l(p, "p") <<", "<< print_m_l(s, "s")<< ", ";
					#endif
				}
				else {
					l(i) = l(i) + s(i);
					s(i) = s(i) * alpha;

					#ifdef VERBOSE
						cout<<"good step, "
							<< print_m_l(l, "l") <<", "<< print_m_l(s, "s")<< ", ";
					#endif
				}
			}
			
			Xopt.ud.add_row(trans(XB.x));
			XB.x = X.x;

			bool change = true;
			for (int i = 0; i < n; ++i)
				if (l(i) == 0 || p(i) == 0)
				{
					change = false;
					#ifdef VERBOSE
						cout<<"no change, ";
				
					#endif

					break;
				}
			if (change)
			{
				#ifdef VERBOSE
					cout <<"step base change, ";
				#endif

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

			if(X.f_calls > Nmax) throw string("Too many f calls!");

			double max_step = abs(s(0));
			for(int i = 1; i < n; i++){
				if(abs(s(i))> max_step) max_step = abs(s(i));
			}


			if(max_step < epsilon){
				#ifdef VERBOSE
					cout<< "max_step: " <<max_step <<", break" << endl;
				#endif

				XB.fit_fun(ff, ud1, ud2);		
				XB.ud = Xopt.ud;		
				break;
			}

			counter++;
		// warunki stopu
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(
	matrix(*ff)(matrix, matrix, matrix), 
	solution(*opt_f)(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix coef,  double epsilon, int Nmax, matrix ud1, matrix ud2),
	matrix x0, matrix coef, double epsilon, int Nmax, matrix innercoef, matrix ud2)
{
	// coef
	// [[c, w],
	//  [dc, 0]]

	#ifdef VERBOSE
		cout<< "Penalty function optimization..." <<endl;
	#endif

	try {
		solution Xopt(x0), X;
		int n = get_dim(Xopt);
		//Tu wpisz kod funkcji
		double c = coef(0,0), dc = coef(1,0), w = coef(0,1);

		#ifdef VERBOSE
			cout<<"[Penalty] Initial state: " << print_m_l(x0, "x0") << ", " << print_m_l(coef, "coef") 
				<< ", c: " << c << ", dc: " << dc << ", w: " << w << ", epsilon: " << epsilon <<endl;
		#endif

		int counter = 0;
		while(Xopt.f_calls < Nmax){
			matrix c_w(2,1); c_w(0,0) = c; c_w(1) = w;

			#ifdef VERBOSE
				cout<< "[Penalty] Iteration " << counter << ", " << print_m_l(Xopt.x, "X") << ", " << print_m_l(c_w, "c_w") << endl;
			#endif

			X = opt_f(ff, Xopt.x, innercoef, epsilon, Nmax, c_w, ud2);

			#ifdef VERBOSE
				cout<< "[Penalty] Inner optimization result " << print_m_l(X, "X") <<endl;
			#endif


			matrix diff = X.x - Xopt.x; double prox = 0;
			for(int i = 0; i < n; i++){
				prox += pow(diff(i),2);
			}
			prox = sqrt(prox);
			
			if(prox < epsilon){
				#ifdef VERBOSE
					cout <<	"[Penalty] Given error tolerance reached - prox: " << prox << ", epsilon: " << epsilon
						 << "; optimization finished with result " << print_m_l(X, "Xopt") << endl;
				#endif

				return X;
			}

			Xopt = X;
			c *= dc;

			counter++;
		}
		throw string("Too many f_calls!");

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
		// //Tu wpisz kod funkcji

		// int n = get_dim(x0);

		// matrix *simplex = new matrix[n+1];
		// simplex[0] = x0;
		// for(int i = 1; i <= n; i++){
		// 	double* temp = new double[n];
		// 	for(int j = 0; j < n; j++) temp[j] = 0;
			
		// 	temp[i-1] = lambda;
		// 	simplex[i] = matrix(x0 + matrix(n, temp));
		// }

		// int max_index =0, min_index=0;
		// while(true){

		// 	// ######## set up of indexes ########
		// 	solution y_min(simplex[0]); y_min.fit_fun(ff, ud1, ud2);
		// 	solution y_max(simplex[0]); y_max.fit_fun(ff, ud1, ud2);
			
		// 	for(int i = 1; i <= n; i++){
		// 		solution X(simplex[i]).fit_fun(ff, ud1, ud2);
		// 		if(X.y < y_min.y){
		// 			y_min.y = X.y;
		// 			min_index = i;
		// 		}
		// 		else if (X.y > y_max.y){
		// 			y_max.y = X.y;
		// 			max_index = i;
		// 		}
		// 	}
		// 	//#############################

		// 	//######## mid-point ########
		// 	matrix sum_filtered(0);
		// 	for(int i = 0; i <= n; i++){
		// 		if(i != max_index) sum_filtered += simplex[i];
		// 	}
		// 	matrix mid_point = sum_filtered / matrix(n);

		// 	//######## reflection ########
		// 	matrix refl_point = simplex[max_index] + alpha * (mid_point - simplex[max_index]);
		// 	solution y_refl(refl_point).fit_fun(ff, ud1, ud2);

		// 	//######## expansion ########
		// 	if(y_refl.y < y_min.y){
		// 		matrix exp_point(mid_point + delta * (refl_point - mid_point));
		// 		solution y_exp(exp_point).fit_fun(ff, ud1, ud2);
		// 		if(y_exp.y < y_refl.y) 
		// 			simplex[max_index] = exp_point;

		// 		else
		// 			simplex[max_index] = refl_point;				
		// 	}
		// 	//######## reflection acceptance ########
		// 	else if(y_refl.y < y_max.y)
		// 		simplex[max_index] = refl_point;

		// 	//######## contraction ########
		// 	else {
		// 		matrix contr_point = mid_point + beta * (simplex[max_index] - mid_point);
		// 		solution y_contr(contr_point).fit_fun(ff, ud1, ud2);
		// 		if(y_contr.y < y_exp.y)
		// 			simplex[max_index] = contr_point;

		// 	//######## reduction ########
		// 		else{
		// 			for(int i = 0; i <=n; i++){
		// 				if(i == min_index) continue;
		// 				else{
		// 					simplex[i] = delta * (simplex[i] + simplex[min_index]);
		// 				}
		// 			}
		// 		}
		// 	}

		// 	if(solution.f_calls > Nmax)
		// 		throw "Too many calls.";

		// 	bool stop_condition = false;
		// 	for(int i = 0; i<=n; i++){
		// 		if(geometric_mean(simplex[i] - simplex[min_index]) < epsilon){
		// 			stop_condition = true;
		// 			break;
		// 		}
		// 	}

		// 	if(stop_condition){
		// 		Xopt.x = simplex[min_index];
		// 		x.fit_fun(ff, ud1, ud2);
		// 		break;
		// 	}

		// }

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
		ofstream SD;
		SD.open("SD.txt");
		solution Xopt;
		Xopt.ud = trans(x0);
		int n = get_len(x0);
		solution X0, X1;
		X0.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;

		while (true)
		{
			d = -X0.grad(gf, ud1, ud2);
			if (h0 < 0)
			{
				P.set_col(X0.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X0.x + h.x * d;
			}
			else
				X1.x = X0.x + h0 * d;

			Xopt.ud.add_row(trans(X1.x));

			if (norm(X0.x - X1.x) < epsilon)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
			X0 = X1;
			SD << X1.x(0) << " " << X1.x(1) << endl;
		}
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
		Xopt.ud = trans(x0);
		int n = get_len(x0);
		solution X0, X1;
		X0.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab{};
		double beta;
		d = -X0.grad(gf);
		while (true)
		{
			if (h0 < 0)
			{
				P.set_col(X0.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X0.x + h.x * d;
			}
			else
				X1.x = X0.x + h0 * d;

			Xopt.ud.add_row(trans(X1.x));

			if (norm(X1.x - X0.x) < epsilon)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}
			if (std::max(solution::f_calls, solution::g_calls) > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
			X1.grad(gf);
			beta = pow(norm(X1.g), 2) / pow(norm(X0.g), 2);
			d = -X1.g + beta * d;
			X0 = X1;
			//            std::cout << "X0.x = " << X0.x << std::endl;
		}
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
		ofstream Newton;
		Newton.open("Newton.txt");
		solution Xopt;
		Xopt.ud = trans(x0);
		int n = get_len(x0);
		solution X0, X1;
		X0.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* ab;
		while (true)
		{
			d = -inv(X0.hess(Hf, ud1, ud1)) * X0.grad(gf, ud1, ud2);
			if (h0 < 0)
			{
				P.set_col(X0.x, 0);
				P.set_col(d, 1);
				ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
				X1.x = X0.x + h.x * d;
			}
			else
				X1.x = X0.x + h0 * d;

			Xopt.ud.add_row(trans(X1.x));

			if (norm(X0.x - X1.x) < epsilon)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 0;
				break;
			}
			if (solution::f_calls > Nmax || solution::g_calls > Nmax)
			{
				Xopt = X1;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				break;
			}
			X0 = X1;
			Newton << X1.x(0) << " " << X1.x(1) << endl;
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{

	#ifdef VERBOSE
		cout<<"Golden ratio optimization..."<<endl;
	#endif

	try
	{
		solution Xopt, X;
		//Tu wpisz kod funkcji
		matrix lower(a), upper(b), x1(a), x2(b);
		double ratio = 2/ (3+ sqrt(5));

		Xopt.x = (1-ratio) * lower + ratio * upper; // x1
		Xopt.fit_fun(ff, ud1, ud2);

		#ifdef VERBOSE
			cout<<"[Golden] Initial state: " << print_m_l(Xopt, "Xopt")<< ", "<< print_m_l(lower, "lower") << ", " << print_m_l(upper, "upper") 
				<< ", epsilon: " << epsilon	<<endl;
		#endif

		int counter = 0;
		while(abs((upper - lower)(0)) > 2 * epsilon){

			X.x = (1 - ratio) * Xopt.x + ratio * upper; // x2
			X.fit_fun(ff, ud1, ud2);
			Xopt.fit_fun(ff, ud1, ud2);
			
			#ifdef VERBOSE
				cout<<"[Golden] Iteration "<< counter<<", " << print_m_l(lower, "lower")<<", "<<print_m_l(upper, "upper")<<", range = "<< abs((upper - lower)(0))<< ", " 
				<< print_m_l(Xopt, "Xopt")<< ", " << print_m_l(X, "X") << ", ";
			#endif

			if(Xopt.y < X.y){

			#ifdef VERBOSE
				cout<< 	"Xopt.y <  X.y";
			#endif

				upper = lower;
				lower = X.x;
			}
			else {
				#ifdef VERBOSE
					cout<< 	"Xopt.y >= X.y";
				#endif
				lower = Xopt.x;
				Xopt = X;
			}
			
			counter++;

			#ifdef VERBOSE
				cout<<endl;
			#endif

			if(X.f_calls > Nmax) throw string("Too many f calls!");
		}	

		Xopt.x = (Xopt.x + X.x)/2;
		Xopt.fit_fun(ff, ud1, ud2);

		#ifdef VERBOSE
				cout<< "[Golden] Optimization finished with result: " << print_m_l(Xopt, "Xopt") << endl << endl;
		#endif

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix coef ,double epsilon, int Nmax, matrix a, matrix ud2)
{
	// matrix a
	// [[a],
	// 	[w]]

	#ifdef VERBOSE
			cout<< "Powell optimization..."<<endl;
	#endif

	try
	{
		solution Xopt(x0);
		//Tu wpisz kod funkcji
		// solution P(x0);
		int n = get_dim(Xopt);
		matrix d(n,n), X(x0), P(x0);
		solution h(matrix(1));

		//for golden method
		const double lower = -5, upper = 5;

		matrix gcoef(2,2);
		gcoef(0,0) = a(0,0);
		gcoef(1,0) = a(1,0);
		// gcoef matrix:
		// [[a, d1],
		//  [w, d2]]

		for(int i = 0; i < n; i++){
			d(i,i) = 1;
		}
		
		#ifdef VERBOSE
			cout<< "[Powell] Initial state: " << print_m_l(X, "X") << ", " << print_m_l(d, "d") << ", "
				<< print_m_l(h, "h") << ", " << "epsilon: " << epsilon << endl;
		#endif

		int counter = 0;
		while(h.f_calls < Nmax){

			#ifdef VERBOSE
				cout << "[Powell] Iteration  " << counter << ", " << print_m_l(X, "X") <<", " <<print_m_l(d, "d") << endl;
			#endif

			P = X;

			for(int i = 0; i < n; i++){
				matrix d_i = trans(get_row(d, i));
				gcoef.set_col(d_i,1);

				#ifdef VERBOSE
					cout <<	"[Powell] Initializing directed optimization with " << print_m_l(P, "x0") << ", " << print_m_l(gcoef, "coef") << endl;
				#endif

				h = golden(ff, lower, upper, epsilon, Nmax, P, gcoef);

				P = P + (h.x * d_i);
				#ifdef VERBOSE
					cout << "[Powell] Shift nr " << i << ": " << print_m_l(h, "h") << ", " << print_m_l(d_i, "d") << ", step " << print_m_l(h.x * d_i, "h*d") << ", result " << print_m_l(P, "P") << endl;
				#endif
			}

			matrix diff = P - X; double prox = 0;
			for(int i = 0; i < n; i++){
				prox += pow(diff(i),2);
			}
			prox = sqrt(prox);

			if(prox < epsilon){
				Xopt = solution(X);
				Xopt.y = h.y;

				#ifdef VERBOSE
					cout <<	"[Powell] Given error tolerance reached - prox: " << prox << ", epsilon: " << epsilon
						 << "; optimization finished with result " << print_m_l(Xopt, "Xopt") << endl;
				#endif
				
			
				return Xopt;
			}

			for(int i = 1; i < n; i++){
				d.set_row(get_row(d, i), i-1);
			}

			d.set_row(trans(P - X), n-1);
			#ifdef VERBOSE
				cout << "[Powell] Updated step vectors: " << print_m_l(d, "d") << endl;
			#endif

			gcoef.set_col(
				trans( 
					get_row(d, n-1) 
				)
			, 1);
	
			// h = Rosen(ff, h.x, s, alpha, beta, epsilon, Nmax, X, coef);

			#ifdef VERBOSE
				cout << "[Powell] Initializing additional directed optimization with " << print_m_l(P, "x0")  << ", " << print_m_l(gcoef, "coef") << endl;
			#endif
			
			h = golden(ff, lower, upper, epsilon, Nmax, P, gcoef);
			X = P + (h.x * trans(get_row(d, n-1)));

			#ifdef VERBOSE
					cout << "[Powell] Additional shift: " << print_m_l(h, "h") << ", " << print_m_l(trans(get_row(d, n-1)), "d") << ", step " << print_m_l(h.x * trans(get_row(d, n-1)), "h*d") << ", result " << print_m_l(X, "X") << endl;
			#endif
			counter++;
		}

		throw string("Too many f calls!");
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
		solution* P = new solution[mi + lambda];
		solution* Pm = new solution[mi];
		matrix IFF(mi, 1), temp(N, 2);
		double r, s, s_IFF;
		double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
		int j_min;
		for (int i = 0; i < mi; ++i)
		{
			P[i].x = matrix(N, 2);
			for (int j = 0; j < N; ++j)
			{
				P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * m2d(rand_mat()) + limits(j, 0);
				P[i].x(j, 1) = sigma0(j);
			}
			if (P[i].fit_fun(ff, ud1, ud2) < epsilon)
			{
				Xopt = P[i];
				Xopt.flag = 0;
				delete[]P;
				delete[]Pm;
				return Xopt;
			}
		}
		while (true)
		{
			s_IFF = 0;
			for (int i = 0; i< mi ; ++i)
			{
				IFF(i) = 1 / m2d(P[i].y);
				s_IFF += IFF(i);
			}
			for (int i = 0; i<lambda ; ++i)
			{
				r = s_IFF * m2d(rand_mat());
				s = 0;
					for (int j = 0; j<mi; ++j)
					{
						s += IFF(j);
						if (r<=s )
						{
							P[mi + i] = P[j] ;
							break;
						}
					}
			}
			for (int i = 0; i< lambda ; ++i)
			{
				r = m2d(rand_mat());
				for (int j = 0; j < N; ++j)
				{
					P[mi + i].x(j, 1) *= exp(tau1 * r +tau *m2d(randn_mat()));
					P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * m2d(randn_mat());
				}
			}
			for (int i = 0; i< lambda; i += 2)
			{
				r = m2d(rand_mat());
				temp = P[mi+i].x ;
				P[mi + i].x = r * P[mi+i].x + (1-r) * P[mi+i+1].x;
				P[mi + i + 1].x = r * P[mi+i+1].x + (1-r) * temp;
			}
			for (int i = 0; i < lambda; ++i)
			{
				if (P[mi+i].fit_fun(ff,ud1, ud2) < epsilon)
				{
		
					Xopt = P[mi+i] ;
					Xopt.flag = 0;
					delete[]P;
					delete[]Pm;
					return Xopt;
				}
			}
			for (int i = 0; i < mi; ++i)
			{
				j_min = 0;
				for (int j = 1; j < mi + lambda ; ++j)
					if (P[j_min].y > P[j].y )
						j_min = j;
				Pm[i] = P[j_min] ;
				P[j_min].y = 1e10;
			}
			for (int i = 0; i < mi; ++i)
				P[i] = Pm[i];
			if (solution::f_calls > Nmax)
			{
				Xopt = P[0];
				Xopt.flag = 1;
				break;
			}
		}
		delete[]P;
		delete[]Pm;
		return Xopt;
	}
	catch(string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}