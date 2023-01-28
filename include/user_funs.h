#pragma once

#include"ode_solver.h"
#include"solution.h"
#include"matrix.h"

matrix df1(double t, matrix Y, matrix ud1, matrix ud2);
matrix df2(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix ff2T(matrix x, matrix ud1, matrix ud2);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix ff1T(matrix x, matrix ud1, matrix ud2);

matrix gf4R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff4R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix Hf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix gf4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff4T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff5_values(matrix x, matrix a, matrix w = NAN);
matrix ff5T(matrix h, matrix x, matrix coef);
matrix ff5(matrix x, matrix a = matrix(1), matrix w = matrix(0.5));

matrix ff5rwP_values(matrix ld);
matrix ff5rwPT(matrix h, matrix ld, matrix coef);