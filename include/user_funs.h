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

matrix ff5a(matrix x, matrix ud1, matrix ud2);
matrix ff5b(matrix x, matrix ud1, matrix ud2);

matrix ff5Ta(matrix h, matrix x, matrix coef);
matrix ff5Tb(matrix h, matrix x, matrix coef);