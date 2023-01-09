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

matrix ff5T(matrix h, matrix x, matrix coef);
matrix ff5(matrix x, matrix a = matrix(1), matrix ud2 = NAN);
