#pragma once

#include"ode_solver.h"


matrix df1(double t, matrix Y, matrix ud1, matrix ud2);
matrix df2(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff2R(matrix x, matrix ud1, matrix ud2);
matrix ff2T(matrix x, matrix ud1, matrix ud2);
matrix ff1R(matrix x, matrix ud1, matrix ud2);
matrix ff1T(matrix x, matrix ud1, matrix ud2);
