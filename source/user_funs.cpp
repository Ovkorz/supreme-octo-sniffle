#include"user_funs.h"
#include<cmath>

double geometric_mean (const matrix &A){
    double sum_of_squares = 0;
    int n = get_dim(A);

    for(int i = 0; i<n; i++){
        sum_of_squares += pow(A[i], 2);
    }

    return sqrt(sum_of_squares);
}