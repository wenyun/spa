// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#pragma once

/* *
 * numerical functions 
 * */

// back substitution using LU decomposition
void lubksb(double* a, int n, int* indx, double* b);
// LU decomposition
int ludcmp(double* a, int n, int* indx, int* d);

/* *
 * vector operation functions
 * */

// x = x + factor * y
void vector_add(double* x, double* y, double factor, int n);
// z = x + factor * y
void vector_add_to_new(double* z,
                       double* x,
                       double* y,
                       double factor,
                       int n);
// x = factor * x
void vector_scale(double* x, double factor, int n);
// x = v
void vector_init(double* x, int n, double v);
// return x' * y
double vector_inner_product(double* x, double* y, int n);
// x = y * y'
void vector_out_product(double* x, double* y, int n);
// x = y
void vector_copy(double* x, double* y, int n);
// return standard deviation for x
double vector_std(double* x, int n);
// normalize x, ||x|| = 1 as a result
void vector_normalize(double* x, int n);
