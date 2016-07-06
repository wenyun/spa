// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#include "spa_util.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#define MAX_DIMENSION 3
#define TINY 1.0e-20

int ludcmp(double *a, int n, int *indx, int *d) {
  int i, imax, j, k;
  double big, dum, sum, temp;
  double vv[MAX_DIMENSION];
  int flag;

  flag = 1;
  imax = 0;
  *d = 1;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++) {
      if ((temp = fabs(a[i * n + j])) > big) 
        big = temp;
    }
    if (big == 0.0) {
      flag = 0;
      return flag;
    }
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i * n + j];
      for (k = 0; k < i; k++) {
        sum -=  a[i * n + k] * a[k * n + j];
      }
      a[i * n + j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i * n + j];
      for (k = 0; k < j; k++)
        sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
      if ( (dum=vv[i] * std::abs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax * n + k];
        a[imax * n + k] = a[j * n + k];
        a[j * n + k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j * n + j] == 0.0) {
      a[j * n + j] = TINY;
      printf("TINYTINY!!!\n");
    }
    if (j != n-1) {
      dum = 1.0 / (a[j * n + j]);
      for (i = j + 1; i < n; i++) {
        a[i * n + j] *= dum;
      }
    }
  }

  return flag;
}

void lubksb(double *a, int n, int *indx, double* b) {
  int i;
  int ii = -1;
  int ip;
  int j;
  double sum;

  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0) {
      for (j = ii; j <= i-1; j++) {
        sum -= a[i*n+j]*b[j];
      }
    } else if (sum) {
      ii = i;
    }
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < n; j++) {
      sum -= a[i * n + j] * b[j];
    }
    b[i] = sum / a[i*n+i];
  }
}

void vector_add(double *x, double *y, double factor, int n) {
  int i;
  for(i = 0; i < n; i++) {
    x[i] += y[i] * factor;
  }
}

void vector_add_to_new(double *z,
                       double *x,
                       double *y,
                       double factor,
                       int n) {
  int i;
  for(i = 0; i < n; i++) {
    z[i] = x[i] + y[i] * factor;
  }
}

void vector_scale(double *x, double factor, int n) {
  int i;
  for(i = 0; i < n; i++) {
    x[i] *= factor;
  }
}

void vector_init(double *x, int n, double v) {
  int i;
  for(i = 0; i < n; i++) {
    x[i] = v;
  }
}

double vector_inner_product(double *x, double *y, int n) {
  double prod = 0;
  int i;
  for(i = 0; i < n; i++) {
    prod += x[i] * y[i];
  }

  return prod;
}

void vector_out_product(double *x, double *y, int n) {
  int i, j;
  for(i = 0; i < n; i++) {
    for(j = 0; j < i; j++) {
      x[i * n + j] = y[i] * y[j];
      x[j * n + i] = x[i * n + j];
    }
    x[i * n + i] = y[i] * y[i];
  }
}

void vector_copy(double *x, double *y, int n) {
  memcpy(x, y, sizeof(double)*n);
}

void vector_normalize(double *x, int n) {
  vector_scale(x, 1 / sqrt(vector_inner_product(x, x, n)), n);
}

double vector_std(double *x, int n) {
  int i;
  double m, s;
  
  m = 0;
  for(i = 0; i < n; i++)
    m += x[i];
  
  m /= n;

  s = 0;
  for(i = 0; i < n; i++) {
    s += pow(x[i] - m, 2);
  }

  s = sqrt(s);
  return s;
}

