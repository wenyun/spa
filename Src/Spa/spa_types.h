// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#pragma once

#include <cstdlib>
#include <cstdio>
#include <string>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

// system macros
#define MAX_DIMENSION 3
#define INF 9e20
// verbose level
#define QUIET 0
#define SHORT 1
#define MEDIUM 2
#define WORDY 3
#define CRAZY 4

// genetics macros
#define MISSING_ALLELE '0'
#define MISSING -1
#define HOMO_MAJOR 0
#define HOMO_MINOR 2
#define HETER 1

// dimension
#define PLANE 2
#define GLOBE 3

// generation
#define SELF 1
#define PARENT 2

// program mode
#define COEF_ONLY 1
#define LOCT_ONLY 2
#define BOTH 3
#define LOCT_GLOBE 4
#define ADMIXED 5

struct spa_parameter {
  
  const char* bfile;
  const char* pfile;
  const char* tfile;
  const char* gfile;
  const char* mfile;
  const char* ilfile;
  const char* olfile;
  const char* imfile;
  const char* omfile;
  
  int generation;
  int dimension;
  int max_iter;
  double step_epsilon;

  int max_sub_iter;
  double epsilon;
  double alpha;
  double beta;

  int verbose;

  int large_step_since;
};

struct snp_info_struct {
  char* chromosome;
  char* snp_id;
  char* morgan;
  int position;
  char snp_minor;
  char snp_major;
};

struct individual_info_struct {
  char* family_id;
  char* individual_id;
  char* paternal_id;
  char* maternal_id;
  char* sex;
  char* phenotype;
};

struct spa_model {
  int n_individual;
  int n_snp;

  snp_info_struct* snp_info;
  individual_info_struct* individual_info;
  
  double** coef_a;
  double* coef_b;
  double** x;

  double* coef_a_space;
  double* x_space;

  double* score;
};

struct spa_data {
  
  int n_individual;
  int n_snp;
  
  snp_info_struct* snp_info;
  individual_info_struct* individual_info;
  
  char** genotype;
  char* genotype_space;
};

void spa_error_exit(const char* msg);
void spa_warning(const char* msg);
void spa_message(const char* msg, 
                 const int level,
                 const spa_parameter *param);
