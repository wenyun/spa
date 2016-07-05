// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#include "spa.h"

#include "spa_io.h"
#include "spa_util.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

using std::string;
using std::map;

// system macro
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define MAX_DIMENSION 3
#define MAX_GENERATION 2
#define MAX_MESSAGE 1024
#define MAX_LINE_INIT 1024

void exit_with_help() {
  printf(
  "Usage: spa [options]\n"
  "options: \n"
  "--bfild bed_prefix : prefix of .bed, .bim and .fam file\n"
  "--pfile ped_prefix : prefix of .ped and .map file\n"
  "--tfile tped_prefix : prefix of .tped and .tfam file\n"
  "--gfile genotype_file : genotype file\n"
  "--mfile 23andme_file : 23andme genotype file\n"
  "--location-input location_file : known locations\n"
  "--model-input model_file : known slope functions\n"
  "--location-output location_file : output file for individual locations\n"
  "--model-output model_file: output file for slope function coefficients" 
    "and SPA score\n"
  "-n generation: number of locations\n"
  "-k dimesion : dimensions of spatial analysis\n"
  "-e epsilon : set tolerance of termination criterion (default 0.01)\n"
  "-r tradeoff : set optimization epsilon tolerance (default 1e-6)."
    "Larger value makes the program run faster but poor accuracy\n"
  "-v verbose : verbose level\n"
  );
  exit(1);
}

void print_version_information() {
  printf(
"@----------------------------------------------------------@\n"
"|         SPA!       |      v1.13      |    4/APR/2012     |\n"
"|----------------------------------------------------------|\n"
"|  (C) 2012 Wen-Yun Yang, GNU General Public License, v2   |\n"
"|----------------------------------------------------------|\n"
"|  For documentation, citation & bug-report instructions:  |\n"
"|             http://genetics.cs.ucla.edu/spa/             |\n"
"@----------------------------------------------------------@\n"
"\n"
);
}

// IO variables
char* line = NULL;
int max_line_len = 1024;
static int* indx;  // less than 10 dimensions Newton's method
static char* error_buffer;

void initialize() {
  indx = Malloc(int, MAX_DIMENSION);
  error_buffer = Malloc(char, MAX_MESSAGE);
  max_line_len = MAX_LINE_INIT;
  line = Malloc(char, max_line_len);
}

void finalize() {
  free(indx);
  free(error_buffer);
  free(line);
}

void set_default_parameter(spa_parameter *param) {
  param->bfile = NULL;
  param->pfile = NULL;
  param->tfile = NULL;
  param->gfile = NULL;
  param->mfile = NULL;
  param->ilfile = NULL;
  param->olfile = NULL;
  param->imfile = NULL;
  param->omfile = NULL;

  // default parameters
  param->generation = 1;
  param->dimension = 2;
  param->max_iter = 1000;
  param->step_epsilon = 0.01;
  param->max_sub_iter = 100;
  param->alpha = 0.01;
  param->beta = 0.5;
  param->epsilon = 1e-6;
  param->verbose = SHORT;

  param->large_step_since = 15;
}

int main (int argc, char** argv)  {
  
  spa_parameter param;
  spa_model model;
  spa_data geno;

  print_version_information();

  if(argc == 1) {
    exit_with_help();
  }

  initialize();
  set_default_parameter(&param);

  parse_input_parameters(argc, argv, &param);

  if(!parameter_sanity_check(&param)) {
    spa_error_exit("Parameter sanity check failed\n");
  }

  // read input file
  if (param.gfile) {
    read_gfile(param.gfile, &geno, &param);
  } else if (param.bfile) {
    read_bedfile(param.bfile, &geno, &param);  
  } else if (param.tfile) {
    read_tpedfile(param.tfile, &geno, &param);
  } else if (param.pfile) {
    read_pedfile(param.pfile, &geno, &param);
  } else if (param.mfile) {
    read_mfile(param.mfile, &geno, &param);
  }
  
  sprintf(line, 
          "%d individual %d snps are read ...", 
          geno.n_individual,
          geno.n_snp);

  spa_message(line, SHORT, &param);

  if (!param.imfile && !param.ilfile) {
    data_sanity_check(&geno);
    spa_message("Unsupervised learning ...", SHORT, &param);
    copy_info(&model, &geno, BOTH);
    allocate_model_for_bootstrap(&model, &param); 
    spa_optimize(&model, &geno, &param, BOTH);   
  } else if(param.imfile) {
    read_model_imfile(param.imfile, &model, &param);
    spa_message("Slope functions are given ...", SHORT, &param);
    
    // make sure model and genotype data have the same number of 
    // SNPs, by the same order and the same minor allele...
    model_sanity_check(&model);
    model_data_consistent(&model, &geno, &param, COEF_ONLY);
    
    switch(param.generation) {
      case SELF: 
        spa_optimize(&model, &geno, &param, LOCT_ONLY); 
        break;
      case PARENT:
        spa_message("Admixed individual localization ...", SHORT, &param);
        spa_optimize(&model, &geno, &param, ADMIXED);
        break;
    }
  } else if(param.ilfile) {
    data_sanity_check(&geno);
    read_location_ilfile(param.ilfile, &model, &param); 
    spa_message("Individual locatioins are given ...", SHORT, &param);
    // make sure the model and genotype data have the same number of
    // individuals, by the same order
    model_data_consistent(&model, &geno, &param, LOCT_ONLY);
    spa_optimize(&model, &geno, &param, COEF_ONLY);  
  }
  
  if(param.olfile) {
    spa_message("Individual location estimated", MEDIUM, &param);
    write_location_olfile(param.olfile, &model, &param);

    if (param.mfile) {
      write_html_location_olfile(param.olfile, &model, &param);
      sprintf(line,
              "CHECK YOUR ANCESTRY ORIGIN BY CLICKING %s.html",
              param.olfile);
      spa_message(line, SHORT, &param);
    }
  }
  if(param.omfile) {
    spa_message("Slope functions estimated", MEDIUM, &param);
    spa_selection(&model, &param);
    write_model_omfile(param.omfile, &model, &param);
  }
  
  free_model(&model);
  free_data(&geno);
  finalize();
  
  return 1;
}

void spa_optimize(spa_model *model,
                  const spa_data *geno,
                  const spa_parameter *param,
                  const int mode) {
  int i, iter;
  double *old_x;
  double step, norm, objective;

  if(mode == BOTH) {
    initialize_random_location(model, geno, param);
    old_x = Malloc(double, geno->n_individual * param->dimension);
  
    for(iter = 0; iter < param->max_iter; iter++) {
      for(i = 0; i < geno->n_snp; i++) {
        spa_sub_optimize(model, geno, param, i, COEF_ONLY);
      }
      
      vector_copy(old_x, model->x_space, geno->n_individual *
                                         param->dimension);
      norm = 
        sqrt(vector_inner_product(old_x, 
                                  old_x, 
                                  geno->n_individual * param->dimension));
  
      for(i = 0; i < geno->n_individual; i++) {
        switch(param->dimension) {
          case PLANE:
            spa_sub_optimize(model, geno, param, i, LOCT_ONLY);
            break;
          case GLOBE:
            spa_sub_optimize(model, geno, param, i, LOCT_GLOBE);
            break;
        }
      }
    
      vector_add(old_x, model->x_space, -1, geno->n_individual *
                                            param->dimension);
      step = sqrt(vector_inner_product(old_x, old_x, geno->n_individual *
                                                     param->dimension));
      step /= norm;

      objective = spa_objective(model, geno, param);
      sprintf(line, 
              "Iter %d: -log likelihood = %.5f, step = %.5f",
              iter+1, 
              objective, 
              step);
      spa_message(line, SHORT, param);
      
      if(step < param->step_epsilon) {
        break;
      }
    }
  
    free(old_x);
  } else if (mode == COEF_ONLY)  {
    for(i = 0; i < geno->n_snp; i++) {
      spa_sub_optimize(model, geno, param, i, COEF_ONLY);
    }
  } else if (mode == LOCT_ONLY) {
    initialize_random_location(model, geno, param);    
    for(i = 0; i < geno->n_individual; i++) {
      switch(param->dimension) {
        case PLANE:
          spa_sub_optimize(model, geno, param, i, LOCT_ONLY);
          break;
        case GLOBE:
          spa_sub_optimize(model, geno, param, i, LOCT_GLOBE);
          break;
      }
    }
  } else if (mode == ADMIXED) {
    initialize_random_location(model, geno, param);
    for(i = 0; i < geno->n_individual; i++) {
      spa_sub_optimize_admixed(model, geno, param, i, 3); // try 3 times
    }
  } else {
    spa_error_exit("You should not see this error message, "
                   "but if you do, please report to wenyun@ucla.edu");
  }
}

void spa_sub_optimize(spa_model *model,
                      const spa_data *geno,
                      const spa_parameter *param,
                      const int i,
                      const int mode) {
  int j;
  int iter;
  double a_grad[MAX_DIMENSION];
  double a_hess[MAX_DIMENSION * MAX_DIMENSION];
  double b_grad;
  double b_hess;
  double gradproj[MAX_DIMENSION];


  double ad[MAX_DIMENSION], bd;
  double atmp[MAX_DIMENSION], htmp[MAX_DIMENSION * MAX_DIMENSION];
  double *pt;
  double pb;

  double f, lambda, t, obj, objt, bound;

  if(mode == COEF_ONLY) {
    // optimize a[i] b[i] 
    pt = model->coef_a[i];
    pb = model->coef_b[i];
    for(iter = 0; iter < param->max_sub_iter; iter++) {  
      vector_init(a_grad, param->dimension, 0); 
      vector_init(a_hess, param->dimension*param->dimension, 0); 

      b_grad = 0;
      b_hess = 0;
      
      // compute gradient and Hessian
      for(j = 0; j < geno->n_individual; j++) {
        f = 1 / (1 + exp(- vector_inner_product(model->coef_a[i],
                                                model->x[j],
                                                param->dimension) 
                         - model->coef_b[i]));

        vector_out_product(htmp, model->x[j], param->dimension);

        switch(get_genotype(geno->genotype, j, i)) {
          case HOMO_MAJOR:
            vector_add(a_grad, model->x[j], 2*f, param->dimension);
            b_grad += 2*f; 
            break;
          case HETER:
            vector_add(a_grad, model->x[j], -1+2*f, param->dimension);
            b_grad += -1+2*f;
            break;
          case HOMO_MINOR:
            vector_add(a_grad, model->x[j], -2+2*f, param->dimension);
            b_grad += -2+2*f;
            break;
          default:
            break;
        }

        vector_add(a_hess, htmp, 2*f*(1-f), param->dimension*param->dimension);
        b_hess += 2*f*(1-f);
      }
      
      // test termination
      vector_copy(ad, a_grad, param->dimension);
      lusolv(a_hess, param->dimension, ad, param);
      vector_scale(ad, -1, param->dimension);

      if(b_hess > 0) {
        bd = - b_grad / b_hess;
      } else {
        bd = - b_grad;
      }
      lambda = (- vector_inner_product(a_grad, ad, param->dimension) 
                - bd * b_grad) 
               / geno->n_individual;
      
      if(lambda < param->epsilon)
        break;

      // line search
      obj = spa_sub_objective(model, geno, param, i, COEF_ONLY);
      
      if(iter > param->large_step_since) {
        t = 1000;
      } else {
        t = 1;
      }

      model->coef_a[i] = atmp;
      vector_add_to_new(atmp, pt, ad, t, param->dimension);
      model->coef_b[i] = pb + bd * t;
      
      objt = spa_sub_objective(model, geno, param, i, COEF_ONLY);
      bound = obj + param->alpha * t *
                    (vector_inner_product(a_grad, ad, param->dimension) + 
                     bd * b_grad);

      while(objt > bound) {
        t = t * param->beta;
        vector_add_to_new(atmp, pt, ad, t, param->dimension);
        model->coef_b[i] = pb + bd * t;
        objt = spa_sub_objective(model, geno, param, i, COEF_ONLY);
        bound = obj + param->alpha * t * 
                      (vector_inner_product(a_grad, ad, param->dimension) + 
                       bd * b_grad);
      }

      vector_copy(pt, atmp, param->dimension);
      model->coef_a[i] = pt;
      pb = model->coef_b[i];

      sprintf(line, 
              "Iter %d: objective = %.10f, gradient norm = %.10f",
              iter,
              obj,
              lambda);
      spa_message(line, WORDY, param);
    }
  } else if (mode == LOCT_ONLY) {
    pt = model->x[i];
  
    for(iter = 0; iter < param->max_sub_iter; iter++) {
      vector_init(a_grad, param->dimension, 0); 
      vector_init(a_hess, param->dimension*param->dimension, 0); 
      
      // compute a_gradient and hessian
      for(j = 0; j < geno->n_snp; j++) {
        f = 1 / (1 + exp(- vector_inner_product(model->coef_a[j],
                                                model->x[i], 
                                                param->dimension)
                         - model->coef_b[j]));
         
        vector_out_product(htmp, model->coef_a[j], param->dimension);

        switch(get_genotype(geno->genotype, i, j)) {
          case HOMO_MAJOR:
            vector_add(a_grad, model->coef_a[j], 2*f, param->dimension);
            break;
          case HETER:
            vector_add(a_grad, model->coef_a[j], -1+2*f, param->dimension);
            break;
          case HOMO_MINOR:
            vector_add(a_grad, model->coef_a[j], -2+2*f, param->dimension);
            break;
        }

        vector_add(a_hess, htmp, 2*f*(1-f), param->dimension*param->dimension);
      }

      // test termination
      vector_copy(ad, a_grad, param->dimension);
      lusolv(a_hess, param->dimension, ad, param);
      vector_scale(ad, -1, param->dimension);

      lambda = sqrt(- vector_inner_product(a_grad, ad, param->dimension))
                 / sqrt((double)geno->n_snp);
      if(lambda < param->epsilon) {
        break;
      }

      // line search
      obj = spa_sub_objective(model, geno, param, i, LOCT_ONLY);
      
      model->x[i] = atmp;
      
      t = 1;
      vector_add_to_new(atmp, pt, ad, t, param->dimension);
      
      objt = spa_sub_objective(model, geno, param, i, LOCT_ONLY);
      bound = obj + param->alpha * t *
                    vector_inner_product(a_grad, ad, param->dimension);

      while(objt > bound) {
        t = t * param->beta;
        vector_add_to_new(atmp, pt, ad, t, param->dimension);

        objt = spa_sub_objective(model, geno, param, i, LOCT_ONLY);
        bound = obj + param->alpha * t *
                      vector_inner_product(a_grad, ad, param->dimension);
      }

      vector_copy(pt, atmp, param->dimension);
      model->x[i] = pt;
  
      sprintf(line,
              "Iter %d: objective = %.10f, gradient norm = %.10f",
              iter,
              obj,
              lambda);
      spa_message(line, WORDY, param);
    }
  } else if(mode == LOCT_GLOBE) {
    pt = model->x[i];
  
    for(iter = 0; iter < param->max_sub_iter; iter++) {
      vector_init(a_grad, param->dimension, 0); 
      vector_init(a_hess, param->dimension*param->dimension, 0); 

      // compute a_gradient and Hessian
      for(j = 0; j < geno->n_snp; j++) {
        f = 1 / (1 + exp(- vector_inner_product(model->coef_a[j],
                                                model->x[i],
                                                param->dimension)
                         - model->coef_b[j]));

        switch(get_genotype(geno->genotype, i, j)) {
          case HOMO_MAJOR:
            vector_add(a_grad, model->coef_a[j], 2*f, param->dimension);
            break;
          case HETER:
            vector_add(a_grad, model->coef_a[j], -1+2*f, param->dimension);
            break;
          case HOMO_MINOR:
            vector_add(a_grad, model->coef_a[j], -2+2*f, param->dimension);
            break;
        }
      }
      
      vector_copy(ad, a_grad, param->dimension);
      vector_scale(ad, -1, param->dimension);

      // line search
      obj = spa_sub_objective(model, geno, param, i, LOCT_ONLY);

      model->x[i] = atmp;

      t = 2 / sqrt(vector_inner_product(ad, ad, param->dimension));
      vector_add_to_new(atmp, pt, ad, t, param->dimension);
      vector_normalize(atmp, param->dimension);
      vector_add_to_new(gradproj, pt, atmp, -1, param->dimension);
      vector_scale(gradproj, 1/t, param->dimension);

      objt = spa_sub_objective(model, geno, param, i, LOCT_ONLY);
      bound = obj - 
              t * vector_inner_product(a_grad,
                                       gradproj,
                                       param->dimension) +
              t / 2 * vector_inner_product(gradproj, 
                                           gradproj,
                                           param->dimension);

      while(objt > bound) {
        t = t * param->beta;
        vector_add_to_new(atmp, pt, ad, t, param->dimension);
        vector_normalize(atmp, param->dimension);
        vector_add_to_new(gradproj, pt, atmp, -1, param->dimension);
        vector_scale(gradproj, 1/t, param->dimension);

        objt = spa_sub_objective(model, geno, param, i, LOCT_ONLY);
        bound = obj - 
                t * vector_inner_product(a_grad,
                                         gradproj,
                                         param->dimension) + 
                t / 2 * vector_inner_product(gradproj,
                                             gradproj,
                                             param->dimension);
      }
      
      // test for termination
      vector_add(pt, atmp, -1, param->dimension);
      lambda = sqrt(vector_inner_product(pt, pt, param->dimension));

      vector_copy(pt, atmp, param->dimension);
      model->x[i] = pt;
        
      sprintf(line,
              "Iter %d: objective = %.10f, step length = %.10f, t = %.10f",
              iter,
              obj,
              lambda,
              t);
      spa_message(line, WORDY, param);

      if(lambda < 1e-3) {
        break;
      }
    }
  } else {
    spa_error_exit("program bug (id 101) detected !!!, "
                   "please report to wenyun@ucla.edu");
  }
}

double spa_objective(const spa_model *model,
                     const spa_data *geno,
                     const spa_parameter *param) {
  int i;
  double objective = 0.0;

  for(i = 0; i < geno->n_individual; i++) {
    objective += spa_sub_objective(model, geno, param, i, LOCT_ONLY);  
  }

  return objective;
}

double spa_sub_objective(const spa_model *model,
                         const spa_data *geno,
                         const spa_parameter *param,
                         int i,
                         int mode) {
  int j;
  double objective = 0.0;
  double f, f1, f2;

  if(mode == COEF_ONLY) { 
    // for a[i] 
    for(j = 0; j < geno->n_individual; j++) {
      f = vector_inner_product(model->coef_a[i],
                               model->x[j],
                               param->dimension) + 
          model->coef_b[i];

      f1 = log(1 + exp(f));
      f2 = log(1 + exp(-f));
      
      switch(get_genotype(geno->genotype, j, i)) {
        case HOMO_MAJOR:
          objective += 2 * f1;
          break;
        case HETER:
          objective += f1 + f2;
          break;
        case HOMO_MINOR:
          objective += 2 * f2;
          break;
        default:
          break;
      }
    }
  } else if(mode == LOCT_ONLY) {
    // for x[i] 
    for(j = 0; j < geno->n_snp; j++) {
      f = vector_inner_product(model->coef_a[j],
                               model->x[i],
                               param->dimension) + 
          model->coef_b[j];
     
      f1 = log(1 + exp(f));
      f2 = log(1 + exp(-f));

      switch(get_genotype(geno->genotype, i, j)) {
        case HOMO_MAJOR:
          objective += 2 * f1;
          break;
        case HETER:
          objective += f1 + f2;
          break;
        case HOMO_MINOR:
          objective += 2 * f2;
          break;
        default:
          break;
      }
    }
  }

  return objective;
}

double spa_sub_objective_admixed(const spa_model *model,
                                 const spa_data *geno,
                                 const spa_parameter *param,
                                 const int i)
{
  int j;
  double objective = 0.0;
  double f, m;
  
  for(j = 0; j < geno->n_snp; j++) {
    f = 1 / (1 + exp( - vector_inner_product(model->coef_a[j],
                                             model->x[i],
                                             param->dimension) 
                      - model->coef_b[j]));

    m = 1 / (1 + exp( - vector_inner_product(model->coef_a[j],
                                             model->x[i] + param->dimension, 
                                             param->dimension) 
                      - model->coef_b[j]));
    
    switch(get_genotype(geno->genotype, i, j)) {
      case HOMO_MAJOR:
        objective += log((1-f)*(1-m));
        break;
      case HETER:
        objective += log(f*(1-m)+m*(1-f));
        break;
      case HOMO_MINOR:
        objective += log(f*m);
        break;
    }
  }

  return (-objective);
}


void spa_sub_optimize_admixed(spa_model *model,
                              const spa_data *geno,
                              const spa_parameter *param,
                              const int i,
                              const int n_trial) {
  int j, trial, iter;
  double grad[MAX_DIMENSION * PARENT];

  double xd[MAX_DIMENSION * PARENT];
  double xtmp[MAX_DIMENSION * PARENT];
  double min_x[MAX_DIMENSION * PARENT];
  double *pt;

  double f, m, lambda, t, obj, objt, bound, min_obj;
  
  pt = model->x[i];
  min_obj = INF;
  for(trial = 0; trial < n_trial; trial++) {
    for(iter = 0; iter < param->max_sub_iter; iter++) {  
      vector_init(grad, param->dimension * param->generation, 0); 

      // compute gradient
      for(j = 0; j < geno->n_snp; j++) {
        f = 1 / (1 + exp(- vector_inner_product(model->coef_a[j], 
                                                model->x[i],
                                                param->dimension) 
                         - model->coef_b[j]));

        m = 1 / (1 + exp(- vector_inner_product(model->coef_a[j], 
                                                model->x[i] + param->dimension,
                                                param->dimension) 
                         - model->coef_b[j]));
        
        switch(get_genotype(geno->genotype, i, j)) {
          case HOMO_MAJOR:
            vector_add(grad, model->coef_a[j], f, param->dimension);
            vector_add(grad + param->dimension, 
                       model->coef_a[j], 
                       m, 
                       param->dimension);
            break;
          case HETER:
            vector_add(grad, 
                       model->coef_a[j],
                       -(1-2*m)*(1-f)*f/(f*(1-m)+m*(1-f)),
                       param->dimension);
            vector_add(grad + param->dimension,
                       model->coef_a[j],
                       -(1-2*f)*(1-m)*m/(f*(1-m)+m*(1-f)),
                       param->dimension);
            break;
          case HOMO_MINOR:
            vector_add(grad, model->coef_a[j], -1+f, param->dimension);
            vector_add(grad + param->dimension,
                       model->coef_a[j],
                       -1+m,
                       param->dimension);
            break;
        }
      }

      // test termination
      vector_copy(xd, grad, param->dimension * param->generation);
      vector_scale(xd, -1, param->dimension * param->generation);

      lambda = 
        sqrt(-vector_inner_product(grad, 
                                   xd, 
                                   param->dimension * param->generation));
      
      if(lambda < param->epsilon)
        break;

      // line search
      obj = spa_sub_objective_admixed(model, geno, param, i);
      model->x[i] = xtmp;
      
      t = 1;
      vector_add_to_new(xtmp, pt, xd, t, param->dimension * param->generation);
      
      objt = spa_sub_objective_admixed(model, geno, param, i);
      bound = obj + param->alpha * 
                    t * 
                    vector_inner_product(grad, 
                                         xd,
                                         param->dimension * param->generation);
      
      while(objt > bound) {
        t = t * param->beta;
        vector_add_to_new(xtmp, pt, xd, t, param->dimension * param->generation);

        objt = spa_sub_objective_admixed(model, geno, param, i);
        bound = obj + param->alpha * 
                      t * 
                      vector_inner_product(grad,
                                           xd,
                                           param->dimension * param->generation);
      }
      vector_copy(pt, xtmp, param->dimension * param->generation);
      model->x[i] = pt;
      
      sprintf(line,
              "Iter %d: objective = %.10f, gradient norm = %.10f",
              iter,
              obj,
              lambda);

      spa_message(line, WORDY, param);
    }
    if (obj < min_obj) {
      min_obj = obj;
      vector_copy(min_x, model->x[i], param->dimension * param->generation);
    }
  }
  vector_copy(model->x[i], min_x, param->dimension * param->generation);
}

void lusolv(double *a, int n, double *b, const spa_parameter *param)
{
  int d;
  int flag;

  flag = ludcmp(a, n, indx, &d);

  if(flag) {
    // if singular, just use gradient
    lubksb(a, n, indx, b);
  } else {
    spa_message("singular hessian! Switch to gradient descent", WORDY, param);
  }
}

void spa_selection(const spa_model *model,
                   const spa_parameter *param) {
  int i, j;
  double *tmp;

  tmp = Malloc(double, model->n_individual);

  for(i = 0; i < model->n_snp; i++) {
    for(j = 0; j < model->n_individual; j++) {
      tmp[j] = 
        1 / (1 + exp(- vector_inner_product(model->coef_a[i], 
                                            model->x[j],
                                            param->dimension) 
                     - model->coef_b[i]));
    }
    model->score[i] = vector_std(tmp, model->n_individual);
  }

  free(tmp);
}

void initialize_random_location(spa_model *model,
                                const spa_data *geno,
                                const spa_parameter *param) {
  int i, j;
  srand(time(NULL) );

  for(i = 0; i < geno->n_individual * 
                 param->dimension * 
                 param->generation; i++) { 
    // random initialization centered at origin
    model->x_space[i] = rand() / (double)RAND_MAX - 0.5;
  }

  if(param->dimension == GLOBE) {
    // normalize to ||x|| = 1
    for(i = 0; i < geno->n_individual; i++) {
      for(j = 0; j < param->generation; j++) {
        vector_normalize(model->x[i] + (param->dimension * j), 
                         param->dimension);
      }
    }
  }
}

void parse_input_parameters(int argc,
                            char** argv,
                            spa_parameter *param) {
  int i;

  for(i = 1; i < argc; i++) {
    if(argv[i][0] != '-') {
      exit_with_help();
    }
    if(++i >= argc) {
      exit_with_help();
    }

    switch(argv[i-1][1]) {
      case '-':
        if(!strcmp(argv[i-1], "--bfile")) {
          param->bfile = argv[i];
        } else if(!strcmp(argv[i-1], "--tfile")) {
          param->tfile = argv[i];
        } else if(!strcmp(argv[i-1], "--pfile")) {
          param->pfile = argv[i];
        } else if(!strcmp(argv[i-1], "--gfile")) {
          param->gfile = argv[i];
        } else if(!strcmp(argv[i-1], "--mfile")) {
          param->mfile = argv[i];
        } else if(!strcmp(argv[i-1], "--location-input")) {
          param->ilfile = argv[i];
        } else if(!strcmp(argv[i-1], "--model-input")) {
          param->imfile = argv[i];
        } else if(!strcmp(argv[i-1], "--location-output")) {
          param->olfile = argv[i];
        } else if(!strcmp(argv[i-1], "--model-output")) {
          param->omfile = argv[i];
        } else {
          sprintf(line, "Unknown options: %s\n", argv[i-1]);
          spa_warning(line);
          exit_with_help();
        }
        break;
      case 'n':
        param->generation = atoi(argv[i]);
        break;
      case 'k':
        param->dimension = atoi(argv[i]);
        break;
      case 'e':
        param->step_epsilon = atof(argv[i]);
        break;
      case 'v':
        param->verbose = atoi(argv[i]);
        break;
      case 'r':
        param->epsilon = atof(argv[i]);
        break;
      default:
        sprintf(line, "Unknown options: -%s\n", argv[i-1]);
        spa_warning(line);
        exit_with_help();
      
    }
  }
}

int parameter_sanity_check(const spa_parameter *param) {
  bool flag = true;

  if ((param->gfile != NULL) +
      (param->pfile != NULL) + 
      (param->tfile != NULL) + 
      (param->bfile != NULL) > 1) {
    spa_warning("More than one genotype/plink files, "
                "make sure they are consistent\n");
  }
  if(param->ilfile && param->imfile) {
    spa_warning("Both location and slope function are given, "
                "nothing to learn...\n");
    flag = false;
  }
  if(param->ilfile && param->olfile) {
    spa_warning("Location file is given, "
                "the output location file would be the same...\n");
    flag = false;
  }
  if(param->imfile && param->omfile) {
    spa_warning("Model file is given, "
                "the output model file would be the same...\n");
    flag = false;
  }
  if(param->generation > MAX_GENERATION) {
    spa_warning("Ancestral inference for more than 2 ancestral locations: "
                "not implemented yet\n");
    flag = false;
  }
  if(param->dimension > MAX_DIMENSION) {
    spa_warning("Dimension larger than 3 is not feasible"
                "geographical coordinates\n");
    flag = false;
  }
  // TODO: remove once method gets extended
  if(param->generation == PARENT && param->dimension == GLOBE) {
    spa_warning("Ancestral inference for admixed individual "
                "in three dimension: not implemented yet\n");
    flag = false;
  }
  if(!param->olfile && !param->omfile) {
    spa_warning("Please specify the output files...\n");
    flag = false;
  }
  if(param->epsilon > 1e-1) {
    spa_warning("The optimization tolerance parameter might be too large\n");
  }

  return flag;
}

void data_sanity_check(const spa_data *geno) {
  int i, j;
  int c0, c1, c2, cm;

  for(i = 0; i < geno->n_snp; i++) {
    c0 = 0;
    c1 = 0;
    c2 = 0;
    cm = 0;
    for(j = 0; j < geno->n_individual; j++) {
      switch(get_genotype(geno->genotype, j, i)) {
        case HOMO_MAJOR:
          c0 ++;
          break;
        case HETER:
          c1 ++;
          break;
        case HOMO_MINOR:
          c2 ++;
          break;
        case MISSING:
          cm ++;
          break;
        default:
          break;
      }
    }

    if(c0 + c1 == 0 || c1 + c2 == 0) {
      sprintf(line, 
              "%d-th SNP is not biallelic. "
              "Please use PLINK software to remove those "
              "non-biallelic SNPs\nHint: use --maf option in PLINK",
              i);
      spa_error_exit(line);
    }
  }
}

void model_sanity_check(const spa_model* model) {

  int i;
  for (i = 0; i < model->n_snp; i++) {
    if (model->snp_info[i].snp_major == MISSING_ALLELE ||
        model->snp_info[i].snp_minor == MISSING_ALLELE) {
      sprintf(line, 
              "Major/Minor SNPs in model can not be missing."
              "Missing is found  in %d-th SNP",
              i);
      spa_error_exit(line);
    }
    if (model->snp_info[i].snp_major == model->snp_info[i].snp_minor) {
      sprintf(line,
              "Major and minor SNPs can not be the same in model."
              "In %d-th SNP, the same SNPs are found",
              i);
      spa_error_exit(line);
    }
  }
}


void model_data_consistent(spa_model *model,
                           spa_data *geno,
                           const spa_parameter *param,
                           const int mode) {
  int i, j;

  map<string, int> order;
  bool* is_valid;

  if (mode == LOCT_ONLY) {

    is_valid = Malloc(bool, model->n_individual);
    check_individual_valid(geno, model, order, is_valid);

    // remove individuals in model but not in genotype data
    j = 0;
    for (i = 0; i < model->n_individual; i++) {
      if (is_valid[i]) {
        model->individual_info[j] = model->individual_info[i];
        model->x[j] = model->x[i];
        j++;
      } else {
        sprintf(line, "%d-th individual (%s) in model "
                      "but not in genotype data, thus removed",
                      i,
                      model->individual_info[i].individual_id);
        spa_message(line, CRAZY, param);
        free_individual_info(model->individual_info + i);
      }
    }
    model->n_individual = order.size();
    copy_info(model, geno, COEF_ONLY);
    allocate_model_coef(model, param);
    reorder_genotype(geno, order, mode, param);
    printf("%d individuals are in common between model and genotype file\n",
           model->n_individual);

  } else if (mode == COEF_ONLY) {
    
    is_valid = Malloc(bool, model->n_snp); 
    check_snp_valid(geno, model, order, is_valid);
    
    // remove SNPs in model but not in genotype data
    j = 0;
    for (i = 0; i < model->n_snp; i++) {
      if (is_valid[i]) {
        model->snp_info[j] = model->snp_info[i];
        model->coef_a[j] = model->coef_a[i];
        model->coef_b[j] = model->coef_b[i];
        model->score[j] = model->score[i];
        j++;
      } else {
        sprintf(line, "%d-th snp (%s) in model but not in genotype data, "
                      "thus removed", i, model->snp_info[i].snp_id);
        spa_message(line, CRAZY, param);
        free_snp_info(model->snp_info + i); 
      }
    }
  
    model->n_snp = order.size();
    copy_info(model, geno, LOCT_ONLY);
    allocate_model_x(model, param);
    reorder_genotype(geno, order, mode, param);

    // check major/minor allele consistency. 
    // Assumption: model should have two non-missing alleles 
    map<int, int> revise_map;
    for (i = 0; i < geno->n_snp; i++) {
      revise_map.clear();
      check_major_minor_allele_consistency(geno->snp_info + i,
                                           model->snp_info + i,
                                           revise_map);
      revise_genotype(geno, i, revise_map);
      geno->snp_info[i].snp_major = model->snp_info[i].snp_major;
      geno->snp_info[i].snp_minor = model->snp_info[i].snp_minor;
    }
    printf("%d SNP are in common between model and genotype file\n",
           model->n_snp);
  } else {
    spa_error_exit("You should not see this error,"
                   "if you do, please let me know (wenyun@ucla.edu)\n");
  }

  free(is_valid);
}

void check_major_minor_allele_consistency(
    const snp_info_struct* geno_snp_info,
    const snp_info_struct* model_snp_info,
    map<int, int>& revise_map) {
  
  if (geno_snp_info->snp_major == model_snp_info->snp_minor) {
    revise_map[HOMO_MAJOR] = HOMO_MINOR;      
  }

  if (geno_snp_info->snp_minor == model_snp_info->snp_major) {
    revise_map[HOMO_MINOR] = HOMO_MAJOR;      
  }
}


void check_individual_valid(const spa_data* geno,
                            const spa_model* model,
                            map<string, int>& order,
                            bool* is_valid) {
  int i, j;
  map<string, int> geno_key;

  for (i = 0; i < geno->n_individual; i++) {
    sprintf(line,
            "%s_%s_%s_%s_%s",
            geno->individual_info[i].family_id,
            geno->individual_info[i].individual_id,
            geno->individual_info[i].paternal_id,
            geno->individual_info[i].maternal_id,
            geno->individual_info[i].sex);
    string key(line);
    geno_key[key] = i;
  }
  
  // remove individuals in model but not in genotype data
  j = 0;
  for (i = 0; i < model->n_individual; i++) {
    sprintf(line,
            "%s_%s_%s_%s_%s",
            model->individual_info[i].family_id,
            model->individual_info[i].individual_id,
            model->individual_info[i].paternal_id,
            model->individual_info[i].maternal_id,
            model->individual_info[i].sex);
    string key(line);
    if (geno_key.count(key) > 0) {
      is_valid[i] = true;
      order[key] = j++;
    } else {
      is_valid[i] = false;
    }
  }
}

void check_snp_valid(const spa_data* geno,
                     const spa_model* model,
                     map<string, int>& order,
                     bool* is_valid) {
  
  int i, j;
  map<string, int> geno_key;

  for (i = 0; i < geno->n_snp; i++) {
    string key(geno->snp_info[i].snp_id);
    geno_key[key] = i;
  }

  j = 0;
  for (i = 0; i < model->n_snp; i++) {
    snp_info_struct& snp_info = model->snp_info[i];
    is_valid[i] = true;  // default
    
    string key(snp_info.snp_id);
    
    // not exist in genotyep data, remove
    if (geno_key.count(key) == 0) {
      is_valid[i] = false;
      continue;
    }

    // might have strand problem, remove
    if ((snp_info.snp_major == 'G' && snp_info.snp_minor == 'C') ||
        (snp_info.snp_major == 'C' && snp_info.snp_minor == 'G') ||
        (snp_info.snp_major == 'A' && snp_info.snp_minor == 'T') ||
        (snp_info.snp_major == 'T' && snp_info.snp_minor == 'A')) {
      is_valid[i] = false;
      continue;
    }

    // inconsistent alleles with genotyep data, remove
    snp_info_struct& g_snp_info = geno->snp_info[geno_key[key]];
    if ((g_snp_info.snp_major != '0' && 
         g_snp_info.snp_major != snp_info.snp_major &&
         g_snp_info.snp_major != snp_info.snp_minor) ||
        (g_snp_info.snp_minor != '0' && 
         g_snp_info.snp_minor != snp_info.snp_major &&
         g_snp_info.snp_minor != snp_info.snp_minor)) {
      is_valid[i] = false;
    }
  
    if (is_valid[i]) {
      order[key] = j++;
    }
  }
}

void revise_genotype(spa_data* geno,
                     const int i,
                     const std::map<int, int>& revise_map) {

  int j;
  int genotype;
  
  if (!revise_map.empty()) {
    for (j = 0; j < geno->n_individual; j++) {
      genotype = get_genotype(geno->genotype, j, i);
      map<int, int>::const_iterator iter = revise_map.find(genotype);
      if (iter != revise_map.end()) {
        set_genotype(geno->genotype, j, i, (char) iter->second);
      }
    }
  }
}

void reorder_genotype(spa_data* geno,
                  const map<string, int>& order,
                  const int mode,
                  const spa_parameter* param) {
  int i, j, k;

  if (mode == LOCT_ONLY) {

    // re-allocate those to replace previous
    int n_individual = (int) order.size();
    individual_info_struct* individual_info =
      Malloc(individual_info_struct, n_individual);
    char** genotype = Malloc(char*, n_individual);

    // reorder genotype data by individual order
    for (i = 0; i < geno->n_individual; i++) {
      sprintf(line,
              "%s_%s_%s_%s_%s",
              geno->individual_info[i].family_id,
              geno->individual_info[i].individual_id,
              geno->individual_info[i].paternal_id,
              geno->individual_info[i].maternal_id,
              geno->individual_info[i].sex);
      string key(line);
      if (order.count(key) > 0) {
        map<string, int>::const_iterator it = order.find(key);
        j = it->second;
        individual_info[j] = geno->individual_info[i];
        genotype[j] = geno->genotype[i];
      } else {
        sprintf(line, "%d-th individual (%s) in genotype data "
                      "but not in model, thus removed",
                      i,
                      geno->individual_info[i].individual_id);
        spa_message(line, CRAZY, param);
        free_individual_info(geno->individual_info + i);
      }
    }
    geno->n_individual = n_individual;
    
    free(geno->individual_info);
    geno->individual_info = individual_info;
    geno->genotype = genotype;

  } else if (mode == COEF_ONLY) {
    
    // keep the old values
    int n_snp = geno->n_snp;
    snp_info_struct* snp_info = geno->snp_info;
    char** genotype = geno->genotype;
    char* genotype_space = geno->genotype_space;

    // reallocate space
    geno->n_snp = (int) order.size();
    geno->snp_info = Malloc(snp_info_struct, geno->n_snp);
    allocate_genotype(geno);
    
    // copy values
    for (i = 0; i < n_snp; i++) {
      string key(snp_info[i].snp_id);
      map<string, int>::const_iterator it = order.find(key);      
      if (it != order.end()) {
        j = it->second;
        geno->snp_info[j] = snp_info[i];
        for (k = 0; k < geno->n_individual; k++) {
          set_genotype(geno->genotype,
                       k,
                       j,
                       get_genotype(genotype, k, i)); 
        }
      } else {
        sprintf(line, "%d-th snp (%s) in genotype data but not in model, "
                      "thus removed", i, snp_info[i].snp_id);
        spa_message(line, CRAZY, param);
        free_snp_info(snp_info + i);
      }
    }

    free(snp_info);
    free(genotype_space);
    free(genotype);
  }
}

void copy_snp_info(snp_info_struct* dest,
                   const snp_info_struct* source) {
  dest->chromosome = Malloc(char, strlen(source->chromosome) + 1);
  dest->snp_id = Malloc(char, strlen(source->snp_id) + 1);
  dest->morgan = Malloc(char, strlen(source->morgan) + 1);

  strcpy(dest->chromosome, source->chromosome);
  strcpy(dest->snp_id, source->snp_id);
  strcpy(dest->morgan, source->morgan);

  dest->position = source->position;
  dest->snp_major = source->snp_major;
  dest->snp_minor = source->snp_minor;
}

void free_snp_info(snp_info_struct* snp_info) {
  free(snp_info->chromosome);
  free(snp_info->snp_id);
  free(snp_info->morgan);
}

void copy_individual_info(individual_info_struct* dest,
                          const individual_info_struct* source) {
  dest->family_id = Malloc(char, strlen(source->family_id) + 1);
  dest->individual_id = Malloc(char, strlen(source->individual_id) + 1);
  dest->paternal_id = Malloc(char, strlen(source->paternal_id) + 1);
  dest->maternal_id = Malloc(char, strlen(source->maternal_id) + 1);
  dest->sex = Malloc(char, strlen(source->sex) + 1);
  dest->phenotype = Malloc(char, strlen(source->phenotype) + 1);

  strcpy(dest->family_id, source->family_id);
  strcpy(dest->individual_id, source->individual_id);
  strcpy(dest->paternal_id, source->paternal_id);
  strcpy(dest->maternal_id, source->maternal_id);
  strcpy(dest->sex, source->sex);
  strcpy(dest->phenotype, source->phenotype);
}

void free_individual_info(individual_info_struct* individual_info) {
  free(individual_info->family_id);
  free(individual_info->individual_id);
  free(individual_info->paternal_id);
  free(individual_info->maternal_id);
  free(individual_info->sex);
  free(individual_info->phenotype);
}

void copy_info(spa_model* model, const spa_data* geno, const int mode) {

  int i;
  if (mode == COEF_ONLY || mode == BOTH) {
    model->n_snp = geno->n_snp;
    model->snp_info = Malloc(snp_info_struct, model->n_snp);
    for (i = 0; i < model->n_snp; i++) {
      copy_snp_info(model->snp_info + i, geno->snp_info + i);
    }
  }
  if (mode == LOCT_ONLY || mode == BOTH) {
    model->n_individual = geno->n_individual;
    model->individual_info =
      Malloc(individual_info_struct, model->n_individual);
    for (i = 0; i < model->n_individual; i++) {
      copy_individual_info(model->individual_info + i,
                           geno->individual_info + i);
    }
  }
}

void allocate_model_x(spa_model *model,
                      const spa_parameter* param) {
  int i;
  model->x_space = Malloc(double, model->n_individual *
                                  param->generation *
                                  param->dimension);
  model->x = Malloc(double*, model->n_individual);
  for(i = 0; i < model->n_individual; i++) {
    model->x[i] = &(model->x_space[i *
                                   param->generation *
                                   param->dimension]);
  }
  vector_init(model->x_space, model->n_individual *
                              param->generation *
                              param->dimension, 0); 
}

void allocate_model_coef(spa_model *model,
                         const spa_parameter* param) {
  int i;
  model->coef_a_space = Malloc(double, model->n_snp * param->dimension);
  model->coef_a = Malloc(double*, model->n_snp);
  for(i = 0; i < model->n_snp; i++) {
    model->coef_a[i] = &(model->coef_a_space[i * param->dimension]);
  }
  model->coef_b = Malloc(double, model->n_snp);
  model->score = Malloc(double, model->n_snp);
  vector_init(model->coef_a_space, model->n_snp * param->dimension, 0);  
  vector_init(model->coef_b, model->n_snp, 0);
  
}

void allocate_model_for_bootstrap(spa_model *model,
                                  const spa_parameter *param) {
  allocate_model_x(model, param);
  allocate_model_coef(model, param);
  
  vector_init(model->coef_a_space, model->n_snp * param->dimension, 0);
  vector_init(model->coef_b, model->n_snp, 0);
}

void free_model(const spa_model *model) {
  int i;

  free(model->coef_a_space);
  free(model->x_space);
  free(model->coef_a);
  free(model->coef_b);
  free(model->x);
  free(model->score);

  for(i = 0; i < model->n_individual; i++) {
    free_individual_info(model->individual_info + i);
  }
  free(model->individual_info);
  
  for(i = 0; i < model->n_snp; i++) {
    free_snp_info(model->snp_info + i);
  }
  free(model->snp_info);
}

void free_data(const spa_data *geno) {
  int i;

  free(geno->genotype);
  free(geno->genotype_space);
  
  for(i = 0; i < geno->n_individual; i++) {
    free_individual_info(geno->individual_info + i);
  }
  free(geno->individual_info);
  
  for(i = 0; i < geno->n_snp; i++) {
    free_snp_info(geno->snp_info + i);
  }
  free(geno->snp_info);
}

void spa_error_exit(const char* msg) {
  fprintf(stderr, "FATAL => %s\n", msg);
  exit(1);
}

void spa_warning(const char* msg) {
  fprintf(stderr, "WARINING => %s\n", msg);
}

void spa_message(const char* msg, 
                 const int level,
                 const spa_parameter *param) {
  if(level <= param->verbose) {
    printf("%s\n", msg);
  }
}
