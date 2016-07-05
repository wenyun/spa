// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#pragma once

#include "spa_types.h"

#include <cstdio>
#include <map>
#include <string>

// initialization functions
void initialize();
void set_default_parameter(spa_parameter* param);
void allocate_model_for_bootstrap(spa_model* model,
                                  const int n_individual,
                                  const int n_snp,
                                  const spa_parameter* param);
void initialize_random_location(spa_model* model,
                                const spa_data* geno,
                                const spa_parameter* param);

// finalization functions
void finalize();
void free_model(const spa_model* model);
void free_data(const spa_data* geno);

// sanity check functions
void parse_input_parameters(int argc, char** argv, spa_parameter* param);
int parameter_sanity_check(const spa_parameter* param);
void model_data_consistent(spa_model* model,
                           spa_data* geno,
                           const spa_parameter* param,
                           const int mode);
void data_sanity_check(const spa_data* geno);
void model_sanity_check(const spa_model* model);

// SPA functions
void spa_selection(const spa_model* model,
                   const spa_parameter* param);
double spa_objective(const spa_model* model,
                     const spa_data* geno,
                     const spa_parameter* param);
void spa_optimize(spa_model* model,
                  const spa_data* geno,
                  const spa_parameter* param,
                  int mode);
double spa_sub_objective(const spa_model* model,
                         const spa_data* geno,
                         const spa_parameter* param,
                         const int index,
                         const int mode);
void spa_sub_optimize(spa_model* model,
                      const spa_data* geno,
                      const spa_parameter* param,
                      const int index,
                      const int mode);
double spa_sub_objective_admixed(const spa_model* model,
                                 const spa_data* geno,
                                 const spa_parameter* param,
                                 const int index);
void spa_sub_optimize_admixed(spa_model* model,
                              const spa_data* geno,
                              const spa_parameter* param,
                              const int index,
                              const int n_trial);

// numerical functions
void lusolv(double* a, int n, double* b, const spa_parameter* param);

// private functions
void copy_info(spa_model* model, const spa_data* geno, const int mode);
void copy_snp_info(snp_info_struct* dest,
                   const snp_info_struct* source);
void copy_individual_info(individual_info_struct* dest,
                          const individual_info_struct* source);
void free_snp_info(snp_info_struct* snp_info);
void free_individual_info(individual_info_struct* individual_info);
void check_individual_valid(const spa_data* geno,
                            const spa_model* model,
                            std::map<std::string, int>& order,
                            bool* is_valid);
void check_snp_valid(const spa_data* geno,
                     const spa_model* model,
                     std::map<std::string, int>& order,
                     bool* is_valid);
void reorder_genotype(spa_data* geno,
                      const std::map<std::string, int>& order,
                      const int mode,
                      const spa_parameter* param);
void check_major_minor_allele_consistency(
    const snp_info_struct* geno_snp_info,
    const snp_info_struct* model_snp_info,
    std::map<int, int>& revise_map); 
void revise_genotype(spa_data* geno,
                     const int snp_i,
                     const std::map<int, int>& revise_map);
void allocate_model_x(spa_model *model,
                      const spa_parameter* param);
void allocate_model_coef(spa_model *model,
                         const spa_parameter* param);
void allocate_model_for_bootstrap(spa_model *model,
                                  const spa_parameter *param);
