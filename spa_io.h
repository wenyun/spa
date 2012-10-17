// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#pragma once

#include "spa_types.h"

#include <cstdio>

// IO functions
FILE* spa_open_file(const char* filename, const char* mode);
void count_file(FILE* fp, int* rows, int* columns);
void read_bedfile(const char* filename,
                  spa_data* geno,
                  const spa_parameter* param);
void read_pedfile(const char* filename,
                  spa_data* geno,
                  const spa_parameter* param);
void read_tpedfile(const char* filename,
                   spa_data* geno,
                   const spa_parameter* param);
void read_gfile(const char* filename,
                spa_data* geno,
                const spa_parameter* param);
void read_mfile(const char* filename,
                spa_data* geno,
                const spa_parameter* param);
void read_location_ilfile(const char* filename, 
                          spa_model* model,
                          const spa_parameter* param);
void read_model_imfile(const char* filename,
                       spa_model* model,
                       const spa_parameter* param);
void write_location_olfile(const char* filename,
                           const spa_model* model,
                           const spa_parameter* param);
void write_html_location_olfile(const char* filename,
                                const spa_model* model,
                                const spa_parameter* param);
void write_model_omfile(const char* filename,
                        const spa_model* model,
                        const spa_parameter* param);
// get genotype 0/1/2 for individual i and SNP j
char get_genotype(char** genotype, const int i, const int j);
// set genotype g = 0/1/2 for individual i and SNP j
void set_genotype(char** genotype, const int i, const int j, const char g);
// flip snp_major and snp_minor if they are wrongly set
void snp_flip(spa_data* geno);
void arrange_genotype_in_number(char* buf, spa_data* geno, int snp_i);


void allocate_genotype(spa_data* geno);
void read_snp_info(snp_info_struct* snp_info, int mode);
void read_individual_info(individual_info_struct* individual_info);
void free_snp_info(snp_info_struct* snp_info);
void free_individual_info(individual_info_struct* individual_info);

