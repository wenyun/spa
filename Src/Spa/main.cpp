#include "spa.h"
#include "spa_io.h"
#include "spa_util.h"

#include <cstdio>

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

