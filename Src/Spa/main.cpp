#include "spa.h"
#include "spa_io.h"
#include "spa_util.h"

#include <boost/program_options.hpp>
#include <string>
#include <iostream>

#define MAX_DIMENSION 3
#define MAX_GENERATION 2

namespace po = boost::program_options;
using namespace std;

void print_version_information() {
  printf(
"@----------------------------------------------------------@\n"
"|         SPA!       |      v1.13      |    6/July/2016     |\n"
"|----------------------------------------------------------|\n"
"|  (C) 2012 Wen-Yun Yang, GNU General Public License, v2   |\n"
"|----------------------------------------------------------|\n"
"|  For documentation, citation & bug-report instructions:  |\n"
"|             http://genetics.cs.ucla.edu/spa/             |\n"
"@----------------------------------------------------------@\n"
"\n"
);
}

constexpr auto kPlinkBedFile = "bfile";
constexpr auto kPlinkTpedFile = "tfile";
constexpr auto kPlinkPedFile = "pfile";
constexpr auto kGenotypeFile = "gfile";
constexpr auto k23andmeFile = "mfile";
constexpr auto kLocationInputFile = "location-input";
constexpr auto kModelInputFile = "model-input";
constexpr auto kLocationOutputFile = "location-output";
constexpr auto kModelOutputFile = "model-output";
constexpr auto kGeneration = "n";
constexpr auto kDimension = "k";
constexpr auto kStepEpsilon = "e";
constexpr auto kVerbose = "v";
constexpr auto kTradeoff = "r";

void parse_input_parameters(const po::variables_map &vm,
                            spa_parameter *param) {

  param->bfile = vm.count(kPlinkBedFile) > 0 ? vm[kPlinkBedFile].as<std::string>().c_str() : NULL;
  param->tfile = vm.count(kPlinkTpedFile) > 0 ? vm[kPlinkTpedFile].as<std::string>().c_str() : NULL;
  param->pfile = vm.count(kPlinkPedFile) > 0 ? vm[kPlinkPedFile].as<std::string>().c_str() : NULL;
  param->gfile = vm.count(kGenotypeFile) > 0 ? vm[kGenotypeFile].as<std::string>().c_str() : NULL;
  param->mfile = vm.count(k23andmeFile) > 0 ? vm[k23andmeFile].as<std::string>().c_str() : NULL;
  param->ilfile = vm.count(kLocationInputFile) > 0 ? vm[kLocationInputFile].as<std::string>().c_str() : NULL;
  param->imfile = vm.count(kModelInputFile) > 0 ? vm[kModelInputFile].as<std::string>().c_str() : NULL;
  param->olfile = vm.count(kLocationOutputFile) > 0 ? vm[kLocationOutputFile].as<std::string>().c_str() : NULL;
  param->omfile = vm.count(kModelOutputFile) > 0 ? vm[kModelOutputFile].as<std::string>().c_str() : NULL;
  param->generation = vm[kGeneration].as<int>();
  param->dimension = vm[kDimension].as<int>();
  param->step_epsilon = vm[kStepEpsilon].as<double>();
  param->epsilon = vm[kTradeoff].as<double>();
  param->verbose = vm[kVerbose].as<int>();
  param->max_iter = 1000;
  param->max_sub_iter = 100;
  param->alpha = 0.01;
  param->beta = 0.5;
  param->large_step_since = 15;
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

int main (int argc, char** argv)  {
  
  spa_parameter param;
  spa_model model;
  spa_data geno;

  print_version_information();

  try {
    
    po::options_description desc("Spatial Ancestry Analysis");
    
    desc.add_options()
    (kPlinkBedFile, po::value<std::string>(), "prefix of .bed, .bim and .fam file")
    (kPlinkTpedFile, po::value<std::string>(), "prefix of .tped and .tfam file")
    (kPlinkPedFile, po::value<std::string>(), "prefix of .ped and .map file")
    (kGenotypeFile, po::value<std::string>(), "genotype file")
    (k23andmeFile, po::value<std::string>(), "23andme genotype file")
    (kLocationInputFile, po::value<std::string>(), "known locations")
    (kModelInputFile, po::value<std::string>(), "known slope functions")
    (kLocationOutputFile, po::value<std::string>(), "output file for individual locations")
    (kModelOutputFile, po::value<std::string>(), "output file for slope function coefficients")
    (kGeneration, po::value<int>()->default_value(1), "number of locations")
    (kDimension, po::value<int>()->default_value(2), "dimentions of spatial analysis")
    (kStepEpsilon, po::value<double>()->default_value(0.01), "set tolerance of termination criterion (default 0.01)")
    (kTradeoff, po::value<double>()->default_value(0.000001), "set optimization epsilon tolerance (default 1e-6). Larger value makes the program run faster but poor accuracy")
    (kVerbose, po::value<int>()->default_value(SHORT), "verbose level")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (argc == 1) {
      std::cout << desc << std::endl;
      return 1;
    }
    
    parse_input_parameters(vm, &param);
    if(!parameter_sanity_check(&param)) {
      spa_error_exit("Parameter sanity check failed\n");
    }
    
  } catch (const std::exception &exception) {
    // Catch the initialization error and bail immediately
    std::cerr << "Error: " << exception.what() << std::endl;
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

  std::cout << geno.n_individual << " individual and " << geno.n_snp << " snps are read ...." << std::endl;
  
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
      std::cout << "CHECK YOUR ANCESTRY ORIGIN BY CLICKING " << param.olfile << ".html" << std::endl;
      
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

