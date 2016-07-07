#define BOOST_TEST_MAIN

#include "spa.h"
#include "spa_types.h"
#include "spa_io.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>
#include <boost/filesystem.hpp>
#include <boost/asio.hpp>
#include <iostream>
#include <string>
#include <fstream>

std::string getHttpData(const std::string& server, const std::string& file)
{
  try
  {
    boost::asio::ip::tcp::iostream s(server, "http");
    s.expires_from_now(boost::posix_time::seconds(60));
    
    if (!s){ throw "Unable to connect: " + s.error().message(); }
    
    // ask for the file
    s << "GET " << file << " HTTP/1.0\r\n";
    s << "Host: " << server << "\r\n";
    s << "Accept: */*\r\n";
    s << "Connection: close\r\n\r\n";
    
    // Check that response is OK.
    std::string http_version;
    s >> http_version;
    unsigned int status_code;
    s >> status_code;
    std::string status_message;
    std::getline(s, status_message);
    if (!s && http_version.substr(0, 5) != "HTTP/"){ throw "Invalid response\n"; }
    if (status_code != 200){ throw "Response returned with status code " + status_code; }
    
    // Process the response headers, which are terminated by a blank line.
    std::string header;
    while (std::getline(s, header) && header != "\r"){}
    
    // Write the remaining data to output.
    std::stringstream ss;
    ss << s.rdbuf();
    return ss.str();
  }
  catch(std::exception& e)
  {
    return e.what();
  }
}

void downloadFileWithUrl(const std::string &server, const std::string &path, const std::string &localFile) {
  
  std::string result = getHttpData(server, path);
  
  std::ofstream of(localFile, std::ios::binary);
  of << result;
  of.close();
}

double standard_deviation(double data[], int n)
{
  float mean=0.0, sum_deviation=0.0;
  int i;
  for(i=0; i<n;++i)
  {
    mean+=data[i];
  }
  mean=mean/n;
  for(i=0; i<n;++i)
    sum_deviation+=(data[i]-mean)*(data[i]-mean);
  return sqrt(sum_deviation/n);
}

BOOST_AUTO_TEST_CASE( SpaTest ) {
  
  boost::filesystem::path tmpDir = boost::filesystem::temp_directory_path();
  const auto prefix = tmpDir / "test_data";
  const auto bedFile = prefix.string() + ".bed";
  downloadFileWithUrl("s3-us-west-2.amazonaws.com", "/wenyun-data/spa_test/test_data.bed", bedFile);
  const auto bimFile = prefix.string() + ".bim";
  downloadFileWithUrl("s3-us-west-2.amazonaws.com", "/wenyun-data/spa_test/test_data.bim", bimFile);
  const auto famFile = prefix.string() + ".fam";
  downloadFileWithUrl("s3-us-west-2.amazonaws.com", "/wenyun-data/spa_test/test_data.fam", famFile);
  
  spa_parameter param;
  spa_model model;
  spa_data geno;

  param.generation = SELF;
  param.dimension = 2;
  param.step_epsilon = 0.01;
  param.epsilon = 0.000001;
  param.verbose = SHORT;
  param.max_iter = 1000;
  param.max_sub_iter = 100;
  param.alpha = 0.01;
  param.beta = 0.5;
  param.large_step_since = 15;
  
  initialize();
  read_bedfile(prefix.c_str(), &geno, &param);
  
  std::cout << geno.n_individual << " individuals and " << geno.n_snp << " snps are read ..." << std::endl;
  
  data_sanity_check(&geno);
  copy_info(&model, &geno, BOTH);
  
  std::cout << "Unsupervised learning ..." << std::endl;
  allocate_model_for_bootstrap(&model, &param);
  spa_optimize(&model, &geno, &param, BOTH);
  
  double standardDeviation = standard_deviation(model.coef_a_space, param.dimension * geno.n_snp);
  std::cout << "standard deviation = " << standardDeviation << std::endl;
  BOOST_CHECK(standardDeviation > 0.1);
  
  double objective = spa_objective(&model, &geno, &param);
  BOOST_CHECK(objective < 1937800);
  
  
  // test prediction
  
  const auto modelFile = prefix.string() + ".model";
  spa_selection(&model, &param);
  write_model_omfile(modelFile.c_str(), &model, &param);
  
  read_model_imfile(modelFile.c_str(), &model, &param);
  model_sanity_check(&model);
  model_data_consistent(&model, &geno, &param, COEF_ONLY);
  spa_optimize(&model, &geno, &param, LOCT_ONLY);
  
  // test prediction of admixed individual
  
  param.generation = PARENT;
  read_model_imfile(modelFile.c_str(), &model, &param);
  model_sanity_check(&model);
  model_data_consistent(&model, &geno, &param, COEF_ONLY);
  spa_optimize(&model, &geno, &param, ADMIXED);
  
  free_model(&model);
  free_data(&geno);
  finalize();
}

