// Copyright 2012-present Wen-Yun Yang. All rights reserved.

#include "spa_io.h"
#include "spa_types.h"

#include <cmath>
#include <cstring>
#include <map>

#define LOC_FILE_HEADER 6
#define MODEL_FILE_HEADER 4
#define MODEL_FILE_TAILER 4 
#define BIM 6
#define SPA_MODEL 66
#define MAP 4
#define TTANDME 8
#define GENOTYPE_PER_BYTE 4
#define MISSING_ALLELE '0'


extern char* line;
extern int max_line_len;

char read_mask_genotype(char byte, char c);
void write_mask_genotype(char* byte, char c, char g);
unsigned char mask[4] = {0x03, /* 00000011 */
                         0x0C, /* 00001100 */
                         0x30, /* 00110000 */
                         0xC0  /* 11000000 */};
const int max_name_len = 256;

static char* readline(FILE* input) {
  int len;
  
  if (fgets(line, max_line_len, input) == NULL) {
    return NULL;
  }

  while(strrchr(line, '\n') == NULL) {
    max_line_len *= 2;
    line = (char*) realloc(line, max_line_len);
    len = (int) strlen(line);
    if(fgets(line + len, max_line_len - len, input) == NULL)
      break;
  }
  return line;
}

FILE* spa_open_file(const char* filename, const char* mode) {
  FILE* fp = fopen(filename, mode);
  if(fp == NULL) {
    sprintf(line, "can't open input file %s\n", filename);
    spa_error_exit(line);
  }
  return fp;
}

void allocate_genotype(spa_data* geno) {
  int i;
  int len = (int) ceil(((double)geno->n_snp) / GENOTYPE_PER_BYTE);
  geno->genotype = Malloc(char*, geno->n_individual);
  geno->genotype_space = Malloc(char, geno->n_individual * len);

  for(i = 0; i < geno->n_individual; i++) {
    geno->genotype[i] = &(geno->genotype_space[i * len]);
  }
}

void read_snp_info(snp_info_struct* snp_info, int mode) {
  char* fid;

  if (mode == TTANDME) {
 
    fid = strtok(line, " \t\n");
    snp_info->snp_id = Malloc(char, strlen(fid) + 1);
    strcpy(snp_info->snp_id, fid);

    fid = strtok(NULL, " \t\n");
    snp_info->chromosome = Malloc(char, strlen(fid) + 2);
    strcpy(snp_info->chromosome, fid);
    
    fid = strtok(NULL, " \t\n");
    snp_info->position = atoi(fid);

    snp_info->morgan = Malloc(char, 2);
    strcpy(snp_info->morgan, "0");

    fid = strtok(NULL, " \t\n");
    if (fid[0] != fid[1]) {
      snp_info->snp_major = fid[0];
      snp_info->snp_minor = fid[1];
    } else if (fid[0] == '-' && fid[1] == '-') {
      snp_info->snp_major = MISSING_ALLELE;
      snp_info->snp_minor = MISSING_ALLELE;
    } else {
      snp_info->snp_major = fid[0];
      snp_info->snp_minor = MISSING_ALLELE;
    }

  } else {

    fid = strtok(line, " \t\n");
    snp_info->chromosome = Malloc(char, strlen(fid) + 2);
    strcpy(snp_info->chromosome, fid);
    
    fid = strtok(NULL, " \t\n");
    snp_info->snp_id = Malloc(char, strlen(fid) + 1);
    strcpy(snp_info->snp_id, fid);
    
    fid = strtok(NULL, " \t\n");
    snp_info->morgan = Malloc(char, strlen(fid) + 1);
    strcpy(snp_info->morgan, fid);

    fid = strtok(NULL, " \t\n");
    snp_info->position = atoi(fid);

    if (mode == BIM) {
      fid = strtok(NULL, " \t\n");
      snp_info->snp_major = fid[0];

      fid = strtok(NULL, " \t\n");
      snp_info->snp_minor = fid[0];

    } else if (mode == SPA_MODEL) {
      fid = strtok(NULL, " \t\n");
      snp_info->snp_minor = fid[0];

      fid = strtok(NULL, " \t\n");
      snp_info->snp_major = fid[0];

    }
  }
}

void read_individual_info(individual_info_struct* individual_info) {
  char* fid;

  fid = strtok(line, " \t\n");
  individual_info->family_id = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->family_id, fid);

  fid = strtok(NULL, " \t\n");
  individual_info->individual_id = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->individual_id, fid);

  fid = strtok(NULL, " \t\n");
  individual_info->paternal_id = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->paternal_id, fid);
  
  fid = strtok(NULL, " \t\n");
  individual_info->maternal_id = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->maternal_id, fid);
  
  fid = strtok(NULL, " \t\n");
  individual_info->sex = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->sex, fid);
  
  fid = strtok(NULL, " \t\n");
  individual_info->phenotype = Malloc(char, strlen(fid) + 1);
  strcpy(individual_info->phenotype, fid);
}

void read_bedfile(const char* filename,
                  spa_data* geno,
                  const spa_parameter* param) {
  int rows, columns, i, j, k, m, l, c;

  char fullname[max_name_len];
  char* buf;
  char header[4];
  FILE* fp;

  // read .bim file
  strcpy(fullname, filename);
  strcat(fullname, ".bim");
  fp = spa_open_file((const char*)fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);

  rewind(fp);
  geno->snp_info = Malloc(snp_info_struct, rows);
  geno->n_snp = rows;
  for(i = 0; i < geno->n_snp; i++) {
    readline(fp);
    read_snp_info(geno->snp_info + i, BIM);
  }

  // read .fam file
  strcpy(fullname, filename);
  strcat(fullname, ".fam");
  fp = spa_open_file((const char*)fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);

  rewind(fp);
  geno->individual_info = Malloc(individual_info_struct, rows);
  geno->n_individual = rows;
  for(i = 0; i < geno->n_individual; i++) {
    readline(fp);
    read_individual_info(geno->individual_info + i);
  }

  // read .bed file
  strcpy(fullname, filename);
  strcat(fullname, ".bed");
  fp = spa_open_file((const char*)fullname, "r");
  sprintf(line, "reading binary file %s", fullname);
  spa_message(line, SHORT, param);

  // read header
  fread(header, 1, 3, fp);
  
  if(header[0] != 0x6C || header[1] != 0x1B) {
    sprintf(line, "%s is an illegal .bed file", fullname);
    spa_error_exit(line);
  }

  allocate_genotype(geno);
  
  if(header[2] == 1) { 
    // snp major
    l = geno->n_individual / GENOTYPE_PER_BYTE;
    buf = Malloc(char, l);
    m = geno->n_individual % GENOTYPE_PER_BYTE;
    
    for(i = 0; i < geno->n_snp; i++) {
      fread(buf, 1, l, fp);

      c = 0;
      for(j = 0; j < l; j++) {
        for(k = 0; k < GENOTYPE_PER_BYTE; k++) {
          set_genotype(geno->genotype, c++, i, read_mask_genotype(buf[j], k));
        }
      }
    
      if(m > 0) {
        fread(buf, 1, 1, fp);
        for(k = 0; k < m; k++) {
          set_genotype(geno->genotype, c++, i, read_mask_genotype(buf[0], k));
        }
      }
    }
    free(buf);
  } else if(header[2] == 0) { 
    // individual major
    fread(geno->genotype_space, 
          1, 
          geno->n_individual * 
            (int)ceil(((double)geno->n_snp) / GENOTYPE_PER_BYTE),
          fp);
  } else {
    sprintf(line, "%s is an illegal .bed file", fullname);
    spa_error_exit(line);
  }

  // flip snp_major and snp_minor if they are wrongly set
  snp_flip(geno);
}

void snp_flip(spa_data* geno) {
  
  int i, j, k, m;
  char sp;

  for(j = 0; j < geno->n_snp; j++) {
    k = 0;
    m = 0;
    for(i = 0; i < geno->n_individual; i++) {
      switch(get_genotype(geno->genotype, i, j)) {
        case HOMO_MAJOR:
          k++;
          break;
        case HOMO_MINOR:
          m++;
          break;
      }
    }

    if(m > k) {
      // swap minor and major allele
      for(i = 0; i < geno->n_individual; i++) {
        if(get_genotype(geno->genotype, i, j) == HOMO_MAJOR) {
          set_genotype(geno->genotype, i, j, HOMO_MINOR);
        } else if(get_genotype(geno->genotype, i, j) == HOMO_MINOR) {
          set_genotype(geno->genotype, i, j, HOMO_MAJOR);
        }
      }

      sp = geno->snp_info[j].snp_major;
      geno->snp_info[j].snp_major = geno->snp_info[j].snp_minor;
      geno->snp_info[j].snp_minor = sp;
    }
  }
}

void read_pedfile(const char* filename,
                  spa_data *geno,
                  const spa_parameter *param) {
  int rows, columns, i, j;

  char fullname[max_name_len];
  char* fid;
  char* buf;
  FILE* fp;

  strcpy(fullname, filename);
  strcat(fullname, ".map");
  fp = spa_open_file((const char*)fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);
  rewind(fp);

  geno->snp_info = Malloc(snp_info_struct, rows);
  geno->n_snp = rows;
  for(i = 0; i < geno->n_snp; i++) {
    readline(fp);
    read_snp_info(geno->snp_info + i, MAP);
  }

  strcpy(fullname, filename);
  strcat(fullname, ".ped");
  fp = spa_open_file((const char*) fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);

  rewind(fp);

  geno->individual_info = Malloc(individual_info_struct, rows);
  geno->n_individual = rows;
  allocate_genotype(geno);
 
  // read the whole file in first. buf is arranged in snp major order
  buf = Malloc(char, geno->n_individual * geno->n_snp * 2);

  for(i = 0; i < geno->n_individual; i++) {
    readline(fp);
    read_individual_info(geno->individual_info + i);
   
    for(j = 0; j < geno->n_snp; j++){  
      fid = strtok(NULL, " \t\n");
      buf[j * (geno->n_individual * 2) + 2 * i] = fid[0];
      fid = strtok(NULL, " \t\n");
      buf[j * (geno->n_individual * 2) + 2 * i + 1] = fid[0];
    }
  }

  for(i = 0; i < geno->n_snp; i++) {
    arrange_genotype_in_number(buf + i * (geno->n_individual * 2), geno, i);
  }
  
  free(buf);
  spa_message("ped file read successfully", SHORT, param);
}

void read_tpedfile(const char* filename,
                   spa_data *geno,
                   const spa_parameter *param) {
  int rows, columns, i, j;

  char fullname[max_name_len];
  char* buf;
  char* fid;
  FILE *fp;

  strcpy(fullname, filename); 
  strcat(fullname, ".tfam");
  fp = spa_open_file((const char*)fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);
  
  rewind(fp);
  geno->individual_info = Malloc(individual_info_struct, rows);
  geno->n_individual = rows;
  for(i = 0; i < geno->n_individual; i++) {
    readline(fp);
    read_individual_info(geno->individual_info + i);
  }

  strcpy(fullname, filename);
  strcat(fullname, ".tped");
  fp = spa_open_file((const char*) fullname, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", fullname, rows, columns);
  spa_message(line, MEDIUM, param);

  rewind(fp);
  
  geno->snp_info = Malloc(snp_info_struct, rows);
  geno->n_snp = rows;
  
  allocate_genotype(geno);

  buf = Malloc(char, geno->n_individual * 2);

  for(i = 0; i < geno->n_snp; i++) {
    readline(fp);
 
    read_snp_info(geno->snp_info + i, MAP);

    for(j = 0; j < geno->n_individual * 2; j++) {
      fid = strtok(NULL, " \t\n");
      buf[j] = fid[0];
    }

    arrange_genotype_in_number(buf, geno, i);
  }

  free(buf);
  spa_message("tped file read successfully", SHORT, param);
}

void arrange_genotype_in_number(char* buf, spa_data* geno, int i) {
   
  int j;
  std::map<char, int> count;
  
  for(j = 0; j < geno->n_individual * 2; j++) {
    if (buf[j] != MISSING_ALLELE) {
      count[buf[j]] ++;
    }
  }
  
  if (count.size() == 0) {
    sprintf(line, "The %d-th SNP is all missing."
                  "Please remove this SNP.", i);
    spa_warning(line);
    geno->snp_info[i].snp_minor = MISSING_ALLELE;
    geno->snp_info[i].snp_major = MISSING_ALLELE;
  } else if (count.size() == 1) {
    std::map<char, int>::const_iterator it = count.begin();
    geno->snp_info[i].snp_minor = MISSING_ALLELE;
    geno->snp_info[i].snp_major = it->first;
  } else {
    if (count.size() > 2) {
      sprintf(line, "The %d-th SNP is not biallelic."
                    "This is supported in current SPA software."
                    "Result might not be accurate", i);
      spa_warning(line);
    }
    
    // get the first two major alleles
    int c[2] = {0, 0};
    char snp[2] = {'0', '0'};
    for (std::map<char, int>::const_iterator it = count.begin();
         it != count.end();
         it++) {

      if (it->second > c[1]) {
        if (it->second > c[0]) {
          snp[1] = snp[0];
          c[1] = c[0];
          snp[0] = it->first;
          c[0] = it->second;
        } else {
          snp[1] = it->first;
          c[1] = it->second;
        }
      } 
    }

    geno->snp_info[i].snp_major = snp[0];
    geno->snp_info[i].snp_minor = snp[1];
  }
 
  for(j = 0; j < geno->n_individual; j++) {
    if (buf[2*j] == MISSING_ALLELE || buf[2*j+1] == MISSING_ALLELE) {
      set_genotype(geno->genotype, j, i, MISSING);
    } else if(buf[2*j] != buf[2*j+1]) {
      set_genotype(geno->genotype, j, i, HETER);
    } else if(buf[2*j] == geno->snp_info[i].snp_major) {
      set_genotype(geno->genotype, j, i, HOMO_MAJOR);
    } else if(buf[2*j] == geno->snp_info[i].snp_minor) {
      set_genotype(geno->genotype, j, i, HOMO_MINOR);
    } else {
      set_genotype(geno->genotype, j, i, MISSING);
    }
  }
}

void read_gfile(const char* filename,
                spa_data *geno,
                const spa_parameter *param) {
  int rows, columns, i, k, g;

  FILE *fp = spa_open_file((const char*)filename, "r");
  
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", filename, rows, columns);
  spa_message(line, MEDIUM, param);

  rewind(fp);
 
  geno->n_individual = rows;
  geno->n_snp = columns;
  geno->individual_info = Malloc(individual_info_struct, rows);
  geno->snp_info = Malloc(snp_info_struct, columns);
  allocate_genotype(geno);
 
  for(i = 0; i < geno->n_individual; i++) {
    sprintf(line, "%d %d %d %d %d %d", i+1, i+1, 0, 0, 0, 0);
    read_individual_info(geno->individual_info + i);
  }

  for(i = 0; i < geno->n_snp; i++) {
    sprintf(line, "%d SNP%d %d %d %d %d", 0, i+1, 0, 0, 1, 2);
    read_snp_info(geno->snp_info + i, SPA_MODEL);
  }

  for(i = 0; i < geno->n_individual; i++) {
    readline(fp);
    for(k = 0; k < geno->n_snp; k++) {
      if (k == 0) {
        g = atoi(strtok(line, " \t\n"));
      } else {
        g = atoi(strtok(NULL, " \t\n"));
      }

      if(g == 0 || g == 1 || g == 2 || g == -1) {
        set_genotype(geno->genotype, i, k, g);
      } else {
        sprintf(line, 
                "genotype %d at %d-th row %d-th column is not allowed. "
                "Only 0/1/2/-1 is allowed\n",
                g,
                i+1,
                k+1);
        spa_error_exit(line);
      }
      set_genotype(geno->genotype, i, k, g);
    }
  }

  spa_message("genotype file read successfully", SHORT, param);
}

void count_file(FILE *fp, int *rows, int *columns) {
  int rw = 0;
  int cl = 0;
  int k;
  char *tok;

  while(readline(fp) != NULL) {
    k = 0;
    if (line[0] != '#') {
      tok = strtok(line, " \t");
      ++k;
      while(1) {
        tok = strtok(NULL, " \t");
        // check '\n' as ' ' may be after the last element
        if(tok == NULL || *tok == '\n') 
          break;
        ++k;
      }
      
      if(k == 0){
        spa_error_exit("empty row\n");
      }
      if(cl == 0) {
        cl = k;
      } else {
        if(cl != k) {
          spa_error_exit("different number of columns\n");
        }
      }
      rw++;
    }
  }
  
  (*rows) = rw;
  (*columns) = cl;
}

void read_mfile(const char* filename,
                spa_data* geno,
                const spa_parameter* param) {
  
  int rows, columns, i;

  FILE *fp;

  fp = spa_open_file(filename, "r");
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", filename, rows, columns);
  spa_message(line, MEDIUM, param);
  
  geno->snp_info = Malloc(snp_info_struct, rows);
  geno->n_snp = rows;
  
  rewind(fp);
  readline(fp);
  while (line[0] == '#') {
    readline(fp);
  }
  for(i = 0; i < geno->n_snp; i++) {
    read_snp_info(geno->snp_info + i, TTANDME);
    readline(fp);
  }
  
  geno->individual_info = Malloc(individual_info_struct, 1);
  geno->n_individual = 1;
  strcpy(line, "1 1 0 0 0 0");
  read_individual_info(geno->individual_info);
  
  allocate_genotype(geno);

  for(i = 0; i < geno->n_snp; i++) {
    if (geno->snp_info[i].snp_major == MISSING_ALLELE &&
        geno->snp_info[i].snp_minor == MISSING_ALLELE) {
      set_genotype(geno->genotype, 0, i, MISSING);
    } else if (geno->snp_info[i].snp_major != MISSING_ALLELE &&
               geno->snp_info[i].snp_minor == MISSING_ALLELE) {
      set_genotype(geno->genotype, 0, i, HOMO_MAJOR);
    } else if (geno->snp_info[i].snp_major != geno->snp_info[i].snp_minor) {
      set_genotype(geno->genotype, 0, i, HETER);
    }
  }

  spa_message("23andme file read successfully", SHORT, param);
}


void read_location_ilfile(const char* filename,
                          spa_model* model,
                          const spa_parameter* param) {
  FILE *fp;
  int rows, columns;
  int i, j;

  fp = spa_open_file(filename, "r");

  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", filename, rows, columns);
  spa_message(line, MEDIUM, param);

  model->n_individual = rows;
  model->individual_info = Malloc(individual_info_struct, rows);
  
  model->x = Malloc(double*, rows);
  model->x_space = Malloc(double, rows * param->dimension);
  for (i = 0; i < rows; i++) {
    model->x[i] = &(model->x_space[i * param->dimension]);
  }
  
  rewind(fp);
  
  if(columns - param->dimension == LOC_FILE_HEADER) {
    for(i = 0; i < model->n_individual; i++) {
      readline(fp);
      read_individual_info(model->individual_info + i);

      for(j = 0; j < param->dimension; j++) {
        model->x[i][j] = atof(strtok(NULL, " \t\n"));
      }
    }
  } else {
    sprintf(line, 
            "dimension inconsistent."
            "File %s should have %d columns but %d columns found\n"
            "Please check SPA manual for details",
            filename,
            LOC_FILE_HEADER + param->dimension,
            columns);
    spa_error_exit(line);
  }  
}

void write_location_olfile(const char* filename,
                           const spa_model* model,
                           const spa_parameter* param) {
  int i, j;
  FILE* fp;

  fp = spa_open_file(filename, "w");

  for(i = 0; i < model->n_individual; i++) {
    fprintf(fp, "%s\t", model->individual_info[i].family_id);
    fprintf(fp, "%s\t", model->individual_info[i].individual_id);
    fprintf(fp, "%s\t", model->individual_info[i].paternal_id);
    fprintf(fp, "%s\t", model->individual_info[i].maternal_id);
    fprintf(fp, "%s\t", model->individual_info[i].sex);
    fprintf(fp, "%s\t", model->individual_info[i].phenotype);
    for(j = 0; j < param->dimension * param->generation - 1; j++) {
      fprintf(fp, "%.10f\t", model->x[i][j]);
    }
    fprintf(fp, "%.10f\n", model->x[i][j]);
  }
  fclose(fp);
}

void write_html_location_olfile(const char* filename,
                                const spa_model* model,
                                const spa_parameter* param) {

  FILE* fp;
  double lat_center;
  double lng_center;

  sprintf(line, "%s.html", filename);
  fp = spa_open_file(line, "w");

  if (param->generation == SELF) {
    
    sprintf(line, 
"var pos = new google.maps.LatLng(%f, %f);\n"
"var infowindow = new google.maps.InfoWindow({\n"
"      map: map,\n"
"      position: pos,\n"
"      content: 'Your ancestry is from here'\n"
"    });",
    model->x[0][0],
    model->x[0][1]
    );
    lat_center = model->x[0][0];
    lng_center = model->x[0][1];
  } else if (param->generation == PARENT) {
    
    sprintf(line, 
"var ppos = new google.maps.LatLng(%f, %f);\n"
"var mpos = new google.maps.LatLng(%f, %f);\n"
"var infowindow = new google.maps.InfoWindow({\n"
"      map: map,\n"
"      position: ppos,\n"
"      content: 'Your first ancestry is from here'\n"
"    });\n"
"var infowindow = new google.maps.InfoWindow({\n"
"      map: map,\n"
"      position: mpos,\n"
"      content: 'Your second ancestry is from here'\n"
"    });",
    model->x[0][0],
    model->x[0][1],
    model->x[0][2],
    model->x[0][3]
    );
    lat_center = model->x[0][0] / 2 + model->x[0][2] / 2;
    lng_center = model->x[0][1] / 2 + model->x[0][3] / 2;
  }
  
  fprintf(fp,
"<!DOCTYPE html>\n"
"  <html>\n"
"    <head>\n"
"      <meta name=\"viewport\" content=\"initial-scale=1.0, user-scalable=no\" />\n"
"      <style type=\"text/css\">\n"
"        html { height: 100%% }\n"
"        body { height: 100%%; margin: 0; padding: 0 }\n"
"      </style>\n"
"      <script type=\"text/javascript\"\n"
"        src=\"http://maps.googleapis.com/maps/api/js?sensor=false&"
"key=AIzaSyDmlI1dvRWM1Ohy_NIvb7on_DbpErFwfN8\">\n"
"      </script>\n"
"      <script type=\"text/javascript\">\n"
"        function initialize() {\n"
"          var mapOptions = {\n"
"            center: new google.maps.LatLng(%f, %f),\n"
"            zoom: 3,\n"
"            mapTypeId: google.maps.MapTypeId.ROADMAP\n"
"          };\n"
"          var map = new google.maps.Map(document.getElementById(\"map_canvas\"),\n"
"              mapOptions);\n"
"%s\n"
"        }\n"
"      </script>\n"
"    </head>\n"
"    <body onload=\"initialize()\">\n"
"      <div id=\"map_canvas\" style=\"width:100%%; height:100%%\"></div>\n"
"    </body>\n"
"  </html>",
  lat_center,
  lng_center,
  line
  );

  fclose(fp);
}

void read_model_imfile(const char* filename, 
                       spa_model *model,
                       const spa_parameter *param) {
  FILE *fp;
  int rows, columns, i, j;

  fp = spa_open_file(filename, "r");
  
  count_file(fp, &rows, &columns);
  sprintf(line, "%s: %d row %d column", filename, rows, columns);
  spa_message(line, MEDIUM, param);

  model->n_snp = rows;
  model->snp_info = Malloc(snp_info_struct, rows);
  
  model->coef_a = Malloc(double*, rows);
  model->coef_a_space = Malloc(double, rows * param->dimension);
  model->score = Malloc(double, rows);
  for (i = 0; i < rows; i++) {
    model->coef_a[i] = &(model->coef_a_space[i * param->dimension]);
  }
  model->coef_b = Malloc(double, rows);

  rewind(fp);

  if(columns - param->dimension == MODEL_FILE_HEADER + MODEL_FILE_TAILER) {
    for(i = 0; i < model->n_snp; i++) {
      readline(fp);
      read_snp_info(model->snp_info + i, SPA_MODEL);

      for(j = 0; j < param->dimension; j++) {  
        model->coef_a[i][j] = atof(strtok(NULL, " \t\n"));
      }
      
      model->coef_b[i] = atof(strtok(NULL, " \t\n"));
      model->score[i] = atof(strtok(NULL, " \t\n"));
    }
  } else {
    sprintf(line,
            "dimension inconsistent\n"
            "File %s should have %d columns, %d found,"
            "Please check SPA manual for detail\n",
            filename,
            param->dimension + MODEL_FILE_HEADER + MODEL_FILE_TAILER,
            columns);
    spa_error_exit(line);
  }
}

void write_model_omfile(const char* filename,
                        const spa_model *model,
                        const spa_parameter *param) {
  int i, j;
  FILE *fp;
  
  fp = spa_open_file(filename, "w");

  for(i = 0; i < model->n_snp; i++) {
    fprintf(fp, "%s\t", model->snp_info[i].chromosome);
    fprintf(fp, "%s\t", model->snp_info[i].snp_id);
    fprintf(fp, "%s\t", model->snp_info[i].morgan);
    fprintf(fp, "%d\t", model->snp_info[i].position);
    fprintf(fp, "%c\t", model->snp_info[i].snp_minor);
    fprintf(fp, "%c\t", model->snp_info[i].snp_major);

    for(j = 0; j < param->dimension; j++) { 
      fprintf(fp, "%.10f\t", model->coef_a[i][j]);
    }
    fprintf(fp, "%.10f\t", model->coef_b[i]);
    fprintf(fp, "%.10f\n", model->score[i]);
  }
  fclose(fp);
}

char get_genotype(char** genotype, int i, int j) {
  int k, m;

  m = j % GENOTYPE_PER_BYTE;
  k = j / GENOTYPE_PER_BYTE;

  return read_mask_genotype(genotype[i][k], m);
}

void set_genotype(char** genotype, int i, int j, char g) {
  int k, m;

  m = j % GENOTYPE_PER_BYTE;
  k = j / GENOTYPE_PER_BYTE;

  write_mask_genotype(&(genotype[i][k]), m, g);
}

char read_mask_genotype(char byte, char c) {
  char genotype = (char)-2;

  switch((unsigned char)(byte & mask[(int)c])) {
    case 0x00:
      genotype = (char) HOMO_MAJOR;
      break;
    case 0x01:
      genotype = (char) MISSING;
      break;
    case 0x02:
      genotype = (char) HETER;
      break;
    case 0x03:
      genotype = (char) HOMO_MINOR;
      break;
    case 0x04:
      genotype = (char) MISSING;
      break;
    case 0x08:
      genotype = (char) HETER;
      break;
    case 0x0C:
      genotype = (char) HOMO_MINOR;
      break;
    case 0x10:
      genotype = (char) MISSING;
      break;
    case 0x20:
      genotype = (char) HETER;
      break;
    case 0x30:
      genotype = (char) HOMO_MINOR;
      break;
    case 0x40:
      genotype = (char) MISSING;
      break;
    case 0x80:
      genotype = (char) HETER;
      break;
    case 0xC0:
      genotype = (char) HOMO_MINOR;
      break;
    default:
      spa_error_exit("read_mask_genotype error\n");
  }

  return genotype;
}

void write_mask_genotype(char* byte, char c, char g) {
  // erase
  switch(c) {
    case 0:
      *byte &= 0xFC;
      break;
    case 1:
      *byte &= 0xF3;
      break;
    case 2:
      *byte &= 0xCF;
      break;
    case 3:
      *byte &= 0x3F;
      break;
    default:
      spa_error_exit("write_mask_genotype error");
  }

  // set new value
  switch(g) {
    case HOMO_MAJOR:
      break;
    case MISSING:
      switch(c) {
        case 0:
          *byte |= 0x01;
          break;
        case 1:
          *byte |= 0x04;
          break;
        case 2:
          *byte |= 0x10;
          break;
        case 3:
          *byte |= 0x40;
          break;
        default:
          spa_error_exit("write_mask_genotype error");
      }
      break;
    case HETER:
      switch(c) {
        case 0:
          *byte |= 0x02;
          break;
        case 1:
          *byte |= 0x08;
          break;
        case 2:
          *byte |= 0x20;
          break;
        case 3:
          *byte |= 0x80;
          break;
        default:
          spa_error_exit("write_mask_genotype error");
      }
      break;
    case HOMO_MINOR:
      switch(c) {
        case 0:
          *byte |= 0x03;
          break;
        case 1:
          *byte |= 0x0C;
          break;
        case 2:
          *byte |= 0x30;
          break;
        case 3:
          *byte |= 0xC0;
          break;
        default:
          spa_error_exit("write_mask_genotype error");
      }
      break;
    default:
      spa_error_exit("write_mask_genotype error");
  }
}

