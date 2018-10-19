#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "model/model.h"

namespace KrigIndic {
class FileParser {

public:
  const int N_CLMNS = 6;        //number of columns for COORD & ZCORNS
  const int N_CLMNS_CELLS = 10; //number of columns for properties

  //read SPECGRID, COORSYS, COORD, ZCORN, ACTNUM, LITO_INDIK
  void readfile(std::string, KrigIndic::Model&);
  //read nx, ny, nz
  void read_SPECGRID(std::ifstream&, KrigIndic::Model&);
  //read pillar grid
  void read_COORD(std::ifstream&, KrigIndic::Model&);
  //read z-values for all points in grid
  void read_ZCORN(std::ifstream&, KrigIndic::Model&);
  //read active-cell flags
  void read_ACTNUM(std::ifstream&, KrigIndic::Model&);
  //read lito_facies indicators flags
  void read_LITO_INDIC(std::ifstream&, KrigIndic::Model&);

  //write data to file
  void writefile(std::string, KrigIndic::Model&);
  //write SPECGRID
  void write_SPECGRID(std::ofstream&, KrigIndic::Model&);
  //write COORDSYS
  void write_COORDSYS(std::ofstream&, KrigIndic::Model&);
  //write pillar grid
  void write_COORD(std::ofstream&, KrigIndic::Model&);
  //write z-values for all points in grid
  void write_ZCORN(std::ofstream&, KrigIndic::Model&);
  //write active-cell flags
  void write_ACTNUM(std::ofstream&, KrigIndic::Model&);
  //write initial lito-facies indicators
  void write_LITO_INDIC(std::ofstream&, KrigIndic::Model&);
  //write generated facies probabilities
  void write_FACIES_PROB(std::ofstream&, KrigIndic::Model&);
  //write generated classification field
  void write_FACIES_CLASS(std::ofstream&, KrigIndic::Model&);
};
}
