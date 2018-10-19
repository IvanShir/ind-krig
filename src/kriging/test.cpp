#include <iostream>
//#include "../FileParser/fileparser.h"
#include "kriging.h"
#include <string>
#include <ctime>

int main() {
  KrigIndic::FileParser parser;
  KrigIndic::Model m;
  const double PI  = 3.141592653589793238463;
  parser.readfile("../../grids/input/BW_3Dmodel.grd", m);
  //parser.readfile("../../grids/input/small_benchmark.grd", m);
  std::cout << "Enter variogram params: " << "\n";
  KrigIndic::Variogram var;
  std::cin >> var.sill >> var.radX >> var.radY >> var.radZ >> var.angXY;
  KrigIndic::Kriging krig;
  clock_t begin = std::clock();
  krig.kriging(m, var);
  clock_t end = std::clock();
  parser.writefile("test.grd", m);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << elapsed_secs << "\n";
}
