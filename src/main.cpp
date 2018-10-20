#include <iostream>
#include <thread>
#include "kriging/kriging.h"
#include "FileParser/fileparser.h"
#include <string>

int main(int argc, char **argv) {
  std::string modelName = "../grids/input/small_benchmark.grd";
  std::string outFile   = "../grids/out/test.grd";
  KrigIndic::Variogram var;
  std::cout << argc << "\n";
  if (argc == 8) {
    modelName = std::string(argv[1]);
    var.sill = std::atoi(argv[2]);
    var.radX = std::atoi(argv[3]);
    var.radY = std::atoi(argv[4]);
    var.radZ = std::atoi(argv[5]);
    var.angXY =std::atoi(argv[6]);
    outFile = std::string(argv[7]);
  } else {
    std::cout << "Enter variogram params: " << "\n";
    std::cin >> var.sill >> var.radX >> var.radY >> var.radZ >> var.angXY;
  }

  KrigIndic::FileParser parser;
  KrigIndic::Model m;
  std::cout << "reading from: " << modelName << "\n";
  parser.readfile(modelName, m);

  KrigIndic::Kriging krig;
  unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
  krig.kriging(m, var, concurentThreadsSupported);
  parser.writefile(outFile, m);
  std::cout << "to file: " << outFile << "\n";
}
