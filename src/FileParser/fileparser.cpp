#include "fileparser.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

namespace KrigIndic {

  //read SPECGRID, COORDSYS, COORD, ZCORN, ACTNUM, LITO_INDIK
  void FileParser::readfile(std::string path, KrigIndic::Model& model) {
    std::cout << "Reading file..." <<"\n";
    std::ifstream infile(path);
    if (infile.fail())
      std::cout << "Error opening file occured: check file name\n";
    else {
      while (!infile.eof()) {
        std::string line;
        getline(infile, line);
        if (line == "SPECGRID")
          read_SPECGRID(infile, model);
        if (line == "COORD")
          read_COORD(infile, model);
        if (line == "ZCORN")
          read_ZCORN(infile, model);
        if (line == "ACTNUM")
          read_ACTNUM(infile, model);
        if (line == "LITO_INDIK")
          read_LITO_INDIC(infile, model);
      }
    }
  }

  void FileParser::read_SPECGRID(std::ifstream& stream, KrigIndic::Model& model) {
    std::string line;
    stream >> model.nx >> model.ny >> model.nz;
  }

  void FileParser::read_COORD(std::ifstream& stream, KrigIndic::Model& model) {
    int n_pairs = (model.nx + 1) * (model.ny + 1); //number of pillars
    for (int i = 0; i < n_pairs; i++) {
      KrigIndic::Pillar p;
      stream >> p.xtop >> p.ytop >> p.ztop >> p.xbottom >> p.ybottom >> p.zbottom;
      model.coord.push_back(p);
    }
  }

  void FileParser::read_ZCORN(std::ifstream& stream, KrigIndic::Model& model) {
    int n_corners = 8; //number of corners in cube
    int n_points = n_corners * model.nx * model.ny * model.nz;
    for (int i = 0; i < n_points; i++) {
      double d;
      stream >> d;
      model.zcorn.push_back(d);
    }
  }

  void FileParser::read_ACTNUM(std::ifstream& stream, KrigIndic::Model& model) {
    int n_cells = model.nx * model.ny * model.nz;
    for (int i = 0; i < n_cells; i++) {
      int d;
      stream >> d;
      model.actnum.push_back(d);
    }
  }

void FileParser::read_LITO_INDIC(std::ifstream& stream, KrigIndic::Model& model) {
  int n_cells = model.nx * model.ny * model.nz;
  for (int i = 0; i < n_cells; i++) {
    double d;
    stream >> d;
    model.lito_indic.push_back(d);
  }
}


//write SPECGRID, COORDSYS, COORD, ZCORN, ACTNUM, LITO_INDIK, facies prob
void FileParser::writefile(std::string path, KrigIndic::Model& model) {
  std::cout << "Writing file..." << "\n";
  std::ofstream outfile(path);
  if (outfile.fail())
    std::cout << "Error writing file occured\n";
  else {
    write_SPECGRID(outfile, model);
    write_COORDSYS(outfile, model);
    write_COORD(outfile, model);
    write_ZCORN(outfile, model);
    write_ACTNUM(outfile, model);
    write_LITO_INDIC(outfile, model);
    if (model.fac_prob.size() > 0) {
      write_FACIES_PROB(outfile, model);
      write_FACIES_CLASS(outfile, model);
    }
    outfile << "/\n";
  }
}

void FileParser::write_SPECGRID(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "SPECGRID\n";
  stream << " " << model.nx << " " << model.ny << " " << model.nz << " 1 F /\n";
  stream << "\n";
}

void FileParser::write_COORDSYS(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "COORDSYS\n";
  stream << " " << "1 " << model.nz << " 'INCOMP  ' /\n";
  stream << "\n";
}

void FileParser::write_COORD(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "COORD\n";
  int n_pairs = (model.nx + 1) * (model.ny + 1);
  for (int i = 0; i < n_pairs; i++) {
    Pillar p = model.coord[i];
    stream << p.xtop << " " << p.ytop << " " << p.ztop << " "
           << p.xbottom << " " << p.ybottom << " " << p.zbottom << "\n";
  }
  stream << "/\n\n";
}

void FileParser::write_ZCORN(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "ZCORN\n";
  int n_corners = 8;
  int n_points = n_corners * model.nx * model.ny * model.nz;
  int counter = 0;
  for (int i = 0; i < n_points; i++){
    stream << model.zcorn[i] << " ";
    if (++counter == N_CLMNS) {
      stream << "\n";
      counter = 0;
    }
  }
  stream << "/\n\n";
}

void FileParser::write_ACTNUM(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "ACTNUM\n";
  int n_cells = model.nx * model.ny * model.nz;
  int counter = 0;
  for (int i = 0; i < n_cells; i++) {
    stream << model.actnum[i] << " ";
    if (++counter == N_CLMNS_CELLS) {
      stream << "\n";
      counter = 0;
    }
  }
  stream << "/\n\n";
}

void FileParser::write_LITO_INDIC(std::ofstream& stream, KrigIndic::Model& model) {
  stream << "LITO_INDIK\n";
  int n_cells = model.NumOfCells();
  int counter = 0;
  for (int i = 0; i < n_cells; i++) {
    stream << model.lito_indic[i] << " ";
    if (++counter == N_CLMNS_CELLS) {
      stream << "\n";
      counter = 0;
    }
  }
  stream << "/\n\n";
}

void FileParser::write_FACIES_PROB(std::ofstream& stream, KrigIndic::Model& model) {
  int counter = 0;
  for (std::map<double, std::vector<double>>::iterator it = model.fac_prob.begin(); it != model.fac_prob.end(); ++it) {
      stream << "LITO_" << it->first << "\n";
      for (int i = 0; i < (it->second).size(); i++) {
        stream << (it->second)[i] << " ";
        if (++counter == N_CLMNS_CELLS) {
          stream << "\n";
          counter = 0;
        }
      }
      stream << "/\n\n";
  }

}

void FileParser::write_FACIES_CLASS(std::ofstream& stream, KrigIndic::Model& model) {
  int n_cells = model.NumOfCells();
  int counter = 0;
  if (model.fac_prob.size() > 0)
  {
    stream << "FACIES_CLASSIFICATION" << "\n";
    for (int i = 0; i < n_cells; i++) {
      stream << model.fac_class[i] << " ";
      if (++counter == N_CLMNS_CELLS) {
        stream << "\n";
        counter = 0;
      }
    }
    stream << "/\n\n";
  }
}
}
