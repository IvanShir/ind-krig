#pragma once
#include <iostream>
#include "FileParser/fileparser.h"
#include "Eigen/Dense"

namespace KrigIndic {

  //Параметры вариограммы.
  //Exponential variogram is used: var = sill(1 - exp(-h))
  //h = sqrt((xdist / radX)^2 + (ydist / radY)^2 + (zdist / radZ)^2)
  //xdist, ydist, zdist - distances btwn coordinates of two points in 3D space
  struct Variogram {
    double sill;
    double radX;  //radius in the main direction
    double radY;  //radius in direction orthogonal to radX in XY
    double radZ;  //vertical radius
    double angXY; //angle of rotation of coordinate axes in XY
  };

  class Kriging {
  public:
    void kriging(KrigIndic::Model&, Variogram, int);
    void lito_classification(KrigIndic::Model&);
    Eigen::VectorXd distances(KrigIndic::Model&, std::vector<KrigIndic::Point>&, int);
    void nearest_cells(KrigIndic::Model& ,Eigen::VectorXd&,
                        std::vector<KrigIndic::Point>&,
                        std::vector<double>&,
                        int, int);
    Eigen::MatrixXd var_matrix(std::vector<KrigIndic::Point>&, Variogram);
    Eigen::VectorXd var_vector(const std::vector<KrigIndic::Point>& v, const KrigIndic::Point p, const KrigIndic::Variogram &vario);
  };
}
