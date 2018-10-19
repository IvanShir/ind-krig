#pragma once

#include <iostream>
#include <vector>
#include <map>

namespace KrigIndic {

  struct Point {
    double x;
    double y;
    double z;
  };

  struct Pillar {
    double xtop;
    double ytop;
    double ztop;
    double xbottom;
    double ybottom;
    double zbottom;
  };

  class Model {
  public:
    //initial data read from file
    int nx, ny, nz;                 //model sizes
    std::vector<Pillar> coord;      //pillars (COORD) of the model
    std::vector<double> zcorn;       //depths (ZCORN) of each point in model
    std::vector<int> actnum;        //active-cell flag (ACTNUM)
    std::vector<Point> corn_coord;  //coordinates of the top left corner of cells
    std::vector<double> lito_indic; //initital facies distrib ution
    //generated data
    int n_facies;                                 //number of facies in model
    std::map<double, std::vector<double>> fac_prob;//facies probabilities
    std::vector<double> fac_class;                //genarated facies classification

    //returns number of cells in model
    int NumOfCells();
    //defines top left corner coordinate of each cell and sets it to "corn_coord"
    //data is used later for calculating distances btwn cells
    void CalcCellsCoordinate();
  };

}
