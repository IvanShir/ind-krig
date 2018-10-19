#include <iostream>
#include <vector>
#include <map>

#include "model.h"

namespace KrigIndic {

  int Model::NumOfCells() {
    return nx * ny * nz;
  }

  //defines top left corner coordinate of each cell and sets it to "corn_coord"
  //data is used later for calculating distances btwn cells
  void Model::CalcCellsCoordinate() {
    int counter = 1;
    int z_counter = 1;
    int n = NumOfCells();
    for (int i = 0; i < n; i++) {
        Point p = {coord[counter - 1].xtop,
                   coord[counter - 1].ytop,
                   zcorn[z_counter - 1] * 2};
        corn_coord.push_back(p);
        if (z_counter % nx == 0)
          z_counter += nx;
        if ((z_counter) % (4 * nx * ny) == 0)
          z_counter += 4 * nx * ny;
        z_counter++;
        counter++;
        if (counter % (nx + 1) == 0)
          counter++;
        if (counter == (ny * (nx + 1) + 1))
          counter = 1;
    }
  }
}
