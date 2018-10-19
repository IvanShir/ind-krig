#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <unordered_map>
#include <math.h>
#include "Eigen/Dense"
#include <ctime>

#include "kriging.h"

using namespace Eigen;

namespace KrigIndic {

MatrixXd Kriging::var_matrix(std::vector<KrigIndic::Point>& v, KrigIndic::Variogram vario) {
  int n = v.size();
  MatrixXd var(n + 1, n + 1);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double xdist = fabs(v[i].x - v[j].x);
      double ydist = fabs(v[i].y - v[j].y);
      double zdist = fabs(v[i].z - v[j].z);
      if (vario.angXY != 0.0) {
        xdist = fabs((v[i].x - v[j].x) * cos(vario.angXY) -
                     (v[i].y - v[j].y) * sin(vario.angXY));
        ydist = fabs((v[i].x - v[j].x) * sin(vario.angXY) +
                     (v[i].y - v[j].y) * cos(vario.angXY));
      }
      double h = sqrt(pow(xdist / vario.radX, 2) + pow(ydist / vario.radY, 2) + pow(zdist / vario.radZ, 2));
      var(i, j) = var(j, i) = vario.sill * (1.0 -  exp(-h));
    }
    var(i, i) = 0.0;
    var(n, i) = var(i, n) = 1.0;
  }
  var(n, n) = 0;
  return var;
}

VectorXd Kriging::var_vector(const std::vector<KrigIndic::Point>& v,
                             const KrigIndic::Point p, const KrigIndic::Variogram &vario) {
  int n = v.size();
  VectorXd var(n + 1);
  for (int i = 0; i < n; i++) {
      double xdist = fabs(v[i].x - p.x);
      double ydist = fabs(v[i].y - p.y);
      double zdist = fabs(v[i].z - p.z);
      if (vario.angXY != 0.0) {
        xdist = fabs((v[i].x - p.x) * cos(vario.angXY) -
                     (v[i].y - p.y) * sin(vario.angXY));
        ydist = fabs((v[i].x - p.x) * sin(vario.angXY) +
                     (v[i].y - p.y) * cos(vario.angXY));
      }
      double h = sqrt(pow(xdist / vario.radX, 2) +
                      pow(ydist / vario.radY, 2) +
                      pow(zdist / vario.radZ, 2));
      var(i) = vario.sill * (1.0 -  exp(-h));
  }
  var(n) = 1.0;
  return var;
}

void Kriging::lito_classification(KrigIndic::Model& model) {
  for (int i = 0; i < model.NumOfCells(); i++) {
    double max = -1;
    double classifier = 0;
    for (std::map<double, std::vector<double>>::iterator it = model.fac_prob.begin();
     it != model.fac_prob.end(); ++it) {
        if ((it->second)[i] > max) {
          max = (it->second)[i];
          classifier = it->first;
        }
    }
  model.fac_class.push_back(classifier);
  }
}

VectorXd Kriging::distances(KrigIndic::Model& model,
                        std::vector<KrigIndic::Point>& initial_lito_points,
                        int p_index) {
  int n = initial_lito_points.size();
  VectorXd dist(n);
  for (int i = 0; i < n; i++) {
      dist(i) = sqrt(pow(initial_lito_points[i].x - model.corn_coord[p_index].x, 2) +
                     pow(initial_lito_points[i].y - model.corn_coord[p_index].y, 2) +
                     pow(initial_lito_points[i].z - model.corn_coord[p_index].z, 2));
  }
  return dist;
}


void Kriging::nearest_cells(KrigIndic::Model& model, VectorXd& distances,
                            std::vector<KrigIndic::Point>& initial_lito_points,
                            std::vector<double>& initial_lito_vals,
                            int point_index, int n_nearest) {

  std::multimap<double, int> sorted_dist;
  for (int i = 0 ; i < initial_lito_points.size(); i++)
    sorted_dist.insert(std::pair<double, int>(distances(i), i));

  int counter = 0;
  std::vector<KrigIndic::Point> upd_initial_lito_points;
  std::vector<double> upd_initial_lito_vals;
  for (std::multimap<double, int>::iterator it = sorted_dist.begin(); counter < n_nearest; ++it) {
    upd_initial_lito_points.push_back(initial_lito_points[(*it).second]);
    upd_initial_lito_vals.push_back(initial_lito_vals[(*it).second]);
    counter++;
  }
  initial_lito_points = upd_initial_lito_points;
  initial_lito_vals = upd_initial_lito_vals;
}


void Kriging::kriging(KrigIndic::Model& model, KrigIndic::Variogram vario) {
  model.CalcCellsCoordinate();
  std::set<double> lito_types;
  std::vector<KrigIndic::Point> initial_lito_points_global;
  std::vector<double> initial_lito_vals;
  int n_nearest = 2000;
  for (int i = 0; i < model.NumOfCells(); i++) {
      if (model.lito_indic[i] != 0) {
        lito_types.insert(model.lito_indic[i]);
        initial_lito_points_global.push_back(model.corn_coord[i]);
        initial_lito_vals.push_back(model.lito_indic[i]);
      }
  }
  std::map<double, double> temp_lito_prob;
  for (std::set<double>::iterator it = lito_types.begin(); it != lito_types.end(); ++it)
    temp_lito_prob.insert(std::pair<double, double>(*it, 0));
  MatrixXd var = var_matrix(initial_lito_points_global, vario);
  MatrixXd varInverse = var.inverse();
  for (int i = 0; i < model.NumOfCells(); i++) {
    clock_t clock0 = std::clock();
    auto initial_lito_points = initial_lito_points_global;

    std::cout<< "Point " << i << "\n";
    if (initial_lito_points.size() >  n_nearest) {
      VectorXd dist = distances(model, initial_lito_points, i);
      nearest_cells(model, dist, initial_lito_points, initial_lito_vals, i, n_nearest);
      var = var_matrix(initial_lito_points, vario);
    }
    clock_t clock1 = std::clock();

    for (int j = 0; j < initial_lito_points.size(); j++) {
      temp_lito_prob[initial_lito_vals[j]] = 0;
    }
    VectorXd vec = var_vector(initial_lito_points, model.corn_coord[i], vario);
    clock_t clock2 = std::clock();
    // VectorXd w = var.householderQr().solve(vec);
    VectorXd w = varInverse * vec;
    clock_t clock3 = std::clock();
    for (int j = 0; j < initial_lito_points.size(); j++)
      temp_lito_prob[initial_lito_vals[j]] += w[j];

    for (std::set<double>::iterator it = lito_types.begin(); it != lito_types.end(); ++it) {
      model.fac_prob[*it].push_back(temp_lito_prob[*it]);
    }

    clock_t clockTot = std::clock();
    std::cout << "var_matrix: " << double(clock1 - clock0) / CLOCKS_PER_SEC << "\n";
    std::cout << "var_vector " << double(clock2 - clock1) / CLOCKS_PER_SEC << "\n";
    std::cout << "varInverse: " << double(clock3 - clock2) / CLOCKS_PER_SEC << "\n";
    std::cout << "Total Time: " << double(clockTot - clock0) / CLOCKS_PER_SEC << "\n";
    std::cout << std::string(30, '*') << "\n";
  }
  lito_classification(model);
}

}
