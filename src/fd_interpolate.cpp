#include "fd_interpolate.h"

void fd_interpolate(const int nx, const int ny, const int nz, const double h,
                    const Eigen::RowVector3d& corner, const Eigen::MatrixXd& P,
                    Eigen::SparseMatrix<double>& W) {
  W.resize(P.rows(), nx * ny * nz);
  typedef Eigen::Triplet<double> T;
  std::vector<T> buf;

  for (int i = 0; i < P.rows(); i++) {
    int x0 = (P(i, 0) - corner(0)) / h;
    int x1 = x0 + 1;
    double dx = (P(i, 0) - corner(0)) / h - x0;

    int y0 = (P(i, 1) - corner(1)) / h;
    int y1 = y0 + 1;
    double dy = (P(i, 1) - corner(1)) / h - y0;

    int z0 = (P(i, 2) - corner(2)) / h;
    int z1 = z0 + 1;
    double dz = (P(i, 2) - corner(2)) / h - z0;
    auto get_grid_index = [x0, x1, y0, y1, z0, z1, nx, ny, nz](bool x, bool y,
                                                               bool z) {
      int grid_index = (x) ? x1 : x0;
      grid_index += (y) ? y1 * nx : y0 * nx;
      grid_index += (z) ? z1 * ny * nx : z0 * ny * nx;
      return grid_index;
    };
    auto get_weight = [dx, dy, dz](bool x, bool y, bool z) {
      double weight = (x) ? dx : (1.0 - dx);
      weight *= (y) ? dy : (1.0 - dy);
      weight *= (z) ? dz : (1.0 - dz);
      return weight;
    };

    // (relative position for each node, grid index, weight)
    buf.push_back(T(i, get_grid_index(0, 0, 0), get_weight(0, 0, 0)));
    buf.push_back(T(i, get_grid_index(1, 0, 0), get_weight(1, 0, 0)));
    buf.push_back(T(i, get_grid_index(0, 1, 0), get_weight(0, 1, 0)));
    buf.push_back(T(i, get_grid_index(1, 1, 0), get_weight(1, 1, 0)));
    buf.push_back(T(i, get_grid_index(0, 0, 1), get_weight(0, 0, 1)));
    buf.push_back(T(i, get_grid_index(1, 0, 1), get_weight(1, 0, 1)));
    buf.push_back(T(i, get_grid_index(0, 1, 1), get_weight(0, 1, 1)));
    buf.push_back(T(i, get_grid_index(1, 1, 1), get_weight(1, 1, 1)));
  }
  W.setFromTriplets(buf.begin(), buf.end());
}
