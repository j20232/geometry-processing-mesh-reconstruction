#include "fd_partial_derivative.h"

void fd_partial_derivative(const int nx, const int ny, const int nz,
                           const double h, const int dir,
                           Eigen::SparseMatrix<double>& D) {
  int nx_D = nx, ny_D = ny, nz_D = nz;
  if (dir == 0) {
    nx_D--;
  } else if (dir == 1) {
    ny_D--;
  } else {
    nz_D--;
  }

  D.resize(nx_D * ny_D * nz_D, nx * ny * nz);
  for (int i = 0; i < nx_D; i++) {
    for (int j = 0; j < ny_D; j++) {
      for (int k = 0; k < nz_D; k++) {
        auto grid_idx = i + nx_D * j + ny_D * nx_D * k;
        int l_prev = i + nx * j + ny * nx * k;
        int l_curr;
        if (dir == 0) {
          l_curr = (i + 1) + nx * j + ny * nx * k;
        } else if (dir == 1) {
          l_curr = i + nx * (j + 1) + ny * nx * k;
        } else {
          l_curr = i + nx * j + ny * nx * (k + 1);
        }
        D.insert(grid_idx, l_prev) = -1;
        D.insert(grid_idx, l_curr) = 1;
      }
    }
  }
}
