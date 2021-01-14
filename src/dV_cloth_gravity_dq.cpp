#include <dV_cloth_gravity_dq.h>
#include <iostream>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    //  M - sparse mass matrix for the entire mesh
//  g - the acceleration due to gravity
//Output:
//  fg - the gradient of the gravitational potential for the entire mesh

    fg.resize(M.rows());
    fg.setZero();

    Eigen::VectorXd ng;
    ng.resize(M.rows());
    for (int i = 0; i < ng.size()/3; i++) {
        ng[3.0 * i] = g[0];
        ng[3.0 * i + 1.0] = g[1];
        ng[3.0 * i + 2.0] = g[2];
    }

    fg = -M * ng;
    
}
