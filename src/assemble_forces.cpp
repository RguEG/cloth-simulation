#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 

    //Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dX - an mx9 matrix which stores the flattened dphi/dX matrices for each tetrahedron.
//       Convert this values back to 3x3 matrices using the following code (NOTE YOU MUST USE THE TEMPORARY VARIABLE tmp_row):
//       Eigen::Matrix<double, 1,9> tmp_row
//       tmp_row = dX.row(ei); //ei is the triangle index. 
//       Eigen::Map<const Eigen::Matrix3d>(tmp_row.data())
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  a0 - the mx1 vector of undeformed triangle areas
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
    Eigen::Vector9d flocal;
    Eigen::Matrix<double, 1, 9> tmp_row;
    Eigen::Matrix3d tmp_matrix;

    f.resize(q.size());
    f.setZero();
    flocal.setZero();
    for (int i = 0; i < F.rows(); i++) {
        tmp_row = dX.row(i);
        tmp_matrix = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());
        dV_membrane_corotational_dq(flocal,q,tmp_matrix,V,F.row(i),a0[i],mu,lambda);
        f[3.0 * F(i, 0)] += flocal[0];
        f[3.0 * F(i, 0) + 1] += flocal[1];
        f[3.0 * F(i, 0) + 2] += flocal[2];
        f[3.0 * F(i, 1)] += flocal[3];
        f[3.0 * F(i, 1) + 1] += flocal[4];
        f[3.0 * F(i, 1) + 2] += flocal[5];
        f[3.0 * F(i, 2)] += flocal[6];
        f[3.0 * F(i, 2) + 1] += flocal[7];
        f[3.0 * F(i, 2) + 2] += flocal[8];

        /*f[3.0 * F(i, 0)] += flocal[0];
        f[3.0 * F(i, 0) + 1] += flocal[3];
        f[3.0 * F(i, 0) + 2] += flocal[6];
        f[3.0 * F(i, 1)] += flocal[1];
        f[3.0 * F(i, 1) + 1] += flocal[4];
        f[3.0 * F(i, 1) + 2] += flocal[7];
        f[3.0 * F(i, 2)] += flocal[2];
        f[3.0 * F(i, 2) + 1] += flocal[5];
        f[3.0 * F(i, 2) + 2] += flocal[8];*/
    }
    f = (-1.0) * f;

    /*std::cout << f;
    std::cout << "\n";
    std::cout << "\n";
    exit(0);*/
 }
