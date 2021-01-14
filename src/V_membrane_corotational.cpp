#include <V_membrane_corotational.h>
#include <iostream>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    //  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  energy- the per-triangle potential energy (the linear model described in the README).

    Eigen::MatrixXd F;
    Eigen::VectorXd n, N;
    Eigen::Matrix34d newX;
    Eigen::Matrix43d newPh;
    Eigen::Vector3d a, b, c, crossA, crossB;
    Eigen::Matrix3d U;
    Eigen::Vector3d S;
    Eigen::Matrix3d W;
    double psi, absn;

    a[0] = q[3.0 * element[0]];
    a[1] = q[3.0 * element[0] + 1.0];
    a[2] = q[3.0 * element[0] + 2.0];
    b[0] = q[3.0 * element[1]];
    b[1] = q[3.0 * element[1] + 1.0];
    b[2] = q[3.0 * element[1] + 2.0];
    c[0] = q[3.0 * element[2]];
    c[1] = q[3.0 * element[2] + 1.0];
    c[2] = q[3.0 * element[2] + 2.0];

    crossA = V.row(element[1]) - V.row(element[0]);
    crossB = V.row(element[2]) - V.row(element[0]);
    N = (crossA).cross(crossB);
    absn = pow(N[0], 2.0) + pow(N[1], 2.0) + pow(N[2], 2.0);
    N = (crossA).cross(crossB) * (1.0 / sqrt(absn));

    crossA = b - a;
    crossB = c - a;
    n = (crossA).cross(crossB);
    absn = pow(n[0], 2.0) + pow(n[1], 2.0) + pow(n[2], 2.0);
    n = (crossA).cross(crossB) * (1.0 / sqrt(absn));


    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < newX.rows(); n++) {
            newX(n, m) = q(3.0 * element[m] + n);
        }
    }

    for (int m = 0; m < 3; m++) {
        newX(m, 3) = n[m];
    }

    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            newPh(m, n) = dX(m, n);
        }
    }

    for (int m = 0; m < 3; m++) {
        newPh(3, m) = N[m];
    }

    F = newX * newPh;

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    S = svd.singularValues();

    psi = mu * (pow((S[0] - 1), 2) + pow((S[1] - 1), 2) + pow((S[2] - 1), 2) + lambda / 2.0 * pow(S[0] + S[1] + S[2] - 3.0 , 2));
    
    energy = area * psi;
}
