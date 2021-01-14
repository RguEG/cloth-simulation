#include <dphi_cloth_triangle_dX.h>
#include <iostream>

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    //Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x3 vertex indices for this tetrahedron
//  X - the 3D position in the underformed space at which to compute the gradient
//Output:
//  dphi - the 3x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX

    Eigen::MatrixXd T, invT, consN;
    Eigen::RowVector2d cons;
    Eigen::RowVector3d consM;
    cons.setZero();
    consM.setZero();
    consN.setZero();
    invT.setZero();
    T.setZero();
    T.resize(3, 2);
    consN.resize(2, 3);

    cons[0] = 1;
    cons[1] = 1;

    for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 3; n++) {
            T(n, m) = V(element[m + 1.0], n) - V(element[0], n);
        }
    }

    consM = cons * (T.transpose()*T).inverse() * T.transpose();
    consN = (T.transpose() * T).inverse() * T.transpose();

    for (int n = 0; n < 3; n++) {
        dphi(0, n) = -consM[n];
    }

    for (int m = 1; m < dphi.rows(); m++) {
        for (int n = 0; n < dphi.cols(); n++) {
            dphi(m, n) = consN(m - 1.0, n);
        }
    }
}