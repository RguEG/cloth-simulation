#include <dV_membrane_corotational_dq.h>
#include <iostream>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).  
    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
    Eigen::MatrixXd F,dF, newn, nNewN, nB;
    Eigen::VectorXd n, N;
    Eigen::Matrix34d newX;
    Eigen::Matrix43d newPh;
    Eigen::RowVector3d a, b, c, crossA, crossB;
    Eigen::MatrixXd dpsi, dpds, dpsiTrans;
    double psi,absn;

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
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient 
    dpds.resize(3, 3);
    dpds.setZero();
    dpds(0, 0) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[0] * 2.0 - 2.0);
    dpds(1, 1) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[1] * 2.0 - 2.0);
    dpds(2, 2) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[2] * 2.0 - 2.0);

    dpsi = U * dpds * W.transpose();

    Eigen::SparseMatrixd B, newN;
    B.resize(9, 9);
    B.setZero();
    newN.resize(9, 3);
    newN.setZero();

    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            B.insert(n, 3.0 * m) = dX(m, n);
            B.insert(n + 3.0, 3.0 * m + 1.0) = dX(m, n);
            B.insert(n + 6.0, 3.0 * m + 2.0) = dX(m, n);
        }
    }
    for (int m = 0; m < 3; m++) {
        newN.insert(m, 0) = N[m];
        newN.insert(m + 3.0, 1) = N[m];
        newN.insert(m + 6.0, 2) = N[m];
    }

    Eigen::MatrixXd newI,deltaX1,deltaX2,IdeltaX1,IdeltaX2;
    newI.resize(3, 3);
    newI.setIdentity();
    
    deltaX1.resize(3, 3);
    deltaX1.setZero();
    deltaX2.resize(3, 3);
    deltaX2.setZero();

    deltaX1(0, 1) = -crossA[2];
    deltaX1(0, 2) = crossA[1];
    deltaX1(1, 0) = crossA[2];
    deltaX1(1, 2) = -crossA[0];
    deltaX1(2, 0) = -crossA[1];
    deltaX1(2, 1) = crossA[0];

    deltaX2(0, 1) = -crossB[2];
    deltaX2(0, 2) = crossB[1];
    deltaX2(1, 0) = crossB[2];
    deltaX2(1, 2) = -crossB[0];
    deltaX2(2, 0) = -crossB[1];
    deltaX2(2, 1) = crossB[0];

    IdeltaX1.resize(3, 9);
    IdeltaX1.setZero();
    IdeltaX2.resize(3, 9);
    IdeltaX2.setZero();

    IdeltaX1(0, 0) = -1;
    IdeltaX1(1, 1) = -1;
    IdeltaX1(2, 2) = -1;
    IdeltaX1(0, 6) = 1;
    IdeltaX1(1, 7) = 1;
    IdeltaX1(2, 8) = 1;

    IdeltaX2(0, 0) = -1;
    IdeltaX2(1, 1) = -1;
    IdeltaX2(2, 2) = -1;
    IdeltaX2(0, 3) = 1;
    IdeltaX2(1, 4) = 1;
    IdeltaX2(2, 5) = 1;
    
    newn = (1.0 / sqrt(absn)) * (newI - n * n.transpose()) * (deltaX1 * IdeltaX1 - deltaX2 * IdeltaX2);

    nB.resize(9, 9);
    nB = Eigen::MatrixXd(B);
    nNewN = Eigen::MatrixXd(newN);
    dpsiTrans = dpsi.transpose();

    dV = area * (nB + nNewN * newn).transpose() * Eigen::Map<const Eigen::VectorXd>(dpsiTrans.data(), dpsiTrans.size());

    /*std::cout << dpsi;
    exit(0);*/

    /*std::cout << (1.0 / sqrt(absn));
    std::cout << "\n";
    std::cout << "\n";
    std::cout << newI;
    std::cout << "\n";
    std::cout << "\n";
    std::cout << deltaX1;
    std::cout << "\n";
    std::cout << "\n";
    std::cout << IdeltaX1;
    std::cout << "\n";
    std::cout << "\n";
    std::cout << deltaX2;
    std::cout << "\n";
    std::cout << "\n";
    std::cout << IdeltaX2;
    std::cout << "\n";
    std::cout << "\n";
    exit(0);*/
}
