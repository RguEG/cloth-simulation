#include <d2V_membrane_corotational_dq2.h>
#include <iostream>

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    
    //  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  H - the per-triangle Hessian of the potential energy (the linear model described in the README).

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here

    Eigen::MatrixXd d2psi;

    Eigen::MatrixXd dF, nB, nNewN;
    Eigen::VectorXd n, N;
    Eigen::Matrix34d newX;
    Eigen::Matrix43d newPh;
    Eigen::Vector3d a, b, c, crossA, crossB;
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
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();
    
    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    //TODO: compute H, the hessian of the corotational energy

    Eigen::Tensor3333d dU;
    Eigen::Tensor333d dS;
    Eigen::Tensor3333d dV;
    Eigen::MatrixXd deltaS,deltaSij,d2psids, dpds, dpsi, flatM, newn, flatMTrans;
    Eigen::VectorXd flatV;
    flatM.resize(3, 3);
    flatMTrans.resize(3, 3);
    deltaS.resize(3, 3);
    deltaSij.resize(3, 1);
    deltaS.setZero();
    d2psids.resize(3, 3);
    d2psids.setZero();
    d2psi.resize(9, 9);
    d2psi.setZero();
    dsvd(dU,dS,dV,F);

    deltaS(0, 0) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[0] * 2.0 - 2.0);
    deltaS(1, 1) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[1] * 2.0 - 2.0);
    deltaS(2, 2) = (lambda * (S[0] * 2.0 + S[1] * 2.0 + S[2] * 2.0 - 6.0)) / 2.0 + mu * (S[2] * 2.0 - 2.0);

    d2psids(0, 0) = lambda + 2 * mu;
    d2psids(0, 1) = lambda;
    d2psids(0, 2) = lambda;
    d2psids(1, 0) = lambda;
    d2psids(1, 1) = lambda + 2 * mu;
    d2psids(1, 2) = lambda;
    d2psids(2, 0) = lambda;
    d2psids(2, 1) = lambda;
    d2psids(2, 2) = lambda + 2 * mu;

    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            deltaSij = d2psids * dS[m][n];
            flatM = dU[m][n] * deltaS * W.transpose() + U * deltaSij.asDiagonal() * W.transpose() + U * deltaS * dV[m][n].transpose();
            flatMTrans = flatM.transpose();
            flatV = Eigen::Map<const Eigen::VectorXd>(flatMTrans.data(), flatMTrans.size());
            d2psi.row(3.0 * m + n) = flatV.transpose();
            /*for (int x = 0; x < 3; x++) {
                for (int y = 0; y < 3; y++) {
                    d2psi(x + 3.0 * m, y + 3.0 * n) = flatM(x, y);
                }
            }*/
        }
    }

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

    Eigen::MatrixXd newI, deltaX1, deltaX2, IdeltaX1, IdeltaX2;
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
    dF = nB + nNewN * newn;

    H = area * dF.transpose() * d2psi * dF;

    /*std::cout << "---------F---------";
    std::cout << "\n";
    std::cout << F;
    std::cout << "\n";
    std::cout << "---------newX---------";
    std::cout << "\n";
    std::cout << newX;
    std::cout << "\n";
    std::cout << "---------newPhi---------";
    std::cout << "\n";
    std::cout << newPh;
    std::cout << "\n";*/


    //-----------------------------------------------
    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
   

    /*std::cout << H;
    exit(0);*/
}
