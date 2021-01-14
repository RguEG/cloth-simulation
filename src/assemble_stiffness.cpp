#include <assemble_stiffness.h>
#include <iostream>
typedef Eigen::Triplet<double> G;
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 

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
//  K- the 3n by 3n sparse stiffness matrix. 

    Eigen::Matrix99d H;
    K.resize(q.size(), q.size());
    K.setZero();
    H.setZero();
    std::vector<G> tripletList;
    tripletList.reserve(9 * V.rows() * V.rows());

    Eigen::Matrix<double, 1, 9> tmp_row;
    Eigen::Matrix3d tmp_matrix;

    for (int i = 0; i < F.rows(); i++) {
        tmp_row = dX.row(i);
        tmp_matrix = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());
        d2V_membrane_corotational_dq2(H,q, tmp_matrix,V,F.row(i),a0[i],mu,lambda);
        for (int m = 0; m < 3; m++) {
            for (int n = 0; n < 3; n++) {
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 0) + n, H(m, n)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 1) + n, H(m, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 2) + n, H(m, n + 6.0)));


                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 0) + n, H(m + 3.0, n)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 1) + n, H(m + 3.0, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 2) + n, H(m + 3.0, n + 6.0)));


                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 0) + n, H(m + 6.0, n)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 1) + n, H(m + 6.0, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 2) + n, H(m + 6.0, n + 6.0)));

                /*tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 0) + n, H(m, n)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 1) + n, H(m + 3.0, n)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 2) + n, H(m + 6.0, n)));


                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 0) + n, H(m, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 1) + n, H(m + 3.0, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 2) + n, H(m + 6.0, n + 3.0)));


                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 0) + n, H(m, n + 6.0)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 1) + n, H(m + 3.0, n + 6.0)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 2) + n, H(m + 6.0, n + 6.0)));*/

            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());

    K = K * (-1.0);

    /*std::cout << H;
    exit(0);*/
    }
