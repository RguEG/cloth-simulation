#include <mass_matrix_mesh.h>
#include <iostream>

typedef Eigen::Triplet<double> G;

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    //Input: 
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions
//  F - the mx3 matrix of triangle-vertex indices
//  density - the density of the cloth material
//  areas - the mx1 vector of undeformed triangle areas
//Output:
//  M - sparse mass matrix for the entire mesh

    Eigen::MatrixXd tmp_M;
    M.resize(q.size(), q.size());
    M.setZero();
    tmp_M.resize(9, 9);
    std::vector<G> tripletList;
    tripletList.reserve(q.size() * q.size());

    for (int i = 0; i < F.rows(); i++) {
        tmp_M.setZero();

        tmp_M(0, 0) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(1, 1) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(2, 2) = 2.0 * density * areas[i] * (1.0 / 12.0);

        tmp_M(0, 3) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(1, 4) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(2, 5) = 2.0 * density * areas[i] * (1.0 / 24.0);

        tmp_M(0, 6) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(1, 7) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(2, 8) = 2.0 * density * areas[i] * (1.0 / 24.0);


        tmp_M(3, 0) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(4, 1) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(5, 2) = 2.0 * density * areas[i] * (1.0 / 24.0);

        tmp_M(3, 3) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(4, 4) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(5, 5) = 2.0 * density * areas[i] * (1.0 / 12.0);

        tmp_M(3, 6) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(4, 7) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(5, 8) = 2.0 * density * areas[i] * (1.0 / 24.0);


        tmp_M(6, 0) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(7, 1) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(8, 2) = 2.0 * density * areas[i] * (1.0 / 24.0);

        tmp_M(6, 3) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(7, 4) = 2.0 * density * areas[i] * (1.0 / 24.0);
        tmp_M(8, 5) = 2.0 * density * areas[i] * (1.0 / 24.0);

        tmp_M(6, 6) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(7, 7) = 2.0 * density * areas[i] * (1.0 / 12.0);
        tmp_M(8, 8) = 2.0 * density * areas[i] * (1.0 / 12.0);
        for (int m = 0; m < 3; m++) {
            for (int n = 0; n < 3; n++) {
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 0) + n, tmp_M(m, n)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 1) + n, tmp_M(m, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 0) + m, 3 * F(i, 2) + n, tmp_M(m, n + 6.0)));

                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 0) + n, tmp_M(m + 3.0, n)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 1) + n, tmp_M(m + 3.0, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 1) + m, 3 * F(i, 2) + n, tmp_M(m + 3.0, n + 6.0)));

                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 0) + n, tmp_M(m + 6.0, n)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 1) + n, tmp_M(m + 6.0, n + 3.0)));
                tripletList.push_back(G(3 * F(i, 2) + m, 3 * F(i, 2) + n, tmp_M(m + 6.0, n + 6.0)));

            }
        }
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}
 