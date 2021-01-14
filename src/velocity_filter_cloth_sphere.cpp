#include <velocity_filter_cloth_sphere.h>
#include <iostream>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {

    //Input:
//  qdot - the 3nx1 generalized velocities of the cloth mesh
//  index - a list of collision vertex indices from the collision detector
//  normals - a list of collision normals from the collision detector
//Output:
//  qdot- the filtered 3nx1 generalized velocities

    double alpha = 0.0;
    Eigen::Vector3d tmp_qdot;

    for (int i = 0; i < indices.size();i++) {

        tmp_qdot[0] = qdot[3.0 * indices[i]];
        tmp_qdot[1] = qdot[3.0 * indices[i] + 1.0];
        tmp_qdot[2] = qdot[3.0 * indices[i] + 2.0];

        if ((-normals.at(i).transpose() * tmp_qdot) >= 0) {
            alpha = 0.0;
        }
        else { 
            alpha = -normals.at(i).transpose() * tmp_qdot; 
            tmp_qdot = tmp_qdot + alpha * normals.at(i);

            qdot[3.0 * indices[i]] -= tmp_qdot[0];
            qdot[3.0 * indices[i] + 1.0] -= tmp_qdot[1];
            qdot[3.0 * indices[i] + 2.0] -= tmp_qdot[2];
        }

    }
}