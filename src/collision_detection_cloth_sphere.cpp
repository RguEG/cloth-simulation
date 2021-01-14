#include <collision_detection_cloth_sphere.h>
#include <iostream>
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();
    //  q - generalized coordinates for the FEM system
//  center - the position of the sphere center in the world space
//  radius - the radius of the sphere in the world space 
//Output:
//  cloth_index - the indices of all vertices currently in contact with the sphere
//  normals - the outward facing contact normals for each contacting vertex. 
    Eigen::Vector3d tmp_q;
    double distance = 0.0;
    double powR = 0.0;


    for (int i = 0; i < q.size() / 3.0; i++) {
        tmp_q[0] = q[3.0 * i];
        tmp_q[1] = q[3.0 * i + 1.0];
        tmp_q[2] = q[3.0 * i + 2.0];

        distance = pow((tmp_q[0] - center[0]), 2.0) + pow((tmp_q[1] - center[1]), 2.0) + pow((tmp_q[2] - center[2]), 2.0);
        powR = pow(radius, 2.0);

        if (distance <= powR) {
            cloth_index.push_back(i);

            tmp_q[0] = center[0] - tmp_q[0];
            tmp_q[1] = center[1] - tmp_q[1];
            tmp_q[2] = center[2] - tmp_q[2];

            normals.push_back(tmp_q);
        }
    }
}