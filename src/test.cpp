#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

int main() {
    Eigen::MatrixXd A(6, 2); //
    A << 1, 0, 
         0, 1,
        -2, 0,
         0, -2, 
         1, 0, 
         0, 1;
    Eigen::MatrixXd smooth_mat = Eigen::MatrixXd::Zero(360, (180 - 2) * 2);
    for (int i = 0; i < 180 - 2; ++i) {
        smooth_mat.block(i * 2, i * 2, 6, 2) = A;
    }
    Eigen::MatrixXd similar_mat = Eigen::MatrixXd::Identity(360, 360);
    Eigen::MatrixXd B(4, 2); //
    B << -1, 0,
         0, -1,
         1, 0,
         0, 1;
    Eigen::MatrixXd compact_mat = Eigen::MatrixXd::Zero(360, (180 - 1) * 2);
    for (int i = 0; i < 180 - 1; ++i) {
        compact_mat.block(i * 2, i * 2, 4, 2) = B;
    }
    // 打印结果
    Eigen::SparseMatrix<double> hessian = (smooth_mat * smooth_mat.transpose() + similar_mat + compact_mat * compact_mat.transpose()).sparseView();
    std::cout << "Result matrix:\n" << smooth_mat << std::endl;
    return 0;
}