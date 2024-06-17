#include "path_smoother.h"
void Path_Smoother::init(const std::vector<waypoint>* const global_path){
    this->_global_path = global_path;
    _projction.init(global_path);
    _ref_path.reserve(100);
    _smooth_ref_path.reserve(100);
}
void Path_Smoother::update_match_index(const location &location){
    _projction.update(location);
    if (_is_first_smooth) {
        _projction.get_match_point_index(_match_point_index);
    }
    else {
        _projction.get_match_point_index(_pre_match_index, _match_point_index);
    }
    // std::cout << "匹配点索引：" << _match_point_index << std::endl;
    _pre_match_index = _match_point_index;
    _is_first_smooth = false;
    _ref_path.clear();
}

void Path_Smoother::get_ref_line(const location &location){
    update_match_index(location);
    int len = this->_global_path->size();
    int start_index;
    //判断能否在匹配点前后得到100个点作为参考线
    if (_match_point_index > 30 && _match_point_index < len -70) {
        start_index = _match_point_index - 30;
    }    
    else if (_match_point_index <= 30) { 
        start_index = 0;
    }
    else {
        start_index = len - 100;
    }    
    for (int i = 0; i < 100; ++i) {
        this->_ref_path.push_back((*_global_path)[start_index + i]);
    }
    if (_pre_start_index == -1) {
        _smooth_points_num = 100;
    }
    else {
        _smooth_points_num = abs(start_index - _pre_start_index) > 100 ? 100 : start_index - _pre_start_index;
    }
    _pre_start_index = start_index;
}

void Path_Smoother::smooth_path() {
    int num = abs(_smooth_points_num); // 优化路径点数
    if (num == 0) {
        return;
    }
    Eigen::MatrixXd A(6, 2);
    A << 1, 0, 
         0, 1,
        -2, 0,
         0, -2, 
         1, 0, 
         0, 1;
    Eigen::MatrixXd smooth_mat = Eigen::MatrixXd::Zero(2 * num, (num - 2) * 2);
    for (int i = 0; i < num - 2; ++i) {
        smooth_mat.block(i * 2, i * 2, 6, 2) = A;
    }
    Eigen::MatrixXd similar_mat = Eigen::MatrixXd::Identity(2 * num, 2 * num);
    Eigen::MatrixXd B(4, 2); //
    B <<-1, 0,
         0, -1,
         1, 0,
         0, 1;
    Eigen::MatrixXd compact_mat = Eigen::MatrixXd::Zero(2 * num, (num - 1) * 2);
    for (int i = 0; i < num - 1; ++i) {
        compact_mat.block(i * 2, i * 2, 4, 2) = B;
    }
    spr_Matrix hessian = (w_smoothing_cost* smooth_mat * smooth_mat.transpose() + w_similar_cost * similar_mat + w_compact_cost * compact_mat * compact_mat.transpose()).sparseView();
    vxd grad(2 * num);
    if (_smooth_points_num > 0) {
        for (int i = 0; i < num; i ++) {
            grad[2 * i] = -_ref_path[(_ref_path.size() - num + i)].x;
            grad[2 * i + 1] = -_ref_path[(_ref_path.size() - num + i)].y;
        }
    }
    else {
        for(int i = 0; i < num; i ++) {
            grad[2 * i] = -_ref_path[i].x;
            grad[2 * i + 1] = -_ref_path[i].y;
        }
    }
    grad = w_similar_cost * grad; // 线性项是从相似代价中得到的
    this->qp_solver.init(2 * num, 0, hessian, grad);
    vxd result;
    this->qp_solver.get_result(result);
    // std::cout << result.size() << "************" << std::endl;
    std::vector<waypoint> old_path = _smooth_ref_path;
    _smooth_ref_path.clear();
    if (_smooth_points_num > 0) {
        for (int i = num; i < old_path.size(); i ++) {
            _smooth_ref_path.push_back(old_path[i]);
        }
        for (int i = 0; i < num; i++) {
            waypoint p;
            p.x = result[2 * i];
            p.y = result[2 * i + 1];
            _smooth_ref_path.push_back(p);
        }
    }
    else {
        for (int i = 0; i < num; i++) {
            waypoint p;
            p.x = result[2 * i];
            p.y = result[2 * i + 1];
            _smooth_ref_path.push_back(p);
        }
        for (int i = 0; i < old_path.size() - num; i ++) {
            _smooth_ref_path.push_back(old_path[i]);
        }
    }

}
bool Path_Smoother::get_final_ref_line(std::vector<waypoint>& final_ref_line) {
    Path my_path(_smooth_ref_path);
    my_path.get_dir_kappa();
    final_ref_line = _smooth_ref_path;
    // for(int i = 0; i < this->_ref_path.size(); i++) {
    //     final_ref_line[i] = this->_ref_path[i];
    // }
    return true;
}