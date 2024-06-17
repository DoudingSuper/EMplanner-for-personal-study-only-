#include "em_planner.h"
void EMplanner::init(const std::vector<waypoint> &global_path, const std::vector<obstacle> &global_obstacles, double T = 0.1) {
    _T = T;
    _path_smoother.init(&global_path);
    // 这里直接将ref_line的地址传给projection进行初始化，但注意此时ref_line是一个空地址
    // _projction.init(&_ref_line);
    _global_obstacles = global_obstacles;
}

bool EMplanner::update_location_info(const location &location) {
    _host_lacation = location;
    _current_t = location.t;
}

bool EMplanner::update_ref_line() {
    _path_smoother.get_ref_line(_host_lacation);
    _path_smoother.smooth_path();
    _path_smoother.get_final_ref_line(_ref_line);
    _projction.init(&_ref_line);
}

bool EMplanner::update_start_point() {
    // 首先将start_point清空
    _start_point.clear();
    _start_point.x = _host_lacation.x;
    _start_point.y = _host_lacation.y;
    _start_point.dir = _host_lacation.dir;
    _start_point.t = _host_lacation.t;
    // // 第一次进行规划
    // if (_is_first_run) {
    //     // 默认初始车辆的速度和加速度为0，故不需要进行外推
    //     _start_point.x = _host_lacation.x;
    //     _start_point.y = _host_lacation.y;
    //     _start_point.dir = _host_lacation.dir;
    //     _start_point.t = _host_lacation.t + _T;
    //     // _is_first_run = false;
    //     return false;
    // }
    // // 若不是第一次规划，也就是说上个周期已经规划出了一条轨迹
    // // 由于定位信息是车辆坐标系下的，需要将其转化到世界坐标系下
    // double vx_gcs = _host_lacation.vx * cos(_host_lacation.dir) - _host_lacation.vy * sin(_host_lacation.dir);
    // double vy_gcs = _host_lacation.vx * sin(_host_lacation.dir) + _host_lacation.vy * cos(_host_lacation.dir);
    // double ax_gcs = _host_lacation.ax * cos(_host_lacation.dir) - _host_lacation.ay * sin(_host_lacation.dir);
    // double ay_gcs = _host_lacation.ax * sin(_host_lacation.dir) + _host_lacation.ay * cos(_host_lacation.dir);
    // double kappa = 0; // 定位信息无法得到kappa，默认为0
    // // 根据时间戳找到对应时刻的参考点
    // waypoint t_match_point;
    // for (int i = 0; i < _pre_trajectory.size() - 1; i++) {
    //     if(_current_t >= _pre_trajectory[i].t && _current_t < _pre_trajectory[i + 1].t) {
    //         t_match_point = _pre_trajectory[i];
    //         Eigen::Vector2f vec = {_host_lacation.x - t_match_point.x,
    //                                _host_lacation.y - t_match_point.y};
    //         // 参考点切向向量
    //         Eigen::Vector2f tor = {cos(t_match_point.dir), sin(t_match_point.dir)};
    //         // 参考点垂向向量
    //         Eigen::Vector2f nor = {-sin(t_match_point.dir), cos(t_match_point.dir)};
    //         // 横向误差
    //         double lat_error = abs(vec.dot(tor));
    //         // 纵向误差
    //         double lon_error = abs(vec.dot(nor));
    //         // 如果控制算法无法跟上，则使用车辆动力学外推作为起点
    //         if (lat_error > 0.5 || lon_error > 2.5) {
    //             Eigen::Vector2d dir = {cos(_host_lacation.dir), sin(_host_lacation.dir)}; // 车辆方向向量、
    //             _start_point.x = vx_gcs * _T + 0.5 * ax_gcs * _T * _T;
    //             _start_point.y = vy_gcs * _T + 0.5 * ay_gcs * _T * _T;
    //             double vx_f = vx_gcs + ax_gcs * _T;
    //             double vy_f = vy_gcs + ay_gcs * _T;
    //             // 车辆绝对速度方向即为轨迹点速度切向方向
    //             _start_point.velocity = sqrt(vx_f * vx_f + vy_f * vy_f);
    //             _start_point.dir = atan2(vy_f, vx_f);
    //             _start_point.acc = sqrt(ax_gcs * ax_gcs + ay_gcs * ay_gcs);
    //             _start_point.k = kappa; // 这一步可以不写，k本就是0
    //             _start_point.t = _host_lacation.t + _T;
    //             return false;
    //         }
    //         // 控制算法可以跟上
    //         else {
    //             for (int j = i; j < _pre_trajectory.size() - 1; j++) {
    //                 if(_current_t + _T >= _pre_trajectory[i].t && _current_t + _T < _pre_trajectory[i + 1].t) {
    //                     _start_point = _pre_trajectory[j];
    //                     _start_point.t = _pre_trajectory[j].t + _T;
    //                     _start_index = j;
    //                     return true;
    //                 }
    //             }
    //         }            
    //     }
    // }
    return false;
}
bool EMplanner::get_splice_points(std::vector<waypoint> &splice_points) {
    // 注意拼接轨迹中不能包括规划起点
    if (_start_index >= 20) {
        for (int i = _start_index - 20; i < _start_index; i++) {
            if (i < _pre_trajectory.size())
                splice_points.push_back(_pre_trajectory[i]);
        } 
    }
    else {
        for (int i = 0; i < _start_index; i++) {
            if (i < _pre_trajectory.size())
                splice_points.push_back(_pre_trajectory[i]);
        }
    }
    return true;

}

bool EMplanner::splice_trajectory() {
    std::vector<waypoint> splice_points;
    // 需要进行轨迹拼接
    if (get_splice_points(splice_points)) {
        std::vector<waypoint> temp = _final_trajectory;
        _final_trajectory.clear();
        for (const auto &p : splice_points) {
            _final_trajectory.push_back(p);
        }
        for (const auto &p : temp) {
            _final_trajectory.push_back(p);
        }
    }
}

bool EMplanner::convert_start_point() {
    // 计算车辆坐标的匹配点和投影点
    _projction.update(_host_lacation);
    int origin_match_index; // 车辆坐标匹配点索引
    _projction.get_match_point_index(origin_match_index);
    waypoint origin_proj_point; // 车辆坐标在参考线上的投影点
    _projction.get_proj_point(origin_proj_point);
    _host_lacation.dir = origin_proj_point.dir; // 由于模拟的定位模块无法得到车辆的dir，用投影点的dir代替
    _converter.init_index_to_s(_ref_line, origin_match_index, origin_proj_point);
    // _converter.cartesian_to_frenet(_host_lacation, origin_proj_point);
    // 计算规划起点的匹配点和投影点
    _projction.update(_start_point);
    int start_match_index;
    _projction.get_match_point_index(start_match_index);
    waypoint start_proj_point;
    _projction.get_proj_point(start_proj_point);
    // 计算规划起点的s坐标
    _converter.cal_s_from_index_to_s(_ref_line, start_match_index, start_proj_point, _start_point.s);
    _converter.cartesian_to_frenet(_start_point, origin_proj_point);
    
}

bool EMplanner::update_obstacle_info() {
    // 首先清空obstacles
    _obstacles.clear();
    for (const auto &obs : _global_obstacles) {
        // 距离向量
        Eigen::Vector2f vec = {obs.x - _host_lacation.x, obs.y - _host_lacation.y};
        // 车辆方向向量
        Eigen::Vector2f tor = {cos(_host_lacation.dir), sin(_host_lacation.dir)};
        Eigen::Vector2f nor = {-sin(_host_lacation.dir), cos(_host_lacation.dir)};
        double lon_dist = vec.dot(tor); // 纵向距离
        double lat_dist = vec.dot(nor); // 横向距离
        if (abs(lon_dist) < 70 && abs(lat_dist) < 10) {
            _obstacles.push_back(obs);
        }
    }
    // std::cout << "障碍物数量：" << _obstacles.size() << std::endl;
    // std::cout << "总障碍物数量：" << _global_obstacles.size() << std::endl;
}

bool EMplanner::convert_obstacles_point() {
    for (auto &obs : _obstacles) {
        // 计算障碍物在ref_line上的匹配点和投影点
        _projction.update(obs);
        int match_index;
        _projction.get_match_point_index(match_index);
        waypoint proj_point;
        _projction.get_proj_point(proj_point);
        // 计算障碍物坐标点s坐标
        _converter.cal_s_from_index_to_s(_ref_line, match_index, proj_point, obs.s);
        // 进行坐标系变换
        _converter.cartesian_to_frenet(obs, proj_point);
    }
}

// 在自然坐标系中进行采样，用于动态规划
bool EMplanner::dp_path_sample(int s_num, int l_num, double delta_s, double delta_l) {
    _convex_space.init(_start_point, s_num, l_num, delta_s, delta_l);
    _convex_space.create_rough_path(_obstacles);

}

void EMplanner::convert_dp_path() {
    _dp_path.clear();
    _convex_space.get_rough_path(_dp_path);
    for (auto &p : _dp_path) {
        // 首先需要将p的投影点找到（根据p的自然坐标）
        waypoint proj_point;
        for (int i = 1; i < _ref_line.size(); i++) {
            if (_ref_line[i].s >= p.s) {
                double ds = p.s - _ref_line[i - 1].s;
                proj_point.x = _ref_line[i - 1].x + ds * cos(_ref_line[i - 1].dir);
                proj_point.y = _ref_line[i - 1].y + ds * sin(_ref_line[i - 1].dir);
                proj_point.k = _ref_line[i - 1].k;
                proj_point.dir = _ref_line[i - 1].dir; // 这里粗略一点表示方向角
                break;
            }
        }
        _converter.frenet_to_cartesian(p, proj_point);
    }
}

void EMplanner::get_qp_bound(std::vector<double> &ub, std::vector<double> &lb) {
    _convex_space.create_convex_space(_obstacles, ub, lb);
}

void EMplanner::create_qp_path(double w_ref, double w_center, double w_dl, double w_ddl, double w_dddl) {
    // sample_on_frenet();
    _qp_path.clear();
    convert_dp_path();
    int n = _dp_path.size(); // 优化目标数量
    Eigen::MatrixXd l(2, 6), A = Eigen::MatrixXd::Zero(2 * (n - 1) + n + 3, 3 * n);
    Eigen::MatrixXd a(3, 1), m_ref = Eigen::MatrixXd::Zero(3 * n, n);
    a << 1, 0, 0;
    Eigen::MatrixXd m_center = m_ref;
    Eigen::MatrixXd b(3, 1), m_dl = Eigen::MatrixXd::Zero(3 * n, n);
    b << 0, 1, 0;
    Eigen::MatrixXd c(3, 1), m_ddl = Eigen::MatrixXd::Zero(3 * n, n);
    c << 0, 0, 1;
    Eigen::MatrixXd d(3, 1), m_dddl = Eigen::MatrixXd::Zero(3 * n, n);
    
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(3 * n);
    std::vector<double> road_ub, road_lb;
    get_qp_bound(road_ub, road_lb);
    Eigen::VectorXd ub = Eigen::VectorXd::Zero(2 * (n - 1) + n + 3);
    Eigen::VectorXd lb = Eigen::VectorXd::Zero(2 * (n - 1) + n + 3);
    double ds = 1;
    l << 1, ds, ds * ds / 3, - 1, 0 , ds * ds / 6,
         0, 1, ds / 2, 0, - 1, ds / 2;
    for (int i = 0; i < n; i++) {
        if (i < n - 2)
            A.block(2 * i, 3 * i, 2, 6) = l;
        A(2 * (n - 1) + i, 3 * i) = 1;
        ub[2 * (n - 1) + i] = road_ub[i];
        lb[2 * (n - 1) + i] = road_lb[i];
        m_ref.block(3 * i, i, 3, 1) = a;
        m_dl.block(3 * i, i, 3, 1) = b;
        m_ddl.block(3 * i, i, 3, 1) = c;
        if (i < n - 1)
            m_dddl.block(3 * (i + 1), i, 3, 1) = a;
        // gradient[3 * i] = -2 * _dp_path[i].l;
        gradient[3 * i] = - w_center * (road_ub[i] + road_lb[i]) - w_ref * 3 ;
    }
    ub[2 * (n - 1)] = std::max(road_ub[0], _host_lacation.l);
    lb[2 * (n - 1)] = std::min(road_lb[0], _host_lacation.l);
    A.block(2 * (n - 1) + n, 0, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
    // 起点约束
    ub[2 * (n - 1) + n] = _host_lacation.l;
    ub[2 * (n - 1) + n + 1] = _host_lacation.d_l;
    ub[2 * (n - 1) + n + 2] = _host_lacation.dd_l;
    lb[2 * (n - 1) + n] = _host_lacation.l;
    lb[2 * (n - 1) + n + 1] = _host_lacation.d_l;
    lb[2 * (n - 1) + n + 2] = _host_lacation.dd_l;
    if (_last_n != n) {
        _LinearMatrix = A.sparseView();
        _Hessian = (w_ref * m_ref * m_ref.transpose() + w_center * m_center * m_center.transpose() +
                                           w_dl * m_dl * m_dl.transpose() + w_ddl * m_ddl * m_ddl.transpose() +
                                           w_dddl * (2 * m_ddl * m_ddl.transpose() - 2 * m_ddl * m_dddl.transpose())).sparseView();
        _last_n = n;
    }
    // for (int i = 0; i < n; i++)
    //     std::cout << road_lb[i] << "  " << road_ub[i] << std::endl;
    Osqp_Solver my_solver;
    my_solver.init(3 * n, 2 * (n - 1) + n + 3, _Hessian, gradient, _LinearMatrix, lb, ub);
    Eigen::VectorXd result;
    my_solver.get_result(result);
    for (int i = 0; i < _dp_path.size(); i++) {
        _dp_path[i].l = result[3 * i];
        _qp_path.push_back(_dp_path[i]);
    }
}

void EMplanner::convert_qp_path(){
    for (auto &p : _qp_path) {
        // 首先需要将p的投影点找到（根据p的自然坐标）
        waypoint proj_point;
        for (int i = 1; i < _ref_line.size(); i++) {
            if (_ref_line[i].s >= p.s) {
                double ds = p.s - _ref_line[i - 1].s;
                proj_point.x = _ref_line[i - 1].x + ds * cos(_ref_line[i - 1].dir);
                proj_point.y = _ref_line[i - 1].y + ds * sin(_ref_line[i - 1].dir);
                proj_point.k = _ref_line[i - 1].k;
                proj_point.dir = _ref_line[i - 1].dir; // 这里粗略一点表示方向角
                break;
            }
        }
        _converter.frenet_to_cartesian(p, proj_point);
    }
}

bool EMplanner::dp_speed_sample(int s_num, int t_num, double delta_s, double delta_t) {
    // 筛选是否有动态障碍物，只有在有动态障碍物的时候才有进行速度动态规划的必要，
    // 否则可以直接进行速度二次规划
    _moving_obstacles.clear();
    for (const auto &obs : _obstacles) {
        if (obs.velocity != 0) {
            _moving_obstacles.push_back(obs);
        }
    }
    if (_moving_obstacles.size() == 0) {
        return false;
    }
    
}

void EMplanner::create_qp_speed(double w_ref, double w_center, double w_dds, double w_ddds) {
}

void EMplanner::create_trajectory() {
    // 
    // 对轨迹进行增密 
    // Path my_path(_qp_path);
    // my_path.dense_path(0.1);
    _final_trajectory.clear();
    _qp_path.front().t = _start_point.t;
    _final_trajectory.push_back(_qp_path.front());
    for (int i = 1; i < _qp_path.size(); i++) {
        _qp_path[i].t = _qp_path[i - 1].t + 0.1;
        _final_trajectory.push_back(_qp_path[i]);
    }
}
void EMplanner::get_final_trajectory(std::vector<waypoint> &trajectory) {
    trajectory.clear();
    trajectory = _final_trajectory;
}

void EMplanner::get_ref_line(std::vector<waypoint> &ref_line) {
    _path_smoother.get_final_ref_line(ref_line);
}

void EMplanner::get_qp_path(std::vector<waypoint> &qp_path) {
    qp_path = _qp_path;
}

void EMplanner::get_dp_path(std::vector<waypoint> &dp_path) {
    dp_path = _dp_path;
}


void EMplanner::run() {
        update_ref_line();
        update_start_point();
        convert_start_point();
        bool need_splice = convert_start_point();
        update_obstacle_info();
        convert_obstacles_point();
        dp_path_sample();
        create_qp_path();
        convert_qp_path();
        dp_speed_sample();
        create_qp_speed();
        create_trajectory();
        if (!_is_first_run) {
            splice_trajectory();
        }
        _pre_trajectory = _final_trajectory;
        _is_first_run = false;
}

