#include "frenetUtils.h"

void readObs(std::vector<obstacle>& obs_arry, const std::string& file_name){
    obstacle point;
    //读取表中的信息
    std::ifstream infile(file_name, std::ios::in);
    if(!infile.is_open()){
        std::cout<<"Can not open the file!!"<< std::endl;
        return;
    }
    std::string line;
    std::stringstream ss;
    while(getline(infile, line)){
        ss<<line;//字符串流输入
        ss>>point.x>>point.y>>point.velocity>>point.acc;//字符串流输出
        obs_arry.push_back(point);
        ss.clear();//将字符串流清空
    }
    infile.close();//结束对表格的占用

}

void readPath(std::vector<waypoint>& path_array, const std::string& file_name){
    waypoint point;
    //读取表中的信息
    std::ifstream infile(file_name, std::ios::in);
    if(!infile.is_open()){
        std::cout<<"Can not open the file!!"<< std::endl;
        return;
    }
    std::string line;
    std::stringstream ss;
    while(getline(infile, line)){
        ss<<line;//字符串流输入
        ss>>point.x>>point.y>>point.velocity>>point.acc;//字符串流输出
        path_array.push_back(point);
        ss.clear();//将字符串流清空
    }
    infile.close();//结束对表格的占用
}

bool Osqp_Solver::init(int n, int m, spr_Matrix hessian, vxd gradient, spr_Matrix linear_matrix, vxd lb, vxd ub) {
    _n = n;
    _m = m;
    _hessian = hessian;
    _gradient = gradient;
    _linear_matrix = linear_matrix;
    _lower_bound = lb;
    _upper_bound = ub;
}
bool Osqp_Solver::get_result(vxd& QP_solution) {
    OsqpEigen::Solver solver;
        // settings
    solver.settings()->setVerbosity(true);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver    
    solver.data()->setNumberOfVariables(_n);   //变量数n
    solver.data()->setNumberOfConstraints(_m); //约束数m
    if (!solver.data()->setHessianMatrix(this->_hessian))
        return false;
    if (!solver.data()->setGradient(this->_gradient))
        return false;
    if(_m > 0) {
        if (!solver.data()->setLinearConstraintsMatrix(this->_linear_matrix))
            return false;
        if (!solver.data()->setLowerBound(_lower_bound))
            return false;
        if (!solver.data()->setUpperBound(_upper_bound))
            return false;
    }
    // instantiate the solver
    if (!solver.initSolver())
        return false;
    
    // solve the QP problem
    if (!solver.solve())
    {
        return false;
    }
    // std::cout << "************";
    QP_solution = solver.getSolution();
        // std::cout << "QPSolution" << std::endl
        //       << QP_solution << std::endl; //输出为m*1的向量
    return true;

}

void Projection::init(const std::vector<waypoint>* const global_path) {
    _global_path = global_path;
}
bool Projection::update(const location& location) {
    _location = location;
    _increase_count = 0;
    _is_get_match_point = false;
    return true; 
}
bool Projection::update(const waypoint& waypoint) {
    _location.x = waypoint.x;
    _location.y = waypoint.y;
    _increase_count = 0;
    _is_get_match_point = false;
    return true; 
}
bool Projection::update(const obstacle& obs) {
    _location.x = obs.x;
    _location.y = obs.y;
    _increase_count = 0;
    _is_get_match_point = false;
    return true; 
}
bool Projection::cal_match_point_index() {
    _is_get_match_point = true;
    int start_index;
    int min_index = 0;
    start_index = 0;
    double min_dis = std::numeric_limits<double>::infinity();// 初始最小距离记为无穷大
    for (int i = start_index; i < (*_global_path).size(); i++) {
        double dis = std::sqrt((_location.x - (*_global_path)[i].x) * (_location.x - (*_global_path)[i].x) +
                                (_location.y - (*_global_path)[i].y) * (_location.y - (*_global_path)[i].y));
        if (dis < min_dis) {
            min_dis = dis;
            min_index = i;
            _increase_count = 0;
        }
        else {
            _increase_count++;
        }
        // 当距离连续增加十次，则认为已经找到距离最近的匹配点
        if (_increase_count > 300) {
            _match_point_index = min_index;
            return true;
        }
    }
    _match_point_index = min_index;
    return true; 
}
bool Projection::cal_match_point_index(const int &pre_match_index) {
    _is_get_match_point = true;
    int start_index;
    int min_index = 0;
    start_index = 0;
    double min_dis = std::numeric_limits<double>::infinity();// 初始最小距离记为无穷大
    start_index = pre_match_index;
    //首先判断前进方向
    Eigen::Vector2f d;// 车辆当前位置与上个匹配点之间的距离向量
    d << (_location.x - (*_global_path)[pre_match_index].x),
            (_location.y - (*_global_path)[pre_match_index].y);
    Eigen::Vector2f tor;
    tor << cos((*_global_path)[pre_match_index].dir / 180 * PI),
            sin((*_global_path)[pre_match_index].dir / 180 * PI);
    //距离向量与前匹配点方向向量同向，方向向前
    if (d.transpose() * tor > 0) {
        // std::cout << "正向前进" << std::endl;
        for (int i = start_index; i < (*_global_path).size(); i++) {
            double dis = std::sqrt((_location.x - (*_global_path)[i].x) * (_location.x - (*_global_path)[i].x) +
                                (_location.y - (*_global_path)[i].y) * (_location.y - (*_global_path)[i].y));
            if (dis < min_dis) {
                min_dis = dis;
                min_index = i;
                _increase_count = 0;
            }
            else {
                _increase_count++;
            }
            // 当距离连续增加十次，则认为已经找到距离最近的匹配点
            if (_increase_count > 300) {
                _match_point_index = min_index;
                return true;
            }
        }
        _match_point_index = min_index;
        return true; 
    }
    else {
        // std::cout << "反向前进" << std::endl;
        for (int i = start_index; i >= 0; i--) {
            double dis = std::sqrt((_location.x - (*_global_path)[i].x) * (_location.x - (*_global_path)[i].x) +
                                (_location.y - (*_global_path)[i].y) * (_location.y - (*_global_path)[i].y));
            if (dis < min_dis) {
                min_dis = dis;
                min_index = i;
                _increase_count = 0;
            }
            else {
                _increase_count++;
            }
            // 当距离连续增加十次，则认为已经找到距离最近的匹配点
            if (_increase_count > 300) {
                _match_point_index = min_index;
                return true;
            }
        }
        _match_point_index = min_index;

        return true; 
    }
}

bool Projection::cal_proj_point() {
    if(!_is_get_match_point){
        cal_match_point_index();
    }
    _match_point = (*_global_path)[_match_point_index];
    Eigen::Vector2f d, tor, ds, match_vec, proj_vec;
    d << _location.x - _match_point.x, _location.y - _match_point.y;
    tor << cos(_match_point.dir),
           sin(_match_point.dir);
    ds = d.transpose() * tor * tor;
    match_vec << _match_point.x, _match_point.y;
    proj_vec = match_vec + ds;
    _proj_point.x = proj_vec[0];
    _proj_point.y = proj_vec[1];
    _proj_point.dir = _match_point.dir + _match_point.k * double(d.transpose() * tor);
    _proj_point.k = _match_point.k;
}

bool Projection::get_match_point_index(int &index) {
    cal_match_point_index();
    index = _match_point_index;
}
bool Projection::get_match_point_index(const int &pre_match_index, int &index) {
    cal_match_point_index(pre_match_index);
    index = _match_point_index;
}

bool Projection::get_proj_point(waypoint &point) {
    cal_proj_point();
    point = _proj_point;
}

Path::Path(std::vector<waypoint> &global_path): _global_path(global_path) {
    
}

bool Path::dense_path(double sample_interval) {
    cal_length(_global_path, _path_len);
    std::vector<waypoint> origin_path(_global_path);
    _global_path.clear();
    _global_path.reserve(_path_len / sample_interval);
    _global_path.push_back(origin_path.front());
    
    for (int i = 1; i < origin_path.size(); i ++) {
        std::vector<waypoint> curve_points = {_global_path.back(),
                                              origin_path.at(i),
                                              i < origin_path.size() - 1 ? origin_path.at(i + 1) : origin_path.at(i)};
        cal_curve(curve_points);
        if (curve_points.at(1).k == 0) {
            resample_on_straight(curve_points, sample_interval);
            // std::cout << "**************" << std::endl;
        }
        else {
            resample_on_curve(curve_points, sample_interval);
            // std::cout << "**************" << std::endl;

        }
    }

}
void Path::cal_length(std::vector<waypoint> &path, double &length) {
    length = 0;
    for(int i = 0; i < path.size() - 1; i ++) {
        Eigen::Vector2f temp_points = {path.at(i).x - path.at(i + 1).x,
                                       path.at(i).y - path.at(i + 1).y};
        length += temp_points.norm();
    }
}
void Path::cal_curve(std::vector<waypoint> &curve_points) {
    if (curve_points.size() != 3) {
        return ;
    }
    waypoint &p0 = curve_points.at(0), &p1 = curve_points.at(1), &p2 = curve_points.at(2);
    // d可以看做p0p1的方向角与p0p2的方向角之差
    double d = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
    // 若角度差小于一定值，则将其视作直线
    if (d < 1e-8) {
        p0.k = 0;
        p0.dir = atan2(p1.y - p0.y, p1.x - p0.x);
        p1.k = 0;
        p1.dir  = p0.dir;
        return ;
    }
    double d_1_0 = p1.x * p1.x - p0.x * p0.x + p1.y * p1.y - p0.y * p0.y;
    double d_2_0 = p2.x * p2.x - p0.x * p0.x + p2.y * p2.y - p0.y * p0.y;
    // cx、cy是由三点确定的圆心坐标
    double cx = (d_1_0 * (p2.y - p0.y) - d_2_0 * (p1.y - p0.y)) / (2 * d);
    double cy = (d_1_0 * (p2.x - p0.x) - d_2_0 * (p1.x - p0.x)) / (-2 * d);
    double R = sqrt((cx - p1.x) * (cx - p1.x) + (cy - p1.y) * (cy - p1.y));
    p0.k = 1 / R;
    p1.k = 1 / R;
    
    // 计算p1方向角

    // 向量p1p2
    Eigen::Vector2f p1p2 = {p2.x - p1.x, p2.y - p1.y};
    // 向量p1c,c为圆心
    Eigen::Vector2f p1c = {cx - p1.x, cy - p1.y};
    // 判断p1c的切线方向
    if ((cy - p1.y) * (p2.x - p1.x) - (cx - p1.x) * (p2.y - p1.y) < 0) {
        p1.dir = atan2(cx - p1.x, -(cy - p1.y));
    }
    else {
        p1.dir = atan2(-(cx - p1.x), cy - p1.y);
    }

    // 计算p0方向角

    // 向量p0p1
    Eigen::Vector2f p0p1 = {p1.x - p0.x, p1.y - p0.y};
    // 向量p0c,c为圆心
    Eigen::Vector2f p0c = {cx - p0.x, cy - p0.y};
    // 判断p0c的切线方向
    if (p0c[1] * p0p1[0] - p0c[0] * p0p1[1] < 0) {
        p0.dir = atan2(p0c[0], -p0c[1]);
    }
    else {
        p0.dir = atan2(-p0c[0], p0c[1]);
    }
}   
void Path::resample_on_curve(const std::vector<waypoint> &points, const double &sample_interval) {
    waypoint p0 = points.at(0), p1 = points.at(1), p2 = points.at(2);
    double d_dir = fabs(p1.dir - p0.dir);
    double dis = d_dir/points.at(1).k;
    double diff_theta = sample_interval * points.at(1).k;
    double R = 1 / points.at(1).k;
    waypoint sample_point = points.at(0);
    int dir; // 用于表示顺时针还是逆时针，顺时针为-1，逆时针为1
    if ((p2.y - p0.y) * (p1.x - p0.x) - (p2.x - p0.x) * (p1.y - p0.y) < 0) {
        dir = -1;
    }
    else {
        dir = 1;
    }
    for (; dis > sample_interval; dis -= sample_interval) {
        Eigen::Vector2f start_vec = {sample_point.x, sample_point.y};
        if (_global_path.size() == _global_path.capacity()) {
            break;
        }
        Eigen::Vector2f vec = {cos(sample_point.dir + dir * diff_theta / 2), sin(points.at(0).dir + dir * diff_theta / 2)};
        Eigen::Vector2f sample_vec = start_vec + 2 * R * sin(diff_theta / 2) * vec;
        sample_point.x = sample_vec[0];
        sample_point.y = sample_vec[1];
        sample_point.dir = sample_point.dir + diff_theta;
        sample_point.k = points.at(1).k;
        _global_path.push_back(sample_point);
    }
}
void Path::resample_on_straight(const std::vector<waypoint> &points, const double &sample_interval) {
    Eigen::Vector2f p0p1 = {points.at(1).x - points.at(0).x, points.at(1).y - points.at(0).y};
    double dis = p0p1.norm();
    waypoint sample_point = points.at(0);
    Eigen::Vector2f vec = {cos(points.at(1).dir), sin(points.at(1).dir)};
    for (; dis > sample_interval; dis -= sample_interval) {
        if (_global_path.size() == _global_path.capacity()) {
            break;
        }
        Eigen::Vector2f start_vec = {sample_point.x, sample_point.y};
        Eigen::Vector2f sample_vec = start_vec + sample_interval * vec;
        sample_point.x = sample_vec[0];
        sample_point.y = sample_vec[1];
        _global_path.push_back(sample_point);
    }
}
/*原理：
输入：path的x,y坐标
输出：path的heading和kappa
heading = arctan（dy / dx），注意，计算方向角时要使用对方向敏感的四象限反正切函数atan2
kappa = d_heading / ds
ds = (dx ^ 2 + dy ^ 2) ^ 0.5
*/
void Path::get_dir_kappa() {
    double dy, dx, ds;
    for (int i = 0; i < _global_path.size(); i++) {
        if (i == 0) {
            dy = _global_path[i + 1].y - _global_path[i].y;
            dx = _global_path[i + 1].x - _global_path[i].x;
            _global_path[i].dir = atan2(dy, dx);
            ds = std::sqrt(dx * dx + dy * dy);
        }
        else if (i < _global_path.size() - 1) {
            double dy0 = _global_path[i].y - _global_path[i - 1].y;
            double dy1 = _global_path[i + 1].y - _global_path[i].y;
            double dx0 = _global_path[i].x - _global_path[i - 1].x;
            double dx1 = _global_path[i + 1].x - _global_path[i].x;
            _global_path[i].dir = (atan2(dy0, dx0) + atan2(dy1, dx1)) / 2;
            ds = (std::sqrt(dx1 * dx1 + dy1 * dy1) + std::sqrt(dx0 * dx0 + dy0 * dy0)) / 2;
        }
        else {
            dy = _global_path[i].y - _global_path[i - 1].y;
            dx = _global_path[i].x - _global_path[i - 1].x;
            _global_path[i].dir = _global_path[i - 1].dir;
            ds = std::sqrt(dx * dx + dy * dy);
        }
        if (i > 0) {
            _global_path[i].k = (_global_path[i].dir - _global_path[i - 1].dir) / ds;
        }
        else {
            _global_path[i].k = 0;
        }
    }
}

void CartesianFrenetConvert::init_index_to_s(std::vector<waypoint> &ref_path, int origin_match_index, waypoint origin_proj_point) {
    _ref_path = ref_path;
    double s = 0;
    for (int i = 1; i < ref_path.size(); i++) {
        waypoint &p1 = ref_path[i - 1], &p2 = ref_path[i];
        p2.s = p1.s + sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
    }
    double s0;
    cal_s_from_index_to_s(ref_path, origin_match_index, origin_proj_point, s0);
    for (int i = 0; i < ref_path.size(); i++) {
        ref_path[i].s -= s0;
    }
}

void CartesianFrenetConvert::cal_s_from_index_to_s(const std::vector<waypoint> &ref_path, int match_point_index, waypoint origin_proj_point, double &s) {
    const waypoint &p1 = ref_path[match_point_index], &p2 = origin_proj_point;
    // 这里直接使用直线距离代表投影点和匹配点之间的距离
    s = p1.s + sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

void CartesianFrenetConvert::cartesian_to_frenet(const double x, const double y, const double theta,
                                                 const double kappa, const double v, const double a,
                                                 const double rx, const double ry, const double rtheta,
                                                 const double rkappa, const double rdkappa, const double rs,
                                                 double &s, double &dot_s, double &ddot_s,
                                                 double &l, double &dot_l, double &ddot_l) {
    s = rs;
    double d_x = x - rx, d_y = y - ry, d_theta = theta - rtheta;
    double C_1 = d_y * cos(rtheta) - d_x * sin(rtheta);
    l = (C_1 > 0 ? 1 : -1) * sqrt(d_x * d_x + d_y * d_y);
    double C_2 = 1 - rkappa * l;
    dot_s = v * cos(d_theta) / C_2;
    dot_l = v * sin(d_theta);
    double d_l = dot_l / dot_s;
    double C_3 = rdkappa * l + rkappa * d_l;
    double C_4 = kappa * C_2 / cos(d_theta) - rkappa;
    ddot_s = (a * cos(d_theta) - dot_s * dot_s * (d_l * C_4 - C_3)) / C_1;
    ddot_l = - C_3 * tan(d_theta) + C_2 / (cos(d_theta) * cos(d_theta)) * C_4;
}
void CartesianFrenetConvert::cartesian_to_frenet(waypoint &point, waypoint &proj_point) {
    cartesian_to_frenet(point.x, point.y, point.dir, point.k, point.velocity, point.acc, 
                        proj_point.x, proj_point.y, proj_point.dir, proj_point.k, proj_point.dk,
                        point.s, point.s, point.dot_s, point.ddot_s, point.l, point.dot_l, point.ddot_l);
}
void CartesianFrenetConvert::cartesian_to_frenet(obstacle &point, waypoint &proj_point) {
    cartesian_to_frenet(point.x, point.y, point.dir, point.k, point.velocity, point.acc, 
                        proj_point.x, proj_point.y, proj_point.dir, proj_point.k, proj_point.dk,
                        point.s, point.s, point.dot_s, point.ddot_s, point.l, point.dot_l, point.ddot_l);
}
void CartesianFrenetConvert::frenet_to_cartesian (const double s, const double dot_s, const double ddot_s,
                                                const double l, const double dot_l, const double ddot_l,
                                                const double rx, const double ry, const double rtheta,
                                                const double rkappa, const double rdkappa,
                                                double &x, double &y, double &theta,
                                                double &kappa, double &v, double &a) {
    double C_1 = 1 - rkappa * l;
    x = rx - l * sin(rtheta);
    y = ry + l * cos(rtheta);
    double d_l = dot_l / dot_s;
    theta = atan(d_l / C_1) + rtheta;
    double C_2 = cos(theta - rtheta);
    double C_3 = rdkappa * l + rkappa * d_l;
    double C_4 = kappa * C_2 / cos(theta - rtheta) - rkappa;
    v = sqrt((dot_s * C_1) * (dot_s * C_1) + (dot_s * d_l) * (dot_s * d_l));
    a = ddot_s * C_1 / C_2 + dot_s * ddot_s / C_2 * (d_l * C_4 - C_3);
    double dd_l = ddot_l / ddot_s;
    kappa = ((dd_l + C_3 * tan(theta - rtheta)) * C_2 * C_2 / C_1 + rkappa) * C_4;
}

void CartesianFrenetConvert::frenet_to_cartesian(waypoint &point, waypoint &proj_point) {
    frenet_to_cartesian(point.s, point.dot_s, point.ddot_s, point.l, point.dot_l, point.ddot_l,
                        proj_point.x, proj_point.y, proj_point.dir, proj_point.k, proj_point.dk,
                        point.x, point.y, point.dir, point.k, point.velocity, point.acc);
}
void CartesianFrenetConvert::frenet_to_cartesian(obstacle &point, waypoint &proj_point) {
    frenet_to_cartesian(point.s, point.dot_s, point.ddot_s, point.l, point.dot_l, point.ddot_l,
                        proj_point.x, proj_point.y, proj_point.dir, proj_point.k, proj_point.dk,
                        point.x, point.y, point.dir, point.k, point.velocity, point.acc);
}

void CreateConvexSpace::init(waypoint &start_point, int s_num, int l_num, double delta_s, double delta_l, double offset) {
    _sps.clear();
    _start_point = start_point;
    _s_num = s_num;
    _l_num = l_num;
    _delta_s = delta_s;
    _delta_l = delta_l;
    _offset = offset;
    int sps_num = 1 + s_num * l_num;
    _edges.resize(sps_num, std::vector<long long>(sps_num, 0));
    _polynomial_sample_points.resize(sps_num, std::vector<std::vector<waypoint>>(sps_num, std::vector<waypoint>()));
    SampleWaypoint sp0;
    sp0.point = start_point;
    sp0.i = 0;
    sp0.j = -1; // 起点的列索引为-1
    _sps.push_back(sp0);
    for (int j = 0; j < s_num; j++) {
        for (int i = 0; i < l_num; i++) {
            SampleWaypoint sp;
            sp.i = i;
            sp.j = j;
            sp.point.s = (j + 1) * delta_s;
            sp.point.l = (l_num / 2 - i) * delta_l;
            _sps.push_back(sp);
        }
    }
}

void CreateConvexSpace::cal_polynomial(const SampleWaypoint &p1, const SampleWaypoint &p2, Eigen::VectorXd &a) {
    Eigen::MatrixXd A(6, 6);
    double s1 = p1.point.s;
    double s2 = s1 * s1;
    double s3 = s2 * s1;
    double s4 = s3 * s1;
    double s5 = s4 * s1;
    double s1_ = p2.point.s;
    double s2_ = s1_ * s1_;
    double s3_ = s2_ * s1_;
    double s4_ = s3_ * s1_;
    double s5_ = s4_ * s1_;
    A << 1, s1, s2, s3, s4, s5, 
         0, 1, 2 * s1, 3 * s2, 4 * s3, 5 * s4,
         0, 0, 2, 6 * s1, 12 * s2, 20 * s3,
         1, s1_, s2_, s3_, s4_, s5_, 
         0, 1, 2 * s1_, 3 * s2_, 4 * s3_, 5 * s4_,
         0, 0, 2, 6 * s1_, 12 * s2_, 20 * s3_;
    Eigen::VectorXd B(6);
    B << p1.point.l,
         p1.point.d_l,
         p1.point.dd_l,
         p2.point.l,
         p2.point.d_l,
         p2.point.dd_l;
    // 计算A的逆矩阵
    Eigen::MatrixXd A_inverse = A.inverse();
    // 得到五次多项式系数向量
    a = A.colPivHouseholderQr().solve(B);
}

void CreateConvexSpace::sample_on_polynomial(const int &i, const int &j, double delta_s) {
    _polynomial_sample_points[i][j].clear();
    Eigen::VectorXd a;
    cal_polynomial(_sps[i], _sps[j], a);
    waypoint p;
    for (double s = _sps[i].point.s + delta_s; s < _sps[j].point.s; s += delta_s) {
        double s1 = s;
        double s2 = s1 * s;
        double s3 = s2 * s;
        double s4 = s3 * s;
        double s5 = s4 * s;
        p.s = s;
        p.l = a[0] + a[1] * s1 + a[2] * s2 + a[3] * s3 + a[4] * s4 + a[5] * s5;
        p.d_l = a[1] + 2 * a[2] * s1 + 3 * a[3] * s2 + 4 * a[4] * s3 + 5 * a[5] * s4;
        p.dd_l = a[2] + 6 * a[3] * s1 + 12 * a[4] * s2 + 20 * a[5] * s3;
        _polynomial_sample_points[i][j].push_back(p);
    }

}

double CreateConvexSpace::cal_obs_dist_cost(const waypoint &p, const std::vector<obstacle> &obstacles) {
    double cost = 0;
    for (const auto &obs : obstacles) {
        double dist = sqrt((p.s - obs.s) * (p.s - obs.s) + (p.l - obs.l) * (p.l - obs.l));
        if (dist < 0.5) {
            cost += MY_MAX; // 距离过近，代价设为极大
        }
        else if (dist >= 0.5 && dist < 2) {
            cost += 1000 * (2 - dist);
        }
        // 若距离大于2m，则代价为0，不进行操作
    }
    return cost;
}


//  w1,  w2,  w3 为平滑代价的三个系数
void CreateConvexSpace::cal_cost(const std::vector<obstacle> &obstacles, double w1, double w2, double w3) {
    // 注意，最后一列节点没有后续节点，这里边界条件不要弄错
    for (int j = -1; j < _s_num - 1; j++) {
        for (int i = 0; i < _l_num; i++) {
            int p1, p2;
            if (j == -1) {
                p1 = 0;
            }
            else {
                p1 = j * _l_num + i + 1;
            }
            for (int k = 0; k < _l_num; k++) {
                p2 = (j + 1) * _l_num + k + 1;
                double smooth_cost = 0, obs_dist_cost = 0, ref_dist_cost = 0;
                // 在第p1和p2点之间的五次多项式上进行采样，并将其存入三维矩阵_polynomial_sample_points
                sample_on_polynomial(p1, p2);
                for (auto point : _polynomial_sample_points[p1][p2]) {
                    // 注意，这里的平滑代价只计算到l的二阶导，如有需要可以提升至三阶导
                    smooth_cost += 10 * point.d_l * point.d_l + 20 * point.dd_l * point.dd_l; 
                    obs_dist_cost += cal_obs_dist_cost(point, obstacles);
                    ref_dist_cost += (point.l - _offset) * (point.l - _offset);
                }
                // 将i和i+j+tem这两个点之间的代价存入二维邻接矩阵_edges
                _edges[p1][p2] = w1 * smooth_cost + w2 * obs_dist_cost + w3 * ref_dist_cost;
            }
            // -1列只有起点一个点
            if (j == -1) {
                break;
            }   
        }
    }
}

void CreateConvexSpace::create_rough_path(const std::vector<obstacle> &obstacles, double w1, double w2, double w3) {
    _rough_path.clear();
    cal_cost(obstacles, w1, w2, w3);
    _sps.at(0).cost = 0;
    for (int i = 0; i < _l_num; i++) {
        _sps.at(i + 1).cost = _sps.at(0).cost + _edges[0][i + 1];
        _sps.at(i + 1).parent = &_sps.at(0);
    }
    // 动态规划确定所有点的最小cost
    for (int j = 1; j < _s_num; j++) {
        for (int i = 0; i < _l_num; i++) {
            double min_cost = MY_MAX;
            int last_i = -1;
            for (int k = 0; k < _l_num; k++) {
                double temp = _sps.at((j - 1) * _l_num + k + 1).cost + _edges[(j - 1) * _l_num + k + 1][j * _l_num + i + 1];
                if (min_cost > temp) {
                    min_cost = temp;
                    last_i = k;
                }
            }
            if (last_i == -1) {
                // 若该点无法到达，则将其代价置为无穷大
                _sps.at(j * _l_num + i + 1).cost = MY_MAX;
                continue;
            }
            _sps.at(j * _l_num + i + 1).cost = min_cost;
            _sps.at(j * _l_num + i + 1).parent = &_sps.at((j - 1) * _l_num + last_i + 1);
        }
    }
    // 在最后一列中挑选cost最小的点作为终点
    double min_cost = MY_MAX;
    int end_i = 0;
    int end_j = _s_num - 1;
    for (int i = 0; i < _l_num; i++) {
        if (min_cost > _sps.at(end_j * _l_num + i + 1).cost) {
            min_cost = _sps.at(end_j * _l_num + i + 1).cost;
            end_i = i;
        }
    }
    SampleWaypoint *end_point = &_sps.at(end_j * _l_num + end_i + 1);
    while (end_point != nullptr) {
        _rough_path.push_back(end_point->point);
        if (end_point->parent == nullptr) {
            break;
        }
        int i = end_point->parent->j * _l_num + end_point->parent->i + 1;
        int j = end_point->j * _l_num + end_point->i + 1;
        if (i < 0) {
            i = 0;
        }
        for (auto it = _polynomial_sample_points[i][j].rbegin(); it != _polynomial_sample_points[i][j].rend(); it++) {
            _rough_path.push_back(*it);
        }
        end_point = end_point->parent; 
    }

    std::reverse(_rough_path.begin(), _rough_path.end());
}

// void CreateConvexSpace::create_rough_path(const std::vector<obstacle> &obstacles, double w1, double w2, double w3) {
//     cal_cost(obstacles, w1, w2, w3);
//     // 这里new一个新内存出来，和原始的地图空间做隔离，其实在本例中并没有必要，
//     // 但在全局地图中直接获取节点信息时，最好进行隔离，这样不会污染原始的地图信息
//     // 其实就是使用深拷贝而不要使用浅拷贝
//     _open_list.push_back(new SampleWaypoint (_sps[0]));
//     int count = 0; //用于计数最后一列已经完成规划的点数
//     // 接下来计算所有达到最后一列采样点的路径cost，得到cost集合，以此为依据选取路径终点
//     SampleWaypoint *cur_point; // 这算不算是个野指针？
//     SampleWaypoint *end_point = new SampleWaypoint; // 用于维护终点
//     end_point->cost = 100000000;
//     int total = 0;
//     while (!_open_list.empty()) { 
//         cur_point = get_least_cost_point();
//         _open_list.remove(cur_point);
//         _close_list.push_back(cur_point);
//         std::cout << cur_point->j << std::endl;
        
//         // 已经到达最后一列
//         if (cur_point->j == _s_num - 1) {
//             if (end_point->cost > cur_point->cost) {
//             std::cout << "已到达" << std::endl;
//                 end_point = cur_point;
//             }
//             count++;
//             if (count == _l_num) {
//                 break;
//             }
//             else {
//                 continue;
//             }
//         }
//         // 还未到达最后一列
//         else {
//             std::vector<SampleWaypoint*> neighbors;
//             get_neighbor_points(cur_point, neighbors);
//             for (auto p : neighbors) {
//                 if (is_in_list(p, _close_list) != nullptr) {
//                     continue;
//                 }
//                 // 
//                 int i;
//                 if (cur_point->j == -1) {
//                     i = 0;
//                 }
//                 else {
//                     i = (cur_point->j) * _l_num + cur_point->i + 1;
//                 }
//                 int j = (p->j) * _l_num + p->i + 1;
//                 double new_cost = cur_point->cost + _edges[i][j];
//                 SampleWaypoint* temp = is_in_list(p, _open_list);
//                 // 如果邻居点已经在openlist中，则更新邻居点地址，并判断是否需要更新信息
//                 if (temp != nullptr) {
//                     if (new_cost < temp->cost) {
//                         temp->parent = cur_point;
//                         temp->cost = new_cost;
//                     }
//                 }
//                 // 否则将邻居点加入openlist
//                 else {
//                     p->cost = new_cost;
//                     p->parent = cur_point;
//                     _open_list.push_back(p);
//                     total++;
//                 }
//             }
//         }
//     }
//     std::cout << "总点数" << total << std::endl;
//     // 将得到的路径粗解存入数组
//     while (end_point != nullptr) {
//         _rough_path.push_back(end_point->point);
//         if (end_point->parent != nullptr) {
//             // 使用多项式采样点进行增密
//             int i = end_point->parent->i + end_point->parent->j * _l_num + 1;
//             int j = end_point->parent->i + end_point->j * _l_num + 1;
//             i = i < 0 ? 0 : i;
//             // 倒序遍历多项式采样点
//             for (auto it = _polynomial_sample_points[i][j].rbegin(); it != _polynomial_sample_points[i][j].rend(); it++) {
//                 _rough_path.push_back(*it);
//             }
//         }
//         end_point = end_point->parent;
//     }
//     // 由于是倒序存入，所以需要将其进行反转。
//     std::reverse(_rough_path.begin(), _rough_path.end());
    
    
// }

SampleWaypoint* CreateConvexSpace::is_in_list(const SampleWaypoint* point, const std::list<SampleWaypoint*> l) {
    for (auto p : l) {
        if (point->i == p->i && point->j == p->j) {    
            return p;
        }
    }
    return nullptr;
}

bool CreateConvexSpace::get_neighbor_points(const SampleWaypoint *p ,std::vector<SampleWaypoint*> &points) {
    int j = p->j + 1;
    for (int i = 0; i < _l_num; i++) {
        int index = j * _l_num + i + 1;
        points.push_back(new SampleWaypoint(_sps[index])); // 这里是一个深拷贝
    }
}

SampleWaypoint* CreateConvexSpace::get_least_cost_point() {
    SampleWaypoint * least_cost_point = _open_list.front();
    for (auto it = _open_list.begin(); it != _open_list.end(); it++) {
        if (least_cost_point->cost > (*it)->cost) {
            least_cost_point = *it;
        }
    }
    return least_cost_point;
}

void CreateConvexSpace::get_rough_path(std::vector<waypoint> &rough_path) {
    rough_path = _rough_path;
}

void CreateConvexSpace::create_convex_space(const std::vector<obstacle> &obstacles, std::vector<double> &ub, std::vector<double> &lb) {
    // 道路上下边界默认为+-5（l坐标）
    ub.resize(_rough_path.size(), 5);
    lb.resize(_rough_path.size(), -5);
    for (const auto &obs : obstacles) {
        // 道路粗解的s坐标肯定在0到60之间
        if (obs.s < 0 || obs.s > 60) {
            continue;
        }
        int left_index = get_left_near_index(obs);
        int right_index = get_right_near_index(obs);
        if (left_index == -1) {
            left_index = 0;
        }
        if (right_index == -1) {
            right_index == _rough_path.size() - 1;
        }
        // 判断凸空间在障碍物的左边or右边
        if (obs.l > _rough_path[(right_index + left_index) / 2].l) {
            // 右边
            for (int i = left_index; i <= right_index; i++) {
                ub[i] = std::min(ub[i], obs.l - OBS_WIDTH / 2);
            }
        }
        else {
            // 右边
            for (int i = left_index; i <= right_index; i++) {
                lb[i] = std::max(lb[i], obs.l - OBS_WIDTH / 2);
            }
        }
    }
}

int CreateConvexSpace::get_left_near_index(const obstacle &obs) {
    for (int i = 1; i < _rough_path.size(); i++) {
        if (_rough_path[i].s > obs.s - OBS_LENGTH / 2) {
            return i - 1;
        }
    }
    return -1;
}

int CreateConvexSpace::get_right_near_index(const obstacle &obs) {
    for (int i = 1; i < _rough_path.size(); i++) {
        if (_rough_path[i].s > obs.s + OBS_LENGTH / 2) {
            return i;
        }
    }
    return -1;
}


