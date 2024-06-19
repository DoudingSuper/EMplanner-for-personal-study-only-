#pragma once
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>
#include "Eigen/Dense"
#include "define.h"
#include "OsqpEigen/OsqpEigen.h"

//定义地图中的设定物信息结构体
struct obstacle{
    double x = 0;
    double y = 0;
    double k = 0;
    double dk = 0;
    double dir = 0;
    // 这里的速度和加速度是在世界坐标系下
    double velocity = 0;
    double acc = 0;
    double t = 0;

    double s = 0;
    double l = 0;
    double dot_s = 0;
    double dot_l = 0;
    double d_l = 0;
    double ddot_s = 0;
    double ddot_l = 0;
    double dd_l = 0;
};
// 车辆定位信息
// 这里单独定义一套定位信息是由于在确定规划起点时，如果控制无法跟上
// 规划，则需要进行运动学或动力学外推，此时需要使用车辆定位信息，包括
// 各项速度和加速度
struct location{
    double x = 0;
    double y = 0;
    double dir = 0;
    // 注意，这里的速度和加速度都是在车辆坐标系下
    double vx = 0;
    double vy = 0;
    double ax = 0;
    double ay = 0;
    double t = 0;

    double s = 0;
    double l = 0;
    double dot_s = 0;
    double dot_l = 0;
    double d_l = 0;
    double ddot_s = 0;
    double ddot_l = 0;
    double dd_l = 0;
};
//定义地图中的路径点信息结构体
struct waypoint{
// 注意，对于轨迹点来说，其将车辆看做一个质点，因此每个轨迹点
// 的速度就是车辆所要跟踪的切向速度，加速度就是车辆的切向加速度，
// 轨迹点中不包括法向速度，因为质点的速度方向就是切向速度方向。
// 但对于实际车辆定位信息来说，其会包含横向速度和纵向速度（切向），
// 这是由于车辆的运动学导致的，所以在进行轨迹规划时需要先进行统一
// 转化到世界坐标系下
    double x = 0;
    double y = 0;
    double k = 0;
    double dk = 0;
    double dir = 0;
    // 这里的速度和加速度是在世界坐标系下
    double velocity = 0;
    double acc = 0;
    double t = 0;

    double s = 0;
    double l = 0;
    double dot_s = 0;
    double dot_l = 0;
    double d_l = 0;
    double ddot_s = 0;
    double ddot_l = 0;
    double dd_l = 0;
    void clear() {
        x = 0;
        y = 0;
        k = 0;
        dk = 0;
        dir = 0;
        velocity = 0;
        acc = 0;
        t = 0;
        s = 0;
        dot_s = 0;
        ddot_s = 0;
        l = 0;
        dot_l = 0;
        ddot_l = 0;
    }
};
// 读取障碍物信息
void readObs(std::vector<obstacle>& obs_array, const std::string& filename);
// 读取路径点信息
void readPath(std::vector<waypoint>& path_array, const std::string& filename);

typedef  Eigen::MatrixXd mxd;
typedef  Eigen::SparseMatrix<double> spr_Matrix;
typedef  Eigen::VectorXd vxd;
class Osqp_Solver{
private:
    // allocate QP problem matrices and vectores
    int _n;
    int _m;
    spr_Matrix _hessian; //P: 二次项系数矩阵 n*n正定矩阵,必须为稀疏矩阵SparseMatrix
    vxd _gradient; // Q: 线性项系数向量,n*1向量
    spr_Matrix _linear_matrix; // A: 线性约束矩阵，m*n矩阵,必须为稀疏矩阵SparseMatrix
    vxd _lower_bound; //L: 约束下界 m*1向量
    vxd _upper_bound; //U: 约束上界 m*1向量

public:
    bool init(int n, int m, spr_Matrix hessian, vxd gradient = vxd(), spr_Matrix linear_matrix = spr_Matrix(), vxd lb = vxd(), vxd ub = vxd());
    bool get_result(vxd& QP_solution);
};

class Projection {
private:
    const std::vector<waypoint>* _global_path;
    int _match_point_index; // 匹配点在全局路径中的索引
    int _increase_count; // 记录匹配点距离连续增大的次数
    bool _is_get_match_point = false;
    location _location;
    waypoint _match_point;
    waypoint _proj_point; // 投影点
    // bool init(const std::vector<waypoint> global_path);
    // 根据车辆位置计算在参考线上的匹配点
    // 这里进行函数重载，获取匹配点有两种情况，对于参考线平滑来说，寻找匹配点不需要从头找起
    bool cal_match_point_index();
    bool cal_match_point_index(const int &pre_match_index);
    // 根据匹配点计算车辆在参考线上的投影点
    bool cal_proj_point();
public:
    void init(const std::vector<waypoint>* const global_path);
    // 每次进行匹配点前先进行定位信息更新
    // 这里进行函数重载，接收到的点可能包括waypoint, location, obstacle三种类型
    bool update(const location& location);
    bool update(const waypoint& waypoint);
    bool update(const obstacle& location);
    // 
    bool get_match_point_index(int& index);
    bool get_match_point_index(const int &pre_match_index, int& index);
    // 
    bool get_proj_point(waypoint& proj_point);
};

class Path{
private:
    std::vector<waypoint> &_global_path;
    double _path_len;
public:
    Path(std::vector<waypoint> &global_path);
    // 对路径进行重新采样增密
    bool dense_path(double sample_interval);
    // 直线插值
    void resample_on_straight(const std::vector<waypoint> &points, const double &sample_interval);
    // 曲线插值
    void resample_on_curve(const std::vector<waypoint> &points, const double &sample_interval);
    // 根据三点计算中间点曲率
    void cal_curve(std::vector<waypoint> &curve_points);
    // 计算路径总长度
    void cal_length(std::vector<waypoint> &path, double &length);
    // 根据路径坐标计算并补充其方向角和曲率信息
    void get_dir_kappa();
};

class CartesianFrenetConvert{
private:
    std::unordered_map<int, double> _index_to_s;
    std::vector<waypoint> _ref_path;
public:
    // 初始化index_to_s, 根据参考线ref_path进行初始化，得到各个路径点对应的s
    void init_index_to_s(std::vector<waypoint> &ref_path, int origin_match_index, waypoint origin_proj_point);
    // 根据匹配点index计算匹配点对应的s，这里的s需要减去host匹配点的s
    void cal_s_from_index_to_s(const std::vector<waypoint> &ref_path, int match_point_index, waypoint origin_proj_point, double &s);
    void cartesian_to_frenet (waypoint &point, waypoint &proj_point);
    void cartesian_to_frenet (obstacle &point, waypoint &proj_point);
    void cartesian_to_frenet (const double x, const double y, const double theta,
                              const double kappa, const double v, const double a,
                              const double rx, const double ry, const double rtheta,
                              const double rkappa, const double rdkappa, const double rs,
                              double &s, double &dot_s, double &ddot_s,
                              double &l, double &dot_l, double &ddot_l);
    void frenet_to_cartesian (waypoint &point, waypoint &proj_point);
    void frenet_to_cartesian (obstacle &point, waypoint &proj_point);
    void frenet_to_cartesian (const double s, const double dot_s, const double ddot_s,
                              const double l, const double dot_l, const double ddot_l,
                              const double rx, const double ry, const double rtheta,
                              const double rkappa, const double rdkappa,
                              double &x, double &y, double &theta,
                              double &kappa, double &v, double &a);
};

class SampleWaypoint {
public:
    int i, j;
    waypoint point;
    double cost;
    SampleWaypoint *parent = nullptr;
    SampleWaypoint () {
        i = 0;
        j = -1;
        cost = 0;
        parent = nullptr;
    }
    SampleWaypoint (const SampleWaypoint &p) {
        i = p.i;
        j = p.j;
        point = p.point;
        cost = p.cost;
        parent = p.parent;
    }
    bool operator>(const SampleWaypoint &p) const {
        return cost > p.cost;
    }
};

class CreateConvexSpace {
    int _s_num, _l_num, _delta_s, _delta_l, _offset;
    waypoint _start_point;
    std::vector<SampleWaypoint> _sps; // 所有初始采样点，构成节点图
    std::vector<std::vector<long long>> _edges; // 初始节点之间的代价
    std::vector<std::vector<std::vector<waypoint>>> _polynomial_sample_points; // 初始节点之间的多项式采样点
    std::list<SampleWaypoint*> _open_list;
    std::list<SampleWaypoint*> _close_list;
    std::vector<waypoint> _rough_path; // 路径粗解
    // 计算
    SampleWaypoint* is_in_list(const SampleWaypoint* point, const std::list<SampleWaypoint*> l);
    bool get_neighbor_points(const SampleWaypoint *p ,std::vector<SampleWaypoint*> &points);
    SampleWaypoint* get_least_cost_point();
    double cal_obs_dist_cost(const waypoint &p, const std::vector<obstacle> &obstacles);
    // 计算节点间的五次多项式系数
    void cal_polynomial(const SampleWaypoint &p1,const SampleWaypoint &p2, Eigen::VectorXd &a);
    // 计算节点之间的路径cost
    void cal_cost(const std::vector<obstacle> &obstacles, double w1, double w2, double w3);
    // 在五次多项式上进行采样，在i和j两点之间进行五次多项式采样
    void sample_on_polynomial(const int &i, const int &j, double delta_s = 1);
public: 
    // 初始化，并确定所有节点s，l坐标
    void init(waypoint &start_point, int s_num, int l_num, double delta_s, double delta_l, double offset);
    // 创建路径粗解
    void create_rough_path(const std::vector<obstacle> &obstacles, double w1 = 10, double w2 = 10, double w3 = 100);
    // 获取路径粗解
    void get_rough_path(std::vector<waypoint> &rough_path);
    // 创建凸空间（获取路径上下边界）
    void create_convex_space(const std::vector<obstacle> &obstacles, std::vector<double> &ub, std::vector<double> &lb);
    // 获取道路上下边界
    // void get_convex_space(std::vector<double> &ub, std::vector<double> &lb);
    // 获取路径粗解中离障碍物左边界最近的点的索引
    int get_left_near_index(const obstacle &obs);
    int get_right_near_index(const obstacle &obs);
};


