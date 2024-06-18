#pragma once
#include "frenetUtils.h"
#include "path_smoother.h"

class EMplanner
{
private:
    /* data */
    Path_Smoother _path_smoother;
    Projection _projction;
    CartesianFrenetConvert _converter;
    // 注意，如果这里_ref_line是一个引用变量，则需要在初始化列表中进行定义，否则会报错
    std::vector<waypoint> _ref_line;
    bool _is_first_run = true;// 是否为第一次运行
    double _T;// 规划周期
    double _current_t;// 当前周期绝对时间
    std::vector<waypoint> _pre_trajectory;// 上个周期规划得到的轨迹，用于与本周期得到的轨迹进行拼接
    location _host_lacation;// 自车定位信息
    waypoint _start_point;// 规划起点信息
    int _start_index;// 规划起点在上个周期规划轨迹中的索引

    // 这里简化感知模块，假设车辆能得到地图上所有障碍物的信息
    std::vector<obstacle> _global_obstacles; // 地图上所有障碍物。
    std::vector<obstacle> _obstacles; // 障碍物信息
    std::vector<obstacle> _moving_obstacles; // 动态障碍物信息
    
    // 这里在SL空间中进行采样并通过二次规划得到规划路径粗解
    // 建立采样空间
    CreateConvexSpace _convex_space;
    std::vector<waypoint> _sample_points;
    std::vector<std::vector<double>> _sample_edges;
    // 最终路径
    std::vector<waypoint> _dp_path;
    std::vector<waypoint> _qp_path;
    Eigen::SparseMatrix<double> _Hessian;
    Eigen::SparseMatrix<double> _LinearMatrix;
    int _last_n;
    // 最终轨迹
    std::vector<waypoint> _final_trajectory;
    // Projection _projection;    
public:
    void init(const std::vector<waypoint> &global_path, const std::vector<obstacle> &global_obstacles, double T);
    // 更新车辆定位信息
    bool update_location_info(const location & location);
    // 更新参考线
    bool update_ref_line();
    // 获取规划起点
    bool update_start_point();
    // 进行轨迹拼接
    bool splice_trajectory();
    // 获取上个周期轨迹中需要拼接的部分
    bool get_splice_points(std::vector<waypoint> &splice_points);
    // 对规划起点进行坐标转化（笛卡尔转自然坐标）
    bool convert_start_point();
    // 获取障碍物信息（筛选一定范围内的障碍物）
    bool update_obstacle_info();
    // 对障碍物坐标进行坐标转化（笛卡尔转自然坐标）
    bool convert_obstacles_point();
    // 在自然坐标系中进行路径点采样，用于动态规划
    bool dp_path_sample(int s_num = 6, int l_num = 9, double delta_s = 10, double delta_l = 1);
    // 在五次多项式连接线上进行采样
    void convert_dp_path();
    // 获取二次规划上下界
    void get_qp_bound(std::vector<double> &ub, std::vector<double> &lb);
    // 使用二次规划进行路径平滑
    /* ***调参注意事项：***
    这里的参数会影响到二次规划能否求解成功，暂时不知道原因
    */
    void create_qp_path(double w_ref = 5, double w_center = 10, double w_dl = 100, double w_ddl = 100, double w_dddl = 1);
    // 获取最终路径
    void convert_qp_path();
    // 在st坐标系中进行采样，用于速度动态规划
    bool dp_speed_sample(int s_num = 6, int t_num = 16, double delta_s = 0.1, double delta_t = 0.5);
    // 进行速度二次规划
    void create_qp_speed(double w_ref = 1, double w_center = 1, double w_dds = 100, double w_ddds = 1000);
    // 合成最终轨迹
    void create_trajectory();
    // 获取最终轨迹
    void get_final_trajectory(std::vector<waypoint> &trajectory);
    // 获取参考线
    void get_ref_line(std::vector<waypoint> &ref_line);
    // 获取qp路径
    void get_qp_path(std::vector<waypoint> &qp_path);
    void get_dp_path(std::vector<waypoint> &dp_path);
    // 规划器运行
    void run();

};

