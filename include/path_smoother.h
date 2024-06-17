# pragma once 

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <Eigen/Dense>
#include "frenetUtils.h"


class Path_Smoother
{
private:
    double w_similar_cost = 10;//相似代价系数
    double w_smoothing_cost = 1000; //平滑代价系数
    double w_compact_cost = 10;//紧凑代价系数
    const std::vector<waypoint>* _global_path; // 全局路径
    Projection _projction; // 用于计算匹配点索引
    waypoint _location; // 车辆定位信息
    int _match_point_index; // 匹配点索引
    int _pre_match_index;// 用于表征上一次平滑时的匹配点，本次寻找参考线时将以这个点作为起点
    int _pre_start_index = -1; //-1表示上一次的参考线起始点还未计算
    int _smooth_points_num; // 需要进行平滑的点数（各次平滑之间存在重叠部分）
    std::vector<waypoint> _ref_path;
    std::vector<waypoint> _smooth_ref_path;
    bool _is_first_smooth = true;//用于表征是否是第一次进行参考线平滑
    Osqp_Solver qp_solver;
public:
    void init(const std::vector<waypoint>* const global_path);//初始化，将全局路径写入
    void update_match_index(const location &location);//更新匹配点，根据匹配点更新参考线
    void get_ref_line(const location &location);//根据匹配点，获取参考线
    void smooth_path();//对参考线做平滑
    bool get_final_ref_line(std::vector<waypoint>& final_ref_line);//获得平滑后的参考线
};

// Path_Smoother::Path_Smoother(/* args */)
// {
// }

// Path_Smoother::~Path_Smoother()
// {
// }
