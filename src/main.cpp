#include <vector>
#include <iostream>
#include "frenetUtils.h"
#include "path_smoother.h"
#include "matplotlibcpp.h"
#include "em_planner.h"
#include <chrono>
namespace plt = matplotlibcpp;
int main(){
    std::vector<waypoint> global_path;//全局路径
    readPath(global_path, "../data/waypoints.txt");
    std::vector<obstacle> obstacle_array;//障碍物信息
    readObs(obstacle_array, "../data/obstacles.txt");
    Path my_path(global_path);
    std::cout << "进行路径点重采样" << std::endl;
    my_path.dense_path(1);
    CartesianFrenetConvert converter;
    // 计算全局路径所有点的frenet坐标（l全为0，只需要计算s即可）
    converter.init_index_to_s(global_path, 0, global_path[0]);
    // 构建道路左右边界
    std::vector<waypoint> road_left_edge, road_right_edge;
    for (auto &p : global_path) {
        waypoint left_p = p, right_p = p;
        left_p.l += 5;
        right_p.l -= 5;
        converter.frenet_to_cartesian(left_p, p);
        converter.frenet_to_cartesian(right_p, p);
        road_left_edge.push_back(left_p);
        road_right_edge.push_back(right_p);
    }
    EMplanner planner;
    planner.init(global_path, obstacle_array, 0.1);
    location host_location;
    host_location.x = 42;
    host_location.y = 436;
    host_location.t = 0;
    std::vector<waypoint> trajectory;
    std::vector<waypoint> ref_line;
    std::vector<waypoint> qp_path;
    std::vector<waypoint> dp_path;
    int control = 0;
    while (true) {
        std::cout << "******周期：" << ++control << "******" << std::endl;
        std::cout << "******车辆位置：" << host_location.x << "   " << host_location.y << "******" << std::endl;
        planner.update_location_info(host_location);
        planner.run();
        planner.get_final_trajectory(trajectory);
        planner.get_ref_line(ref_line);
        planner.get_qp_path(qp_path);
        planner.get_dp_path(dp_path);
        plt::cla();
        plt::plotTrajectory(global_path);
        plt::plotTrajectory(road_left_edge, "g");
        plt::plotTrajectory(road_right_edge, "g");
        plt::plotTrajectory(ref_line, "y");
        plt::plotTrajectory(dp_path, "blue");
        plt::plotTrajectory(trajectory, "purple");
        plt::plotTrajectory(obstacle_array, ".r");
        plt::plot(std::vector<double> {host_location.x}, std::vector<double> {host_location.y}, "vc");
        plt::xlim(host_location.x - 40,host_location.x + 80);
        plt::ylim(host_location.y - 30,host_location.y + 80);
        plt::pause(0.1);
        host_location.x = qp_path.at(2).x;
        host_location.y = qp_path.at(2).y;
        host_location.l = qp_path.at(2).l;
        host_location.d_l = qp_path.at(2).d_l;
        host_location.dd_l = qp_path.at(2).dd_l;
        host_location.t = qp_path.at(2).t;
    }
    // // // plt::plotTrajectory(global_path, "y");
    // std::cout << "总路径点数" << global_path.size() << std::endl;
    // // my_path.get_dir_kappa();
    // // plt::plotKappa(global_path, "g");

    // location location1, location2;
    // location1.x = 65;
    // location1.y = 455;
    // location2.x = 250;
    // location2.y = 510;
    // EMplanner planner;
    // planner.init(global_path, obstacle_array, 0.1);
    // planner.update_location_info(location1);
    // planner.update_ref_line();
    // planner.update_start_point();
    // planner.convert_start_point();
    // planner.splice_trajectory();
    // planner.update_obstacle_info();
    // planner.convert_obstacles_point();
    // planner.dp_path_sample();
    // std::vector<waypoint> rough_path;
    // planner.get_dp_path(rough_path);
    // std::cout << "路径大小：" << rough_path.size() << std::endl;
    // planner.create_qp_path();
    // std::vector<waypoint> final_path;
    // planner.get_qp_path(final_path);
    // Path_Smoother path_smoother;
    // std::vector<waypoint> final_ref_line1, final_ref_line2;
    // path_smoother.init(&global_path);
    // auto t1 = std::chrono::steady_clock::now();
    // path_smoother.get_ref_line(location1);
    // std::cout << "进行参考线平滑" << std::endl;
    // path_smoother.smooth_path();

    // path_smoother.get_final_ref_line(final_ref_line1);
    // auto t2 = std::chrono::steady_clock::now();
    // path_smoother.get_ref_line(location2);
    // path_smoother.smooth_path();
    // path_smoother.get_final_ref_line(final_ref_line2);
    // auto t3 = std::chrono::steady_clock::now();
    // auto diff1 = std::chrono::duration <double, std::milli> (t2 - t1).count();
    // auto diff2 = std::chrono::duration <double, std::milli> (t3 - t2).count();
    // std::cout << "第一次运行时间"<< diff1 << std::endl << "第二次运行时间：" << diff2 << std::endl;
    // plt::cla();
    // plt::plotTrajectory(global_path);
    // plt::plotTrajectory(road_left_edge, "g");
    // plt::plotTrajectory(road_right_edge, "g");
    // plt::plotTrajectory(final_ref_line1, "y");
    // // plt::plotTrajectory(final_ref_line2, "r");
    // plt::plotTrajectory(rough_path, "black");
    // plt::plotTrajectory(final_path, "purple");
    // plt::plotTrajectory(obstacle_array, ".r");
    // plt::plot(std::vector<double> {location1.x}, std::vector<double> {location1.y}, "vc");

    // plt::pause(-1);

    // Projection proj;
    // proj.init(&global_path);
    // location.x = 123;
    // location.y = 500;
    // proj.update(location);
    // proj.cal_match_point_index();
    // int match_index;
    // proj.get_match_point_index(match_index);
    // std::cout << match_index << std::endl;
    // location.x = 50;
    // location.y = 100;
    // proj.update(location);
    // proj.cal_match_point_index();
    // proj.get_match_point_index(match_index);
    // std::cout << match_index << std::endl;
    // // while(true){
    // //     plt::cla();
    // //     plt::plotTrajectory(global_path);
    // //     plt::plotTrajectory(obstacle_array, "r");
    // //     plt::title("Visualize");
    // //     plt::grid(true);
    // //     plt::pause(0.1);
    // // }
        
    return 0;
}