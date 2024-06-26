#include <vector>
#include <iostream>
#include "frenetUtils.h"
#include "path_smoother.h"
#include "matplotlibcpp.h"
#include "em_planner.h"
#include <chrono>
#include <thread>
#include <atomic>
#include <fstream>
namespace plt = matplotlibcpp;
std::atomic<bool> running(true); //原子变量可以保证其每一次读写操作不可分割，提供了一种轻量级的线程互斥功能
void signalHandler(int a) {
    running = false;
}
int main(){
    std::ofstream file("timeWithParallel.csv");
    signal(SIGINT, signalHandler);
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
    EMplanner planner1, planner2;
    planner1.init(global_path, obstacle_array, 0.1);
    planner2.init(global_path, obstacle_array, 0.1);
    location host_location;
    host_location.x = 42;
    host_location.y = 436;
    host_location.t = 0;
    std::vector<waypoint> trajectory1, trajectory2;
    std::vector<waypoint> ref_line;
    std::vector<waypoint> qp_path;
    std::vector<waypoint> dp_path1, dp_path2;
    int control = 0;
    while (control < 160) {
        auto start_t = std::chrono::steady_clock::now();
        std::cout << "******周期：" << ++control << "******" << std::endl;
        std::cout << "******车辆位置：" << host_location.x << "   " << host_location.y << "******" << std::endl;
        planner1.update_location_info(host_location);
        planner2.update_location_info(host_location);
        std::thread t1(&EMplanner::run, &planner1, 2);
        std::thread t2(&EMplanner::run, &planner2, -2);
        t1.join();
        t2.join(); // 阻塞主线程，使其等待子线程结束
        // planner1.run();
        planner1.get_final_trajectory(trajectory1);
        planner2.get_final_trajectory(trajectory2);
        planner1.get_ref_line(ref_line);
        planner2.get_qp_path(qp_path);
        planner1.get_dp_path(dp_path1);
        planner2.get_dp_path(dp_path2);
        plt::cla();
        plt::plotTrajectory(global_path);
        plt::plotTrajectory(road_left_edge, "g");
        plt::plotTrajectory(road_right_edge, "g");
        plt::plotTrajectory(ref_line, "y");
        plt::plotTrajectory(obstacle_array, ".r");
        // 不知道为什么在t1这个线程中最后得到的轨迹在开头会多出五个点，非常抽象
        std::vector<waypoint> new_trajectory1(trajectory1.begin() + 5, trajectory1.end() - 3);
        std::vector<waypoint> new_trajectory2(trajectory2.begin(), trajectory2.end() - 3);
        plt::plotTrajectory(dp_path1, "pink");
        plt::plotTrajectory(dp_path2, "pink");
        plt::plotTrajectory(new_trajectory1, "purple");
        plt::plotTrajectory(new_trajectory2, "blue");
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
        auto end_t = std::chrono::steady_clock::now();
        auto diff = std::chrono::duration <double, std::milli> (end_t - start_t).count();
        std::cout << "总用时" << diff << std::endl;
        if (!file.is_open()) {
            std::cerr << "Error opening file" << std::endl;
            continue;
        }
        file << diff;
        file << "\n"; // 分行符
    }
    file.close();
    // // // plt::plotTrajectory(global_path, "y");
    // std::cout << "总路径点数" << global_path.size() << std::endl;
    // // my_path.get_dir_kappa();
    // // plt::plotKappa(global_path, "g");

    // location location1, location2;
    // location1.x = 65;
    // location1.y = 455;
    // location2.x = 250;
    // location2.y = 510;
    // EMplanner1 planner1;
    // planner1.init(global_path, obstacle_array, 0.1);
    // planner1.update_location_info(location1);
    // planner1.update_ref_line();
    // planner1.update_start_point();
    // planner1.convert_start_point();
    // planner1.splice_trajectory();
    // planner1.update_obstacle_info();
    // planner1.convert_obstacles_point();
    // planner1.dp_path_sample();
    // std::vector<waypoint> rough_path;
    // planner1.get_dp_path(rough_path);
    // std::cout << "路径大小：" << rough_path.size() << std::endl;
    // planner1.create_qp_path();
    // std::vector<waypoint> final_path;
    // planner1.get_qp_path(final_path);
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