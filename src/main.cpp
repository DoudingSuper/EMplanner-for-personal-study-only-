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
#include <mpi.h>

namespace plt = matplotlibcpp;
std::atomic<bool> running(true); //原子变量可以保证其每一次读写操作不可分割，提供了一种轻量级的线程互斥功能
void signalHandler(int a) {
    running = false;
}
int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int world_rank; // 进程标识
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size; // 进程数量
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size != 2) {
        std::cerr << "This program requires exactly 2 MPI processes." << std::endl;
        MPI_Finalize();
        return 1;
    }
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
    // 定义 MPI 数据类型
    MPI_Datatype waypoint_type;
    int block_lengths[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[15] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                              MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint offsets[15];
    offsets[0] = offsetof(waypoint, x);
    offsets[1] = offsetof(waypoint, y);
    offsets[2] = offsetof(waypoint, k);
    offsets[3] = offsetof(waypoint, dk);
    offsets[4] = offsetof(waypoint, dir);
    offsets[5] = offsetof(waypoint, velocity);
    offsets[6] = offsetof(waypoint, acc);
    offsets[7] = offsetof(waypoint, t);
    offsets[8] = offsetof(waypoint, s);
    offsets[9] = offsetof(waypoint, l);
    offsets[10] = offsetof(waypoint, dot_s);
    offsets[11] = offsetof(waypoint, dot_l);
    offsets[12] = offsetof(waypoint, d_l);
    offsets[13] = offsetof(waypoint, ddot_s);
    offsets[14] = offsetof(waypoint, ddot_l);

    MPI_Type_create_struct(15, block_lengths, offsets, types, &waypoint_type);
    MPI_Type_commit(&waypoint_type);
    std::ofstream file("timeWithMPI&Parallel.csv");
    std::vector<waypoint> trajectory1, trajectory2;
    std::vector<waypoint> ref_line;
    std::vector<waypoint> qp_path;
    std::vector<waypoint> dp_path1, dp_path2;
    int size1, size2, size_ref, size_qp,size_dp1, size_dp2;
    location host_location;
    host_location.x = 42;
    host_location.y = 436;
    host_location.t = 0;
    if (world_rank == 0) {
        int control = 0;
        EMplanner planner1, planner2;
        planner1.init(global_path, obstacle_array, 0.1);
        planner2.init(global_path, obstacle_array, 0.1);
        while (control < 160) {
            auto start_t = std::chrono::steady_clock::now();
            std::cout << "******周期：" << ++control << "******" << std::endl;
            MPI_Send(&control, 1,  MPI_INT, 1, 12, MPI_COMM_WORLD);
            std::cout << "******车辆位置：" << host_location.x << "   " << host_location.y << "******" << std::endl;
            planner1.update_location_info(host_location);
            planner2.update_location_info(host_location);
            std::thread t1(&EMplanner::run, &planner1, 2);
            std::thread t2(&EMplanner::run, &planner2, -2);
            t1.join();
            t2.join();
            planner1.get_final_trajectory(trajectory1);
            planner1.get_ref_line(ref_line);
            planner1.get_dp_path(dp_path1);
            size1 = trajectory1.size();
            size_ref = ref_line.size();
            size_dp1 = dp_path1.size();
            MPI_Send(&size1, 1,  MPI_INT, 1, 6, MPI_COMM_WORLD);
            MPI_Send(&size_ref, 1,  MPI_INT, 1, 8, MPI_COMM_WORLD);
            MPI_Send(&size_dp1, 1,  MPI_INT, 1, 10, MPI_COMM_WORLD);
            MPI_Send(trajectory1.data(), trajectory1.size(), waypoint_type, 1, 0, MPI_COMM_WORLD);
            MPI_Send(ref_line.data(), ref_line.size(), waypoint_type, 1, 2, MPI_COMM_WORLD);
            MPI_Send(dp_path1.data(), dp_path1.size(), waypoint_type, 1, 4, MPI_COMM_WORLD);
            // planner2.run(2);
            // planner1.run();
            planner2.get_final_trajectory(trajectory2);
            planner2.get_qp_path(qp_path);
            planner2.get_dp_path(dp_path2);
            size2 = trajectory2.size();
            size_qp = qp_path.size();
            size_dp2 = dp_path2.size();
            MPI_Send(&size2, 1,  MPI_INT, 1, 7, MPI_COMM_WORLD);
            MPI_Send(&size_qp, 1,  MPI_INT, 1, 9, MPI_COMM_WORLD);
            MPI_Send(&size_dp2, 1,  MPI_INT, 1, 11, MPI_COMM_WORLD);
            MPI_Send(trajectory2.data(), trajectory2.size(), waypoint_type, 1, 1, MPI_COMM_WORLD);
            MPI_Send(qp_path.data(), qp_path.size(), waypoint_type, 1, 3, MPI_COMM_WORLD);
            MPI_Send(dp_path2.data(), dp_path2.size(), waypoint_type, 1, 5, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (qp_path.size() > 2) {
                host_location.x = qp_path.at(2).x;
                host_location.y = qp_path.at(2).y;
                host_location.l = qp_path.at(2).l;
                host_location.d_l = qp_path.at(2).d_l;
                host_location.dd_l = qp_path.at(2).dd_l;
                host_location.t = qp_path.at(2).t;
            }
            auto end_t = std::chrono::steady_clock::now();
            auto diff = std::chrono::duration <double, std::milli> (end_t - start_t).count();
            std::cout << "总用时" << diff << std::endl;
            if (!file.is_open()) {
                std::cerr << "Error opening file" << std::endl;
                // continue;
            }
            file << diff;
            file << "\n"; // 分行符
        }
    }
    if (world_rank == 1 ) {
        int control = 0;
        while (control < 160) {
            plt::cla();
            plt::plotTrajectory(global_path);
            plt::plotTrajectory(road_left_edge, "g");
            plt::plotTrajectory(road_right_edge, "g");
            plt::plotTrajectory(obstacle_array, ".r");
            MPI_Recv(&control, 1, MPI_INT, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&size1, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&size_ref, 1, MPI_INT, 0, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&size_dp1, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            trajectory1.resize(size1);
            ref_line.resize(size_ref);
            dp_path1.resize(size_dp1);
            MPI_Recv(trajectory1.data(), size1 * sizeof(waypoint), waypoint_type, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(ref_line.data(), size_ref, waypoint_type, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(dp_path1.data(), size_dp1, waypoint_type, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 不知道为什么在t1这个线程中最后得到的轨迹在开头会多出五个点，非常抽象
            std::vector<waypoint> new_trajectory1(trajectory1.begin(), trajectory1.end() - 5);
            std::vector<waypoint> new_ref_line(ref_line.begin(), ref_line.end() - 10);
            std::vector<waypoint> new_dp_path1(dp_path1.begin(), dp_path1.end() - 5);
            plt::plotTrajectory(new_ref_line, "y");
            plt::plotTrajectory(new_trajectory1, "purple");
            plt::plotTrajectory(new_dp_path1, "pink");
            MPI_Recv(&size2, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&size_qp, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&size_dp2, 1, MPI_INT, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            trajectory2.resize(size2);
            qp_path.resize(size_qp);
            dp_path2.resize(size_dp2);
            MPI_Recv(trajectory2.data(), size2, waypoint_type, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(qp_path.data(), size_qp, waypoint_type, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(dp_path2.data(), size_dp2, waypoint_type, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "进程1************" << std::endl;
            std::vector<waypoint> new_trajectory2(trajectory2.begin(), trajectory2.end() - 5);
            std::vector<waypoint> new_dp_path2(dp_path2.begin(), dp_path2.end() - 5);
            plt::plotTrajectory(new_trajectory2, "blue");
            plt::plotTrajectory(new_dp_path2, "pink");
            if (qp_path.size() > 2) {
                host_location.x = qp_path.at(2).x;
                host_location.y = qp_path.at(2).y;
                host_location.l = qp_path.at(2).l;
                host_location.d_l = qp_path.at(2).d_l;
                host_location.dd_l = qp_path.at(2).dd_l;
                host_location.t = qp_path.at(2).t;
            }
            plt::plot(std::vector<double> {host_location.x}, std::vector<double> {host_location.y}, "vc");
            plt::xlim(host_location.x - 40,host_location.x + 80);
            plt::ylim(host_location.y - 30,host_location.y + 80);
            plt::pause(0.1);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    file.close();
    MPI_Finalize();
    return 0;
}