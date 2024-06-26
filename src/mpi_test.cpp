#include <mpi.h>
#include <iostream>
#include <cstring>

int main(int argc, char** argv) {
    // 初始化 MPI 环境
    MPI_Init(&argc, &argv);

    // 获取所有进程中的当前进程的 rank（ID）
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // 获取所有进程的总数
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // 发送和接收的缓冲区
    const int MAX_MESSAGE_SIZE = 100;
    char message[MAX_MESSAGE_SIZE];

    if (world_rank == 0) {
        // 如果是进程 0，则发送消息
        std::strcpy(message, "Hello, World");
        for (int i = 1; i < world_size; i++) {
            MPI_Send(message, std::strlen(message) + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
        std::cout << "Process 0 sent message: " << message << std::endl;
    } else {
        // 如果是其他进程，则接收消息
        MPI_Recv(message, MAX_MESSAGE_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "Process " << world_rank << " received message: " << message << std::endl;
    }

    // 结束 MPI 环境
    MPI_Finalize();

    return 0;
}