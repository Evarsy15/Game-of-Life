#include <stdio.h>
#include <stdlib.h>
#ifdef BUILD_WITH_OPENMPI
#include <mpi.h>
#endif

#include "cell_board.h"
#include "aux.h"

using namespace Nix;

#ifdef BUILD_WITH_OPENMPI
int main(int argc, char* argv[]) {
    Timer openmpi_timer;
    openmpi_timer.Start();

    // Initialize MPI Environment.
    MPI_Init(&argc, &argv);

    // Get MPI Process Informations.
    int num_proc, proc_id;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    // Get Program Arguments by standard input stream.
    // ※ Note : By default, standard input of an MPI program is 
    //           only processed in Process ID 0. 
    uint board_size, target_gen, ghost_size;
    if (proc_id == 0) {
        std::cin >> board_size;
        std::cin >> target_gen;
        std::cin >> ghost_size;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    
    char **init_cell_board;
    if (proc_id == 0) {
        init_cell_board = new char*[board_size];
        for (int i = 0; i < board_size; i++)
            init_cell_board[i] = new char[board_size];
    }

    if (proc_id == 0) {
        for (int i = 0; i < board_size; i++) {
            for (int j = 0; j < board_size; j++) {
                std::cin >> init_cell_board[i][j];
            }
        }
    }
    
    // if (proc_id == 0) {
    //     for (int i = 0; i < board_size; i++) {
    //         for (int j = 0; j < board_size; j++)
    //             std::cout << init_cell_board[i][j];
    //         std::cout << "\n";
    //     }
    // }

    // Spread Program Arguments (except initial cell-board status)
    MPI_Bcast(&board_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&target_gen, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ghost_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

#if defined(CARTESIAN_PARTITION)
#ifdef ROWWISE_PARTITION
#error CARTESIAN_PARTITION and ROWWISE_PARTITION cannot be defined at once.
#endif
    
    // Construct 2D-Mesh Topology
    int cart_dims[2], cart_coords[2];
    int cart_periods[2] = {0, 0}; // 2D-Mesh Topology
    compute_cart_dims(cart_dims, num_proc);
    
    MPI_Comm cartesian_topology;
    MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dims, cart_periods, 0, &cartesian_topology);
    MPI_Cart_coords(cartesian_topology, proc_id, 2, cart_coords);
//  compute_cart_coords(cart_coords, proc_id, cart_dims);
    
    // Spread initial cell-board status
    uint cell_base[2], cell_size[2];
    compute_cart_cell_info(cell_base, cell_size, cart_dims, cart_coords, board_size);
    #ifdef DEBUG
    printf("Process ID %d : Dim(%d, %d), Coord(%d, %d), Base(%d, %d), Size(%d, %d)\n",
            proc_id, cart_dims[0], cart_dims[1], cart_coords[0], cart_coords[1],
            cell_base[0], cell_base[1], cell_size[0], cell_size[1]);
    #endif
    
    char **my_init_cell_board = new char*[cell_size[0]];
    for (uint i = 0; i < cell_size[0]; i++)
        my_init_cell_board[i] = new char[cell_size[1]];
    
    if (proc_id == 0) {
        for (int i = 1; i < num_proc; i++) {
            int dest_cart_coords[2];
            uint dest_cell_base[2], dest_cell_size[2];
            MPI_Cart_coords(cartesian_topology, i, 2, dest_cart_coords);
            // compute_cart_coords(dest_cart_coords, i, cart_dims);
            compute_cart_cell_info(dest_cell_base, dest_cell_size, cart_dims, 
                                   dest_cart_coords, board_size);
            
            uint row_base = dest_cell_base[0]; uint col_base = dest_cell_base[1];
            uint row_cnt  = dest_cell_size[0]; uint col_cnt  = dest_cell_size[1];
            #ifdef DEBUG
            printf("Process ID 0 : Sending init_cell_board[%d][%d] ~ init_cell_board[%d][%d] to Process %d...\n",
                    row_base, col_base, row_base+row_cnt-1, col_base+col_cnt-1, i);
            #endif
            for (int j = 0; j < row_cnt; j++) {
                MPI_Send(&init_cell_board[row_base+j][col_base], col_cnt,
                          MPI_CHAR, i, j, cartesian_topology);//, &req[i-1][j]);
            }
            #ifdef DEBUG
            printf("Process ID 0 : Sent to Process %d\n", i);
            #endif
        }
        for (int i = 0; i < cell_size[0]; i++)
            for (int j = 0; j < cell_size[1]; j++)
                my_init_cell_board[i][j] = init_cell_board[i][j];
    } else {
        uint row_base = cell_base[0]; uint col_base = cell_base[1];
        uint row_cnt  = cell_size[0]; uint col_cnt  = cell_size[1];
        #ifdef DEBUG
        printf("Process ID %d : Receiving init_cell_board[%d][%d] ~ init_cell_board[%d][%d] from Process 0...\n",
                proc_id, row_base, col_base, row_base+row_cnt-1, col_base+col_cnt-1);
        #endif
        for (int j = 0; j < row_cnt; j++) {
            MPI_Recv(&my_init_cell_board[j][0], col_cnt, MPI_CHAR, 
                     0, j, cartesian_topology, MPI_STATUS_IGNORE);
        }
        #ifdef DEBUG
        printf("Process ID %d : Received from Process 0\n", proc_id);
        #endif
    }

    if (proc_id == 0) {
        for (int i = 0; i < board_size; i++)
            delete[] init_cell_board[i];
        delete[] init_cell_board;
    }

#ifdef DEBUG
    std::string text_str = "input_partition" + std::to_string(proc_id) + ".txt";
    FILE *input_partition = fopen(text_str.c_str(), "w");
    fprintf(input_partition, "Process ID : %d\n", proc_id);
    fprintf(input_partition, "Cartesian Coordinate : (%d, %d)\n", cart_coords[0], cart_coords[1]);
    fprintf(input_partition, "Cell Base : (%d, %d)\n", cell_base[0], cell_base[1]);
    fprintf(input_partition, "Cell Size : (%d, %d)\n", cell_size[0], cell_size[1]);
    fprintf(input_partition, "Board Input :\n");
    for (int i = 0; i < cell_size[0]; i++){
        for (int j = 0; j < cell_size[1]; j++)
            fprintf(input_partition, "%c", my_init_cell_board[i][j]);
        fprintf(input_partition, "\n");
    }
    fclose(input_partition);
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    // Configure Ghost-size to minimum cell-tile size.
    uint row_min = board_size / cart_dims[0];
    uint col_min = board_size / cart_dims[1];
    ghost_size = std::min(ghost_size, std::min(row_min, col_min));

    // Run Game-of-Life Game
    uint global_board_size[2] = {board_size, board_size};
    CartesianCellBoard my_cell_board(proc_id, cell_base, cell_size,
                                     global_board_size, ghost_size,
                                     my_init_cell_board, &cartesian_topology);
    for (int i = 0; i < target_gen; i += ghost_size) {
        // std::cout << "Process ID " << proc_id << " : Generation " << i << "\n";
        for (int j = 0; j < ghost_size; j++) {
            if (i+j >= target_gen)
                break;
            my_cell_board.cycle();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Gather results
    bool **final_cell_board;
    if (proc_id == 0) {
        final_cell_board = new bool*[board_size];
        for (uint i = 0; i < board_size; i++)
            final_cell_board[i] = new bool[board_size];
    }

    bool **my_final_cell_board = my_cell_board.get_board_status(target_gen);

#ifdef DEBUG
    std::string output_str = "output_partition" + std::to_string(proc_id) + ".txt";
    FILE *output_partition = fopen(output_str.c_str(), "w");
    fprintf(output_partition, "Process ID : %d\n", proc_id);
    fprintf(output_partition, "Cartesian Coordinate : (%d, %d)\n", cart_coords[0], cart_coords[1]);
    fprintf(output_partition, "Cell Base : (%d, %d)\n", cell_base[0], cell_base[1]);
    fprintf(output_partition, "Cell Size : (%d, %d)\n", cell_size[0], cell_size[1]);
    fprintf(output_partition, "Board Output :\n");
    int k = ghost_size;
    for (int i = 0; i < cell_size[0]; i++) {
        for (int j = 0; j < cell_size[1]; j++) {
            char c = my_final_cell_board[i+k][j+k] ? '#' : '.';
            fprintf(output_partition, "%c", c);
        }
        fprintf(output_partition, "\n");
    }
    fclose(output_partition);
#endif
    
    if (proc_id == 0) {
        for (int i = 1; i < num_proc; i++) {
            int dest_cart_coords[2];
            uint dest_cell_base[2], dest_cell_size[2];
            MPI_Cart_coords(cartesian_topology, i, 2, dest_cart_coords);
            compute_cart_cell_info(dest_cell_base, dest_cell_size, cart_dims, 
                                   dest_cart_coords, board_size);
            
            uint row_base = dest_cell_base[0]; uint col_base = dest_cell_base[1];
            uint row_cnt  = dest_cell_size[0]; uint col_cnt  = dest_cell_size[1];
            #ifdef DEBUG
            printf("Process 0 : Receiving final_cell_board[%d][%d] ~ final_cell_board[%d][%d] from Process %d...\n",
                    row_base, col_base, row_base+row_cnt-1, col_base+col_cnt-1, i);
            #endif
            for (int j = 0; j < row_cnt; j++) {
                MPI_Recv(&final_cell_board[row_base+j][col_base], col_cnt, MPI_C_BOOL,
                         i, board_size + j, cartesian_topology, MPI_STATUS_IGNORE);
            }
            #ifdef DEBUG
            printf("Process 0 : Received from Process %d\n", i);
            #endif
        }
        for (int i = 0; i < cell_size[0]; i++) {
            for (int j = 0; j < cell_size[1]; j++)
                final_cell_board[i][j] = my_final_cell_board[i+ghost_size][j+ghost_size];
        }
    } else {
        uint row_base = cell_base[0]; uint col_base = cell_base[1];
        uint row_cnt  = cell_size[0]; uint col_cnt  = cell_size[1];
        #ifdef DEBUG
        printf("Process %d : Sending final_cell_board[%d][%d] ~ final_cell_board[%d][%d] from Process 0...\n",
                proc_id, row_base, col_base, row_base+row_cnt-1, col_base+col_cnt-1);
        #endif
        for (int j = 0; j < row_cnt; j++) {
            MPI_Send(&my_final_cell_board[ghost_size+j][ghost_size], col_cnt, MPI_C_BOOL, 
                     0, board_size + j, cartesian_topology);
        }
        #ifdef DEBUG
        printf("Process %d : Sent to Process 0\n", proc_id);
        #endif
    }

    if (proc_id == 0) {
        for (int i = 0; i < board_size; i++) {
            for (int j = 0; j < board_size; j++) {
                char c = final_cell_board[i][j] ? '#' : '.';
                std::cout << c;
            }
            std::cout << "\n";
        }
    }

#elif defined(ROWWISE_PARTITION)
#ifdef CARTESIAN_PARTITION
#error CARTESIAN_PARTITION and ROWWISE_PARTITION cannot be defined at once.
#endif
    // Construct 1D Mesh Topology

#else
// #error Define CARTESIAN_PARTITION or ROWWISE_PARTITION.
#endif

    MPI_Finalize();

    openmpi_timer.Stop();
    if (proc_id == 0)
        openmpi_timer.Print("Running Time");
}
#else
int main(int argc, char* argv[]) {
    Timer single_timer;
    single_timer.Start();

    uint board_size, target_gen, ghost_size;
    #ifdef DEBUG
    single_input_prompt("Enter the size of (square) board", board_size);
    single_input_prompt("Enter the final generation to display", target_gen);
    single_input_prompt("Enter the size of ghost cell", ghost_size);
    #else
    std::cin >> board_size;
    std::cin >> target_gen;
    std::cin >> ghost_size;
    #endif
    char **init_cell_board = new char*[board_size];
    for (uint i = 0; i < board_size; i++)
        init_cell_board[i] = new char[board_size];
    #ifdef DEBUG
    cell_board_input_prompt("Enter the initial cell-board status (# : Alive, . : Dead)",
                           init_cell_board, board_size);
    #else
    for (uint i = 0; i < board_size; i++) {
        scanf("%s", init_cell_board[i]);
    }
    #endif

    CellBoard cell_board(board_size, board_size, init_cell_board);
    for (int i = 0; i < target_gen; i++)
        cell_board.cycle();
    cell_board.print();

    for (uint i = 0; i < board_size; i++)
        delete[] init_cell_board[i];
    delete[] init_cell_board;

    single_timer.Stop();
    single_timer.Print("Running Time");
}
#endif