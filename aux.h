#ifndef _AUX_H_
#define _AUX_H_

#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

using uint = unsigned;

namespace Nix {

template <typename T>
void single_input_prompt(std::string prompt_msg, T& input) {
    std::cout << prompt_msg << " : ";
    std::cin >> input;
}

void cell_board_input_prompt(std::string prompt_msg, char** cell_board, uint board_size) {
    std::cout << prompt_msg << " : \n";
    for (uint i = 0; i < board_size; i++) {
        scanf("%s", cell_board[i]);
    }
}

#ifdef CARTESIAN_PARTITION
void compute_cart_dims(int cart_dims[], uint num_proc) {
    int divisor = static_cast<int>(std::sqrt((double)num_proc));
    while (divisor > 0) {
        if (num_proc % divisor == 0)
            break;
        divisor--;
    }
    cart_dims[0] = num_proc / divisor;
    cart_dims[1] = divisor;
}

void compute_cart_coords(int cart_coords[], int proc_id, int cart_dims[]) {
    cart_coords[0] = proc_id / cart_dims[0];
    cart_coords[1] = proc_id % cart_dims[0];
}

void compute_cart_cell_info(uint cell_base[], uint cell_size[], 
                       int cart_dims[], int cart_coords[], uint board_size) {
    uint rdim = cart_dims[0], cdim = cart_dims[1];
    uint rid = cart_coords[0], cid = cart_coords[1];

    cell_base[0] = (board_size * rid) / rdim;
    cell_base[1] = (board_size * cid) / cdim;

    cell_size[0] = ((board_size * (rid+1)) / rdim) - cell_base[0];
    cell_size[1] = ((board_size * (cid+1)) / cdim) - cell_base[1];
}
#endif

#ifdef ROWWISE_PARTITION

#endif

struct Timer {
    enum TimerStatus {
        RUNNING,
        PAUSED,
        STOPPED,
        NUM_TIMER_STATUS
    };

    using clock = std::chrono::high_resolution_clock;
    using time_point_t = std::chrono::high_resolution_clock::time_point;
    using duration_t = std::chrono::duration<double>;
    
    TimerStatus stat = STOPPED;
    time_point_t start, end;
    double duration;

    void Start() {
        // assert(stat == STOPPED);

        stat = RUNNING;
        duration = 0.0;
        start = clock::now();
    }

    void Pause() {
        end = clock::now();

        // assert(stat == RUNNING);
        stat = PAUSED;

        duration_t delta = end - start;
        duration += static_cast<double>(delta.count());
    }

    void Continue() {
        // assert(stat == PAUSED);

        stat = RUNNING;
        start = clock::now();
    }

    void Stop() {
        end = clock::now();

        if (stat == RUNNING) {
            duration_t delta = end - start;
            duration += static_cast<double>(delta.count());
        }
        
        stat = STOPPED;
    }

    void Print(const std::string info) {
        std::cout << info << " : " << duration << "s \n";
    }
};

} // namespace Nix

#endif