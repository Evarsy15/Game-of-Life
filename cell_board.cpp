
#include "cell_board.h"

namespace Nix {

CellBoard::CellBoard(uint nrow, uint ncol, char **init_cell_board) {
    m_cell_size[0] = nrow;
    m_cell_size[1] = ncol;

    for (int k = 0; k < 2; k++) {
        m_cell_board[k] = new bool*[nrow];
        for (uint i = 0; i < nrow; i++)
            m_cell_board[k][i] = new bool[ncol];
    }
    
    for (uint i = 0; i < nrow; i++) {
        for (uint j = 0; j < ncol; j++) {
            if (init_cell_board[i][j] == '#')
                m_cell_board[0][i][j] = true;
            else if (init_cell_board[i][j] == '.')
                m_cell_board[0][i][j] = false;
            else {
                std::cout << "Warning: Unsupported Initial State " 
                          << init_cell_board[i][j] << "at (" << i << ", " << j << ") "  
                          << "Detected.\n";
                std::cout << "Initialize as dead cell.\n";
                    
                m_cell_board[0][i][j] = false;
            }
            m_cell_board[1][i][j] = false;
        }
    }
    
//  m_alive_cells.clear();
    m_cycle = 0;
}

CellBoard::~CellBoard() {
    for (int k = 0; k < 2; k++) {
        for (uint i = 0; i < m_cell_size[0]; i++)
            delete[] m_cell_board[k][i];
        delete[] m_cell_board[k];
    }
}

void CellBoard::cycle() {
    uint nrow = m_cell_size[0];
    uint ncol = m_cell_size[1];
    uint _cur = m_cycle % 2;
    uint _new = (_cur + 1) % 2;

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            uint alive_neighbor_cell = 0;
            for (int p = -1; p <= 1; p++) {
                for (int q = -1; q <= 1; q++) {
                    // Out-of-bound
                    if (i+p < 0 || i+p >= nrow || j+q < 0 || j+q >= ncol)
                        continue;
                    // Self
                    if (p == 0 && q == 0)
                        continue;
                    // Count alive neighbor cells
                    if (m_cell_board[_cur][i+p][j+q])
                        alive_neighbor_cell++;
                }
            }
            if (alive_neighbor_cell < 2 || alive_neighbor_cell > 3)
                m_cell_board[_new][i][j] = false;
            else if (alive_neighbor_cell == 2)
                m_cell_board[_new][i][j] = m_cell_board[_cur][i][j];
            else if (alive_neighbor_cell == 3)
                m_cell_board[_new][i][j] = true;
        }
    }
    
    m_cycle++;
}

void CellBoard::print() {
    uint nrow = m_cell_size[0];
    uint ncol = m_cell_size[1];
    uint cur = m_cycle % 2;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            char c = m_cell_board[cur][i][j] ? '#' : '.';
            std::cout << c;
        }
        std::cout << "\n";
    }
}

#ifdef BUILD_WITH_OPENMPI

#ifdef CARTESIAN_PARTITION

CartesianCellBoard::CartesianCellBoard(uint proc_id,
                                       uint cell_base[],
                                       uint cell_size[],
                                       uint board_size[],
                                       uint ghost_size,
                                       char **init_cell_board,
                                       MPI_Comm *cartesian_topology) {
    m_proc_id = proc_id;

    m_cell_base[0] = cell_base[0];
    m_cell_base[1] = cell_base[1];
    m_cell_size[0] = cell_size[0];
    m_cell_size[1] = cell_size[1];
    m_board_size[0] = board_size[0];
    m_board_size[1] = board_size[1];
    m_ghost_size = ghost_size;
    
    uint nrow = m_cell_size[0];
    uint ncol = m_cell_size[1];
    uint k = m_ghost_size;

    for (uint m = 0; m < 2; m++) {
        m_cell_board[m] = new bool*[nrow+2*k];
        for (uint i = 0; i < nrow+2*k; i++)
            m_cell_board[m][i] = new bool[ncol+2*k];
    }

    for (uint i = 0; i < nrow+2*k; i++) {
        for (uint j = 0; j < ncol+2*k; j++) {
            if (is_ghost(i, j))
                m_cell_board[0][i][j] = false;
            else if (init_cell_board[i-k][j-k] == '#')
                m_cell_board[0][i][j] = true;
            else if (init_cell_board[i-k][j-k] == '.')
                m_cell_board[0][i][j] = false;
            else {
                std::cout << "Warning: Unsupported Initial State " 
                          << init_cell_board[i-k][j-k] << "at (" << i-k << ", " << j-k << ") "  
                          << "Detected.\n";
                std::cout << "Initialize as dead cell.\n";
                    
                m_cell_board[0][i][j] = false;
            }
            m_cell_board[1][i][j] = false;
        }
    }

    m_cycle = 0;
    m_cartesian_topology = cartesian_topology;

    other_left_ghost  = new bool[k * nrow];
    other_right_ghost = new bool[k * nrow];
    other_up_ghost    = new bool[k * (ncol+2*k)];
    other_down_ghost  = new bool[k * (ncol+2*k)];
    my_left_ghost  = new bool[k * nrow];
    my_right_ghost = new bool[k * nrow];
    my_up_ghost    = new bool[k * (ncol+2*k)];
    my_down_ghost  = new bool[k * (ncol+2*k)];

    for (int i = 0; i < k*nrow; i++) {
            my_left_ghost[i] = 0;
            my_right_ghost[i] = 0;
    }
    for (int i = 0; i < k*(ncol+2*k); i++) {
        my_up_ghost[i] = 0;
        my_down_ghost[i] = 0;
    }
}

CartesianCellBoard::~CartesianCellBoard() {
    for (int k = 0; k < 2; k++) {
        for (uint i = 0; i < m_cell_size[0] + 2 * m_ghost_size; i++)
            delete[] m_cell_board[k][i];
        delete[] m_cell_board[k];
    }

    delete[] other_left_ghost;
    delete[] other_right_ghost;
    delete[] other_up_ghost;
    delete[] other_down_ghost;
    delete[] my_left_ghost;
    delete[] my_right_ghost;
    delete[] my_up_ghost;
    delete[] my_down_ghost;
}

void CartesianCellBoard::cycle() {
    // std::cout << "Process ID " << m_proc_id << " : cycle() = " << m_cycle << "\n";

    int nrow = m_cell_size[0];
    int ncol = m_cell_size[1];
    int k = m_ghost_size;
    uint _cur = m_cycle % 2;
    uint _new = (_cur + 1) % 2;
    
    if (is_update_cycle()) {
        // Tag as current cycle, to synchronize well.
        int tag = m_cycle;

        // Send-Receive Horizontally First
        for (int i = k; i < nrow+k; i++) {
            for (int j = ncol; j < ncol+k; j++)
                other_left_ghost[(i-k)*k+(j-ncol)] = m_cell_board[_cur][i][j];
            for (int j = k; j < 2*k; j++)
                other_right_ghost[(i-k)*k+(j-k)] = m_cell_board[_cur][i][j];
        }
        int left_proc, right_proc;
        MPI_Cart_shift(*m_cartesian_topology, 1, 1, &left_proc, &right_proc);
        // std::cout << "Process ID " << m_proc_id << " : (LeftProc, RightProc) = (" << left_proc << ", " << right_proc << ")\n";
        MPI_Sendrecv(other_left_ghost, k*nrow, MPI_C_BOOL, right_proc, tag, // Send my right-border to right-proc's left-ghost
                     my_left_ghost,    k*nrow, MPI_C_BOOL, left_proc,  tag, // Receive my left-ghost from left-proc's right-border
                     *m_cartesian_topology, MPI_STATUS_IGNORE);
        MPI_Sendrecv(other_right_ghost, k*nrow, MPI_C_BOOL, left_proc,  tag,
                     my_right_ghost,    k*nrow, MPI_C_BOOL, right_proc, tag,
                     *m_cartesian_topology, MPI_STATUS_IGNORE);
        for (int i = k; i < nrow+k; i++) {
            for (int j = 0; j < k; j++)
                m_cell_board[_cur][i][j] = my_left_ghost[(i-k)*k+j];
            for (int j = ncol+k; j < ncol+2*k; j++)
                m_cell_board[_cur][i][j] = my_right_ghost[(i-k)*k+(j-ncol-k)];
        }

        // Send-Receive Vertically
        for (int j = 0; j < ncol+2*k; j++) {
            for (int i = k; i < 2*k; i++)
                other_down_ghost[(i-k)*(ncol+2*k)+j] = m_cell_board[_cur][i][j];
            for (int i = nrow; i < nrow+k; i++)
                other_up_ghost[(i-nrow)*(ncol+2*k)+j] = m_cell_board[_cur][i][j];
        }
        int up_proc, down_proc;
        MPI_Cart_shift(*m_cartesian_topology, 0, 1, &up_proc, &down_proc);
        // std::cout << "Process ID " << m_proc_id << " : (UpProc, DownProc) = (" << up_proc << ", " << down_proc << ")\n";
        MPI_Sendrecv(other_up_ghost, k*(ncol+2*k), MPI_C_BOOL, down_proc, tag,
                     my_up_ghost,    k*(ncol+2*k), MPI_C_BOOL, up_proc,   tag, 
                     *m_cartesian_topology, MPI_STATUS_IGNORE);
        MPI_Sendrecv(other_down_ghost, k*(ncol+2*k), MPI_C_BOOL, up_proc,   tag,
                     my_down_ghost,    k*(ncol+2*k), MPI_C_BOOL, down_proc, tag,
                     *m_cartesian_topology, MPI_STATUS_IGNORE);
        for (int j = 0; j < ncol+2*k; j++) {
            for (int i = 0; i < k; i++)
                m_cell_board[_cur][i][j] = my_up_ghost[i*(ncol+2*k)+j];
            for (int i = nrow+k; i < nrow+2*k; i++)
                m_cell_board[_cur][i][j] = my_down_ghost[(i-nrow-k)*(ncol+2*k)+j];
        }
    }

    for (int i = 0; i < nrow + 2*k; i++) {
        for (int j = 0; j < ncol + 2*k; j++) {
            // Do not include out-of-border (in global sense)
            if (is_out_of_border(i, j))
                continue;
            
            uint alive_neighbor_cell = 0;
            for (int p = -1; p <= 1; p++) {
                for (int q = -1; q <= 1; q++) {
                    // Out-of-bound
                    if (!is_valid(i+p, j+q))
                        continue;
                    // Self
                    if (p == 0 && q == 0)
                        continue;
                    // Count alive neighbor cells
                    if (m_cell_board[_cur][i+p][j+q])
                        alive_neighbor_cell++;
                }
            }
            if (alive_neighbor_cell < 2 || alive_neighbor_cell > 3)
                m_cell_board[_new][i][j] = false;
            else if (alive_neighbor_cell == 2)
                m_cell_board[_new][i][j] = m_cell_board[_cur][i][j];
            else if (alive_neighbor_cell == 3)
                m_cell_board[_new][i][j] = true;
        }
    }

    m_cycle++;
}

#endif

#ifdef ROWWISE_PARTITION
#endif

#endif
} // namespace Nix