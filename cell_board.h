#ifndef _CELL_BOARD_H_
#define _CELL_BOARD_H_

#include <iostream>
#include <vector>
#ifdef BUILD_WITH_OPENMPI
#include <mpi.h>
#endif

using uint = unsigned;
// typedef struct {
//     int row;
//     int col;
// } Coord;

namespace Nix {

class CellBoard {
public:
    CellBoard() = delete;
    CellBoard(uint nrow, uint ncol, char **init_cell_board);
    virtual ~CellBoard();

    virtual void cycle();
    virtual void print();

// private:
protected:
    uint m_cell_size[2];    // Row, Col
    bool **m_cell_board[2]; // Visible Board Status
//  std::vector<Coord> m_alive_cells; // Alive cells outside the board
    uint m_cycle;

};

#ifdef BUILD_WITH_OPENMPI

#ifdef CARTESIAN_PARTITION
class CartesianCellBoard /* : public CellBoard */ {
public:
    CartesianCellBoard() = delete;
    CartesianCellBoard(uint proc_id,
                       uint cell_base[], 
                       uint cell_size[],
                       uint board_size[],
                       uint ghost_size,
                       char **init_cell_board,
                       MPI_Comm *cartesian_topology);
    ~CartesianCellBoard();

    void cycle();
    void print() {}
    bool **get_board_status() { return m_cell_board[(m_cycle % 2)]; }
    bool **get_board_status(int gen) { return m_cell_board[gen % 2]; }

    // void get_horizontal_ghosts(bool *left_ghost, bool *right_ghost);
    // void get_vertical_ghosts(bool *up_ghost, bool *down_ghost);
    
private:
    uint m_cell_base[2];
    uint m_cell_size[2];    // Row, Col
    uint m_board_size[2];   // Global Board Size
    uint m_ghost_size;
    bool **m_cell_board[2]; // Visible Board Status
    uint m_cycle;

    uint m_proc_id;
    MPI_Comm *m_cartesian_topology;

    bool is_valid(int i, int j) {
        return (i >= 0 && i < m_cell_size[0] + 2 * m_ghost_size)
            && (j >= 0 && j < m_cell_size[1] + 2 * m_ghost_size); 
    }
    bool is_ghost(uint i, uint j) {
        // assert(is_valid(i, j));
        return (i < m_ghost_size || i >= (m_cell_size[0] + m_ghost_size))
            || (j < m_ghost_size || j >= (m_cell_size[1] + m_ghost_size));
    }
    bool is_out_of_border(int i, int j) {
        int i_glob = i - m_ghost_size + m_cell_base[0];
        int j_glob = j - m_ghost_size + m_cell_base[1];
        return (i_glob < 0 || i_glob >= m_board_size[0])
            || (j_glob < 0 || j_glob >= m_board_size[1]);
    }
    bool is_update_cycle() {
        // return false;
        return (m_ghost_size == 1 || m_cycle % m_ghost_size == 0);
    }

    bool *other_left_ghost;
    bool *other_right_ghost;
    bool *other_up_ghost;
    bool *other_down_ghost;
    bool *my_left_ghost;
    bool *my_right_ghost;
    bool *my_up_ghost;
    bool *my_down_ghost;    
};
#endif

#ifdef ROWWISE_PARTITION
class RowwiseCellBoard : public CellBoard {
public:

private:
    MPI_Comm *linear_topology;
};
#endif

#endif

} // namespace Nix

#endif 