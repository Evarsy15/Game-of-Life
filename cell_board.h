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
    ~CellBoard();

    virtual void cycle();
    virtual void print();

protected:
    uint m_cell_size[2];    // Row, Col
    bool **m_cell_board[2]; // Visible Board Status
//  std::vector<Coord> m_alive_cells; // Alive cells outside the board
    uint m_cycle;

};

#ifdef BUILD_WITH_OPENMPI

#ifdef CARTESIAN_PARTITION
class CartesianCellBoard : public CellBoard {
public:
    CartesianCellBoard() = delete;
    CartesianCellBoard(uint cell_size[], char **init_cell_board)
        : CellBoard(cell_size[0], cell_size[1], init_cell_board) { }
    ~CartesianCellBoard() { }

    void cycle();
    void print() { }
    
private:
    MPI_Comm *cartesian_topology;
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