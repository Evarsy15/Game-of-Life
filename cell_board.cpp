
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
                          << init_cell_board[i][j] << "at (" << i-1 << ", " << j-1 << ") "  
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
            delete m_cell_board[k][i];
        delete m_cell_board[k];
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
void CartesianCellBoard::cycle() {

}

#endif

#ifdef ROWWISE_PARTITION
#endif
#endif
} // namespace Nix