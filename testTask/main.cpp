#include <iostream>
#include <vector>
#include <random>
#include <time.h>

/*
You are given a locked container represented as a two-dimensional grid of boolean values (true = locked, false = unlocked). 
Your task is to write an algorithm that fully unlocks the box, i.e., 
transforms the entire matrix into all false.

Implement the function:
bool openBox(uint32_t y, uint32_t x);
This function should:
    - Use the SecureBox public API (toggle, isLocked, getState).
    - Strategically toggle cells to reach a state where all elements are false.
    - Return true if the box remains locked, false if successfully unlocked.
You are not allowed to path or modify the SecureBox class.

Evaluation Criteria:
    - Functional correctness
    - Computational efficiency
    - Code quality, structure, and comments
    - Algorithmic insight and clarity
*/

class SecureBox
{
private:
    std::vector<std::vector<bool>> box;

public:
    //================================================================================
    // Constructor: SecureBox
    // Description: Initializes the secure box with a given size and 
    //              shuffles its state using a pseudo-random number generator 
    //              seeded with current time.
    //================================================================================
    SecureBox(uint32_t y, uint32_t x): ySize(y), xSize(x)
    {
        rng.seed(time(0));
        box.resize(y);
        for (auto& it : box)
            it.resize(x);
        shuffle();
    }

    //================================================================================
    // Method: toggle
    // Description: Toggles the state at position (x, y) and also all cells in the
    //              same row above and the same column to the left of it.
    //================================================================================
    void toggle(uint32_t y, uint32_t x)
    {
        box[y][x] = !box[y][x];
        for (uint32_t i = 0; i < xSize; i++)
            box[y][i] = !box[y][i];
        for (uint32_t i = 0; i < ySize; i++)
            box[i][x] = !box[i][x];
    }

    //================================================================================
    // Method: isLocked
    // Description: Returns true if any cell 
    //              in the box is true (locked); false otherwise.
    //================================================================================
    bool isLocked()
    {
        for (uint32_t x = 0; x < xSize; x++)
            for (uint32_t y = 0; y < ySize; y++)
                if (box[y][x])
                    return true;

        return false;
    }

    //================================================================================
    // Method: getState
    // Description: Returns a copy of the current state of the box.
    //================================================================================
    std::vector<std::vector<bool>> getState()
    {
        return box;
    }

private:
    std::mt19937_64 rng;
    uint32_t ySize, xSize;

    //================================================================================
    // Method: shuffle
    // Description: Randomly toggles cells in the box to 
    // create an initial locked state.
    //================================================================================
    void shuffle()
    {
        for (uint32_t t = rng() % 1000; t > 0; t--)
            toggle(rng() % ySize, rng() % xSize);
    }
};


/**
 * @brief Attempts to solve the SecureBox by toggling one cell at a time.
 * 
 * This function iterates over each cell of the box and toggles it. If after toggling
 * the box becomes unlocked (i.e. isLocked() returns false), then it returns true.
 * Otherwise, it reverts the toggle and continues searching.
 *
 * @param box The SecureBox object.
 * @param y The number of rows.
 * @param x The number of columns.
 * @return true if toggling a single cell can unlock the box, false otherwise.
 */
bool oneToggleSolve(SecureBox& box, uint32_t y, uint32_t x)
{
    if (box.isLocked())
        for (uint32_t i = 0; i < y; i++)
            for (uint32_t j = 0; j < x; j++)
            {
                box.toggle(i, j);
                if (!box.isLocked())
                    return true;
                box.toggle(i, j);
            }
    return false;
}

/**
 * @brief Calculates the parity (0 or 1) for each row and column of the box's state.
 * 
 * For the given box state (where true equals 1, false equals 0), this function calculates
 * the XOR (parity) for all elements in each row and each column.
 * Then it computes the toggle matrix based on the idea:
 * T[i][j] = rowParity[i] XOR colParity[j] XOR matrix[i][j].
 *
 * @param box The SecureBox object.
 * @param y The number of rows.
 * @param x The number of columns.
 * @return A 2D vector representing the computed toggle pattern.
 */
std::vector<std::vector<bool>> oddRowColumnSums(SecureBox& box, uint32_t y, uint32_t x)
{
    // Get the current state of the box
    auto matrix = box.getState();

    // Vectors to store parity for rows and columns
    std::vector<bool> rowParity(y, false);
    std::vector<bool> colParity(x, false);

    // Calculate the parity (XOR) for each row
    for (uint32_t i = 0; i < y; i++)
        for (uint32_t j = 0; j < x; j++)
            rowParity[i] = rowParity[i] ^ matrix[i][j];

    // Calculate the parity (XOR) for each column
    for (uint32_t j = 0; j < x; j++)
        for (uint32_t i = 0; i < y; i++)
            colParity[j] = colParity[j] ^ matrix[i][j];

    // Now compute the toggle matrix.
    // Idea: Suppose T[i][j] indicates if toggling cell (i,j) is needed.
    // We set T[i][j] = rowParity[i] XOR colParity[j] XOR matrix[i][j].
    std::vector<std::vector<bool>> result(y, std::vector<bool>(x, false));
    for (uint32_t i = 0; i < y; i++)
        for (uint32_t j = 0; j < x; j++)
            result[i][j] = rowParity[i] ^ colParity[j] ^ matrix[i][j];

    return result;
}

/**
 * @brief Solves a system of linear equations over GF(2) using Gaussian elimination.
 * 
 * The system is represented as A * x = b, where A is an N x N matrix (elements 0 or 1),
 * and b is a vector of size N (elements 0 or 1). The solution vector is filled with 0s and 1s.
 *
 * @param A The coefficient matrix (size N x N) over GF(2).
 * @param b The right-hand side vector (size N) over GF(2).
 * @param solution The vector that will store the solution.
 * @return true if a solution is found, false if the system is unsolvable.
 */
bool solveGF2(std::vector<std::vector<int>> A, std::vector<int> b, std::vector<int>& solution)
{
    int N = A.size();
    // Forward elimination (transform to row echelon form)
    int row = 0;
    for (int col = 0; col < N && row < N; ++col, row++)
    {
        // Find a row with a 1 in the current column
        int pivot = row;
        while (pivot < N && A[pivot][col] == 0)
            pivot++;
        if (pivot == N)
            continue; // If the entire column is 0, go to the next column

        // Swap the current row with the pivot row
        std::swap(A[row], A[pivot]);
        std::swap(b[row], b[pivot]);

        // Eliminate the current column elements in all other rows
        for (int i = 0; i < N; i++)
            if (i != row && A[i][col] == 1)
            {
                // Add (XOR) the current row to row i over GF(2)
                for (int j = col; j < N; j++)
                    A[i][j] ^= A[row][j];

                b[i] ^= b[row];
            }
    }

    // Check for consistency
    for (int i = row; i < N; i++)
        if (b[i] != 0)
            return false;

    // Back substitution: form the solution vector
    solution.assign(N, 0);
    for (int i = 0; i < N; i++)
    {
        // Find the leading column in row i
        int leading = -1;
        for (int j = 0; j < N; j++)
            if (A[i][j] == 1)
            {
                leading = j;
                break;
            }
        if (leading == -1) continue; // Free equation, skip
        solution[leading] = b[i];
    }
    return true;
}

/**
 * @brief Builds the toggle coefficient matrix for a box of size y x x.
 *
 * The matrix A (of size N x N, where N = y * x) is built such that
 * A[k][l] = 1 if toggling cell l (i.e., (i, j) with l = i * x + j)
 * affects cell k (i.e., (a, b) with k = a * x + b).
 * In other words, A[k][l] = 1 if (a == i || b == j).
 * All operations are considered modulo 2.
 *
 * @param y The number of rows.
 * @param x The number of columns.
 * @return The toggle coefficient matrix.
 */
std::vector<std::vector<int>> buildToggleMatrix(uint32_t y, uint32_t x)
{
    int N = y * x;
    std::vector<std::vector<int>> A(N, std::vector<int>(N, 0));
    for (uint32_t a = 0; a < y; a++)
        for (uint32_t b = 0; b < x; b++)
            for (uint32_t i = 0,
                          eq = a * x + b; // Equation number for cell (a,b)
                 i < y;
                 i++)
                for (uint32_t j = 0; j < x; j++)
                    // If toggling cell (i,j) affects cell (a,b)
                    if (a == i || b == j)
                        A[eq][i * x + j] = 1;

    return A;
}

/**
 * @brief Flattens a 2D boolean matrix into a 1D vector of 0s and 1s.
 *
 * Converts the boolean state matrix into a vector where true is 1 and false is 0.
 *
 * @param state The 2D boolean matrix.
 * @return A 1D vector of integers (0 or 1) representing the state.
 */
std::vector<int> flattenState(const std::vector<std::vector<bool>>& state)
{
    std::vector<int> b;
    for (const auto& row : state)
        for (bool cell : row)
            b.push_back(cell ? 1 : 0);

    return b;
}

//================================================================================
// Function: openBox
// Description: Your task is to implement this function to unlock the SecureBox.
//              Use only the public methods of SecureBox (toggle, getState, isLocked).
//              You must determine the correct sequence of toggle operations to make
//              all values in the box 'false'. The function should return false if
//              the box is successfully unlocked, or true if any cell remains locked.
//================================================================================
bool openBox(uint32_t y, uint32_t x)
{
    if (y <= 0 || x <= 0) return false;
    // Initialize the box (the SecureBox constructor shuffles the state)
    SecureBox box(y, x);

    // Special case for boxes with a single row or column
    if (y == 1 || x == 1)
    {
        if (box.isLocked())
            box.toggle(0, 0);
        return false;
    }

    // Try one-toggle solution first
    if (!oneToggleSolve(box, y, x))
    {
        // For boxes with even (y + x), use oddRowColumnSums strategy
        if ((y + x) % 2 == 0)
        {
            const auto computedOddToggles = oddRowColumnSums(box, y, x);
            for (uint32_t i = 0; i < y; ++i)
                for (uint32_t j = 0; j < x; ++j)
                    if (computedOddToggles[i][j])
                        box.toggle(i, j);
        }
        else
        {
            // For boxes with odd (y + x), use GF(2) system solving
            // Get the initial state of the box
            auto initState = box.getState();
            // Build the coefficient matrix A (size N x N, where N = y * x)
            auto A = buildToggleMatrix(y, x);
            // Flatten the initial state into a vector (0s and 1s)
            auto b = flattenState(initState);

            // Solve the system A * sol = b modulo 2
            std::vector<int> sol;
            if (solveGF2(A, b, sol))
                // Apply toggles according to the solution:
                // For each variable, if sol[index] is 1, toggle the corresponding cell.
                for (uint32_t i = 0; i < y; i++)
                    for (uint32_t j = 0; j < x; j++)
                    {
                        int index = i * x + j;
                        if (sol[index])
                            box.toggle(i, j);
                    }
        }
    }

    return box.isLocked();
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: testTask <rows> <columns>\n";
        return 1;
    }
    uint32_t y = std::atol(argv[1]);
    uint32_t x = std::atol(argv[2]);

    bool state = openBox(y, x);
    if (state)
        std::cout << "BOX: LOCKED!" << std::endl;
    else
        std::cout << "BOX: OPENED!" << std::endl;
    /**/
    return state;
}
