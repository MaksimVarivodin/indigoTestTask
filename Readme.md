
# SecureBox Solver

This project implements algorithms to solve the SecureBox puzzle. The objective of the puzzle is to unlock a box represented as a two-dimensional grid of boolean values (where `true` indicates a locked cell and `false` indicates an unlocked cell). The box is unlocked when all cells become `false`.

## Overview

The project features several methods to solve the SecureBox:

- **One-Toggle Solve:** Tries to unlock the box by toggling each cell individually.
- **Odd Row-Column Sums:** Computes a toggle matrix based on the parity (XOR) of row and column sums. This method is used when the sum of the dimensions is even.
- **GF(2) System Solver:** For cases where the sum of the dimensions is odd, the system of equations is solved over GF(2) using Gaussian elimination. This method constructs a coefficient matrix, flattens the box state, and finds the appropriate toggle operations by solving a linear system modulo 2.

## Files

- **main.cpp:**  
  Contains the definition and implementation of the `SecureBox` class and its public API (`toggle`, `isLocked`, `getState`, etc.).

- **Solver Functions:**  
  Implementations for:
  - `oneToggleSolve`
  - `oddRowColumnSums`
  - `buildToggleMatrix`
  - `flattenState`
  - `solveGF2`
  - `openBox` (the main function to solve/unlock the box)
## Building

To compile the project, use your preferred C++ compiler. For example, with `g++`:

```bash
g++ -std=c++17 -O2 -o SecureBoxSolver main.cpp SecureBox.cpp Matrix.cpp
```

Make sure to include all source files and headers in the compilation process.

## Running

The project executable can be run from the command line. The program accepts two arguments: the number of rows and the number of columns of the SecureBox. For example:

```bash
./testTask 6 5
```

The solver will attempt to unlock the box and print the result. If the box is unlocked successfully, it will output:

```
BOX: OPENED!
```

Otherwise, it prints an error or that the box remains locked.

```
BOX: LOCKED!
```

## Notes

- The project includes multiple strategies to solve the SecureBox puzzle. The GF(2) solver is utilized when the sum of the dimensions is odd to handle cases where a simple parity-based method may fail.


## License

This project is provided under the MIT License.