# Free-Majorana-Model – Ground State Energy via NPA

Code for calculating the ground state energy for the free-Majorana-model with all-to-all interactions. 

## Physics
The Hamiltonian at hand takes following form:
$$H = \sum_{i<j} i  c_i c_j$$.

In the code the Hamiltonian is constructed from the moment matrix (`Solve.jl`). Different Hamiltonians can be constructed by extracting different terms form this moment matrix.

## Method
The ground state energy is solved for using NPA hierarchy.

Formalizing the model: `JuMP`
Solving the problem : `SCS` (splitting cone solver)

## Scripts

| File | Description |
|------|-------------|
| `Algebra.jl` | Defines the algebra of the Majoranas and words and operators constructed form these Majoranas. |
| `Moment-Matrix.jl` | Containts functions to build the moment-matrix and its entries. |
| `Operators.jl` | Defines Majoranas, Words and Operators as objects. |
| `SDP.jl` | Builds the SDP. |
| `Solve.jl` | Here you define the system size and degree of the NPA. The Hamiltonian is constructed (i.e. the operator to minimize). |

To solve the problem, you run the `Solve.jl` file. In order for it to work, all other files must be included, i.e. run the file in an environment containing the other files.
## Dependencies
- `JuMP.jl`
- `Combinatorics.jl`
- `SCS.jl`
- `MathOptInterface.jl`
