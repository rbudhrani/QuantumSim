# QuantumSim

![Generic badge](https://img.shields.io/badge/Version-0.10-blue.svg)
[![Generic badge](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

QuantumSim is a simulator for quantum circuits.

## Setting up

Clone the repository
```zsh 
git clone https://github.com/rbudhrani/QuantumSim
```

Run the tests from the root folder
```zsh 
cd <path-to-repository>/QuantumSim && make check
```
___
## Usage
The `QubitLayer` class in `src/QubitLayer.cpp` defines the functions for the quantum gates and the arrays that store the amplitudes for the qubits. The table below lists the functions that can be used by a `QubitLayer` object. `precision` is a `typedef` for `double`. The supported gates are:
| Quantum Gate                  | Function                                                     |
| ------------------------------|--------------------------------------------------------------|
| Pauli X                       | `applyPauliX(int target)`                                    |
| Pauli Y                       | `applyPauliY(int target)`                                    |
| Pauli Z                       | `applyPauliZ(int target)`                                    |
| Hadamard                      | `applyHadamard(int target)`                                  |
| Rotation X                    | `applyRx(int target, precision theta)`                       |
| Rotation Y                    | `applyRy(int target, precision theta)`                       |
| Rotation Z                    | `applyRz(int target, precision theta)`                       |
| Controlled NOT                | `applyCnot(int control, int target)`                         |
| Toffoli                       | `applyToffoli(int control1, int control2, int target)`       |
| Multiple controlled CNOT      | `applyMcnot(int *controls, int numControls, int target)`     |
| Controlled Z                  | `applyCz(int control, int target)`                           |
| Multiple controlled Z         | `applyMcz(int *controls, int numControls, int target)`       | 
___
## Example

Running a quantum circuit is done by writing the circuit in `src/main.cpp`.
```cpp 
#include <iostream>
#include <complex>
#include "QubitLayer.hpp"
#include "qAlgorithms.hpp"

int main(){
    // create a QubitLayer object
    QubitLayer q;
    // apply a hadamard gate
    q.applyHadamard(0);
    // apply cnot on the 2 qubits
    q.applyCnot(0, 1);
    // print qubits
    q.printQubits();
}
```
This circuit creates an EPR pair.

By default when a `QubitLayer` object is created, it is initialised with 2 qubits. To setup a larger quantum circuit, a `QubitLayer` object can be constructed using arguments.
```cpp 
#include <iostream>
#include <complex>
#include "QubitLayer.hpp"
#include "qAlgorithms.hpp"

int main(){
    // create a QubitLayer object with 4 qubits
    QubitLayer q(4);
    // apply hadamard gates
    q.applyHadamard(0);
    q.applyHadamard(2);
    // apply cnot gates
    q.applyCnot(0, 1);
    q.applyCnot(2, 3);
    // print qubits
    q.printQubits();
}
```
___
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.