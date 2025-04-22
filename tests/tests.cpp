#include <iostream>
#include <cassert>
#include "../src/QubitLayer.hpp"
#include "tests.hpp"

// list of quantum gates
enum Gates
{
    X,
    Y,
    Z,
    rx,
    ry,
    rz,
    cnot,
    cphase,
    toffoli,
    mcnot,
    Hadamard,
    mcphase
};

bool testGate(Gates gate)
{
    // initalise state to |111> so all gates can change the state
    qubitLayer testInput[8] = {zeroComplex, zeroComplex, zeroComplex, zeroComplex, zeroComplex, zeroComplex, zeroComplex, {1, 0}};
    QubitLayer q = QubitLayer(3, testInput);
    // define array of control qubits for multiple control qubit gates
    int ctrlQubits[2]{0, 1};
    TesterGate testerGate;
    switch (gate)
    {
    case X:
        q.applyPauliX(0);
        testerGate = pauliXtester;
        break;
    case Y:
        q.applyPauliY(0);
        testerGate = pauliYtester;
        break;
    case Z:
        q.applyPauliZ(0);
        testerGate = pauliZtester;
        break;
    case Hadamard:
        q.applyHadamard(0);
        testerGate = hadamardTester;
        break;
    case rx:
        q.applyRx(0, pi);
        testerGate = rxTester;
        break;
    case ry:
        q.applyRy(0, pi);
        testerGate = ryTester;
        break;
    case rz:
        q.applyRz(0, pi);
        testerGate = rzTester;
        break;
    case cnot:
        q.applyCnot(0, 1);
        testerGate = cnotTester;
        break;
    case cphase:
        q.applyCz(0, 1);
        testerGate = czTester;
        break;
    case toffoli:
        q.applyToffoli(0, 1, 2);
        testerGate = toffoliTester;
        break;
    case mcnot:
        q.applyMcnot(ctrlQubits, 2, 2);
        testerGate = mcnotTester;
        break;
    case mcphase:
        q.applyMcphase(ctrlQubits, 2, 2);
        testerGate = mcphaseTester;
        break;
    default:
        q.applyPauliX(0);
        break;
    }
    // define variable to store the result of the test
    bool testResult = true;
    // iterate over the states to check if the gates work
    for (unsigned long long int i = 0; i < q.getNumStates(); i++)
        testResult = *(q.getQubitLayerOdd() + i) == testerGate.outputState[i] && testResult;
    std::cout << testerGate.gateName << (testResult ? " \033[32;32m[PASSED]\033[m" : " \033[31;31m[FAILED]\033[m") << std::endl;
    return testResult;
}

int main(int argc, char *argv[])
{
    // define variable to store result of the tests
    bool testResult = true;
    std::cout << "\033[34;34m===========Test Results===========\033[m" << std::endl;
    for (int gate = X; gate <= mcphase; gate++)
        testResult = testGate(static_cast<Gates>(gate)) && testResult;
    return testResult ? EXIT_SUCCESS : EXIT_FAILURE;
}