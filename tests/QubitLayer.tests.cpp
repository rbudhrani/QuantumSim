#include <iostream>
#include <cassert>
#include <complex>
#include <stdexcept>
#include "../src/QubitLayer.hpp"

precision cosHalfPi = cos(pi / 2);
precision sinHalfPi = sin(pi / 2);

// Custom assertion macro
#define CUSTOM_ASSERT(condition, message)                                             \
    do                                                                                \
    {                                                                                 \
        if (!(condition))                                                             \
        {                                                                             \
            std::cerr << "Assertion failed: " << (message) << "\n"                    \
                      << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl; \
            throw std::runtime_error(message);                                        \
        }                                                                             \
    } while (false)

// Helper function to run a test and print the result
void runTest(const std::string &testName, const std::function<void()> &testFunction)
{
    try
    {
        testFunction();
        std::cout << testName << "   \033[32;32m[PASSED]\033[m" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << testName << "   \033[31;31m[FAILED]\033[m " << e.what() << std::endl;
    }
}

// Helper function to assert the state of a QubitLayer
void assertQubitLayerState(QubitLayer &q, const std::complex<precision> *expectedState, unsigned long long int numStates)
{
    for (unsigned long long int i = 0; i < numStates; ++i)
    {
        CUSTOM_ASSERT(*(q.getQubitLayerOdd() + i) == expectedState[i], "State mismatch at index " + std::to_string(i));
    }
}

// Test: Pauli-X gate flips qubit
void testPauliX()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyPauliX(0);

    // Assert
    std::complex<precision> expectedState[2] = {{0, 0}, {1, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Pauli-Y gate flips and applies phase on state |0>
void testPauliYState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyPauliY(0);

    // Assert
    std::complex<precision> expectedState[2] = {{0, 0}, {0, 1}};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Pauli-Y gate flips and applies phase on state |1>
void testPauliYState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyPauliY(0);

    // Assert
    std::complex<precision> expectedState[2] = {{0, -1}, {0, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Pauli-Z gate does nothing to state |0>
void testPauliZState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyPauliZ(0);

    // Assert
    std::complex<precision> expectedState[2] = {{1, 0}, {0, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Pauli-Z gate applies phase to state |1>
void testPauliZState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyPauliZ(0);

    // Assert
    std::complex<precision> expectedState[2] = {{0, 0}, {-1, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Hadamard gate bring state |0> to superposition
void testHadamardState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyHadamard(0);

    // Assert
    std::complex<precision> expectedState[2] = {hadamardCoef, hadamardCoef};
    assertQubitLayerState(q, expectedState, 2);
}

// Test: Hadamard gate bring state |1> to superposition
void testHadamardState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyHadamard(0);

    // Assert
    std::complex<precision> expectedState[2] = {hadamardCoef, -hadamardCoef};
    assertQubitLayerState(q, expectedState, 2);
}

void testRxState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyRx(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{cosHalfPi, 0}, {0, -sinHalfPi}};
    assertQubitLayerState(q, expectedState, 2);
}

void testRxState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyRx(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{0, -sinHalfPi}, {cosHalfPi, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

void testRyState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyRy(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{cosHalfPi, 0}, {sinHalfPi, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

void testRyState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyRy(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{-sinHalfPi, 0}, {cosHalfPi, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

void testRzState0()
{
    // Arrange
    QubitLayer q(1);

    // Act
    q.applyRz(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{cosHalfPi, -sinHalfPi}, {0, 0}};
    assertQubitLayerState(q, expectedState, 2);
}

void testRzState1()
{
    // Arrange
    std::complex<precision> input[2] = {{0, 0}, {1, 0}};
    QubitLayer q(1, input);

    // Act
    q.applyRz(0, pi);

    // Assert
    std::complex<precision> expectedState[2] = {{0, 0}, {cosHalfPi, sinHalfPi}};
    assertQubitLayerState(q, expectedState, 2);
}

void testCnotState00()
{
    // Arrange
    QubitLayer q(2);

    // Act
    q.applyCnot(0, 1);

    // Assert
    std::complex<precision> expectedState[4] = {{1, 0}, {0, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expectedState, 4);
}

void testCnotState01()
{
    // Arrange
    std::complex<precision> input[4] = {{0, 0}, {0, 0}, {1, 0}, {0, 0}};
    QubitLayer q(2, input);

    // Act
    q.applyCnot(0, 1);

    // Assert
    std::complex<precision> expectedState[4] = {{0, 0}, {0, 0}, {1, 0}, {0, 0}};
    assertQubitLayerState(q, expectedState, 4);
}

void testCnotState10()
{
    // Arrange
    std::complex<precision> input[4] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}};
    QubitLayer q(2, input);

    // Act
    q.applyCnot(0, 1);

    // Assert
    std::complex<precision> expectedState[4] = {{0, 0}, {0, 0}, {0, 0}, {1, 0}};
    assertQubitLayerState(q, expectedState, 4);
}

void testCnotState11()
{
    // Arrange
    std::complex<precision> input[4] = {{0, 0}, {0, 0}, {0, 0}, {1, 0}};
    QubitLayer q(2, input);

    // Act
    q.applyCnot(0, 1);

    // Assert
    std::complex<precision> expectedState[4] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expectedState, 4);
}

void testToffoliState000()
{
    // Arrange
    QubitLayer q(3);

    // Act
    q.applyToffoli(0, 1, 2);

    // Assert
    std::complex<precision> expectedState[8] = {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expectedState, 8);
}

void testToffoliState100()
{
    // Arrange
    std::complex<precision> input[8] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}};
    QubitLayer q(3, input);

    // Act
    q.applyToffoli(0, 1, 2);

    // Assert
    assertQubitLayerState(q, input, 8);
}

void testToffoliState010()
{
    // Arrange
    std::complex<precision> input[8] = {{0, 0}, {0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    QubitLayer q(3, input);

    // Act
    q.applyToffoli(0, 1, 2);

    // Assert
    assertQubitLayerState(q, input, 8);
}

void testToffoliState110()
{
    // Arrange
    std::complex<precision> input[8] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}, {0, 0}};
    QubitLayer q(3, input);

    // Act
    q.applyToffoli(0, 1, 2);

    // Assert
    std::complex<precision> expected[8] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {1, 0}};
    assertQubitLayerState(q, expected, 8);
}

void testMcnotState0000()
{
    // Arrange
    std::complex<precision> input[16] = {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    QubitLayer q(4, input);
    int controls[3] = {3, 2, 1};

    // Act
    q.applyMcnot(controls, 3, 0);

    // Assert
    std::complex<precision> expected[16] = {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expected, 16);
}

void testMcnotState0100()
{
    // Arrange
    std::complex<precision> input[16] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    QubitLayer q(4, input);
    int controls[3] = {3, 2, 1};

    // Act
    q.applyMcnot(controls, 3, 0);

    // Assert
    std::complex<precision> expected[16] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expected, 16);
}

void testMcnotState1110()
{
    // Arrange
    std::complex<precision> input[16] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    QubitLayer q(4, input);
    int controls[3] = {3, 2, 1};

    // Act
    q.applyMcnot(controls, 3, 0);

    // Assert
    std::complex<precision> expected[16] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    assertQubitLayerState(q, expected, 16);
}

// Main function to run all tests
int main()
{
    std::cout << "\033[34;34m===========Test Results===========\033[m" << std::endl;

    runTest("Pauli-X Gate flips state                         ", testPauliX);
    runTest("Pauli-Y Gate flips and applies phase on state |0>", testPauliYState0);
    runTest("Pauli-Y Gate flips and applies phase on state |1>", testPauliYState1);
    runTest("Pauli-Z Gate does nothing to state |0>           ", testPauliZState0);
    runTest("Pauli-Z Gate applies phase to state |1>          ", testPauliZState1);
    runTest("Hadamard Gate puts |0> into superposition        ", testHadamardState0);
    runTest("Hadamard Gate puts |1> into superposition        ", testHadamardState1);
    runTest("Rx Gate rotates |0> by π about X                 ", testRxState0);
    runTest("Rx Gate rotates |1> by π about X                 ", testRxState1);
    runTest("Ry Gate rotates |0> by π about Y                 ", testRyState0);
    runTest("Ry Gate rotates |1> by π about Y                 ", testRyState1);
    runTest("Rz Gate rotates |0> by π about Z                 ", testRzState0);
    runTest("Rz Gate rotates |1> by π about Z                 ", testRzState1);
    runTest("CNOT Gate does nothing to state |00>             ", testCnotState00);
    runTest("CNOT Gate does nothing to state |01>             ", testCnotState01);
    runTest("CNOT Gate |10> -> |11>                           ", testCnotState10);
    runTest("CNOT Gate |11> -> |10>                           ", testCnotState11);
    runTest("Toffoli Gate |000> -> |000>                      ", testToffoliState000);
    runTest("Toffoli Gate |100> -> |100>                      ", testToffoliState100);
    runTest("Toffoli Gate |010> -> |010>                      ", testToffoliState010);
    runTest("Toffoli Gate |110> -> |111>                      ", testToffoliState110);
    runTest("MCNOT Gate |0000> -> |0000>                      ", testMcnotState0000);
    runTest("MCNOT Gate |0100> -> |0100>                      ", testMcnotState0100);

    std::cout << "\033[34;34m==================================\033[m" << std::endl;

    return 0;
}