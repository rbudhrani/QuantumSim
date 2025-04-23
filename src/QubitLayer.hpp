#ifndef QUBITLAYER_H
#define QUBITLAYER_H
#include <bitset>
#include "definitions.hpp"

struct qProb
{
    std::bitset<maxQubits> state;
    precision prob;
};

/**
 * @class QubitLayer
 * @brief Represents a QubitLayer, which is a pair of qubit states, the input and the ouptut.
 */
class QubitLayer
{
public:
    /**
     * @brief Default constructor.
     * Initialises the qOdd layer to the state |0...0> and the qEven layer to the state |0...01> unless an input is provided.
     *
     * @param[in] numQubits number of qubits.
     * @param[in] qL optional pointer to the array to initialise qEven to.
     */
    QubitLayer(unsigned int numQubits, qubitLayer *qL = nullptr);

    /**
     * @brief Default destructor.
     * Deletes the memory allocated for qEven and qOdd.
     */
    ~QubitLayer();

    /**
     * @brief Applies the Pauli X gate to a specific qubit.
     *
     * The Pauli X gate flips the state of the target qubit (|0> ↔ |1>).
     *
     * @param[in] target The index of the qubit to apply the gate to.
     */
    void applyPauliX(int target);

    /**
     * @brief Applies the Pauli Y gate to a specific qubit.
     *
     * The Pauli Y gate flips the state of the target qubit and adds a phase (|0> → i|1>, |1> → -i|0>).
     *
     * @param[in] target The index of the qubit to apply the gate to.
     */
    void applyPauliY(int target);
    /**
     * @brief Applies the Pauli Z gate to a specific qubit.
     *
     * The Pauli Z gate adds a phase to |1> and does nothing to |0> (|1> → -|1>).
     *
     * @param[in] target The index of the qubit to apply the gate to.
     */

    void applyPauliZ(int target);

    /**
     * @brief Applies the Hadamard gate to a specific qubit.
     *
     * The Hadamard gate puts a qubit into superposition (|0> → 1/√2(|0> + |1>), |1> → 1/√2(|0> - |1>)).
     *
     * @param[in] target The index of the qubit to apply the gate to.
     */
    void applyHadamard(int target);

    /**
     * @brief Applies the Rx rotation gate to a specific qubit.
     *
     * Rotates the state of the target qubit around the X-axis by the specified angle.
     *
     * @param[in] target The index of the qubit to apply the gate to.
     * @param[in] theta The rotation angle in radians.
     */
    void applyRx(int target, precision theta);

    /**
     * @brief Applies the Ry rotation gate to a specific qubit.
     *
     * Rotates the state of the target qubit around the Y-axis by the specified angle.
     *
     * @param[in] target The index of the qubit to apply the gate to.
     * @param[in] theta The rotation angle in radians.
     */
    void applyRy(int target, precision theta);

    /**
     * @brief Applies the Rz rotation gate to a specific qubit.
     *
     * Rotates the state of the target qubit around the Z-axis by the specified angle.
     *
     * @param[in] target The index of the qubit to apply the gate to.
     * @param[in] theta The rotation angle in radians.
     */
    void applyRz(int target, precision theta);

    /**
     * @brief Applies the CNOT gate on a qubit given the control qubit.
     *
     * @param[in] control The control qubit.
     * @param[in] target The target qubit.
     */
    void applyCnot(int control, int target);

    /**
     * @brief Applies the Toffoli gate on a qubit given the control qubits.
     *
     * @param[in] control1 The first control qubit.
     * @param[in] control2 The second control qubit.
     * @param[in] target The target qubit.
     */
    void applyToffoli(int control1, int control2, int target);

    /**
     * @brief Applies the multi controlled CNOT gate on a qubit given the control qubits.
     *
     * @param[in] controls Pointer to the array of control qubits.
     * @param[in] numControls The number of control qubits.
     * @param[in] target The target qubit.
     */
    void applyMcnot(int *controls, int numControls, int target);

    /**
     * @brief Applies the controlled phase/CZ gate on a qubit given the control qubit.
     *
     * @param[in] control The control qubit.
     * @param[in] target The target qubit.
     */
    void applyCz(int control, int target);

    /**
     * @brief Applies the multi controlled phase/CZ gate on a qubit given the control qubits.
     *
     * @param[in] controls Pointer to the array of control qubits.
     * @param[in] numControls The number of control qubits.
     * @param[in] target The target qubit.
     */
    void applyMcphase(int *controls, int numControls, int target);

    /**
     * @brief Gets the maximum amplitude of the QubitLayer.
     *
     * @return A `qProb` structure containing the state and its probability.
     */
    qProb getMaxAmplitude();

    /**
     * @brief Prints the measurement result, which is the state with the highest probability.
     */
    void printMeasurement();

    /**
     * @brief Prints the state of the qubits in the odd layer and even layer.
     */
    void printQubits();

    /**
     * @brief Gets the even-parity qubit layer.
     *
     * @return A pointer to the even qubit layer.
     */
    qubitLayer *getQubitLayerEven();

    /**
     * @brief Gets the odd-parity qubit layer.
     *
     * @return A pointer to the odd qubit layer.
     */
    qubitLayer *getQubitLayerOdd();

    /**
     * @brief Gets the total number of states in the QubitLayer.
     *
     * @return The total number of states (2^numQubits).
     */
    unsigned long long int getNumStates();

    /**
     * @brief Gets the number of qubits in the QubitLayer.
     *
     * @return The number of qubits.
     */
    unsigned int getNumQubits();

private:
    /**
     * @brief Checks if the control qubits are high for a given state.
     *
     * @param[in] controls Pointer to the array of control qubits.
     * @param[in] numControls The number of control qubits.
     * @param[in] state The current state of the qubits.
     * @return True if the control qubits are all high, false otherwise.
     */
    bool checkControls(int *controls, int numControls, std::bitset<maxQubits> state);

    /**
     * @brief Checks if a state has a probability amplitude that is 0.
     *
     * @param[in] state The index of the state to check.
     * @return True if the state has a probability amplitude that is 0, false otherwise.
     */
    bool checkZeroState(unsigned long long int state);

    /**
     * @brief Updates the QubitLayer based on the parity. If the parity is true then the even
     * layer is reset to 0, else the odd layer is reset to 0. Then the parity is inverted to
     * prepare for the next gate.
     */
    void updateLayer();

    /**
     * @brief Toggles the parity of the QubitLayer, so the instance know if the odd layer is
     * the input or output. Same for the even layer.
     */
    void toggleParity();
    unsigned int numQubits_;
    unsigned long long int numStates_;
    qubitLayer *qEven_;
    qubitLayer *qOdd_;
    bool parity = true;
};

#endif