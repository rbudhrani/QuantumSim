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
     * @brief Applies Pauli X gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     */
    void applyPauliX(int target);
    /**
     * @brief Applies Pauli Y gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     */
    void applyPauliY(int target);
    /**
     * @brief Applies Pauli Z gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     */
    void applyPauliZ(int target);
    /**
     * @brief Applies Hadamard gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     */
    void applyHadamard(int target);
    /**
     * @brief Applies Rx gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     * @param[in] theta angle to rotate qubit along X axis.
     */
    void applyRx(int target, precision theta);
    /**
     * @brief Applies Ry gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     * @param[in] theta angle to rotate qubit along Y axis.
     */
    void applyRy(int target, precision theta);
    /**
     * @brief Applies Rz gate on a qubit.
     *
     * @param[in] target qubit to apply the gate to.
     * @param[in] theta angle to rotate qubit along Z axis.
     */
    void applyRz(int target, precision theta);
    /**
     * @brief Applies CNOT gate on a qubit given the control qubit.
     *
     * @param[in] control control qubit.
     * @param[in] target target qubit.
     */
    void applyCnot(int control, int target);
    /**
     * @brief Applies Toffoli gate on a qubit given the control qubits.
     *
     * @param[in] control1 control qubit 1.
     * @param[in] control1 control qubit 2.
     * @param[in] target target qubit.
     */
    void applyToffoli(int control1, int control2, int target);
    /**
     * @brief Applies MCNOT gate on a qubit given the control qubits.
     *
     * @param[in] controls pointer to array of control qubits.
     * @param[in] numControls number of control qubits.
     * @param[in] target target qubit.
     */
    void applyMcnot(int *controls, int numControls, int target);
    /**
     * @brief Applies CZ gate on a qubit given the control qubit.
     *
     * @param[in] control control qubit.
     * @param[in] target target qubit.
     */
    void applyCz(int control, int target);
    /**
     * @brief Applies MCPHASE gate on a qubit given the control qubit.
     *
     * @param[in] controls pointer to array of control qubits.
     * @param[in] numControls number of control qubits.
     * @param[in] target target qubit.
     */
    void applyMcphase(int *controls, int numControls, int target);
    /**
     * @brief Gets the maximum amplitude of the QubitLayer.
     */
    qProb getMaxAmplitude();
    void printMeasurement();
    void printQubits();
    qubitLayer *getQubitLayerEven();
    qubitLayer *getQubitLayerOdd();
    unsigned long long int getNumStates();
    unsigned int getNumQubits();

private:
    bool checkControls(int *controls, int numControls, std::bitset<maxQubits> state);
    bool checkZeroState(int qubit);
    void updateLayer();
    void toggleParity();
    unsigned int _numQubits;
    unsigned long long int _numStates;
    bool _parity = true;
    qubitLayer *_qEven;
    qubitLayer *_qOdd;
};

#endif