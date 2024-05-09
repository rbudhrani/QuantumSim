#ifndef QUBITLAYER_H
#define QUBITLAYER_H
#include <bitset>
#include "definitions.hpp"

struct qProb
{
    std::bitset<maxQubits> state;
    precision prob;
};

class QubitLayer
{
public:
    QubitLayer(unsigned int numQubits, qubitLayer *qL = nullptr);
    ~QubitLayer();
    void applyPauliX(int target);
    void applyPauliY(int target);
    void applyPauliZ(int target);
    void applyHadamard(int target);
    void applyRx(int target, precision theta);
    void applyRy(int target, precision theta);
    void applyRz(int target, precision theta);
    void applyCnot(int control, int target);
    void applyToffoli(int control1, int control2, int target);
    void applyMcnot(int *controls, int numControls, int target);
    void applyCz(int control, int target);
    void applyMcphase(int *controls, int numControls, int target);
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
    unsigned int _execTime;
    qubitLayer *_qEven;
    qubitLayer *_qOdd;
};

#endif