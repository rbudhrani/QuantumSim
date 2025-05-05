#include <complex>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "QubitLayer.hpp"
#include "TanglrBitset.hpp"

QubitLayer::QubitLayer(unsigned int numQubits, qubitLayer *qL)
{
    // calculate the number of states
    numStates_ = 1;
    for (unsigned int i = 0; i < numQubits; i++)
        numStates_ *= 2;
    // allocate memory for arrays
    qEven_ = new qubitLayer[numStates_];
    qOdd_ = new qubitLayer[numStates_];
    // if input is provided then use that to fill the input qubit state
    if (!(qL == nullptr))
        for (unsigned long long int row = 0; row < numStates_; row++)
            qEven_[row] = qL[row];
    else
        qEven_[0] = {1, 0};
}

QubitLayer::~QubitLayer()
{
    delete[] qEven_;
    delete[] qOdd_;
}

void QubitLayer::updateLayer()
{
    parity ? std::fill(qEven_, qEven_ + numStates_, zeroComplex) : std::fill(qOdd_, qOdd_ + numStates_, zeroComplex);
    toggleParity();
}

bool QubitLayer::checkZeroState(unsigned long long int state)
{
    return parity ? qEven_[state].imag() || qEven_[state].real() : qOdd_[state].imag() || qOdd_[state].real();
}

void QubitLayer::applyPauliX(int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip the target qubit in the state
            state.flip(target);
            updateState(state.to_ullong(), i, 1);
        }
    updateLayer();
}

void QubitLayer::applyPauliY(int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip qubit
            state.flip(target);
            // add phase of -i if bit was 1 else i
            updateState(state.to_ullong(), i, state.test(target) ? complexImg : -complexImg);
        }
    updateLayer();
}

void QubitLayer::applyPauliZ(int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // add phase of -1 if bit is 1
            updateState(i, i, state.test(target) ? -1 : 1);
        }
    updateLayer();
}

void QubitLayer::applyHadamard(int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            updateState(i, i, state.test(target) ? -hadamardCoef : hadamardCoef, false);
            state.flip(target);
            updateState(state.to_ulong(), i, hadamardCoef, false);
        }
    updateLayer();
}

void QubitLayer::applyRx(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
    // map |0> to cosTheta*|0> and |1> to cosTheta*|1>
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            updateState(i, i, cosTheta, false);
            std::bitset<maxQubits> state = i;
            state.flip(target);
            updateState(state.to_ullong(), i, -complexImg * sinTheta, false);
        }
    updateLayer();
}

void QubitLayer::applyRy(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
    // map |1> to cosTheta*|1> and |0> to cosTheta*|0>
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            updateState(i, i, cosTheta, false);
            std::bitset<maxQubits> state = i;
            state.flip(target);
            updateState(state.to_ullong(), i, state.test(target) ? sinTheta : -sinTheta, false);
        }
    updateLayer();
}

void QubitLayer::applyRz(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // action if bit is 1 (i.e. set)
            if (state.test(target))
                parity ? qOdd_[i] += cosTheta * qEven_[i] + complexImg * sinTheta * qEven_[i] : qEven_[i] += cosTheta * qOdd_[i] + complexImg * sinTheta * qOdd_[i];
            // action if bit is 0 (i.e. not set)
            else
                parity ? qOdd_[i] += cosTheta * qEven_[i] - complexImg * sinTheta * qEven_[i] : qEven_[i] += cosTheta * qOdd_[i] - complexImg * sinTheta * qOdd_[i];
        }
    updateLayer();
}

bool QubitLayer::checkControls(int *controls, int numControls, std::bitset<maxQubits> state)
{
    int finalControl{0};
    for (int i = 0; i < numControls; i++)
        finalControl += state.test(controls[i]);
    return (numControls == finalControl);
}

void QubitLayer::applyCnot(int control, int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
    {
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bit is 1 (i.e. set)
            if (state.test(control))
                updateState(state.flip(target).to_ullong(), i, 1);
            else
                updateState(i, i, 1);
        }
    }
    updateLayer();
}

void QubitLayer::applyToffoli(int control1, int control2, int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bits are 1 (i.e. set)
            if (state[control1] && state[control2])
                updateState(state.flip(target).to_ullong(), i, 1);
            else
                updateState(i, i, 1);
        }
    updateLayer();
}

void QubitLayer::applyMcnot(int *controls, int numControls, int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bit(s) is 1 (i.e. set)
            if (checkControls(controls, numControls, state))
                updateState(state.flip(target).to_ullong(), i, 1);
            else
                updateState(i, i, 1);
        }
    updateLayer();
}

void QubitLayer::applyCz(int control, int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // add phase to target qubit if control bit and target bits are 1 (i.e. set)
            updateState(i, i, (state.test(control) && state.test(target)) ? -1 : 1);
        }
    updateLayer();
}

void QubitLayer::applyMcphase(int *controls, int numControls, int target)
{
    for (unsigned long long int i = 0; i < numStates_; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            updateState(i, i, (checkControls(controls, numControls, state) && state.test(target)) ? -1 : 1);
        }
    updateLayer();
}

qProb QubitLayer::getMaxAmplitude()
{
    std::bitset<maxQubits> state;
    qProb result;
    precision currentProb{0};
    precision previousProb{0};
    for (unsigned long long int i = 0; i < numStates_; i++)
    {
        state = i;
        currentProb = parity ? abs(qEven_[i]) * abs(qEven_[i]) : abs(qOdd_[i]) * abs(qOdd_[i]);
        if (currentProb > previousProb)
        {
            result.state = state;
            result.prob = currentProb;
            previousProb = currentProb;
        }
    }
    return result;
}

void QubitLayer::toggleParity()
{
    parity = !parity;
}

void QubitLayer::printMeasurement()
{
    qProb q = getMaxAmplitude();
    std::cout << "Measurement outcome:        |" << q.state << ">" << std::endl;
    std::cout << "Probability of outcome:     " << q.prob << std::endl;
}

void QubitLayer::printQubits()
{
    std::cout << "Amplitude, "
              << "State \n";
    for (unsigned long long int i = 0; i < numStates_; i++)
    {
        TanglrBitset<maxQubits> binaryRep = i;
        std::string state = binaryRep.toString();
        std::cout << qEven_[i] << " " << qOdd_[i] << " ";
        std::cout << "|" << state << ">\n";
    }
}

void QubitLayer::updateState(unsigned long long int targetIndex, unsigned long long int sourceIndex, std::complex<precision> scalingValue, bool assignmentOperation)
{
    if (assignmentOperation)
        if (parity)
            qOdd_[targetIndex] = scalingValue * qEven_[sourceIndex];
        else
            qEven_[targetIndex] = scalingValue * qOdd_[sourceIndex];
    else if (parity)
        qOdd_[targetIndex] += scalingValue * qEven_[sourceIndex];
    else
        qEven_[targetIndex] += scalingValue * qOdd_[sourceIndex];
}

qubitLayer *QubitLayer::getQubitLayerEven() { return qEven_; }

qubitLayer *QubitLayer::getQubitLayerOdd() { return qOdd_; }

unsigned long long int QubitLayer::getNumStates() { return numStates_; }

unsigned int QubitLayer::getNumQubits() { return numQubits_; }