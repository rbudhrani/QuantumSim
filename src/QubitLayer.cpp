#include <complex>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "QubitLayer.hpp"

QubitLayer::QubitLayer(unsigned int numQubits, qubitLayer *qL)
{
    // calculate the number of states
    _numStates = 1;
    for (unsigned int i = 0; i < numQubits; i++)
        _numStates *= 2;
    // allocate memory for arrays
    _qEven = new qubitLayer[_numStates];
    _qOdd = new qubitLayer[_numStates];
    // if input is provided then use that to fill the input qubit state
    if (!(qL == nullptr))
        for (unsigned int row = 0; row < _numStates; row++)
            _qEven[row] = qL[row];
    else
        // the default state is |000...0>
        _qEven[0] = {1, 0};
}

QubitLayer::~QubitLayer()
{
    delete[] _qEven;
    delete[] _qOdd;
}

void QubitLayer::updateLayer()
{
    _parity ? std::fill(_qEven, _qEven + _numStates, zeroComplex) : std::fill(_qOdd, _qOdd + _numStates, zeroComplex);
    toggleParity();
}

bool QubitLayer::checkZeroState(int qubit)
{
    return _parity ? _qEven[qubit].imag() || _qEven[qubit].real() : _qOdd[qubit].imag() || _qOdd[qubit].real();
}

void QubitLayer::applyPauliX(int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            state.flip(target);
            _parity ? _qOdd[state.to_ulong()] = _qEven[i] : _qEven[state.to_ulong()] = _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyPauliY(int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip qubit
            state.flip(target);
            // add phase of -i if bit was 1 (i.e. set) and flip it
            if (!state.test(target))
                _parity ? _qOdd[state.to_ulong()] = -complexImg * _qEven[i] : _qEven[state.to_ulong()] = -complexImg * _qOdd[i];
            // add phase of i if bit was 0 (i.e. set) and flip it
            else
                _parity ? _qOdd[state.to_ulong()] = complexImg * _qEven[i] : _qEven[state.to_ulong()] = complexImg * _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyPauliZ(int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // add phase if bit is 1 (i.e. it is set)
            if (state.test(target))
                _parity ? _qOdd[i] = -_qEven[i] : _qEven[i] = -_qOdd[i];
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyHadamard(int target)
{
// map |1> to -hadamardCoef*|1> and |0> to hadamardCoef*|0>
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            if (state.test(target))
                _parity ? _qOdd[i] -= hadamardCoef * _qEven[i] : _qEven[i] -= hadamardCoef * _qOdd[i];
            else
                _parity ? _qOdd[i] += hadamardCoef * _qEven[i] : _qEven[i] += hadamardCoef * _qOdd[i];
            if (!isParallel)
            {
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] += hadamardCoef * _qEven[i] : _qEven[state.to_ulong()] += hadamardCoef * _qOdd[i];
            }
        }
    if (isParallel)
    {
#pragma omp barrier
// map |0> to hadamardCoef*|1> and |1> to hadamardCoef*|0>
#pragma omp parallel for shared(qOdd_, qEven_)
        for (unsigned long long int i = 0; i < _numStates; i++)
            if (checkZeroState(i))
            {
                std::bitset<maxQubits> state = i;
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] += hadamardCoef * _qEven[i] : _qEven[state.to_ulong()] += hadamardCoef * _qOdd[i];
            }
    }
    updateLayer();
}

void QubitLayer::applyRx(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
// map |0> to cosTheta*|0> and |1> to cosTheta*|1>
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            _parity ? _qOdd[i] += cosTheta * _qEven[i] : _qEven[i] += cosTheta * _qOdd[i];
            if (!isParallel)
            {
                std::bitset<maxQubits> state = i;
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] += -complexImg * sinTheta * _qEven[i] : _qEven[state.to_ulong()] += -complexImg * sinTheta * _qOdd[i];
            }
        }
    if (isParallel)
    {
#pragma omp barrier
// map |1> to -isinTheta*|0> and |0> to -isineTheta*|1>
#pragma omp parallel for shared(qOdd_, qEven_)
        for (unsigned long long int i = 0; i < _numStates; i++)
            if (checkZeroState(i))
            {
                std::bitset<maxQubits> state = i;
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] += -complexImg * sinTheta * _qEven[i] : _qEven[state.to_ulong()] += -complexImg * sinTheta * _qOdd[i];
            }
    }
    updateLayer();
}

void QubitLayer::applyRy(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
// map |1> to cosTheta*|1> and |0> to cosTheta*|0>
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            _parity ? _qOdd[i] += cosTheta * _qEven[i] : _qEven[i] += cosTheta * _qOdd[i];
            if (!isParallel)
            {
                std::bitset<maxQubits> state = i;
                state.flip(target);
                // action if bit is 1 (i.e. set)
                if (state.test(target))
                    _parity ? _qOdd[state.to_ulong()] += sinTheta * _qEven[i] : _qEven[state.to_ulong()] += sinTheta * _qOdd[i];
                // action if bit is 0 (i.e. not set)
                else
                    _parity ? _qOdd[state.to_ulong()] -= sinTheta * _qEven[i] : _qEven[state.to_ulong()] -= sinTheta * _qOdd[i];
            }
        }
    if (isParallel)
    {
#pragma omp barrier
// map |0> to sinTheta*|1> and |1> to -sinTheta*|0>
#pragma omp parallel for shared(qOdd_, qEven_)
        for (unsigned long long int i = 0; i < _numStates; i++)
            if (checkZeroState(i))
            {
                std::bitset<maxQubits> state = i;
                state.flip(target);
                // action if bit is 1 (i.e. set)
                if (state.test(target))
                    _parity ? _qOdd[state.to_ulong()] += sinTheta * _qEven[i] : _qEven[state.to_ulong()] += sinTheta * _qOdd[i];
                // action if bit is 0 (i.e. not set)
                else
                    _parity ? _qOdd[state.to_ulong()] -= sinTheta * _qEven[i] : _qEven[state.to_ulong()] -= sinTheta * _qOdd[i];
            }
    }
    updateLayer();
}

void QubitLayer::applyRz(int target, precision theta)
{
    // compute the sine and cosine of the rotation angle
    precision cosTheta = cos(theta / 2);
    precision sinTheta = sin(theta / 2);
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // action if bit is 1 (i.e. set)
            if (state.test(target))
                _parity ? _qOdd[i] += cosTheta * _qEven[i] + complexImg * sinTheta * _qEven[i] : _qEven[i] += cosTheta * _qOdd[i] + complexImg * sinTheta * _qOdd[i];
            // action if bit is 0 (i.e. not set)
            else
                _parity ? _qOdd[i] += cosTheta * _qEven[i] - complexImg * sinTheta * _qEven[i] : _qEven[i] += cosTheta * _qOdd[i] - complexImg * sinTheta * _qOdd[i];
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
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
    {
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bit is 1 (i.e. set)
            if (state.test(control))
            {
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] = _qEven[i] : _qEven[state.to_ulong()] = _qOdd[i];
            }
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    }
    updateLayer();
}

void QubitLayer::applyToffoli(int control1, int control2, int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bits are 1 (i.e. set)
            if (state.test(control1) && state.test(control2))
            {
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] = _qEven[i] : _qEven[state.to_ulong()] = _qOdd[i];
            }
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyMcnot(int *controls, int numControls, int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // flip target qubit if control bit(s) is 1 (i.e. set)
            if (checkControls(controls, numControls, state))
            {
                state.flip(target);
                _parity ? _qOdd[state.to_ulong()] = _qEven[i] : _qEven[state.to_ulong()] = _qOdd[i];
            }
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyCz(int control, int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // add phase to target qubit if control bit and target bits are 1 (i.e. set)
            if (state.test(control) && state.test(target))
                _parity ? _qOdd[i] = -_qEven[i] : _qEven[i] = -_qOdd[i];
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    updateLayer();
}

void QubitLayer::applyMcphase(int *controls, int numControls, int target)
{
#pragma omp parallel for shared(qOdd_, qEven_)
    for (unsigned long long int i = 0; i < _numStates; i++)
        if (checkZeroState(i))
        {
            std::bitset<maxQubits> state = i;
            // add phase to target qubit if control bit(s) and target bit is 1 (i.e. set)
            if (checkControls(controls, numControls, state) && state.test(target))
                _parity ? _qOdd[i] = -_qEven[i] : _qEven[i] = -_qOdd[i];
            else
                _parity ? _qOdd[i] = _qEven[i] : _qEven[i] = _qOdd[i];
        }
    updateLayer();
}

qProb QubitLayer::getMaxAmplitude()
{
    std::bitset<maxQubits> state;
    qProb result;
    precision currentProb{0};
    precision previousProb{0};
    for (unsigned long long int i = 0; i < _numStates; i++)
    {
        state = i;
        currentProb = _parity ? abs(_qEven[i]) * abs(_qEven[i]) : abs(_qOdd[i]) * abs(_qOdd[i]);
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
    _parity = !_parity;
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
    for (unsigned long long int i = 0; i < _numStates; i++)
    {
        std::bitset<maxQubits> binaryRep = i;
        std::string state = binaryRep.to_string();
        std::cout << _qEven[i] << " " << _qOdd[i] << " ";
        std::cout << "|" << state << ">\n";
    }
}

qubitLayer *QubitLayer::getQubitLayerEven()
{
    return _qEven;
}

qubitLayer *QubitLayer::getQubitLayerOdd()
{
    return _qOdd;
}

unsigned long long int QubitLayer::getNumStates() { return _numStates; }

unsigned int QubitLayer::getNumQubits() { return _numQubits; }