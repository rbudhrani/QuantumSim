#include <iostream>
#include "qAlgorithms.hpp"

QubitLayer grover(unsigned int numQubits, int unsigned dSolution, int unsigned numReps)
{
    QubitLayer q(numQubits);
    unsigned long long int numStates = q.getNumStates();
    if (dSolution > numStates)
    {
        std::cout << "\033[31;31m[Error]\033[m" << std::endl;
        std::cout << "Number of states: " << numStates << std::endl;
        std::cout << "Your solution:    " << dSolution << std::endl;
        exit(EXIT_FAILURE);
    }
    if (numReps == 0)
        numReps = static_cast<int>(std::sqrt(numStates) / 4 * pi);
    std::bitset<maxQubits> bSolution = dSolution;
    int ctrlQubits[numQubits - 1];
    for (int i = 0; i < (numQubits - 1); i++)
        ctrlQubits[i] = i;
    // initiliase qubits to a superposition of all states
    for (int i = 0; i < numQubits; i++)
        q.applyHadamard(i);
    for (int rep = 0; rep < numReps; rep++)
    {
        // oracle to tag solution
        for (int i = 0; i < numQubits; i++)
            if (!bSolution.test(i))
                q.applyPauliX(i);
        q.applyMcphase(ctrlQubits, numQubits - 1, numQubits - 1);
        for (int i = 0; i < numQubits; i++)
            if (!bSolution.test(i))
                q.applyPauliX(i);
        // grover diffusion operator (inversion about mean and amplitude amplification)
        for (int i = 0; i < numQubits; i++)
            q.applyHadamard(i);
        for (int i = 0; i < numQubits; i++)
            q.applyPauliX(i);
        q.applyMcphase(ctrlQubits, numQubits - 1, numQubits - 1);
        for (int i = 0; i < numQubits; i++)
            q.applyPauliX(i);
        for (int i = 0; i < numQubits; i++)
            q.applyHadamard(i);
    }
    return q;
}

/*QubitLayer repCode3(unsigned int numQubits, unsigned int errorLocation, pauliError errorType){
    QubitLayer q(numQubits);
    if (numQubits != 3){
        std::cout<<"\033[31;31m[Error]\033[m"<<std::endl;
        std::cout<<"Number of qubits defined:                            "<<numQubits<<std::endl;
        std::cout<<"Number of qubits used by 3 qubit repetition code:    "<<3<<std::endl;
        exit(EXIT_FAILURE);
    }
    //initialise state to anything you want
    q.applyRx(0, pi/3);
    std::cout<<"Initial state:"<<std::endl;
    q.printQubits();
    //"copy" state
    q.applyCnot(0, 1);
    q.applyCnot(0, 2);
    //inject error
    switch(errorType){
        case errorX: q.applyPauliX(errorLocation); break;
        case errorY: q.applyPauliY(errorLocation); break;
        case errorZ: q.applyPauliZ(errorLocation); break;
    }
    //undo the "copy"
    q.applyCnot(0, 1);
    q.applyCnot(0, 2);
    q.applyToffoli(2, 1, 0);
    std::cout<<"Final state:"<<std::endl;
    return q;
}

QubitLayer genEPR(unsigned int numQubits, unsigned int q1, unsigned int q2){
    if (numQubits < 2){
        std::cout<<"\033[31;31m[Error]\033[m"<<std::endl;
        std::cout<<"Number of qubits defined:           "<<numQubits<<std::endl;
        std::cout<<"Minimum number of qubits needed:    "<<2<<std::endl;
        exit(EXIT_FAILURE);
    }
    QubitLayer q(numQubits);
    //apply hadamard to the first qubit
    q.applyHadamard(q2);
    //apply cnot on the 2 qubits
    q.applyCnot(q2, q1);
    return q;
}*/