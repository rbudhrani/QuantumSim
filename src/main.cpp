#include <iostream>
#include <complex>
#include <chrono>
#include "QubitLayer.hpp"
#include "../examples/qAlgorithms.hpp"

int main(int argc, char *argv[])
{
    auto start = std::chrono::steady_clock::now();
    QubitLayer q = grover(4, 0);
    auto stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    std::cout << "Execution time: " << duration << " µs" << std::endl;
    q.printQubits();
}