#include <bitset>
#include <iostream>
#include "TanglrBitset.hpp"

template <unsigned long long int N>
TanglrBitset<N>::TanglrBitset(unsigned long long int value) : bits_(value)
{
}

template <unsigned long long int N>
bool TanglrBitset<N>::operator[](unsigned long long int index) const { return bits_[N - 1 - index]; }

template <unsigned long long int N>
std::string TanglrBitset<N>::toString() const
{
    std::string reversedString;
    for (unsigned long long int i = 0; i < N; ++i)
        reversedString = ((*this)[i] ? '1' : '0') + reversedString;
    return reversedString;
}