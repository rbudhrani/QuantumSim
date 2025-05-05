#ifndef TANGLR_BITSET_HPP
#define TANGLR_BITSET_HPP

#include <bitset>
#include <string>

/**
 * @class TanglrBitset
 * @brief A wrapper around std::bitset that provides reversed bit access and representation.
 *
 * The TanglrBitset class allows you to work with a bitset where the least significant bit (LSB)
 * is treated as the leftmost bit, and the most significant bit (MSB) is treated as the rightmost bit.
 * It provides functionality to initialize the bitset, access bits in reversed order, and generate
 * a reversed string representation of the bitset.
 *
 * @tparam N The size of the bitset (number of bits).
 */
template <unsigned long long int N>
class TanglrBitset
{
public:
    /**
     * @brief Constructs a TanglrBitset from an unsigned long long integer.
     *
     * Initializes the bitset with the binary representation of the given value.
     * The value is truncated if it exceeds the size of the bitset (N bits).
     *
     * @param[in] value The unsigned long long integer to initialize the bitset.
     */
    TanglrBitset(unsigned long long int value);

    /**
     * @brief Accesses a bit in reversed order.
     *
     * This operator allows you to access bits in reversed order, where index 0 corresponds
     * to the least significant bit (LSB) and index N-1 corresponds to the most significant bit (MSB).
     *
     * @param[in] index The index of the bit to access (0-based, reversed order).
     * @return The value of the bit at the specified index (true for 1, false for 0).
     */
    bool operator[](unsigned long long int index) const;

    /**
     * @brief Generates a string representation of the bitset in reversed order.
     *
     * The string representation treats the least significant bit (LSB) as the leftmost character
     * and the most significant bit (MSB) as the rightmost character. This is the reverse of the
     * default order used by std::bitset.
     *
     * @return A string representing the bitset in reversed order, where the LSB is the first character
     *         and the MSB is the last character.
     *
     * @example
     * TanglrBitset<8> bits(179); // 179 in binary is 10110011
     * std::string reversedBits = bits.toString();
     * // reversedBits will be "11001101"
     */
    std::string toString() const;

private:
    /**
     * @brief The underlying std::bitset that stores the bits.
     *
     * The bitset uses the default bit order, where the least significant bit (LSB) is the rightmost bit,
     * and the most significant bit (MSB) is the leftmost bit.
     */
    std::bitset<N> bits_;
};

#include "TanglrBitset.tpp"

#endif // TANGLR_BITSET_HPP