#include <iostream>
#include <cassert>
#include "../src/TanglrBitset.hpp"

/**
 * @brief Helper function to run a test and print the result.
 */
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

/**
 * @brief Test: operator[] retrieves bits in reversed order.
 */
void testOperatorAccess()
{
    TanglrBitset<8> bits(178); // 178 in binary is 10110010
    assert(bits[0] == 1);      // LSB (leftmost in reversed order)
    assert(bits[7] == 0);      // MSB (rightmost in reversed order)
    assert(bits[1] == 0);
    assert(bits[6] == 1);
}

/**
 * @brief Test: toString generates the correct reversed string.
 */
void testToString()
{
    // Test case 1: 8-bit bitset
    TanglrBitset<8> bits1(178);             // 178 in binary is 10110010
    assert(bits1.toString() == "01001101"); // Reversed: 01001101

    // Test case 2: 4-bit bitset
    TanglrBitset<4> bits2(6);           // 6 in binary is 0110
    assert(bits2.toString() == "0110"); // Reversed: 0110

    // Test case 3: All bits set
    TanglrBitset<8> bits3(255);             // 255 in binary is 11111111
    assert(bits3.toString() == "11111111"); // Reversed: 11111111

    // Test case 4: All bits unset
    TanglrBitset<8> bits4(0);               // 0 in binary is 00000000
    assert(bits4.toString() == "00000000"); // Reversed: 00000000

    // Test case 5: Larger bitset
    TanglrBitset<16> bits5(43690);                  // 43690 in binary is 1010101010101010
    assert(bits5.toString() == "0101010101010101"); // Reversed: 0101010101010101
}

/**
 * @brief Main function to run all tests.
 */
int main()
{
    std::cout << "\033[34;34m=========== TanglrBitset Tests ===========\033[m" << std::endl;

    runTest("operator[] retrieves bits in reversed order       ", testOperatorAccess);
    runTest("toString generates the correct reversed string    ", testToString);

    std::cout << "\033[34;34mAll tests completed!\033[m" << std::endl;
    return 0;
}