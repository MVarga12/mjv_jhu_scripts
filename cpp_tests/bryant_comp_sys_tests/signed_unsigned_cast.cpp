#include <iostream>
#include <limits>

int main() {
    // Casting between signed and unsigned numbers (Ch. 2.2.4, pp. 65-69)
    int x { -12345 };
    unsigned y { 4294954951 };
    std::cout << static_cast<unsigned>(x) << '\n';
    std::cout << static_cast<int>(y) << '\n';
}
