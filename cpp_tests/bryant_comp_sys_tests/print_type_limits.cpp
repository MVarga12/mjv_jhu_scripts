#include <iostream>
#include <limits>
#include <iomanip>

std::ostream& bon(std::ostream& os)
{
    return os << "\e[1m";
}

std::ostream& boff(std::ostream& os)
{
    return os << "\e[0m";
}

int main() {
    std::cout << std::setw(20) << bon << "type" << boff << '\n';
    std::cout << std::setw(25) << "char " << std::numeric_limits<char>::max() << '\n';
    std::cout << std::setw(25) << "unsigned char " << std::numeric_limits<unsigned char>::max() << '\n';
    std::cout << std::setw(25) << "short int " << std::numeric_limits<short>::max() << '\n';
    std::cout << std::setw(25) << "signed short int " << std::numeric_limits<signed short>::max() << '\n';
    std::cout << std::setw(25) << "int " << std::numeric_limits<int>::max() << '\n';
    std::cout << std::setw(25) << "signed int " << std::numeric_limits<signed int>::max() << '\n';
    std::cout << std::setw(25) << "unsigned int " << std::numeric_limits<unsigned int>::max() << '\n';
    std::cout << std::setw(25) << "long int " << std::numeric_limits<long int>::max() << '\n';
    std::cout << std::setw(25) << "signed long int " << std::numeric_limits<signed long int>::max() << '\n';
    std::cout << std::setw(25) << "unsigned long int " << std::numeric_limits<unsigned long int>::max() << '\n';
    std::cout << std::setw(25) << "signed long long int " << std::numeric_limits<signed long long int>::max() << '\n';
    std::cout << std::setw(25) << "unsigned long long int " << std::numeric_limits<unsigned long long int>::max() << '\n';
    std::cout << std::setw(25) << "float " << std::numeric_limits<float>::max() << '\n';
    std::cout << std::setw(25) << "double " << std::numeric_limits<double>::max() << '\n';
    std::cout << std::setw(25) << "long double " << std::numeric_limits<long double>::max() << '\n';
}
