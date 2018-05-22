/*
 * Created: 16-May-2018 10:30:22 AM EDT
 * Modified: 16-May-2018 10:41:19 AM EDT
 * Created by: Matthew Varga
 * Purpose: testing function to see if all vector elements are the same
 * if using C++14, could use std::adjacent_find, but we're not
 */

#include <iostream>
#include <vector>

template<class X>
bool allSame(const std::vector<X>& v) {
    X prevType;
    auto elemItr {v.begin()};
    bool isSame {true};

    while (isSame) {
        prevType = *elemItr;
        ++elemItr;
        if (elemItr == v.end())
            break;
        std::cout << *elemItr << ' ' << prevType << '\n';
        isSame = (*elemItr == prevType) ? true : false;
    }
    return isSame;
}

int main() {
    std::vector<int> test {1, 2, 3, 4};
    std::vector<int> test2 {1, 1, 1, 1};

    std::cout << allSame(test) << '\n';
    std::cout << allSame(test2) << '\n';
}

