/*
 * Created: 25-Apr-2018 10:32:11 AM EDT
 * Modified: 25-Apr-2018 01:03:08 PM EDT
 * Created by: Matthew Varga
 * Purpose: testing an easier way to determine state changes
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

template <class T>
void outp(T &vec) {
    for (auto &elem : vec)
        std::cout << ' ' << elem;
    std::cout << '\n';
}

int main() {
    std::vector<std::string> reactants = {"m2muh~U", "m2muh~P"};
    std::vector<char> states{};

    {
        for (auto iface : reactants) {
            size_t pos{0};
            while ((pos = iface.find('~')) != std::string::npos) {
                // using the [0] index below converts the size 1 string into a char
                std::cout << iface.substr(0, pos) << '\n';
                states.push_back(iface.substr(pos + 1, iface.size())[0]);
                iface.erase(0, pos + 1);
            }
        }
        outp(states);
        std::cout << std::endl << std::boolalpha << bool{states[0] == states[1]};
    }
}
