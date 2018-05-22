/*
 * Created: 27-Apr-2018 08:27:49 AM EDT
 * Modified: 27-Apr-2018 09:26:39 AM EDT
 * Created by: Matthew Varga
 * Purpose: testing std::multimap to use for cooperativity and
 * omega reactions
 */

#include <iostream>
#include <map>

int main()
{
    // std::multimap<X key, Y value>
    std::multimap<std::string, int> reqd_bonds;
    reqd_bonds.insert({ "m2muh", 1 });
    reqd_bonds.insert({ "td1", 2 });
    reqd_bonds.insert({ "m2muh", 3 });

    // print entire map
    for (auto& elem : reqd_bonds)
        std::cout << '{' << elem.first << ", " << elem.second << "} ";
    std::cout << '\n';

    // this finds the first key-pair which have the stated key
    auto pos = reqd_bonds.find("m2muh");
    std::cout << "std::map::find: ";
    if (pos != reqd_bonds.end()) {
        std::cout << '{' << pos->first << ", " << pos->second << "} ";
    }
    std::cout << '\n';

    // this finds all key pairs which have the stated key
    auto range = reqd_bonds.equal_range("m2muh");
    std::cout << "std::map::equal_range: ";
    for (auto i{ range.first }; i != range.second; i++)
        std::cout << "{" << i->first << ", " << i->second << "} ";
    std::cout << '\n';
}
