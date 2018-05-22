/*
 * Created:  8-May-2018 01:35:25 PM EDT
 * Modified:  9-May-2018 02:47:32 PM EDT
 * Created by: Matthew Varga
 * Purpose: just testing tolower
 */

#include <fstream>
#include <iostream>
#include <sstream>

int main()
{
    { // testing if it works if there's a numeral
        std::string test{ "M2MUH" };
        for (auto& ch : test) {
            ch = std::tolower(ch);
        }
        std::cout << test << '\n';
    }

    { // testing tolower on string, with a loop
        std::string test{ "1#123" };
        for (auto& c : test.erase(test.find('#'), test.size()))
            c = std::tolower(c);
        std::cout << test << '\n';
    }

    {
        // if the string doesn't have #, it throws a std::out_of_range error
        std::string test{ "123" };
        try {
            std::cout << test.erase(test.find('#'), std::string::npos) << '\n';
            std::cout << "hey\n";
        } catch (std::out_of_range) {
            std::cout << "I had an error, but I'm going to continue\n";
        }
        std::cout << test << '\n';
    }

    { // testing tellg() and seekg()
        std::ifstream test("in.dat");
        auto pos{ test.tellg() };

        while (!test.eof()) {
            std::string tmp;
            test >> tmp;
            std::cout << tmp << ' ';
        }
        std::cout << '\n';

        test.clear();
        test.seekg(pos);
        while (!test.eof()) {
            std::string tmp;
            test >> tmp;
            std::cout << tmp << ' ';
        }
    }
}
