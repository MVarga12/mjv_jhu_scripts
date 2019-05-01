#include <fstream>
#include <iostream>
#include <algorithm>

int main() {
    bool val1 { false };

    std::ifstream file{ "test.txt" };
    std::string line { };

    while (!file.eof()) {
        file >> line;
        std::transform(line.begin(), line.end(), line.end(), ::tolower);
        if (line == "0" || line == "false")
            val1 = false;
        else if (line == "1" || line == "true")
            val1 = true;
        else {
            std::cerr << "FATAL ERROR: Cannot read boolean value.\n";
            exit(1);
        }
        std::cout << val1 << '\n';
    }
}
