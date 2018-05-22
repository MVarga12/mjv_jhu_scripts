/*
 * Created: 28-Apr-2018 03:41:01 PM EDT
 * Modified: 28-Apr-2018 03:54:45 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <string>

int main()
{
    std::string test_str{ "this is a test string" };

    std::cout << "array-buffer method\n";
    // using the current method, make a buffer and fill it
    {
        char* buffer = new char[100];
        int i{ 0 };
        for (auto it{ test_str.begin() }; it != test_str.end(); ++it) {
            if (isalnum(*it))
                buffer[i] = *it;
            else if (isspace(*it)) {
                buffer[i] = '\0';
                std::cout << buffer << '\n';
                std::fill(buffer, buffer + sizeof(buffer), '\0');
                i = -1;
            }
            ++i;
        }
        delete[] buffer;
    }

    // using std::string::append
    {
        std::cout << "\nstring-buffer method\n";
        std::string buffer{ "" };
        int i{ 0 };
        for (auto it{ test_str.begin() }; it != test_str.end(); ++it) {
            if (isalnum(*it))
                buffer.append(1, *it);
            else if (isspace(*it)) {
                std::cout << buffer << '\n';
                buffer = "";
                i = -1;
            }
            ++i;
        }
    }
}
