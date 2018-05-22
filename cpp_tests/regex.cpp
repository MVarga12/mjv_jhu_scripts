/*
 * Created: 10-May-2018 09:42:19 AM EDT
 * Modified: 10-May-2018 04:26:20 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <regex>

int main()
{
    std::vector<std::string> test_strings{ "0false", "False", "fAlse", "FALSE", "0" };
    std::vector<std::string> pi_strs{ "M_PI", "pi", "m_pi", "PI", "PU" };
    std::vector<std::string> null_strs{ "null", "NULL", "0NULL", "0", "O" };
    std::regex pi_match("[mM]*_*[pP][iI]");
    std::regex match_str("[0]|[fF][aA][lL][sS][eE]");
    std::regex null_match("[0]|[oO]|[nN][uU][lL][lL]");

    for (const auto& candidate : test_strings)
        std::cout << candidate << ' ' << std::regex_match(candidate, match_str) << '\n';
    for (const auto& candidate : pi_strs)
        std::cout << candidate << ' ' << std::regex_match(candidate, pi_match) << '\n';
    for (const auto& candidate : null_strs)
        std::cout << candidate << ' ' << std::regex_match(candidate, null_match) << '\n';
}
