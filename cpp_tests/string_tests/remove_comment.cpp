/*
 * Created: 10-May-2018 02:20:20 PM EDT
 * Modified: 10-May-2018 02:25:05 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <vector>

void remove_comment(std::string& line)
{
    std::cout << line.size() << " : " << line << '\n';
    auto commentPos{ line.find('#') };
    if (commentPos != std::string::npos)
        line.erase(commentPos, std::string::npos);
    std::cout << line.size() << " : " << line << '\n';
}

int main()
{
    std::vector<std::string> test_str{ "123 # 456", "123456" };
    for (auto& str : test_str)
        remove_comment(str);
}
