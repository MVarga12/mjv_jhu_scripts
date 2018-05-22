/*
 * Created:  4-May-2018 02:51:35 PM EDT
 * Modified:  4-May-2018 02:53:57 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <vector>

int main()
{
    std::vector<int> info{ 3, 1, 2, 4 };
    std::cout << "unsorted:";
    for (auto& elem : info)
        std::cout << ' ' << elem;

    std::sort(info.begin(), info.end());
    std::cout << "\nsorted:";
    for (auto& elem : info)
        std::cout << ' ' << elem;
}
