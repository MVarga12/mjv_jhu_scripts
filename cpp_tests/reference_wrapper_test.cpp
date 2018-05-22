/*
 * Created: 18-May-2018 07:57:41 AM EDT
 * Modified: 18-May-2018 08:28:32 AM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <functional>
#include <iostream>
#include <vector>

struct A {
    int x{ 0 };
    int y{ 0 };
    int z{ 0 };

    friend std::ostream& operator<<(std::ostream& os, const A& a);
    friend void operator*=(A& a, int val);

    A() = default;
    A(int _x, int _y, int _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const A& a) { return os << '(' << a.x << ", " << a.y << ", " << a.z << ')'; }
void operator*=(A& a, int val)
{
    a.x = a.x * val;
    a.y = a.y * val;
    a.z = a.z * val;
}

int main()
{
    std::vector<A> test{ A(1, 2, 3), A(4, 5, 6), A(7, 8, 9) };
    std::vector<std::reference_wrapper<A>> test_ref1{*test.begin()};
    std::vector<std::reference_wrapper<A>> test_ref2{*(test.begin()+1)};

    std::cout << "Orig Vector: " << '\n';
    for (const auto& elem : test)
        std::cout << elem << '\n';

    std::cout << "References1: " << '\n';
    for (const auto& elem : test_ref1)
        std::cout << elem << '\n';

    std::cout << "References2: " << '\n';
    for (const auto& elem : test_ref2)
        std::cout << elem << '\n';

    for (auto& elem : test)
        elem *= 2;

    std::cout << "Orig Vector *2: " << '\n';
    for (const auto& elem : test)
        std::cout << elem << '\n';

    std::cout << "References1 *2: " << '\n';
    for (const auto& elem : test_ref1)
        std::cout << elem << '\n';

    std::cout << "References2 *2: " << '\n';
    for (const auto& elem : test_ref2)
        std::cout << elem << '\n';

    test.erase(test.begin() + 1);

    std::cout << "Orig Vector erase: " << '\n';
    for (const auto& elem : test)
        std::cout << elem << '\n';

    std::cout << "References1 erase: " << '\n';
    for (const auto& elem : test_ref1)
        std::cout << elem << '\n';

    std::cout << "References2 erase: " << '\n';
    for (const auto& elem : test_ref2)
        std::cout << elem << '\n';

    test.insert(test.begin() + 1, A(8, 9 , 10));

    std::cout << "Orig Vector insert: " << '\n';
    for (const auto& elem : test)
        std::cout << elem << '\n';

    std::cout << "References1 insert: " << '\n';
    for (const auto& elem : test_ref1)
        std::cout << elem << '\n';

    std::cout << "References2 insert: " << '\n';
    for (const auto& elem : test_ref2)
        std::cout << elem << '\n';
}
