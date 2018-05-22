/*
 * Created: 11-May-2018 09:04:34 AM EDT
 * Modified: 11-May-2018 01:33:43 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <vector>

struct A {
    std::string name;

    A() = default;
    explicit A(std::string _name)
        : name(_name)
    {
    }
};

struct B {
    std::string name;

    B() = default;
    explicit B(std::string _name)
        : name(_name)
    {
    }
};

template <class X, class Y>
bool haveSameName(const X& item1, const Y& item2)
{
    return item1.name == item2.name;
}

int main()
{
    std::vector<A> testA{ A("a1"), A("a2"), A("a1"), A("a4") };
    std::vector<B> testB{ B("a1"), B("a2"), B("a3"), B("a4") };

    for (const auto& oneB : testB) {
        auto oneANameItr
            = std::find_if(testA.begin(), testA.end(), [&](const A& oneA) -> bool { return oneA.name == oneB.name; });
        std::cout << oneB.name << ' ' << oneANameItr->name << '\n';
    }

    { // what if there are multiples?
        auto oneVal{ "a1" };
        auto oneNameItr
            = std::find_if(testA.begin(), testA.end(), [&](const A& oneA) -> bool { return oneA.name == oneVal; });
        std::cout << std::distance(testA.begin(), oneNameItr) << '\n';
        while (oneNameItr != testA.end()) {
            ++oneNameItr;
            oneNameItr
                = std::find_if(oneNameItr, testA.end(), [&](const A& oneA) -> bool { return oneA.name == oneVal; });
            // this is a ternary operator (condition) ? (do this if true) : (do this if false)
            std::cout << (oneNameItr == testA.end() ? "No More Values"
                                                    : std::to_string(std::distance(testA.begin(), oneNameItr)))
                      << '\n';
        }
    }

    // testing using boolean returned functions as predicate
    for (const auto& oneB : testB) {
        auto oneANameItr
            = std::find_if(testA.begin(), testA.end(), [&](const A& oneA) -> bool { return haveSameName(oneB, oneA); });
        std::cout << oneB.name << ' ' << oneANameItr->name << '\n';
    }
}
