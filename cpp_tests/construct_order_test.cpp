/*
 * Created: 17-May-2018 01:10:13 PM EDT
 * Modified: 17-May-2018 01:14:49 PM EDT
 * Created by: Matthew Varga
 * Purpose: testing timing of constructor actions
 */

#include <iostream>

class A {
    private:
        int a{1};

    public:
        int get_a() const {return a;};

        A() = default;
        A(int _a)
            :a(_a)
        {
            std::cout << "In constructor: " << get_a() << '\n';
        }
};

int main() {
    A a1;
    A a2(3);

    std::cout << "Default outside: " << a1.get_a() << '\n';
    std::cout << "Constructed outside: " << a2.get_a() << '\n';
}
