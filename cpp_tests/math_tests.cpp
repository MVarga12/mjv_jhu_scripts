/*
 * Created: 19-May-2018 02:52:09 PM EDT
 * Modified: 19-May-2018 03:01:08 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>

struct A {
    int x {0};
    A operator/=(int val) {return {this->x/val};}
    A operator/(int val) {return {this->x/val};}

    void blah() {*this = *this /x;}
    A()=default;
    A(int _x)
        :x(_x)
    {}
};

int main() {

    A test(4);
    test.blah();
    std::cout << test.x << '\n';
}
