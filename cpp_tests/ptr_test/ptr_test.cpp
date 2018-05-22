#include <iostream>
#include <memory>

class B {};

class A {
public:
    int a;
    std::unique_ptr<B> Bptr;
};

int main() {
}
