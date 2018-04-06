#include <iostream>
#include <memory>
#include <vector>

class Base {
public:
    int x;

    virtual void display() = 0;
    Base() {}
    Base(int x1)
        : x(x1) {}
};

class Derived_a : public Base {
public:
    int a;

    void display();
    Derived_a(int a1)
        : a(a1) {}
};

void Derived_a::display() {
    std::cout << x << ' ' << a << std::flush;
}

class Derived_b : public Base {
public:
    int b;

    void display();
    Derived_b(int b1)
        : b(b1) {}
};

void Derived_b::display() {
    std::cout << x << ' ' << b << std::endl;
}

int main() {
    // this is a est comment
    std::vector<std::unique_ptr<Base>> vec;
    std::unique_ptr<Derived_a> a(new Derived_a(10));
    std::unique_ptr<Derived_b> b(new Derived_b(20));
    a->x = 5;
    b->x = 10;
    vec.push_back(std::move(a));
    vec.push_back(std::move(b));
    vec.push_back(std::unique_ptr<Derived_a>(new Derived_a(10)));

    vec.at(0)->display();
    vec.at(1)->display();
    vec.at(2)->display();
}
