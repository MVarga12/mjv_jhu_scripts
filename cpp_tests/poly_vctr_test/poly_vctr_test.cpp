#include <vector>
#include <memory>
#include <iostream>

class Base {
    public:
        int x;

        virtual void display() = 0;
        virtual void new_back_rxn(std::vector<std::unique_ptr<Base> > &v) = 0;

        Base() {}
        Base(int x1)
            :x(x1)
        {}
};

class Derived_a : public Base {
    public:
        int a;
        class Derived_b *conj_b;

        void display();
        void new_back_rxn(std::vector<std::unique_ptr<Base> > &v);

        Derived_a(int a1)
            : a(a1)
        {}
};

void Derived_a::display() {
    std::cout << x << ' ' << a << std::endl;
}

void Derived_a::new_back_rxn(std::vector<std::unique_ptr<Base> > &v) {

}

class Derived_b : public Base {
    public:
        int b;
        class Derived_a *conj_a;

        void display();
        void new_back_rxn(std::vector<std::unique_ptr<Base> > &v) = 0;

        Derived_b(int b1)
            : b(b1)
        {}
};

void Derived_b::display() {
    std::cout << x << ' ' << b << std::endl;
}

class ReactionList {
    std::string name;
    std::vector<std::unique_ptr<Base> > rxns;

    ReactionList() {}
    ReactionList(std::string &n)
        :name(n)
    {}
}


int main() {
    {
        std::vector<std::shared_ptr<Base> > vec;
        std::shared_ptr<Derived_a> a (new Derived_a(10));
        std::shared_ptr<Derived_b> b (new Derived_b(20));
        a->x = 5;
        b->x = 10;
        vec.push_back(std::move(a));
        vec.push_back(std::move(b));
        vec.push_back(std::shared_ptr<Derived_a> (new Derived_a(10)));

        vec.at(0)->display();
        vec.at(1)->display();
        vec.at(2)->display();
    }
    { // vec of references
        std::vector<std::unique_ptr<std::vector<std::unique_ptr<Base> > > > Rxnlist;
    }
}
