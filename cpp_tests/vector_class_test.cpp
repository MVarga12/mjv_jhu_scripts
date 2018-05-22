/*
 * Created: 18-May-2018 03:51:38 PM EDT
 * Modified: 18-May-2018 04:02:12 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <cmath>

std::ostream& bon(std::ostream& os) { 
    return os << "\e[1m";
}

std::ostream& boff(std::ostream& os) { 
    return os << "\e[0m";
}

class Coord {
public:
    double x;
    double y;
    double z;

    void display() const;
    Coord() = default;
    Coord(double _x, double _y, double _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
    }
};

void Coord::display() const {
    std::cout << x << ' '  << y << ' ' << z << '\n';
}

struct Vector : public Coord {
    /*!
     * \brief Holds a vector with the origin as the start point
     */

    double magnitude{ 0.0 };

    void display() const;
    Vector() = default;
    Vector(Coord& _coord)
    {
        x = _coord.x;
        y = _coord.y;
        z = _coord.z;

        magnitude = sqrt(x * x + y * y + z * z);
    }
};

void Vector::display() const {
    std::cout << x << ' ' << y << ' ' << z << " mag: |" << magnitude << "|\n";
}

int main() {
    Coord test{1, 2, 3};
    Vector test_v{test};

    test.display();
    test_v.display();
}
