/*
 * Created: Thu 19 Apr 2018 11:02:47 EDT
 * Modified: Thu 19 Apr 2018 11:09:15 EDT
 * Created by: Matthew Varga
 * Purpose: test public/private class members
 */

#include <iostream>

class Coord {
private:
    double x;
    double y;
    double z;

public:
    void set_crds(int nx, int ny, int nz);
    void display();

    Coord() {
    }
    Coord(double x1, double y1, double z1)
        : x(x1)
        , y(y1)
        , z(z1) {
    }
};

void Coord::set_crds(int nx, int ny, int nz) {
    x = nx;
    y = ny;
    z = nz;
}

void Coord::display() {
    std::cout << '(' << x << ", " << y << ", " << z << ")\n";
}

int main() {
    Coord new_crd{Coord(0, 0, 0)};

    new_crd.set_crds(1, 0, 0);
    new_crd.x = 1;
    new_crd.display();
}
