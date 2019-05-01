#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

struct Coord {
    double x { 0 };
    double y { 0 };
    double z { 0 };

    explicit Coord() = default;
    Coord(double _x, double _y, double _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
    }
};

struct Molecule {
    std::string name;
    Coord comCoord;
    std::vector<Coord> iCoords;

    explicit Molecule() = default;
    Molecule(std::string _name, Coord _comCoord)
        : name(_name)
        , comCoord(_comCoord)
    {
    }
};

int main()
{
    std::ofstream outFile { "test.pdb" };
    std::string classification { "TEST PDB FILE" };

    std::vector<Molecule> moleculeList;
    moleculeList.emplace_back("clat", Coord { 0.0, 0.0, 0.0 });
    moleculeList.emplace_back("clat", Coord { 1.0, 1.0, 1.0 });
    moleculeList.emplace_back("clat", Coord { -1, -1, -1 });

    // print header
    outFile << std::left << std::setw(6) << "TITLE" << ' ' << std::setw(70) << classification << '\n';

    int i { 0 };
    for (auto mol : moleculeList) {
        ++i;
        outFile << std::right << "ATOM  " << std::setw(5) << i << ' ' << std::setw(4) << mol.name.substr(0, 4) << ' '
                << std::setw(3) << mol.name.substr(0, 3) << ' ' << std::setw(4) << i << ' '
                << std::setw(8) << mol.comCoord.x << std::setw(8) << mol.comCoord.y << std::setw(8) << mol.comCoord.z
                << std::setw(6) << 0.00 << std::setw(6) << 0.00 << std::left << std::setw(2) << "CL" << '\n';
    }

    outFile << "CONECT" << std::setw(5) << 1 << std::setw(5) << 2 << '\n';
    outFile << "END";
}
