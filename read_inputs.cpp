#include <array>
#include <vector>
#include <limits>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

/* BNGL parser which reads a formatted file including BNGL reaction
 * and creates a reaction for use in some molecular dynamics style simulation package
 *
 * FILE FORMAT:
 * Num reactions #integer
 * -----
 *  BNGL reaction
 *  double #onrate
 *  double #offrate
 *  double, double, double #coordinates for separation of the two molecules after association
 *  double, double, double #angles between the two molecules after association
 *  int # 2D or 3D, 1 or 0
 */

class Interface {
    public:
    std::string name;
    int index;

    Interface() {}
    Interface(std::string name1, int index1)
        :name(name1), index(index1)
    {}
};

class Molecule {
    public:
        int protype;
        std::string name;
        std::vector<Interface> Interfaces;

        Molecule() {}
        Molecule(int protype1, std::string name1, std::vector<Interface> Int1)
            :protype(protype1), name(name1), Interfaces(Int1)
        {}
};

class Reaction {
    public:
        bool onmem;
        double onrate;
        double offrate;

        std::array<double, 3> sigma; // TODO: this should be a Coord class to make it easier?
        std::array<double, 3> angles;
        std::vector<int> Reactants;
        std::vector<int> Products;
        std::vector<std::string> prodspecies; // strings showing the BNGL products, for comparison


        void display() {
            std::cout << "Reactants:";
            for (auto &specie : Reactants) 
                std::cout << ' ' << specie;
            std::cout << std::endl;
            std::cout << "Products:";
            for (auto &specie : Products) 
                std::cout << ' ' << specie;
            std::cout << std::endl;
            std::cout << "Onrate: " << onrate << std::endl;
            std::cout << "Offrate: " << offrate << std::endl;
            std:: cout << "Sigma:";
            for (auto &elem : sigma) 
                std::cout << ' ' << elem;
            std::cout << std::endl;
            std:: cout << "Angles:";
            for (auto &elem : angles) 
                std::cout << ' ' << elem;
            std::cout << std::endl;
            std::cout << "MemStat: " << onmem << std::endl;
        }
        Reaction() {}
};

bool parse_reaction(int &totspecies, std::string reaction, std::vector<Molecule> &wholep, std::vector<Reaction> &forRxns) {
    struct parsedMol {
        std::string name;
        std::vector<std::string> ifaces;
        int *used; // this is the interface used to react

        void display() {
            std::cout << "Mol: " << name << std::endl;
            std::cout << '\t' << "ifaces:";
            for (auto &iface : ifaces) 
                std::cout << ' ' << iface;
            std::cout << std::endl;
        }
    };
    struct parsedRxn {
        std::vector<std::string> reactants;
        std::vector<std::string> ifaces;
        std::vector<int> bonds;

        void display() {
            std::cout << "mol:";
            for (auto &elem : reactants)
                std::cout << ' ' << elem;
            std::cout << std::endl << "reactants:";
            for (auto &elem : ifaces)
                std::cout << ' ' << elem;
            std::cout << std::endl << "bonds:";
            for (auto &elem : bonds)
                std::cout << ' ' << elem;
            std::cout << std::endl;
        }
    };

    bool revflag {false};
    std::string prods;
    std::string reacts;
    { // break into reactant and product side using <-> or -> 
        if (reaction.find("<->") != std::string::npos) {
            std::cout << "Reversible reaction." << std::endl;
            size_t position = reaction.find("<->");
            reacts = reaction.substr(0, position);
            prods = reaction.substr(position + 4, std::string::npos); // +4 is to remove the delimiter and space
            revflag = true;
        }
        else if (reaction.find("->") != std::string::npos) {
            std::cout << "Irreversible reaction." << std::endl;
            size_t position = reaction.find("->");
            reacts = reaction.substr(0, position);
            prods = reaction.substr(position + 3, std::string::npos); // +3 is to remove delimiter and space
        }
        else {
            std::cerr << "Error: Missing reaction arrow, not a valid reaction. Exiting..." << std::endl;
            std::cout << reaction << std::endl;
            exit(1);
        }
    }

    std::vector<std::string> reactvec;
    std::vector<std::string> prodvec;
    { // break into species based on '+'
        size_t position {0};
        while ((position = reacts.find('+')) != std::string::npos) {
            reactvec.push_back(reacts.substr(0, position));
            reacts.erase(0, position+1);
        }
        reactvec.push_back(reacts.substr(0, std::string::npos));
        for (auto &specie : reactvec) {
            std::remove_if(specie.begin(), specie.end(), isspace);
        }
        
        while ((position = prods.find('+')) != std::string::npos) {
            prodvec.push_back(prods.substr(0, position));
            prods.erase(0, position+1);
        }
        prodvec.push_back(prods.substr(0, std::string::npos));
        for (auto &specie : prodvec) {
            std::remove_if(specie.begin(), specie.end(), isspace);
        }
    }

    parsedRxn tmprxn;
    { // Parse the BNGL
        // Parse reactants
        std::vector<parsedMol> reactmols;
        for (auto &specie : reactvec) {
            int i{0};
            char buffer[100];
            parsedMol tmpmol;

            // TODO: need to account for ! in the reactants
            for (std::string::iterator it = specie.begin(); it != specie.end(); it++) {
                if (isalnum(*it)) {
                    buffer[i] = (*it);
                }
                else if ((*it) == '(') {
                    buffer[i] = '\0';
                    std::string mol(buffer);
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    tmpmol.name = mol;
                    i = -1;
                }
                else if ((*it) == ',') {
                    buffer[i] = '\0';
                    std::string iface(buffer);
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    tmpmol.ifaces.push_back(iface);
                    i = -1;
                }
                else if ((*it) == ')') {
                    buffer[i] = '\0';
                    std::string iface(buffer);
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    tmpmol.ifaces.push_back(iface);
                    i = -1;
                }
                i++;
            }
            reactmols.push_back(tmpmol);
        }

        // Parse products
        std::vector<parsedMol> prodmols;
        for (auto &specie : prodvec) {
            int i{0};
            char buffer[100];
            buffer[0] = '\0';
            parsedMol tmpmol;

            // TODO: need to account for ! in the reactants
            for (std::string::iterator it = specie.begin(); it != specie.end(); it++) {
                std::string tmp(buffer);
                if (isalnum(*it)) {
                    buffer[i] = (*it);
                }
                else if ((*it) == '(') {
                    buffer[i] = '\0';
                    std::string mol(buffer);
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    tmprxn.reactants.push_back(mol);
                    i = -1;
                }
                else if ((*it) == '!') {
                    buffer[i] = '\0';
                    std::string iface(buffer);
                    std::fill(buffer, buffer +sizeof(buffer), '\0');
                    tmprxn.ifaces.push_back(iface);
                    it++;
                    tmprxn.bonds.push_back((*it) - '0'); // this weirdness is to convert char to int literally, not to ASCII
                    i = -1;
                }
                else if ((*it) == ',') { // flush unused ifaces
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    i = -1;
                }
                else if ((*it) == '.') { // flush when discovering new bound specie
                    std::fill(buffer, buffer + sizeof(buffer), '\0');
                    i = -1;
                }
                i++;
            }
        }
    }

    
    Reaction newrxn;
    { // Determine mol and iface from wholep which corresponds with parsed mol
        // Add reactant iface(s)
        for (unsigned i{0}; i < tmprxn.reactants.size(); i++) {
            for (unsigned j{0}; j < wholep.size(); j++) {
                if (tmprxn.reactants[i] == wholep[j].name) {
                    for (unsigned k{0}; k < wholep[j].Interfaces.size(); k++) {
                        if (tmprxn.ifaces[i] == wholep[j].Interfaces[k].name)
                            newrxn.Reactants.push_back(wholep[j].Interfaces[k].index);
                    }
                }
            }
        }
        // Add product iface(s)
        for (auto &Rxn : forRxns) {
            for (auto &specie : prodvec) {
                std::vector<std::string>::iterator pos = std::find(Rxn.prodspecies.begin(), Rxn.prodspecies.end(), specie); 
                if (pos != Rxn.prodspecies.end()) {
                    std::cout << "Reaction with duplicate products found." << std::endl;
                    //do something
                }
                else {
                    for (unsigned l{0}; l < prodvec.size(); l++) {
                        newrxn.Products.push_back(totspecies);
                        totspecies++;
                    }
                }
            }
        }
        if (forRxns.size() == 0) {
            for (unsigned l{0}; l < prodvec.size(); l++) {
                newrxn.Products.push_back(totspecies);
                totspecies++;
            }
        }
    }
    forRxns.push_back(newrxn);
    return revflag;
}

int main() {
    std::string filename = "rxn2.inp";
    std::ifstream infile;
    infile.open(filename);

    // TMP VARS //
    std::vector<Reaction> Rxnlist;
    std::vector<Reaction> forRxns;
    std::vector<Reaction> backRxns;
    std::vector<Molecule> wholep;
    int totspecies{5}; // should start out as number of interfaces-1
    for (int i{0}; i < 1; i++) {  
        std::string name = "MolecX";
        std::vector<Interface> ifaces;
        Molecule tmpmol(0, name, ifaces);
        Interface iface1("Site1", 0);
        Interface iface2("Site2", 1);
        Interface iface3("Site3", 2);
        tmpmol.Interfaces.push_back(iface1);
        tmpmol.Interfaces.push_back(iface2);
        tmpmol.Interfaces.push_back(iface3);
        wholep.push_back(tmpmol);
    }
    for (int i{0}; i < 1; i++) {  
        std::string name = "Mol1";
        std::vector<Interface> ifaces;
        Molecule tmpmol(0, name, ifaces);
        Interface iface1("Site11", 3);
        Interface iface2("Site12", 4);
        tmpmol.Interfaces.push_back(iface1);
        tmpmol.Interfaces.push_back(iface2);
        wholep.push_back(tmpmol);
    }
    // END TMP VARS //
    
    if (!infile) {
        std::cerr << "Reaction file cannot be opened. Exiting..." << std::endl;
        exit(1);
    }
    else {
        int numrxns = 2;
        int numrxns2;

        infile.ignore(500, '\n');
        infile >> numrxns2;
        if (numrxns2 != numrxns) {
            std::cout << "Error, wrong number of reactions!" << std::endl;
            exit(1);
        }
        infile.ignore(500, '\n');

        for (int i{0}; i < numrxns; i++) {
            infile.ignore(500, '\n');

            // Parse reaction
            std::string reaction;
            getline(infile, reaction);
            bool revflag = parse_reaction(totspecies, reaction, wholep, forRxns);
            
            { // Read onrate and offrate
                if (infile >> forRxns[i].onrate) {
                    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                else {
                    std::cerr << "Error: cannot read onrate. Exiting..." << std::endl;
                    exit(1);
                }

                if (infile >> forRxns[i].offrate) {
                    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
                else {
                    std::cerr << "Error: cannot read offrate. Exiting..." << std::endl;
                }
            }

            { // Read sigma
                std::string sigma;
                getline(infile, sigma);
                size_t pos = sigma.find('#');
                sigma = sigma.substr(0, pos);
                for (unsigned j{0}; j < 3; j++) {
                    pos = sigma.find(',');
                    double token = std::stod(sigma.substr(0, pos));
                    forRxns[i].sigma[j] = token;
                    sigma.erase(0, pos + 1);
                }
            }

            { // Read angles
                std::string angles;
                getline(infile, angles);
                size_t pos = angles.find('#');
                angles = angles.substr(0, pos);
                for (unsigned j{0}; j < 3; j++) {
                    pos = angles.find(',');
                    double token = std::stod(angles.substr(0, pos));
                    forRxns[i].angles[j] = token;
                    angles.erase(0, pos + 1);
                }
            }

            { // Read membrane status
                std::string memstat;
                getline(infile, memstat);
                size_t pos = memstat.find('#');
                int memstat2 = std::stoi(memstat.substr(0, pos));
                if (memstat2 == 0)
                    forRxns[i].onmem = false;
                else if (memstat2 == 1)
                    forRxns[i].onmem = true;
                else {
                    std::cerr << "Error, cannot read membrane status. Exiting..." << std::endl;
                    exit(1);
                }
            }

            if (revflag == true) { // create back reaction
                Reaction tmprxn;
                tmprxn.Products = forRxns[i].Reactants;
                tmprxn.Reactants = forRxns[i].Products;
                tmprxn.onrate = forRxns[i].offrate;
                tmprxn.offrate = forRxns[i].onrate;
                tmprxn.sigma = forRxns[i].sigma;
                tmprxn.angles = forRxns[i].angles;
                tmprxn.onmem = forRxns[i].onmem;
                backRxns.push_back(tmprxn);
            }
        }
    }

    Rxnlist = forRxns;
    Rxnlist.insert(Rxnlist.end(), backRxns.begin(), backRxns.end());

    for (auto &Rxn : Rxnlist) 
        Rxn.display();
}
