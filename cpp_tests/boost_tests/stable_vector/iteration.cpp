/*
 * Created: 22-May-2018 08:55:12 AM EDT
 * Modified: 22-May-2018 09:30:52 AM EDT
 * Created by: Matthew Varga
 * Purpose: Simple code to test iteration persistence after element deletion in boost::vector
 */

#include <boost/container/stable_vector.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

int main()
{
    using boost::container::stable_vector;

    // stable_vector of strings
    stable_vector<std::string> test{ "cow", "pig", "puppy", "bird", "elephant" };

    // iterate over and print entire container
    std::cout << "Before deletion, entire container:\n";
    for (auto vecItr = test.begin(); vecItr != test.end(); ++vecItr)
        std::cout << *vecItr << '\n';

    // find iterators corresponding to "pig" and "elephant"
    stable_vector<stable_vector<std::string>::iterator> singleItrs;
    singleItrs.emplace_back(
        std::find_if(test.begin(), test.end(), [&](std::string animal) -> bool { return std::move(animal) == "pig"; }));
    singleItrs.emplace_back(std::find_if(
        test.begin(), test.end(), [&](std::string animal) -> bool { return std::move(animal) == "elephant"; }));

    // print them
    std::cout << "\nBefore deletion, single iterators:";
    for (auto& itr : singleItrs)
        std::cout << ' ' << *itr;

    std::cout << '\n' << std::setw(15) << std::setfill('-') << ' ' << std::setfill(' ') << '\n';
    // erase the first element
    test.erase(test.begin());

    // iterate over and print the entire container, after deletion
    std::cout << "After deletion, entire container:\n";
    for (auto vecItr = test.begin(); vecItr != test.end(); ++vecItr)
        std::cout << *vecItr << '\n';

    std::cout << "\nAfter deletion, single iterators:";
    for (auto& itr : singleItrs)
        std::cout << ' ' << *itr;
}
