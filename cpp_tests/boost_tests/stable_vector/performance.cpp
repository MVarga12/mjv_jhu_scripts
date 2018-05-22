/*
 * Created: 18-May-2018 08:45:53 AM EDT
 * Modified: 22-May-2018 08:45:39 AM EDT
 * Created by: Matthew Varga
 * Purpose: testing stable_vector performance
 */

#include <boost/container/stable_vector.hpp>
#include <chrono>
#include <iostream>
#include <vector>

struct A {
    std::vector<double> data;

    A() {
        data.assign(1000, 0);
    }
};

int main()
{
    // get initial time
    std::vector<A> vec_t;

    std::chrono::steady_clock::time_point start_v1 = std::chrono::steady_clock::now();
    for (int i{ 0 }; i < 10000; ++i)
        vec_t.emplace_back();

    std::chrono::steady_clock::time_point end_v1 = std::chrono::steady_clock::now();
    std::cout << "Adding to std::vector took "
              << std::chrono::duration<double, std::milli>(end_v1 - start_v1).count() << " ms\n";

    using boost::container::stable_vector;
    stable_vector<A> vec_st;

    auto start_v2 = std::chrono::steady_clock::now();
    for (int i{ 0 }; i < 10000; i++)
        vec_st.emplace_back();

    auto end_v2 = std::chrono::steady_clock::now();
    std::cout << "Adding to boost::stable_vector took "
              << std::chrono::duration<double, std::milli>(end_v2 - start_v2).count() << " ms\n";

    start_v1 = std::chrono::steady_clock::now();
    for (int i{ 0 }; i < 10000; ++i)
        ;

    end_v1 = std::chrono::steady_clock::now();
    std::cout << "Iterating over std::vector took "
              << std::chrono::duration<double, std::milli>(end_v1 - start_v1).count() << " ms\n";

    start_v2 = std::chrono::steady_clock::now();
    for (int i{ 0 }; i < 10000; i++)
        ;
    end_v2 = std::chrono::steady_clock::now();
    std::cout << "Iterating over boost::stable_vector took "
              << std::chrono::duration<double, std::milli>(end_v2 - start_v2).count() << " ms\n";

    start_v1 = std::chrono::steady_clock::now();
    vec_t.erase(vec_t.begin(), vec_t.end());
    end_v1 = std::chrono::steady_clock::now();
    std::cout << "Erasing std::vector took "
              << std::chrono::duration<double, std::milli>(end_v1 - start_v1).count() << " ms\n";

    start_v2 = std::chrono::steady_clock::now();
    vec_st.erase(vec_st.begin(), vec_st.end());
    end_v2 = std::chrono::steady_clock::now();
    std::cout << "Erasing boost::stable_vector took "
              << std::chrono::duration<double, std::milli>(end_v2 - start_v2).count() << " ms\n";
}
