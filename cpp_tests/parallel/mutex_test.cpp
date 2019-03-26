#include <iostream>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>

// So this mutex object must be like a lock so all instances of addMutex.lock()/unlock()
// will talk to each other. But if we have subtractMutex, it probably won't talk to addMutex
std::mutex addMutex {};

// This would be safe, since we're using a mutex lock
void add_to_vector_safe(std::vector<int>& vec, int scal)
{
    addMutex.lock(); // lock and tell others to wait
    vec.push_back(scal); // push back value
    addMutex.unlock(); // unlock and return
}

// this would be unsafe, since there's no condition to wait
// if other threads try to access it at the same time
void add_to_vector_unsafe(std::vector<int>& vec, int scal)
{
    vec.push_back(scal);
}

int main() {
    std::vector<int> vec { };

    std::chrono::steady_clock::time_point start { std::chrono::steady_clock::now() };
    for (int i { 0 }; i < 10000; i += 2) {
        std::thread t1(&add_to_vector_safe, std::ref(vec), i); // start first thread
        std::thread t2(&add_to_vector_safe, std::ref(vec), i+1); // start second thread

        // join them together (should do whenever you're done with them, I guess?
        t1.join();
        t2.join();
    }
    std::cout << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';

    std::vector<int> vec2 { };
    std::chrono::steady_clock::time_point start2 { std::chrono::steady_clock::now() };
    for (int i { 0 }; i < 10000; ++i) {
        vec2.push_back(i);
    }
    std::cout << std::chrono::duration<double>(std::chrono::steady_clock::now() - start2).count() << '\n';
}
