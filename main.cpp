
#include <atomic>
#include <random>
#include <memory>
#include <functional>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>
#include <iostream>
#include "assert.h"
#include "sort_variants.h"

template <class Task>
auto timeit(Task task) -> size_t{

    auto then = std::chrono::high_resolution_clock::now();
    task();
    auto now = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(now - then).count();
} 

int main(){

    // const size_t SZ = size_t{1} << 23;

    auto random_device_1    = std::bind(std::uniform_int_distribution<uint32_t>(0, 2048), std::mt19937{});
    auto random_device_2    = std::bind(std::uniform_int_distribution<uint32_t>(0, 256), std::mt19937{});
    size_t iter             = 0u; 

    while (true){
        const size_t SZ = random_device_1();
        std::vector<uint32_t> vec(SZ);

        if (random_device_1() % 2 == 0){
            std::generate(vec.begin(), vec.end(), std::ref(random_device_1));
        } else{
            std::generate(vec.begin(), vec.end(), std::ref(random_device_2));
        }

        std::vector<uint32_t> vec2 = vec;

        std::sort(vec.begin(), vec.end());
        dg::sort_variants::quicksort::quicksort(vec2.data(), std::next(vec2.data(), vec2.size()));

        if (vec != vec2){
            std::cout << "mayday" << std::endl;
            std::abort();
        }

        iter++;

        if (iter % 10000 == 0){
            std::cout << iter << std::endl;
        }
    }
}