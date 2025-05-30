
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


    const size_t SZ = size_t{1} << 23;

    std::vector<uint32_t> vec1(SZ);

    auto random_device = std::bind(std::uniform_int_distribution<uint32_t>{}, std::mt19937{});
    std::generate(vec1.begin(), vec1.end(), random_device);
    std::vector<uint32_t> vec2 = vec1;

    auto task1 = [&]{
        std::sort(vec1.begin(), vec1.end());
    };

    auto task2 = [&]{
        dg::sort_variants::quicksort::quicksort(vec2.data(), std::next(vec2.data(), vec2.size()));
    };

    std::cout << timeit(task1) << "<ms>" << "<std_qs>" << std::endl;
    std::cout << timeit(task2) << "<ms>" << "<qs>" << std::endl;

    assert(vec1 == vec2);

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

        std::vector<uint32_t> vec3 = vec;

        std::sort(vec.begin(), vec.end());
        dg::sort_variants::quicksort::quicksort(vec3.data(), std::next(vec3.data(), vec3.size()));

        if (vec != vec3){
            std::cout << "mayday" << std::endl;
            std::abort();
        }

        iter++;

        if (iter % 10000 == 0){
            std::cout << iter << std::endl;
        }
    }
}