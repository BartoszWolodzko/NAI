#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <map>
#include <cmath>
#include <chrono>

std::random_device rd;
std::mt19937 mt_generator(rd());

/**
 * domain - generate domain points. Throws exception when all the points were returned
 */
auto brute_force_or_random_probe = [](auto f, auto domain, int how_many=1000) {
    auto best = domain();
    auto best_value = f(best);
    try{
        for(int i=0;i<how_many;i++){
            auto point = domain();
            auto value = f(point);
            if(value < best_value){
                best = point;
                best_value = value;
            }
        }
    } catch(std::exception& e){
        return std::make_pair(best, best_value);
    }
    return std::make_pair(best, best_value);
};

using domain_t = std::vector<double>;
using my_function_t = std::function<double(domain_t)>;

int main() {
    std::map<std::string, my_function_t> domains;
    domains["sphere"] = [](domain_t x) {
        double sum = 0;
        for(auto& xi : x){
            sum += xi*xi;
        }
        return sum;
    };
    domains["rosenbrock"] = [](domain_t x) {
        double sum = 0;
        for(int i=0;i<x.size()-1;i++){
            sum += 100*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (1-x[i])*(1-x[i]);
        }
        return sum;
    };
    domains["ackley"] = [](domain_t x) {
        double sum1 = 0;
        double sum2 = 0;
        for(auto& xi : x){
            sum1 += xi*xi;
            sum2 += std::cos(2*M_PI*xi);
        }
        return -20*std::exp(-0.2*std::sqrt(sum1/x.size())) - std::exp(sum2/x.size()) + 20 + std::exp(1);
    };
    domains["gomez_and_levy"] = [](domain_t x) {
        double sum = 0;
        for(int i=0;i<x.size()-1;i++){
            sum += std::pow(std::sin(3*M_PI*x[i]),2) + std::pow(x[i]-1,4)*(1+std::pow(std::sin(3*M_PI*x[i+1]),2)) + std::pow(x[i+1]-1,2)*(1+std::pow(std::sin(2*M_PI*x[i+1]),2));
        }
        return sum;
    };
    /*auto rosenbrock_disk_f = [](auto x){
        return (1 - x[0]) * (1 - x[0]) + 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };*/

    auto domain = [](){
        std::uniform_real_distribution<double> dist(-2.048, 2.048);
        return domain_t{dist(mt_generator), dist(mt_generator)};
    };

    for (int i = 0; i < 20; ++i) {
        std::cout << "Pomiar "<<i+1 << " Ilosc probow "<< 100000*(i+1) << std::endl;
        auto begin = std::chrono::high_resolution_clock::now();
        auto result = brute_force_or_random_probe(domains["gomez_and_levy"], domain, 100000*(i+1));
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cout << "Best point: " << result.first[0] << ", " << result.first[1] << std::endl<<std::endl;
    }
    return 0;
}
/*
Pomiar 1 Ilosc probow 100000
Time measured: 0.067 seconds.
Best point: 0.997718, 0.98911

Pomiar 2 Ilosc probow 200000
Time measured: 0.132 seconds.
Best point: 0.998483, 1.01244

Pomiar 3 Ilosc probow 300000
Time measured: 0.198 seconds.
Best point: 0.999384, 0.999796

Pomiar 4 Ilosc probow 400000
Time measured: 0.273 seconds.
Best point: 1.00004, 1.01014

Pomiar 5 Ilosc probow 500000
Time measured: 0.327 seconds.
Best point: 0.999517, 1.00732

Pomiar 6 Ilosc probow 600000
Time measured: 0.391 seconds.
Best point: 0.998748, 0.994603

Pomiar 7 Ilosc probow 700000
Time measured: 0.459 seconds.
Best point: 1.0004, 0.994566

Pomiar 8 Ilosc probow 800000
Time measured: 0.524 seconds.
Best point: 0.999589, 1.00605

Pomiar 9 Ilosc probow 900000
Time measured: 0.587 seconds.
Best point: 1.00059, 0.992003

Pomiar 10 Ilosc probow 1000000
Time measured: 0.651 seconds.
Best point: 1.00046, 1.00354

Pomiar 11 Ilosc probow 1100000
Time measured: 0.714 seconds.
Best point: 1.00012, 0.998024

Pomiar 12 Ilosc probow 1200000
Time measured: 0.780 seconds.
Best point: 0.999636, 0.991534

Pomiar 13 Ilosc probow 1300000
Time measured: 0.843 seconds.
Best point: 0.999436, 1.00427

Pomiar 14 Ilosc probow 1400000
Time measured: 0.911 seconds.
Best point: 1.00011, 0.993225

Pomiar 15 Ilosc probow 1500000
Time measured: 0.974 seconds.
Best point: 0.999616, 1.00826

Pomiar 16 Ilosc probow 1600000
Time measured: 1.039 seconds.
Best point: 1.00035, 0.995099

Pomiar 17 Ilosc probow 1700000
Time measured: 1.102 seconds.
Best point: 1.00024, 1.0068

Pomiar 18 Ilosc probow 1800000
Time measured: 1.172 seconds.
Best point: 0.999905, 0.996703

Pomiar 19 Ilosc probow 1900000
Time measured: 1.249 seconds.
Best point: 1.00037, 0.997238

Pomiar 20 Ilosc probow 2000000
Time measured: 1.305 seconds.
Best point: 0.99988, 0.995512
 */
