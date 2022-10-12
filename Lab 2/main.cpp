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
        auto begin = std::chrono::high_resolution_clock::now();
        auto result = brute_force_or_random_probe(domains["gomez_and_levy"], domain, std::pow(4,i));
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        std::cout << "Best point: " << result.first[0] << ", " << result.first[1] << std::endl;
        std::cout << i << std::endl;
    }




    return 0;
}

