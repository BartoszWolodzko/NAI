#include <iostream>
#include <vector>
#include <functional>
#include <random>

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

int main() {
    auto rosenbrock_disk_f = [](auto x){
        return (1 - x[0]) * (1 - x[0]) + 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
    };



    double min = -2.05;
    double max = 2.05;
    std::uniform_real_distribution<double> dist(min, max);

    auto rosenbrock_disk_domain = [&](){
        domain_t p(2);
        p[0] = dist(mt_generator);
        p[1] = dist(mt_generator);
        return p;
    };

    auto best_point = brute_force_or_random_probe(rosenbrock_disk_f, rosenbrock_disk_domain);
    std::cout << "Best point: " << best_point.first[0] << ", " << best_point.first[1] << std::endl;

    return 0;
}

