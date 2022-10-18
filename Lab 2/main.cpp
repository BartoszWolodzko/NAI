#include <iostream>
#include <vector>
#include <functional>
#include <random>

using domain_t = std::vector<double>;
std::random_device rd;
std::mt19937 mt_generator(rd());

domain_t hill_climbing(const std::function<double(domain_t)> &f, domain_t minimal_d, domain_t maximal_d, int max_iterations){
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    domain_t current_d = minimal_d;
    for(int i = 0; i < max_iterations; i++){
        domain_t new_d = current_d;
        for(int j = 0; j < new_d.size(); j++){
            new_d[j] = current_d[j] + distribution(mt_generator) * (maximal_d[j] - minimal_d[j]);
        }
        if(f(new_d) < f(current_d)){
            current_d = new_d;
        }
    }
    return current_d;
}
int main() {
    //rosenbrock hill climbing
    auto rosenbrock = [](domain_t p) {
        double sum = 0;
        for (int i = 0; i < p.size() - 1; i++) {
            sum += 100 * (p[i + 1] - p[i] * p[i]) * (p[i + 1] - p[i] * p[i]) + (p[i] - 1) * (p[i] - 1);
        }
        return sum;
    };
    domain_t minimal_d = {-2.048, -2.048};
    domain_t maximal_d = {2.048, 2.048};
    auto rosenbrock_hill_climbing = [&rosenbrock, &minimal_d, &maximal_d](int max_iterations) {
        return hill_climbing(rosenbrock, minimal_d, maximal_d, max_iterations);
    };

    std::cout << "Rosenbrock function" << std::endl;
    for (int i = 0; i < 10; i++) {
        auto p = rosenbrock_hill_climbing(1000);
        std::cout << "x = " << p[0] << ", y = " << p[1] << ", f(x, y) = " << rosenbrock(p) << std::endl;
    }

    //ackley hill climbing
    auto ackley = [](domain_t p) {
        double sum1 = 0;
        double sum2 = 0;
        for (int i = 0; i < p.size(); i++) {
            sum1 += p[i] * p[i];
            sum2 += cos(2 * M_PI * p[i]);
        }
        return -20 * exp(-0.2 * sqrt(sum1 / p.size())) - exp(sum2 / p.size()) + 20 + exp(1);
    };
    minimal_d = {-32.768, -32.768};
    maximal_d = {32.768, 32.768};
    auto ackley_hill_climbing = [&ackley, &minimal_d, &maximal_d](int max_iterations) {
        return hill_climbing(ackley, minimal_d, maximal_d, max_iterations);
    };

    std::cout << "Ackley function" << std::endl;
    for (int i = 0; i < 10; i++) {
        auto p = ackley_hill_climbing(1000);
        std::cout << "x = " << p[0] << ", y = " << p[1] << ", f(x, y) = " << ackley(p) << std::endl;
    }

    //griewank hill climbing
    auto griewank = [](domain_t p) {
        double sum = 0;
        double product = 1;
        for (int i = 0; i < p.size(); i++) {
            sum += p[i] * p[i];
            product *= cos(p[i] / sqrt(i + 1));
        }
        return sum / 4000 - product + 1;
    };
    minimal_d = {-600, -600};
    maximal_d = {600, 600};
    auto griewank_hill_climbing = [&griewank, &minimal_d, &maximal_d](int max_iterations) {
        return hill_climbing(griewank, minimal_d,maximal_d, max_iterations);
    };

    std::cout << "Griewank function" << std::endl;
    for (int i = 0; i < 10; i++) {
        auto p = griewank_hill_climbing(1000);
        std::cout << "x = " << p[0] << ", y = " << p[1] << ", f(x, y) = " << griewank(p) << std::endl;
    }
    return 0;
}