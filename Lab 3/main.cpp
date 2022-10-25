#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <cmath>

std::random_device rd;
std::mt19937 mt_generator(rd());

auto brute_force = [](
        const std::function<double(std::pair<double, double>)> &f,
        const std::function<std::pair<double, double>(int, int)> &domain,
        int iterations, int a, int b) {

    auto current_p = domain(a, b);
    auto best_point = current_p;
    for (int i = 0; i < iterations; ++i) {
        if (f(current_p) < f(best_point)) {
            best_point = current_p;
        }
        current_p = domain(a, b);
    }
    return best_point;
};

auto hill_climbing = [](
        const std::function<double(std::pair<double, double> pair)> &f,
        const std::function<std::pair<double, double>(int, int)> &start_point,
        const std::function<std::vector<std::pair<double, double>>(std::pair<double, double>, int, int)> &get_close_points,
int max_iterations, int a, int b) {

std::pair<double, double> best_p = start_point(a, b);
for (int iteration = 0; iteration < max_iterations; iteration++) {
auto close_points = get_close_points(best_p, a, b);
auto best_neighbour = *min_element(
        close_points.begin(),
        close_points.end(),
        [f](auto a, auto b) {
            return f(a) > f(b);
        }
);
if (f(best_neighbour) < f(best_p))
best_p = best_neighbour;
}
return best_p;
};

auto simulated_annealing = [](
        const std::function<double(std::pair<double, double> pair)> &f,
        const std::function<std::pair<double, double>(int, int)> &domain,
        int iterations, int a, int b) {

    std::vector<std::pair<double, double>> pairsVector;
    std::uniform_real_distribution<double> uk(0, 1);
    double ukValue = uk(mt_generator);
    auto best_point = domain(a, b);
    pairsVector.push_back(best_point);
    for (int i = 0; i < iterations; ++i) {
        auto tk = domain(a, b);
        if (f(tk) <= f(best_point)) {
            best_point = tk;
            pairsVector.push_back(best_point);
        } else {
            if (ukValue < exp(-(abs(f(tk) - f(best_point)) / (1 / log(i))))) {
                best_point = tk;
                pairsVector.push_back(best_point);
            }
        }
    }
    return best_point;
};

auto domain_generator = [](int a, int b) {
    std::uniform_real_distribution<> dis(a, b);
    return std::pair<double, double>(dis(mt_generator), dis(mt_generator));
};

auto close_points_generator = [](std::pair<double, double> p, int a, int b) -> std::vector<std::pair<double, double>> {
    std::uniform_real_distribution<> dis(a, b);
    return {std::pair<double, double>(dis(mt_generator), dis(mt_generator))};
};

void execute(const std::function<double(std::pair<double, double>)> &f,
             const std::function<std::pair<double, double>(int, int)> &domain,
             const std::function<std::vector<std::pair<double, double>>(std::pair<double, double>, int, int)> &get_close_points,
int a, int b) {

int iterations = 1000000;

auto brute_force_result = brute_force(f, domain, iterations, a, b);
auto hill_climbing_result = hill_climbing(f, domain, get_close_points, iterations, a, b);
auto simulated_annealing_result = simulated_annealing(f, domain, iterations, a, b);

std::cout << "brute f(" << brute_force_result.first << ", " << brute_force_result.second << ") = " << f(brute_force_result)
<< std::endl;
std::cout << "hill f(" << hill_climbing_result.first << ", " << hill_climbing_result.second << ") = " << f(hill_climbing_result)
<< std::endl;
std::cout << "simulated f(" << simulated_annealing_result.first << ", " << simulated_annealing_result.second << ") = "
<< f(simulated_annealing_result) << std::endl;
std::cout << std::endl;
}


int main() {
    using namespace std;
    auto ackley_f = [](pair<double, double> pair) {
        return -20.0 * exp(-0.2 * sqrt(0.5 * (pow(pair.first, 2) + pow(pair.second, 2)))) -
               exp(0.5 * (cos(2 * M_PI * pair.first) + cos(2 * M_PI * pair.second))) + exp(1) + 20;
    };

    auto himmelblau_f = [](pair<double, double> pair) {
        return pow(pow(pair.first, 2) + pair.second - 11, 2) + pow(pair.first + pow(pair.second, 2) - 7, 2);
    };

    auto holderTable_f = [](pair<double, double> pair) {
        return -abs(sin(pair.first) * cos(pair.second) *
                    exp(abs(1 - (sqrt(pow(pair.first, 2) + pow(pair.second, 2)) / M_PI))));
    };

    cout << "ackley" << endl<<endl;
    execute(ackley_f, domain_generator, close_points_generator, -5, 5);

    cout << "himmelblau" << endl<<endl;
    execute(himmelblau_f, domain_generator, close_points_generator, -5, 5);

    cout << "holderTable" << endl<<endl;
    execute(holderTable_f, domain_generator, close_points_generator, -10, 10);

    return 0;
}