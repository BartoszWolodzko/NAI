#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <cmath>

std::random_device rd;
std::mt19937 mt_generator(rd());

auto domain_generator_simulated_annealing = [](auto p) {
    std::normal_distribution<double> n(0.0, 0.3);
    for (auto& e : p) {
        e = e + n(mt_generator);
    }
    return p;
};

auto domain_generator = [](int a, int b) {
    std::uniform_real_distribution<> dis(a, b);
    return std::pair<double, double>(dis(mt_generator), dis(mt_generator));
};

auto get_close_points_random = [](auto p0) -> std::vector<double> {
    std::uniform_real_distribution<double> distr(-1, 1);
    return {{p0.first + distr(mt_generator), p0.second+distr(mt_generator)}};
};

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
        const std::function<double(std::pair<double, double>)> &f,
        const std::function<std::pair<double, double>(int, int)> &domain,
        int iterations, int a, int b) {

    auto current_p = domain(a, b);
    auto best_point = current_p;
    for (int i = 0; i < iterations; ++i) {
        auto new_p = domain(a, b);
        if (f(new_p) < f(current_p)) {
            current_p = new_p;
        }
        if (f(current_p) < f(best_point)) {
            best_point = current_p;
        }
    }
    return best_point;
};

auto simulated_annealing = [](
        const std::function<double(std::pair<double, double> pair)> &f,
        const std::function<std::pair<double, double>(int, int)> &domain,
        int K, int a, int b) {

    std::pair<double, double> current_p = domain(a, b);
    std::pair<double, double> best_p = current_p;


    std::uniform_real_distribution<double> uk(0, 1);

    for (int i = 0; i < K; ++i) {

        std::vector<double> p;
        p.push_back(current_p.first);
        p.push_back(current_p.second);
        auto tmp_new_p = domain_generator_simulated_annealing(p);
        std::pair<double, double> new_p = std::make_pair(tmp_new_p[0], tmp_new_p[1]);

        //sorry for casting so much times doing this fast

        if (f(current_p) < f(best_p)) {
            best_p = current_p;
        } else {
            if (uk(mt_generator) < exp(-(f(new_p) - f(current_p)) / K)) {
                current_p = new_p;
            }
        }
    }
    return best_p;
};

void execute(const std::function<double(std::pair<double, double>)> &f, const std::function<std::vector<double>(std::pair<double, double>)> &get_close_points,int a, int b) {
    int iterations = 1000000;

    auto brute_force_result = brute_force(f, domain_generator, iterations, a, b);
    auto hill_climbing_result = hill_climbing(f, domain_generator, iterations, a, b);
    auto simulated_annealing_result = simulated_annealing(f, domain_generator, iterations, a, b);
    using namespace std;
    cout << "brute f(" << brute_force_result.first << ", " << brute_force_result.second << ") = " << f(brute_force_result)<< endl;
    cout << "hill f(" << hill_climbing_result.first << ", " << hill_climbing_result.second << ") = " << f(hill_climbing_result)<< endl;
    cout << "simulated f(" << simulated_annealing_result.first << ", " << simulated_annealing_result.second << ") = "<< f(simulated_annealing_result) << endl;
    cout << endl;
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
    execute(ackley_f, get_close_points_random, -5, 5);

    cout << "himmelblau" << endl<<endl;
    execute(himmelblau_f, get_close_points_random, -5, 5);

    cout << "holderTable" << endl<<endl;
    execute(holderTable_f, get_close_points_random, -10, 10);

    return 0;
}