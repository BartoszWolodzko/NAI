#include <iostream>
#include <vector>
#include <functional>
#include <random>

using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
const int chromosome_size_const = 150;

auto ackley = [](double x, double y) -> double {
    using namespace std;
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + M_E +
           20.0;
};

const double ackley_domain_min = -32.768;
const double ackley_domain_max = 32.768;

std::random_device rd;
std::mt19937 mt_generator(rd());

std::pair<double,double> decode_chromosome(const chromosome_t &chromosome) {
    using namespace std;
    double x = 0.0, y = 0.0;
    for (int i = 0; i < chromosome_size_const/2; ++i) {
        x += chromosome[i] * pow(2, i);
        y += chromosome[i + chromosome_size_const/2] * pow(2, i);
    }

    x = ackley_domain_min + (ackley_domain_max - ackley_domain_min) * x / (pow(2, chromosome_size_const/2) - 1);
    y = ackley_domain_min + (ackley_domain_max - ackley_domain_min) * y / (pow(2, chromosome_size_const/2) - 1);

    return {x, y};
}

double ackley_fitness(const chromosome_t &chromosome) {
    auto [x, y] = decode_chromosome(chromosome);
    return 100 - ackley(x, y);
}

std::vector<double> fintess_function(population_t pop) {
    std::vector<double> fitness;
    for (auto chromosome: pop) {
        fitness.push_back(ackley_fitness(chromosome));
    }
    return fitness;
}

chromosome_t mutation(chromosome_t parents, double p_mutation) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0, 1.0);
    for (auto &gene: parents) {
        if (uniform(mt_generator) < p_mutation) {
            gene = 1 - gene;
        }
    }
    return parents;
}

std::vector<chromosome_t> crossover(std::vector<chromosome_t> parents) {
    using namespace std;
    uniform_int_distribution<int> uniform(0, chromosome_size_const - 1);
    int crossover_point = uniform(mt_generator);
    for (int i = crossover_point; i < chromosome_size_const; ++i) {
        swap(parents[0][i], parents[1][i]);
    }
    return parents;
}

std::vector<int> selection(std::vector<double> fitnesses) {
    using namespace std;
    discrete_distribution<int> discrete(fitnesses.begin(), fitnesses.end());
    vector<int> parents_indexes;
    for (int i = 0; i < fitnesses.size(); ++i) {
        parents_indexes.push_back(discrete(mt_generator));
    }
    return parents_indexes;
}

double min_fitness(const std::vector<double> &fitnesses) {
    double min = fitnesses[0];
    for (auto fitness: fitnesses) {
        if (fitness < min) min = fitness;
    }
    return min;
}

double avg_fitness(const std::vector<double> &fitnesses) {
    double sum = 0.0;
    for (auto fitness: fitnesses) {
        sum += fitness;
    }
    return sum / fitnesses.size();
}

double max_fitness(const std::vector<double> &fitnesses) {
    double max = fitnesses[0];
    for (auto fitness: fitnesses) {
        if (fitness > max) max = fitness;
    }
    return max;
}

double standard_deviation(const std::vector<double> &fitnesses) {
    double avg = avg_fitness(fitnesses);
    double sum = 0.0;
    for (auto fitness: fitnesses) {
        sum += (fitness - avg) * (fitness - avg);
    }
    return sqrt(sum / fitnesses.size());
}

bool term_condition(const std::vector<double> &fitnesses, double termination_value) {
    return standard_deviation(fitnesses) < termination_value;
}

population_t initial_population(int population_size) {
    using namespace std;
    uniform_int_distribution<int> uniform(0, 1);
    population_t population;
    for (int i = 0; i < population_size; ++i) {
        chromosome_t chromosome;
        for (int j = 0; j < chromosome_size_const; ++j) {
            chromosome.push_back(uniform(mt_generator));
        }
        population.push_back(chromosome);
    }
    return population;
}

void print_best_chromosome(population_t population) {
    auto fitnesses = fintess_function(population);
    int best_chromosome_index = 0;
    for (int i = 0; i < fitnesses.size(); ++i) {
        if (fitnesses[i] > fitnesses[best_chromosome_index]) {
            best_chromosome_index = i;
        }
    }
    auto [x, y] = decode_chromosome(population[best_chromosome_index]);
    std::cout << "x: " << x << " y: " << y << " fitness: " << fitnesses[best_chromosome_index] << std::endl;
}

void print_stats(const std::vector<double> &fitnesses) {
    std::cout << "min: " << min_fitness(fitnesses) << " avg: " << avg_fitness(fitnesses) << " max: "
              << max_fitness(fitnesses) << " std: " << standard_deviation(fitnesses) << std::endl;
}

auto genetic_algorithm = [](
        auto initial_population, auto fitness, auto term_condition,
        auto selection, double p_crossover,
        auto crossover, double p_mutation, auto mutation, double termination_value) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0, 1.0);
    auto population = initial_population;
    vector<double> population_fit = fitness(population);
    while (!term_condition(population_fit, termination_value)) {
        auto parents_indexes = selection(population_fit);
        decltype(population) new_population;
        for (int i = 0; i < parents_indexes.size(); i += 2) {
            decltype(initial_population) offspring = {population[i], population[i + 1]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring);
            }
            for (auto chromosome: offspring) new_population.push_back(chromosome);
        }
        for (auto &chromosome: new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }
        population = new_population;
        population_fit = fitness(population);

        print_stats(population_fit);
        print_best_chromosome(population);

    }
    return population;
};

int main() {
    using namespace std;

    auto result = genetic_algorithm(initial_population(100),
                                    fintess_function, term_condition,
                                    selection, 0.5,
                                    crossover,
                                    0.01, mutation, 0.5);

    for (chromosome_t chromosome: result) {
        cout << "[";
        for (int p: chromosome) {
            cout << p;
        }
        cout << "] ";
    }
    cout << endl;
    return 0;
}