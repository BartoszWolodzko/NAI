#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

using namespace std;

std::random_device rd;
std::mt19937 mt_generator(rd());
using chromosome = std::vector<bool>;
const int chromosome_size_const = 100 + (23136 % 10 * 2);
const double double_precision = 10000000000000;

vector<bool> addOne(vector<bool> binary) {
    vector<bool> result = binary;
    for (int i = result.size() - 1; i >= 0; i--) {
        if (result[i] == 0) {
            result[i] = 1;
            break;
        }
        else {
            result[i] = 0;
        }
    }
    return result;
}

vector<bool> integerToBinary(int integer, int size) {
    if (integer > 0){
        vector<bool> binary(size, 0);
        for (int i = size - 1; i >= 0; i--) {
            if (integer % 2 == 1) {
                binary[i] = 1;
            }
            integer /= 2;
        }
        return binary;
    }
    else {
        integer = abs(integer);
        vector<bool> binary(size, 1);
        for (int i = size - 1; i >= 0; i--) {
            if (integer % 2 == 1) {
                binary[i] = 0;
            }
            integer /= 2;
        }
        return addOne(binary);
    }
}

long long binaryToInteger(vector<bool> binary) {
    long long integer = 0;
    if (binary[0] == 1){
        integer = -pow(2, binary.size()-1);
    }
    for (int i = 1; i < binary.size(); i++) {
        if (binary[i] == 1) {
            integer += pow(2, binary.size() - i-1);
        }
    }

    return integer;
}

chromosome encode_chromosome(std::pair<double, double> p) {
    chromosome result;
    int x = p.first * double_precision;
    int y = p.second * double_precision;
    vector<bool> x_binary = integerToBinary(x, chromosome_size_const / 2);
    vector<bool> y_binary = integerToBinary(y, chromosome_size_const / 2);
    result.insert(result.end(), x_binary.begin(), x_binary.end());
    result.insert(result.end(), y_binary.begin(), y_binary.end());
    return result;
}

std::pair<double, double> decode_chromosome(chromosome c) {
    vector<bool> x_binary(c.begin(), c.begin() + chromosome_size_const / 2);
    vector<bool> y_binary(c.begin() + chromosome_size_const / 2, c.end());
    long long x = binaryToInteger(x_binary);
    long long y = binaryToInteger(y_binary);
    return std::make_pair(x / double_precision, y / double_precision);
}

auto ackley = [](double x, double y) -> double {
    using namespace std;
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + M_E + 20.0;
};

const double ackley_domain_min = -32.768;
const double ackley_domain_max = 32.768;

auto fitness = [](chromosome c) -> double {
    auto p = decode_chromosome(c);
    return 100-ackley(p.first, p.second);
};

std::vector<chromosome> crossover(std::vector<chromosome> chromosomes) {
    using namespace std;
    uniform_int_distribution<int> uniform(0, chromosome_size_const - 1);
    int crossover_point = uniform(mt_generator);
    for (int i = crossover_point; i < chromosome_size_const; ++i) {
        swap(chromosomes[0][i], chromosomes[1][i]);
    }
    return chromosomes;
}

std::vector<chromosome> two_point_crossover(chromosome c1, chromosome c2) {
    using namespace std;
    uniform_int_distribution<int> uniform(0, chromosome_size_const - 1);
    int crossover_point1 = uniform(mt_generator);
    int crossover_point2 = uniform(mt_generator);
    if (crossover_point1 > crossover_point2) {
        swap(crossover_point1, crossover_point2);
    }
    for (int i = crossover_point1; i < crossover_point2; ++i) {
        swap(c1[i], c2[i]);
    }
    return { c1, c2 };
}

chromosome uniform_mutation(chromosome c, double mutation_rate) {
    for (int i = 0; i < chromosome_size_const; i++) {
        if (std::uniform_real_distribution<double> (0, 1)(mt_generator) < mutation_rate) {
            c[i] = !c[i];
        }
    }
    return c;
}

chromosome multi_point_mutation(chromosome c, int number_of_points_to_mutate){
    //select a random mutation point
    for (int i = 0; i < number_of_points_to_mutate; i++) {
        int mutation_point = std::uniform_int_distribution<int> (0, chromosome_size_const - 1)(mt_generator);
        //perform mutation
        c[mutation_point] = !c[mutation_point];
    }
    return c;
}

int roulette_wheel_selection(std::vector<double> population_fit) {
    double total_fitness = 0;

    for (auto f : population_fit) {
        total_fitness += f;
    }

    uniform_real_distribution<double> distribution(0, total_fitness);
    double random = distribution(mt_generator);
    double sum = 0;

    for (int i = 0; i < population_fit.size(); i++) {
        sum += population_fit[i];
        if (sum >= random) {
            return i;
        }
    }
    return 0;
}


vector<chromosome> initial_population(int population_size) {
    using namespace std;
    uniform_int_distribution<int> uniform(0, 1);
    vector<chromosome> population;
    for (int i = 0; i < population_size; ++i) {
        chromosome chromosome;
        for (int j = 0; j < chromosome_size_const; ++j) {
            chromosome.push_back(uniform(mt_generator));
        }
        population.push_back(chromosome);
    }
    return population;
}

vector<double> fintess_function(vector<chromosome> population) {
    vector<double> fitnesses;
    for (auto c : population) {
        fitnesses.push_back(fitness(c));
    }
    return fitnesses;
}

void print_best_chromosome(vector<chromosome> population) {
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

// select chromosome for crossover
int tournament_selection(std::vector<double> population_fit, int tournament_size) {
    std::vector<int> tournament;
    for (int i = 0; i < tournament_size; i++) {
        int random_index = std::uniform_int_distribution<int> (0, population_fit.size() - 1)(mt_generator);
        tournament.push_back(random_index);
    }
    int best_index = 0;
    for (int i = 0; i < tournament.size(); i++) {
        if (population_fit[tournament[i]] > population_fit[tournament[best_index]]) {
            best_index = i;
        }
    }
    return tournament[best_index];
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

void print_stats(const vector<chromosome> &population, const std::vector<double> &fitnesses) {
    std::cout << endl;
    std::cout << "Population size: " << population.size() << std::endl;
    std::cout << "Min fitness: " << min_fitness(fitnesses) << std::endl;
    std::cout << "Avg fitness: " << avg_fitness(fitnesses) << std::endl;
    std::cout << "Max fitness: " << max_fitness(fitnesses) << std::endl;
    std::cout << "Standard deviation: " << standard_deviation(fitnesses) << std::endl;
    print_best_chromosome(population);
    std::cout << endl;
}

auto genetic_algorithm = [](
        auto initial_population, auto fitness, auto term_condition,
        auto selection, double p_crossover,
        auto crossover, double p_mutation, auto mutation, double termination_value,bool show_progress,
        int number_of_generations) {
    using namespace std;

    uniform_real_distribution<double> uniform(0.0, 1.0);

    //initial population
    vector<chromosome> population = initial_population;
    //evaluate fitness
    vector<double> population_fit = fitness(population);

    int generation = 0;
    while (!term_condition(population_fit, termination_value) && generation < number_of_generations) {
        if(show_progress){
            cout << "----------------------------------"<<endl;
            cout << "generation: " << generation << endl;
            print_stats(population, population_fit);
        }

        vector<int> parents_indexes(population.size());
        vector<chromosome> new_population(population.size());
        // select parents
        transform(population_fit.begin(), population_fit.end(),
                  parents_indexes.begin(),
                  [&](auto e) { return selection(population_fit,6); });

        if(false){
            cout << "parents fitness: ";
            // calculate parents fitness
            vector<double> parents_fitness;
            for (auto i : parents_indexes) {
                parents_fitness.push_back(population_fit[i]);
            }
            // print parents fitness
            vector<chromosome> parents;
            for (auto i : parents_indexes) {
                parents.push_back(population[i]);
            }
            print_stats(parents, parents_fitness);
        }

        // perform crossover operations on parents
        for (int i = 0; i < parents_indexes.size() - 1; i += 2) {
            //take two parents
            vector<chromosome> offspring = {population[parents_indexes[i]], population[parents_indexes[i + 1]]};
            //perform crossover
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring[0], offspring[1]);
            }
            // replace parents with offspring
            new_population[i] = offspring[0];
            new_population[i + 1] = offspring[1];
        }

        // perform mutation operations on new population
        for (auto &chromosome : new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }

        // replace old population with new population
        population = new_population;
        // calculate fitness of new population
        population_fit = fitness(population);

        generation++;
        if (show_progress){
            cout << "----------------------------------"<<endl;
        }
    }
    return population;
};



int main(int argc, char **argv) {
    using namespace std;

    // 100 10000 0.01 0.7 1 0.01

    int population_size = (argc > 1) ? std::stoi(argv[1]) : 5000;
    int number_of_generations = (argc > 2) ? std::stoi(argv[2]) : 100000000;
    double mutation_rate = (argc > 3) ? std::stod(argv[3]) : 0.02;
    double crossover_rate = (argc > 4) ? std::stod(argv[4]) : 0.8;
    bool show_progress = (argc > 5) ? std::stoi(argv[5]) : true;
    double standard_deviation_value = (argc > 6) ? std::stod(argv[6]) : 0.001;

    auto result = genetic_algorithm(
            initial_population(population_size),
            fintess_function,
            term_condition,
            tournament_selection,
            crossover_rate,
            two_point_crossover,
            mutation_rate,
            uniform_mutation,
            standard_deviation_value,
            show_progress,
            number_of_generations);

    print_best_chromosome(result);
    return 0;
}