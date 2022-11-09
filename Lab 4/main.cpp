#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

using namespace std;

std::random_device rd;
std::mt19937 mt_generator(rd());

// Add 1 to a binary number
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

//convert integer to binary u2
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

//convert binary u2 to integer
int binaryToInteger(vector<bool> binary) {
    long long int integer = 0;
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

using chromosome = std::vector<bool>;
const int chromosome_size_const = 100 + (23136 % 10 * 2);
const double double_precision = 1000000;

//encode the chromosome from a pair of double into u2
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

//decode the chromosome into a pair of double
std::pair<double, double> decode_chromosome(chromosome c) {
    vector<bool> x_binary(c.begin(), c.begin() + chromosome_size_const / 2);
    vector<bool> y_binary(c.begin() + chromosome_size_const / 2, c.end());
    int x = binaryToInteger(x_binary);
    int y = binaryToInteger(y_binary);
    return std::make_pair(x / double_precision, y / double_precision);
}


//ackley function
auto ackley = [](double x, double y) -> double {
    using namespace std;
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + M_E + 20.0;
};

//ackley domain
const double ackley_domain_min = -32.768;
const double ackley_domain_max = 32.768;

//fitness function
auto fitness = [](chromosome c) -> double {
    auto p = decode_chromosome(c);
    return 100-ackley(p.first, p.second);
};

//generate a random chromosome
chromosome generate_random_chromosome() {
    chromosome c;
    for (int i = 0; i < chromosome_size_const; i++) {
        c.push_back(mt_generator() % 2);
    }
    return c;
}

//generate a random population
std::vector<chromosome> generate_random_population(int population_size) {
    std::vector<chromosome> population;
    for (int i = 0; i < population_size; i++) {
        population.push_back(generate_random_chromosome());
    }
    return population;
}

//generate a random chromosome with encoding
chromosome generate_random_chromosome_encoding() {
    chromosome c;
    uniform_real_distribution<double> distribution(ackley_domain_min, ackley_domain_max);
    double x = distribution(mt_generator);
    double y = distribution(mt_generator);
    c = encode_chromosome(make_pair(x, y));
    return c;
}

//generate a random population with encoding
std::vector<chromosome> generate_random_population_encoding(int population_size) {
    std::vector<chromosome> population;
    for (int i = 0; i < population_size; i++) {
        population.push_back(generate_random_chromosome_encoding());
    }
    return population;
}

//generate chromosome with encoding
chromosome generate_chromosome_encoding(double x, double y) {
    chromosome c;
    c = encode_chromosome(make_pair(x, y));
    return c;
}

int main() {

    //generate a random population with encoding
    std::vector<chromosome> population = generate_random_population_encoding(10);
    for (auto c : population) {
        auto p = decode_chromosome(c);
        cout << "x: " << p.first << " y: " << p.second << " fitness: "<<fitness(c)<< endl;
    }

    //generate chromosome with encoding
    chromosome c = generate_chromosome_encoding(0, 0);
    cout << endl;
    auto p = decode_chromosome(c);
    cout << "x: " << p.first << " y: " << p.second << " fitness: " << fitness(c) << endl;
    return 0;
}
