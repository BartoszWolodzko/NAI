#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

std::random_device rd;
std::mt19937 mt_generator(rd());

const int chromosome_size = 100+(23136%10*2);
struct chromosome {

    std::vector<bool> genes_x = std::vector<bool>(chromosome_size/4);
    std::vector<bool> genes_x_fraction = std::vector<bool>(chromosome_size/4);
    std::vector<bool> genes_y = std::vector<bool>(chromosome_size/4);
    std::vector<bool> genes_y_fraction = std::vector<bool>(chromosome_size/4);
    double fitness;
};

//decodes the chromosome into a pair of double
std::pair<double, double> decode_chromosome(chromosome c) {
    using namespace std;
    double x = 0;
    double y = 0;
    for (int i = 0; i < chromosome_size/4; i++) {
        x += c.genes_x[i] * pow(2, chromosome_size/4 - i - 1);
        y += c.genes_y[i] * pow(2, chromosome_size/4 - i - 1);
    }
    for (int i = 0; i < chromosome_size/4; i++) {
        x += c.genes_x_fraction[i] * pow(2, -i - 1);
        y += c.genes_y_fraction[i] * pow(2, -i - 1);
    }
    return make_pair(x, y);
}

//encodes a chromosome
chromosome encode_chromosome(double x, double y) {
    chromosome c;
    //split the double into integer and fraction
    int x_int = floor(x);
    int y_int = floor(y);
    double x_fraction = x - x_int;
    double y_fraction = y - y_int;
    //change fraction to integer
    x_fraction *= pow(2, chromosome_size/4);
    y_fraction *= pow(2, chromosome_size/4);
    //convert the integer part to binary
    for (int i = 0; i < chromosome_size/4; i++) {
        c.genes_x[i] = x_int % 2;
        x_int /= 2;
        c.genes_y[i] = y_int % 2;
        y_int /= 2;
    }

    //convert the fraction part to binary
    for (int i = 0; i < chromosome_size/4; i++) {
        c.genes_x_fraction[i] = (int)x_fraction % 2;
        x_fraction /= 2;
        c.genes_y_fraction[i] = (int)y_fraction % 2;
        y_fraction /= 2;
    }
    reverse(c.genes_x.begin(), c.genes_x.end());
    reverse(c.genes_x_fraction.begin(), c.genes_x_fraction.end());
    reverse(c.genes_y.begin(), c.genes_y.end());
    reverse(c.genes_y_fraction.begin(), c.genes_y_fraction.end());
    return c;
}

//ackley function
auto ackley = [](double x, double y) -> double {
    using namespace std;
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + M_E + 20.0;
};

//ackley domain
double ackley_domain_min = -32.768;
double ackley_domain_max = 32.768;


//fitness function
double fitness(chromosome c , double(*f)(double, double)) {
    auto decoded = decode_chromosome(c);
    return 100-f(decoded.first, decoded.second);
    //return 100-ackley(decoded.first, decoded.second);
}

//generates a random chromosome
chromosome generate_chromosome(double min = ackley_domain_min, double max = ackley_domain_max, double(*f)(double, double) = ackley) {
    using namespace std;
    chromosome c;
    uniform_real_distribution<double> distribution(min, max);
    auto encoded = encode_chromosome(distribution(mt_generator), distribution(mt_generator));
    c.genes_x = encoded.genes_x;
    c.genes_x_fraction = encoded.genes_x_fraction;
    c.genes_y = encoded.genes_y;
    c.genes_y_fraction = encoded.genes_y_fraction;
    c.fitness = fitness(encoded, f);
    return c;
}

chromosome generate_chromosome_set(double x, double y, double(*f)(double, double) = ackley) {
    chromosome c;
    auto encoded = encode_chromosome(x, y);
    c.genes_x = encoded.genes_x;
    c.genes_x_fraction = encoded.genes_x_fraction;
    c.genes_y = encoded.genes_y;
    c.genes_y_fraction = encoded.genes_y_fraction;
    c.fitness = fitness(encoded,f);
    return c;
}

std::vector<chromosome> generate_population(int size, double min = ackley_domain_min, double max = ackley_domain_max, double(*f)(double, double) = ackley) {
    using namespace std;
    vector<chromosome> population(size);
    for (int i = 0; i < size; i++) {
        population[i] = generate_chromosome(min, max, f);
    }
    return population;
}

int main(){
    using namespace std;
    vector<chromosome> population = generate_population(100);
    for (int i = 0; i < population.size(); i++) {
        cout << "Chromosome " << i << ": ";
        cout << "x: " << decode_chromosome(population[i]).first;
        cout << "   y: " << decode_chromosome(population[i]).second;
        cout << "   fitness: " << population[i].fitness << endl;
    }

    //test
    cout << "Perfect Test assert true: ";
    chromosome c = generate_chromosome_set(0,0);
    cout << "x: " << decode_chromosome(c).first;
    cout << "   y: " << decode_chromosome(c).second;
    cout << "   fitness: " << c.fitness << endl;

    return 0;
}

