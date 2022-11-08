#include <iostream>
#include <vector>
#include <functional>
#include <random>

std::random_device rd;
std::mt19937 mt_generator(rd());

using namespace std;

using chromosome_t = std::vector<int>;          //osobnik

auto ackley = [](double x, double y) -> double {
    using namespace std;
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x * x + y * y))) - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + M_E + 20.0;
};

auto decode = [](const chromosome_t chromosome){
    double x = 0.0;
    double y = 0.0;
    for (int i = 0; i < chromosome.size()/2; i++) {
        x += x + chromosome[i];
    }
    for (int i = chromosome.size()/2; i < chromosome.size(); i++) {
        y += y + chromosome[i];
    }
    return std::make_pair(x, y);
};

//fitness function
auto fitness = [](const chromosome_t& chromosome) -> double {
    auto [x, y] = decode(chromosome);
    return 1.0/(1.0 + std::abs(ackley(x, y)));
};

int main(){
    std::uniform_int_distribution<int> dist(0, 1);
    for (int i = 0; i < 10; i++) {
        chromosome_t chromosome(100+(23136%10*2));
        for (int j = 0; j < chromosome.size(); j++) {
            chromosome[j] = dist(mt_generator);
        }

        //nimam pojecia co ja robie ale dziala
        //miło by było gdyby ktos mi to wytłumaczył
        //w sumie nie moja specka

        auto xy = decode(chromosome);
        std::cout << "x: " << xy.first << "      y: " << xy.second << "      function: "
        << ackley(xy.first, xy.second) << "      fitness: " << fitness(chromosome) << std::endl;
    }

    //test ackley
    std::cout << "ackley(0,0): " << ackley(0, 0) << std::endl;
    std::cout << "ackley(1,1): " << ackley(1, 1) << std::endl;

    //test optimal point
    {
        chromosome_t chromosome(100+(23136%10*2));
        cout << "chromosome size: " << chromosome.size() << endl;
        for (int j = 0; j < chromosome.size(); j++) {
            chromosome[j] = 0;
        }
        auto xy = decode(chromosome);
        std::cout << "x: " << xy.first << "      y: " << xy.second << "      function: "
                  << ackley(xy.first, xy.second) << "      fitness: " << fitness(chromosome) << std::endl;
    }

    {
        chromosome_t chromosome(100+(23136%10*2));
        cout << "chromosome size: " << chromosome.size() << endl;
        for (int j = 0; j < chromosome.size(); j++) {
            chromosome[j] = 1;
        }
        auto xy = decode(chromosome);
        std::cout << "x: " << xy.first << "      y: " << xy.second << "      function: "
                  << ackley(xy.first, xy.second) << "      fitness: " << fitness(chromosome) << std::endl;
    }





    return 0;
}