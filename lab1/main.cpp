#include <any>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

using mojamapa_t = std::map<std::string, std::vector<double>>;
using mojafunkcja_t = std::function<std::double_t (std::vector<double>)>;

void wypisz(mojamapa_t mapa, mojafunkcja_t fun) {
  using namespace std;
  for (auto kv : mapa) {
    auto [k, v] = kv;
    cout << "klucz: " << k << "; wartosc " << fun(v) << endl;
  }
}

int main(int argc, char **argv) {
  using namespace std;
  map<string, mojafunkcja_t> formatery;
  formatery["sin"] = [](vector<double> x) { return sin(x[0]); };
  formatery["add"] = [](vector<double> x) { return x[0] + x[1]; };
  formatery["mod"] = [](vector<double> x) { return (int)x[0] % (int)x[1]; };
  formatery[""] = [](vector<double> x) {
      cout << "sin x - oblicza sinus od x\n add x y - oblicza sume\n mod x y - oblicza x%y";
      return 0;
  };

  try {
    vector<string> argumenty(argv, argv + argc);
    string klucz = argumenty[1];

    vector<double> wartosci;
    for (int i = 2; i < argc; i++) {
      wartosci.push_back(stod(argumenty[i]));
    }
    mojamapa_t mojamapa = {{klucz, wartosci}};
    wypisz(mojamapa, formatery[klucz]);
    } catch (exception &e) {
    cout << "Blad: " << e.what() << endl;
      return 1;
    }
    return 0;
}
