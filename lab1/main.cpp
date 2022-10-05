#include <any>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

using mojafunkcja_t = std::function<std::double_t (std::vector<double>)>;

void wypisz(std::vector<double> wartosci, std::string klucz) {
  using namespace std;
    map<string, mojafunkcja_t> formatery;
    formatery["sin"] = [](vector<double> x) { return sin(x[0]); };
    formatery["add"] = [](vector<double> x) { return x[0] + x[1]; };
    formatery["mod"] = [](vector<double> x) { return (int)x[0] % (int)x[1]; };
  cout << "Wynik "<< klucz <<": " << formatery[klucz](wartosci) << endl;
}

int main(int argc, char **argv) {
  using namespace std;
  try {
    vector<string> argumenty(argv, argv + argc);
    string klucz = argumenty[1];
    vector<double> wartosci;
    for (int i = 2; i < argc; i++) {
      wartosci.push_back(stod(argumenty[i]));
    }
    wypisz(wartosci, klucz);
    } catch (exception &e) {
      if (argc == 1) {
        cout << "sin x - oblicza sinus od x\n add x y - oblicza sume\n mod x y - oblicza x%y" << endl;
      } else {
        cout << "Niepoprawne argumenty" << endl;
        cout << "Blad: " << e.what() << endl;
      }
      return 1;
    }
    return 0;
}
