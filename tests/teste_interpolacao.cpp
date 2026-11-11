#include <iostream>
#include "interpolacao.hpp"

using namespace std;

int main() {
  try {
    vector<Ponto> pontos = {
    {1.0, 2.0},
    {2.0, 4.0},
    {3.0, 6.0}
    };
    double x = 2.5;

    double resultadoLagrange = lagrange(pontos, x);
    double resultadoNewton = newton(pontos, x);

    cout << "Valor de x: " << x << endl;
    cout << "\nMetodo de Lagrange: " << resultadoLagrange << endl;
    cout << "Metodo de Newton: " << resultadoNewton << endl;
  } catch (const exception& erro) {
    cout << "Erro: " << erro.what() << endl;
  }
  return 0;
}
