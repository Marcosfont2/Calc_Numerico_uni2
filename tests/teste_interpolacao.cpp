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
    double resultadoGauss = interpolacaoGauss(pontos, x);

    cout << "Teste de Interpolacao" << endl;
    cout << "Valor de x: " << x << endl;
    cout << "\nMetodo de Lagrange: " << resultadoLagrange << endl;
    cout << "Metodo de Newton: " << resultadoNewton << endl;
    cout << "Metodo por sistema linear (Gauss): " << resultadoGauss << endl;
  } catch (const exception& erro) {
    cout << "Erro: " << erro.what() << endl;
  }
  cout << "\nTeste de validacao dos pontos" << endl;
  vector<Ponto> pontosInvalidos = {
    {1.0, 2.0},
    {1.0, 4.0}  
  };
  try {
    validarPontos(pontosInvalidos);
    cout << "Erro: pontos repetidos foram aceitos" << endl;
  } catch(const invalid_argument& erro) {
    cout << "Validacao funcionando: " << erro.what() << endl;
  };
  
  return 0;
}
