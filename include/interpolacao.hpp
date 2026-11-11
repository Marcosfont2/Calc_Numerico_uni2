#ifndef INTERPOLACAO_HPP
#define INTERPOLACAO_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

struct Ponto {
double x;
double y;
};

// Verifica se existem pontos e se os valores de x sao distintos
void validarPontos(const vector<Ponto>& pontos) {
  if (pontos.empty()){
    throw invalid_argument("A lista de pontos nao pode estar vazia.");
  }
  for (int i = 0; i < pontos.size(); i++) {
    for (int j = i + 1; j < pontos.size(); j++) {
      if (pontos[i].x == pontos[j].x) {
        throw invalid_argument("Existem pontos com valores de x repetidos");
      }
    }
  }
}

// Método de Lagrange
double lagrange(const vector<Ponto>& pontos, double x) {
  validarPontos(pontos);
  double resultado = 0.0;
  int n = pontos.size();
  for (int i = 0; i < n; i++) {
    // Calcula cada termo da formula de Lagrange
    double termo  = pontos[i].y;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        termo *= (x - pontos[j].x) / (pontos[i].x - pontos[j].x);
      }
    }
    // Soma o termo calculado ao resultado final
    resultado += termo;
  }
  return resultado;
}

// Método de Newton
double newton(const vector<Ponto>& pontos, double x) {
  validarPontos(pontos);
  int n = pontos.size();
  // Armazena os valores de y, que sao transformados nos coeficientes
  vector<double> coeficientes(n);
  for (int i = 0; i < n; i++) {
    coeficientes[i] = pontos[i].y;
  }
  // Calcula os coeficientes por meio das diferenas divididas
  for (int i = 1; i < n; i++) {
    for (int j = n - 1; j >= i; j--) {
      coeficientes[j] = (coeficientes[j] - coeficientes[j - 1]) / (pontos[j].x - pontos[j - i].x);
    }
  }
  // Avalia o polinomio utilizando a forma de Newton
  double resultado = coeficientes[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    resultado = coeficientes[i] + (x - pontos[i].x) * resultado;
  }
  return resultado;
}

#endif
