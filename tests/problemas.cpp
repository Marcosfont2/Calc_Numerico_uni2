#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "interpolacao.hpp"

using namespace std;

vector<Ponto> gerarPontosExponencial(int grau) {
    vector<Ponto> pontos;
    int quantidade = grau + 1;
    double inicio = 1.0;
    double fim = 3.0;
    double intervalo = (fim - inicio) / grau;
    for (int i = 0; i < quantidade; i++) {
        double x = inicio + i * intervalo;
        double y = exp(x);
        pontos.push_back({x, y});
    }
    return pontos;
}

void executarProblema2_2() {
    int graus[] = {3, 4, 5, 10};
    double xTeste = 2.25;
    double valorReal = exp(xTeste);
    cout << fixed << setprecision(10) << endl;
    cout << "Problema 2.2" << endl;
    cout << "Valor de x: " << endl;
    cout << "Valor real de e^x: " << valorReal << endl;
    for (int grau : graus) {
        vector <Ponto> pontos = gerarPontosExponencial(grau);
        double resultadoGauss = interpolacaoGauss(pontos, xTeste);
        double resultadoLagrange = lagrange(pontos, xTeste);
        double resultadoNewton = newton(pontos, xTeste);
        cout << "\nGrau do polinomio: " << grau << endl;
        cout << "Sistema linear: " << resultadoGauss << " | Erro: " << fabs(valorReal - resultadoGauss) << endl;
        cout << "Lagrange: " << resultadoLagrange << " | Erro: " << fabs(valorReal - resultadoLagrange) << endl;
        cout << "Newton: " << resultadoNewton << " | Erro: " << fabs(valorReal - resultadoNewton) << endl;
    }
}

int main() {
    try {
        executarProblema2_2();
    } catch (const exception& erro) {
        cout << "Erro: " << erro.what() << endl;
    }
    return 0;
}