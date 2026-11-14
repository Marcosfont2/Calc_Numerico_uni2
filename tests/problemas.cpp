#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "interpolacao.hpp"

using namespace std;

// gera os pontos de interpolacao da funcao exponencial
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

// gera os pontos de interpolacao da funcao racional
vector<Ponto> gerarPontosRacional(int grau) {
    vector<Ponto> pontos;
    int quantidade = grau + 1;
    double inicio = 0.0;
    double fim = 4.0;
    double intervalo = (fim - inicio) / grau;
    for (int i = 0; i < quantidade; i++) {
        double x = inicio + i * intervalo;
        double y = 1.0 / (1.0 + 25.0 * pow(x, 2));
        pontos.push_back({x, y});
    }
    return pontos;
}

// dados historicos da concentracao de C02 (grafico)
vector<Ponto> obterDadosCO2() {
    vector <Ponto> dados = {
        {1960, 315},
        {1970, 327},
        {1980, 337},
        {1990, 335},
        {2000, 370},
        {2010, 388},
        {2020, 411}
    };
    return dados;
}

// executa o problema 2.2 comparando os 3 metodos
void problema2_2() {
    int graus[] = {3, 4, 5, 10};
    double xTeste = 2.25;
    double valorReal = exp(xTeste);
    cout << fixed << setprecision(10) << endl;
    cout << "Problema 2.2" << endl;
    cout << "Valor de x: " << xTeste << endl;
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

// executa o probleam 2.3 com o metodo do sistema linear
void problema2_3() {
    int graus[] = {4, 10, 15};
    double xTeste = 2.5;
    double valorReal = 1.0 / (1.0 + 25.0 * pow(xTeste, 2));
    cout << "\nProblema 2.3" << endl;
    cout << "Valor de x: " << xTeste << endl;
    cout << "Valor real da funcao: " << valorReal << endl;
    cout << fixed << setprecision(10);
    for (int grau : graus) {
        vector<Ponto> pontos = gerarPontosRacional(grau);
        double resultado = interpolacaoGauss(pontos, xTeste);
        cout << "\nGrau do polinomio: " << grau << endl;
        cout << "Interpolacao por sistema linear: " << resultado << endl;
        cout << "Erro absoluto: " << fabs(valorReal - resultado) << endl;
    }
}

// executa o item 2.4.a, estimando a concetracao de C02 em 2005 com o polinomio interpolador de grau 3
void problema2_4a() {
    vector<Ponto> dados = obterDadosCO2();
    vector<Ponto> pontos = {
        dados[2], // 1980
        dados[3], // 1990
        dados[4], // 2000
        dados[5]  // 2010
    };
    double ano = 2005.0;
    double estimativa = lagrange(pontos, ano);
    cout << "\nProblema 2.4.a" << endl;
    cout << "Estimativa para 2005: " << estimativa << "ppm" << endl;
}

// executa o item 2.4.b, estimando a concetracao de CO2 em 2008 com polimonios de grau 2 e 3
void problema2_4b(){
    vector<Ponto> dados = obterDadosCO2();
    vector<Ponto> pontos2 = {
        dados[3], // 1990
        dados[4], // 2000
        dados[5]  // 2010
    };
    vector<Ponto> pontos3 = {
        dados[2], // 1980
        dados[3], // 1990
        dados[4], // 2000
        dados[5]  // 2010
    };
    double ano = 2008.0;
    double valorObservado = 381.0;
    double estimativa2 = lagrange(pontos2, ano);
    double estimativa3 = lagrange(pontos3, ano);
    cout << "\nProblema 2.4.b" << endl;
    cout << "Valor observado em 2008: " << valorObservado << "ppm" << endl;
    cout << "\nEstimativa com grau 2: " << estimativa2 << "ppm" << endl;
    cout << "Diferenca para o observado: " << fabs(valorObservado - estimativa2) << "ppm" << endl;
    cout << "\nEstimativa para o grau 3: " << estimativa3 << "ppm" << endl;
    cout << "Diferenca para o observado: " << fabs(valorObservado - estimativa3) << "ppm" << endl;
}

// executa o item 2.4.c, estimando o ano em que a concentracao de C02 atingiu 350 ppm com polinomio de grau 2
void problema2_4c() {
    vector<Ponto> dados = obterDadosCO2();
    // tres pontos proximos a 350
    vector<Ponto> pontos = {
        dados[2], // 1980
        dados[3], // 1990
        dados[4]  // 2000
    };
    double valorAlvo = 350.0;
    // intervalo em que a concetracao passou de 350
    double inicio = 1990.0;
    double fim = 2000.0;
    // metodo da bissecao
    for (int i = 0; i < 100; i++) {
        double meio = (inicio + fim) / 2.0;
        double valorInicio = lagrange(pontos, inicio) - valorAlvo;
        double valorMeio = lagrange(pontos, meio) - valorAlvo;
        if (valorInicio * valorMeio <= 0) {
            fim = meio;
        } else {
            inicio = meio;
        }
    }
    double anoEstimado = (inicio + fim) / 2.0;
    cout << "\nProblema 2.4.c" << endl;
    cout << "Concentrada desejada: " << valorAlvo << "ppm" << endl;
    cout << "Ano estimado: " << anoEstimado << endl;
}

int main() {
    try {
        problema2_2();
        problema2_3();
        problema2_4a();
        problema2_4b();
        problema2_4c();
    } catch (const exception& erro) {
        cout << "Erro: " << erro.what() << endl;
    }
    return 0;
}