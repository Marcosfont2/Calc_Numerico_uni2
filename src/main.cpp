#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include "matriz_market.hpp"
#include "sistemas_lineares.hpp"

// Função para gerar e salvar o vetor b aleatório entre -50 e 50
void gerarESalvarVetorB(int n, const std::string& caminho_saida) {
    std::ofstream arquivo(caminho_saida);
    if (!arquivo.is_open()) {
        std::cerr << "Não foi possível criar o arquivo do vetor b em: " << caminho_saida << std::endl;
        return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-50.0, 50.0);

    // Grava a dimensão na primeira linha e os elementos nas seguintes
    arquivo << n << "\n";
    for (int i = 0; i < n; ++i) {
        double val = dis(gen);
        arquivo << val << "\n";
    }

    std::cout << "O Vetor b de dimensão " << n << " gerado e salvo em " << caminho_saida << std::endl;
}

int main() {

    try {
        // Caminho para a matriz Harwell-Boeing (.rsa) com o exemplo bcsstk22
        std::string caminho_matriz = "data/bcsstk22.rsa";
        
        std::cout << "[INFO] Carregando matriz do arquivo: " << caminho_matriz << "..." << std::endl;
        std::vector<std::vector<double>> A = MatrixMarketReader::carregarDensO(caminho_matriz);

        int n = A.size();
        std::cout << "[INFO] Dimensão da matriz A obtida: " << n << "x" << n << std::endl;

        // 2. Geração do vetor b aleatório correspondente à dimensão da matriz (Intervalo: -50 a 50)
        std::string caminho_b = "data/vetor_b_bcsstk22.txt";
        gerarESalvarVetorB(n, caminho_b);

        // O próximo passo aqui é carregar esse vetor b do arquivo e chamar seus métodos lineares:
        // Ex: std::vector<double> b = carregarVetorB(caminho_b);
        // Ex: std::vector<double> x = SistemasLineares::eliminacaoGauss(A, b);

    } catch (const std::exception& e) {
        std::cerr << "[EXCEÇÃO] " << e.what() << std::endl;
    }

    return 0;
}