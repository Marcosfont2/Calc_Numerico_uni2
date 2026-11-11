#include "sistemas_lineares.hpp"
#include "matriz_market.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using namespace linear_systems;
using Clock = std::chrono::steady_clock;

// Entrada: dimensão positiva e caminho do TXT de b. Saída: arquivo com a
// dimensão na primeira linha e n valores aleatórios em [-50, 50] nas seguintes.
// Preserva a geração do grupo e lança exceção se não conseguir salvar o vetor.
void gerarESalvarVetorB(int n, const std::string& caminho_saida) {
    if (n <= 0) throw std::invalid_argument("A dimensão do vetor b deve ser positiva.");
    std::ofstream arquivo(caminho_saida);
    if (!arquivo) throw std::runtime_error("Não foi possível criar o vetor b: " + caminho_saida);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-50.0, 50.0);
    arquivo << std::setprecision(std::numeric_limits<double>::max_digits10);
    arquivo << n << '\n';
    for (int i = 0; i < n; ++i) arquivo << dis(gen) << '\n';
    arquivo.close();
    if (!arquivo) throw std::runtime_error("Falha ao salvar o vetor b: " + caminho_saida);
    std::cout << "Vetor b de dimensão " << n << " gerado e salvo em " << caminho_saida << '\n';
}

// Entrada: instante inicial. Saída: segundos transcorridos até esta chamada.
double elapsedSeconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// Entrada: caminho de saída e diagnóstico do Jacobi. Saída: histórico CSV;
// iterações que falharam numericamente não possuem uma nova aproximação válida.
void saveHistory(const std::string& path, const JacobiResult& result) {
    std::ofstream file(path);
    file << "iteration,absolute_distance,relative_distance\n" << std::setprecision(17);
    for (const auto& step : result.history) {
        file << step.iteration << ',' << step.absoluteDistance << ',' << step.relativeDistance << '\n';
    }
    file.close();
    if (!file) throw std::runtime_error("Não foi possível salvar o histórico: " + path);
}

// Entrada: sistema original, resultado e pasta de saída. Saída: arquivos da
// solução e do resíduo; escreve normas no resumo e no terminal. Resíduo não
// representável é registrado como indisponível, sem substituir por zeros.
void saveSolution(const Matrix& a, const Vector& b, const Vector& x,
                  const std::string& directory, std::ostream& summary) {
    saveVector(directory + "/solution.txt", x);
    Vector residual;
    double norm, normalized;
    try {
        residual = calculateResidual(a, b, x);
        norm = infinityNorm(residual);
        normalized = normalizedResidual(a, b, x, residual);
    } catch (const std::runtime_error& error) {
        summary << "Resíduo indisponível: " << error.what() << '\n';
        std::cout << "  Resíduo indisponível: " << error.what() << '\n';
        return;
    }
    // Escrita fora do bloco acima: falha no disco não é um erro numérico.
    saveVector(directory + "/residual.txt", residual);
    summary << "Norma do resíduo (máximo): " << norm << '\n'
            << "Resíduo normalizado: " << normalized << '\n';
    std::cout << "  ||r||inf = " << norm << "; normalizado = " << normalized << '\n';
}

// Entrada: estado do Jacobi. Saída: descrição em português para os resultados.
std::string describeStatus(JacobiStatus status) {
    if (status == JacobiStatus::Converged) return "Ambos os critérios de parada atendidos";
    if (status == JacobiStatus::IterationLimit) return "Limite de iterações atingido";
    return "Falha numérica; preservada a última aproximação finita";
}

// Entrada: matriz identificada, sistema e parâmetros. Saída: execução sequencial
// dos três métodos e arquivos separados. Uma falha de um método não impede os
// demais; retorna false se ocorrer falha de entrada, escrita ou método direto.
bool runMethods(const std::string& name, const Matrix& a, const Vector& b,
                double tolerance, int maxIterations) {
    bool success = true;
    for (const std::string method : {"gauss", "lu", "jacobi"}) {
        const std::string directory = "results/" + name + "/" + method;
        std::filesystem::create_directories(directory);
        // Remove saídas anteriores para não apresentar resultados antigos em
        // uma execução que falhe ou que não consiga calcular seu resíduo.
        for (const std::string file : {"solution.txt", "residual.txt", "history.csv"}) {
            std::filesystem::remove(directory + "/" + file);
        }
        std::ofstream summary(directory + "/summary.txt");
        if (!summary) throw std::runtime_error("Não foi possível criar o resumo em " + directory);
        summary << std::setprecision(17) << "Matriz: " << name << "\nMétodo: " << method
                << "\nOrdem: " << a.size() << "\nVetor b: data/vetor_b_" << name << ".txt\n";
        std::cout << "\n" << name << " / " << method << std::endl;
        const auto start = Clock::now();
        try {
            Vector x;
            double seconds;
            if (method == "gauss") {
                const auto result = gaussianElimination(a, b);
                seconds = elapsedSeconds(start);
                x = result.x;
                summary << "Estado: solução calculada sem pivoteamento\n";
                if (std::isfinite(result.determinant) && result.determinant != 0.0) {
                    summary << "Determinante: " << result.determinant << '\n';
                } else {
                    summary << "Determinante: produto não representável em double\n";
                }
            } else if (method == "lu") {
                const auto factors = factorLU(a);
                x = solveLU(factors, b);
                seconds = elapsedSeconds(start);
                summary << "Estado: solução calculada por A=LU, sem pivoteamento\n";
            } else {
                const auto result = gaussJacobi(a, b, tolerance, maxIterations);
                seconds = elapsedSeconds(start);
                x = result.x;
                const auto status = describeStatus(result.status);
                summary << "Estado: " << status << "\nTolerância: " << tolerance
                        << "\nLimite de iterações: " << maxIterations
                        << "\nIterações tentadas: " << result.iterations
                        << "\nDistância absoluta (última iteração válida): " << result.absoluteDistance
                        << "\nDistância relativa (última iteração válida): " << result.relativeDistance
                        << "\nCritério estrito das linhas: " << (result.rowCriterion ? "atendido" : "não atendido")
                        << "\nCritério estrito das colunas: " << (result.columnCriterion ? "atendido" : "não atendido")
                        << "\nAproximação inicial: b_i/a_ii\n"
                        << "As distâncias não garantem erro da solução inferior à tolerância.\n";
                std::cout << "  " << status << "; iterações: " << result.iterations << '\n';
                if (!result.rowCriterion && !result.columnCriterion) {
                    std::cout << "  Sem garantia de convergência pelos critérios das linhas e colunas.\n";
                }
                saveHistory(directory + "/history.csv", result);
            }
            summary << "Tempo do método (segundos): " << seconds << '\n';
            std::cout << "  Tempo: " << seconds << " s\n";
            saveSolution(a, b, x, directory, summary);
        } catch (const std::exception& error) {
            success = false;
            summary << "Falha: " << error.what() << '\n';
            std::cerr << "  Falha: " << error.what() << '\n';
        }
        summary.close();
        if (!summary) throw std::runtime_error("Falha ao escrever resumo em " + directory);
    }
    return success;
}

// Entrada: argumentos opcionais --tolerance e --max-iterations.
// Saída: processa todas as matrizes; retorna 0 para execução concluída e 1 em
// falhas de entrada/escrita ou dos métodos diretos. Não convergência de Jacobi
// é um resultado numérico registrado, sem interromper as demais execuções.
int main(int argc, char* argv[]) {
    double tolerance = 1e-6;
    int maxIterations = 10000;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help") {
                std::cout << "Uso: ./build/trabalho [--tolerance 1e-6] "
                             "[--max-iterations 10000]\n"
                             "Execute a partir da raiz do projeto.\n";
                return 0;
            }
            if (i + 1 == argc) throw std::invalid_argument("Falta o valor de " + option);
            const std::string value = argv[++i];
            std::size_t parsed;
            if (option == "--tolerance") {
                tolerance = std::stod(value, &parsed);
                if (!std::isfinite(tolerance) || tolerance <= 0) throw std::invalid_argument("Tolerância inválida.");
            } else if (option == "--max-iterations") {
                maxIterations = std::stoi(value, &parsed);
                if (maxIterations <= 0) throw std::invalid_argument("Limite de iterações inválido.");
            } else {
                throw std::invalid_argument("Opção desconhecida: " + option);
            }
            if (parsed != value.size()) throw std::invalid_argument("Valor inválido para " + option);
        }

        std::cout << std::scientific << std::setprecision(6);
        bool success = true;
        for (const std::string name : {"bcsstk22", "bcsstk23", "bcsstk24"}) {
            try {
                const Matrix a = MatrixMarketReader::carregarDensO("data/" + name + ".rsa");
                const std::string path = "data/vetor_b_" + name + ".txt";
                const int n = static_cast<int>(a.size());
                if (!std::filesystem::exists(path)) gerarESalvarVetorB(n, path);
                const Vector b = loadVector(path, n);
                if (!runMethods(name, a, b, tolerance, maxIterations)) success = false;
            } catch (const std::exception& error) {
                success = false;
                std::cerr << name << ": " << error.what() << '\n';
            }
        }
        std::cout << "\nArquivos de saída disponíveis em results/.\n";
        return success ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "Erro: " << error.what() << '\n';
        return 1;
    }
}
