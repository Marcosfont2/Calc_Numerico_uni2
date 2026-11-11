#include "sistemas_lineares.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace linear_systems {

// Entrada: caminho e valores finitos. Saída: TXT com um valor por linha.
void saveVector(const std::string& path, const Vector& values) {
    infinityNorm(values); // Valida antes de abrir e eventualmente substituir o arquivo.
    std::ofstream file(path);
    if (!file) throw std::runtime_error("Não foi possível criar o vetor: " + path);
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (double value : values) file << value << '\n';
    file.close();
    if (!file) throw std::runtime_error("Falha ao salvar o vetor: " + path);
}

// Entrada: TXT com dimensão na primeira linha e tamanho esperado.
// Saída: vetor finito com exatamente essa quantidade de componentes.
// Cada palavra deve representar um número completo; texto extra é recusado.
Vector loadVector(const std::string& path, int expectedSize) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Não foi possível abrir o vetor: " + path);
    if (expectedSize <= 0) throw std::invalid_argument("Tamanho esperado deve ser positivo.");
    std::string header, extraHeader;
    if (!std::getline(file, header)) throw std::runtime_error("Falta a dimensão do vetor: " + path);
    std::istringstream dimensionLine(header);
    int dimension;
    if (!(dimensionLine >> dimension) || dimensionLine >> extraHeader || dimension != expectedSize) {
        throw std::runtime_error("Dimensão inválida no cabeçalho do vetor: " + path);
    }
    Vector values;
    std::string token, extra;
    while (file >> token) {
        std::istringstream number(token);
        double value;
        if (!(number >> value) || number >> extra || !std::isfinite(value)) {
            throw std::runtime_error("Valor inválido no vetor: " + path);
        }
        values.push_back(value);
        if (values.size() > static_cast<std::size_t>(expectedSize)) {
            throw std::runtime_error("O vetor tem mais componentes que a matriz: " + path);
        }
    }
    if (file.bad() || values.size() != static_cast<std::size_t>(expectedSize)) {
        throw std::runtime_error("Não foi possível ler o vetor com o tamanho esperado: " + path);
    }
    return values;
}

} // namespace linear_systems
