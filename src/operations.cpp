#include "sistemas_lineares.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace linear_systems {

// Entrada: A e b. Saída: nenhuma; garante as pré-condições dos métodos.
void validateSystem(const Matrix& a, const Vector& b) {
    if (a.empty() || b.size() != a.size()) {
        throw std::invalid_argument("O sistema deve ser não vazio e ter dimensões compatíveis.");
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != a.size()) {
            throw std::invalid_argument("A matriz deve ser quadrada e ter linhas de mesmo tamanho.");
        }
        if (!std::isfinite(b[i])) throw std::invalid_argument("O vetor contém valor não finito.");
        for (double value : a[i]) {
            if (!std::isfinite(value)) throw std::invalid_argument("A matriz contém valor não finito.");
        }
    }
}

// Entrada: A e x compatíveis. Saída: produto matriz-vetor calculado por linhas.
Vector multiply(const Matrix& a, const Vector& x) {
    validateSystem(a, x);
    Vector product(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a.size(); ++j) product[i] += a[i][j] * x[j];
        if (!std::isfinite(product[i])) {
            throw std::runtime_error("Produto matriz-vetor não representável em double.");
        }
    }
    return product;
}

// Entrada: A, b e x originais. Saída: vetor resíduo, com convenção r=b-A*x.
Vector calculateResidual(const Matrix& a, const Vector& b, const Vector& x) {
    validateSystem(a, b);
    Vector residual = multiply(a, x);
    for (std::size_t i = 0; i < b.size(); ++i) {
        residual[i] = b[i] - residual[i];
        if (!std::isfinite(residual[i])) throw std::runtime_error("Resíduo não representável.");
    }
    return residual;
}

// Entrada: vetor finito. Saída: norma do máximo, usando o módulo de cada valor.
double infinityNorm(const Vector& values) {
    double norm = 0.0;
    for (double value : values) {
        if (!std::isfinite(value)) throw std::invalid_argument("Norma recebeu valor não finito.");
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

// Entrada: sistema, aproximação e resíduo. Saída: resíduo normalizado.
// O denominador usa long double para evitar overflow intermediário desnecessário.
double normalizedResidual(const Matrix& a, const Vector& b,
                          const Vector& x, const Vector& residual) {
    validateSystem(a, b);
    if (x.size() != a.size() || residual.size() != a.size()) {
        throw std::invalid_argument("Dimensões inválidas no resíduo normalizado.");
    }
    long double matrixNorm = 0.0L;
    for (const Vector& row : a) {
        long double sum = 0.0L;
        for (double value : row) sum += std::abs(static_cast<long double>(value));
        matrixNorm = std::max(matrixNorm, sum);
    }
    const long double denominator = matrixNorm * infinityNorm(x) + infinityNorm(b);
    const double numerator = infinityNorm(residual);
    if (denominator == 0.0L) {
        return numerator == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return static_cast<double>(numerator / denominator);
}

// Entrada: aproximações consecutivas. Saída: máximo dos módulos das diferenças.
double absoluteDistance(const Vector& current, const Vector& previous) {
    if (current.size() != previous.size()) throw std::invalid_argument("Vetores incompatíveis.");
    double distance = 0.0;
    for (std::size_t i = 0; i < current.size(); ++i) {
        if (!std::isfinite(current[i]) || !std::isfinite(previous[i])) {
            throw std::invalid_argument("Distância recebeu valor não finito.");
        }
        distance = std::max(distance, std::abs(current[i] - previous[i]));
    }
    return distance;
}

// Entrada: distância absoluta e vetor atual. Saída: distância relativa,
// com tratamento explícito dos casos 0/0 e distância positiva dividida por zero.
double relativeDistance(double absolute, const Vector& current) {
    const double denominator = infinityNorm(current);
    if (denominator == 0.0) {
        return absolute == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return absolute / denominator;
}

// Entrada: matriz quadrada. Saída: verificação estrita por linhas, com módulos.
bool hasRowDominance(const Matrix& a) {
    validateSystem(a, Vector(a.size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < a.size(); ++j) {
            if (i != j) sum += std::abs(a[i][j]);
        }
        if (!(std::abs(a[i][i]) > sum)) return false;
    }
    return true;
}

// Entrada: matriz quadrada. Saída: verificação estrita por colunas, com módulos.
bool hasColumnDominance(const Matrix& a) {
    validateSystem(a, Vector(a.size(), 0.0));
    for (std::size_t j = 0; j < a.size(); ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (i != j) sum += std::abs(a[i][j]);
        }
        if (!(std::abs(a[j][j]) > sum)) return false;
    }
    return true;
}

// Entrada: L triangular inferior e b. Saída: solução de L*y=b por substituição.
Vector forwardSubstitution(const Matrix& lower, const Vector& b) {
    validateSystem(lower, b);
    Vector y(b.size(), 0.0);
    for (std::size_t i = 0; i < b.size(); ++i) {
        if (lower[i][i] == 0.0) throw std::runtime_error("Diagonal nula na substituição direta.");
        double sum = 0.0;
        for (std::size_t j = 0; j < i; ++j) sum += lower[i][j] * y[j];
        y[i] = (b[i] - sum) / lower[i][i];
        if (!std::isfinite(y[i])) throw std::runtime_error("Falha numérica na substituição direta.");
    }
    return y;
}

// Entrada: U triangular superior e b. Saída: solução por retrosubstituição.
Vector backSubstitution(const Matrix& upper, const Vector& b) {
    validateSystem(upper, b);
    Vector x(b.size(), 0.0);
    for (int i = static_cast<int>(b.size()) - 1; i >= 0; --i) {
        if (upper[i][i] == 0.0) throw std::runtime_error("Diagonal nula na retrosubstituição.");
        double sum = 0.0;
        for (std::size_t j = i + 1; j < b.size(); ++j) sum += upper[i][j] * x[j];
        x[i] = (b[i] - sum) / upper[i][i];
        if (!std::isfinite(x[i])) throw std::runtime_error("Falha numérica na retrosubstituição.");
    }
    return x;
}

} // namespace linear_systems
