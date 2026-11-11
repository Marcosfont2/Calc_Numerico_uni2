#include "sistemas_lineares.hpp"

#include <cmath>
#include <stdexcept>

namespace linear_systems {

// Entrada: sistema A*x=b. Saída: solução e determinante por Gauss sem trocas
// de linhas. Cópias de A e b preservam os originais para o cálculo do resíduo.
GaussResult gaussianElimination(const Matrix& a, const Vector& b) {
    validateSystem(a, b);
    Matrix upper = a;
    Vector rhs = b;
    const int n = static_cast<int>(a.size());

    for (int k = 0; k < n; ++k) {
        if (upper[k][k] == 0.0) {
            throw std::runtime_error("Pivô nulo em Gauss, etapa " + std::to_string(k + 1)
                                     + "; a versão sem pivoteamento não pode prosseguir.");
        }
        if (!std::isfinite(upper[k][k])) throw std::runtime_error("Pivô não finito em Gauss.");
        for (int i = k + 1; i < n; ++i) {
            // Um elemento já nulo dispensa uma operação de eliminação.
            if (upper[i][k] == 0.0) continue;
            const double multiplier = upper[i][k] / upper[k][k];
            if (!std::isfinite(multiplier)) throw std::runtime_error("Multiplicador não finito em Gauss.");
            upper[i][k] = 0.0;
            for (int j = k + 1; j < n; ++j) upper[i][j] -= multiplier * upper[k][j];
            rhs[i] -= multiplier * rhs[k];
            if (!std::isfinite(rhs[i])) throw std::runtime_error("Termo independente não finito em Gauss.");
        }
    }

    Vector x = backSubstitution(upper, rhs);
    double determinant = 1.0;
    for (int i = 0; i < n; ++i) determinant *= upper[i][i];
    return {x, determinant};
}

} // namespace linear_systems
