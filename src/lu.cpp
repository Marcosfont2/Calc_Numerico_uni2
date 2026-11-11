#include "sistemas_lineares.hpp"

#include <cmath>
#include <stdexcept>

namespace linear_systems {

// Entrada: matriz quadrada A. Saída: A=L*U, com diagonal unitária em L.
// Os multiplicadores de Gauss são guardados em L; não há pivoteamento.
LUFactors factorLU(const Matrix& a) {
    const int n = static_cast<int>(a.size());
    validateSystem(a, Vector(n, 0.0));
    Matrix lower(n, Vector(n, 0.0));
    Matrix upper = a;
    for (int i = 0; i < n; ++i) lower[i][i] = 1.0;

    for (int k = 0; k < n; ++k) {
        if (upper[k][k] == 0.0) {
            throw std::runtime_error("Pivô nulo em LU, etapa " + std::to_string(k + 1)
                                     + "; a versão sem pivoteamento não pode prosseguir.");
        }
        if (!std::isfinite(upper[k][k])) throw std::runtime_error("Pivô não finito em LU.");
        for (int i = k + 1; i < n; ++i) {
            if (upper[i][k] == 0.0) continue;
            lower[i][k] = upper[i][k] / upper[k][k];
            if (!std::isfinite(lower[i][k])) throw std::runtime_error("Multiplicador não finito em LU.");
            upper[i][k] = 0.0;
            for (int j = k + 1; j < n; ++j) upper[i][j] -= lower[i][k] * upper[k][j];
        }
    }
    validateSystem(upper, Vector(n, 0.0));
    return {lower, upper};
}

// Entrada: fatores já calculados e b. Saída: solução obtida primeiro de Ly=b
// e depois de Ux=y. L e U não são modificados e podem ser reutilizados.
Vector solveLU(const LUFactors& factors, const Vector& b) {
    const Vector y = forwardSubstitution(factors.lower, b);
    return backSubstitution(factors.upper, y);
}

} // namespace linear_systems
