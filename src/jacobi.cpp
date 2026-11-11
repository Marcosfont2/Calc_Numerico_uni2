#include "sistemas_lineares.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace linear_systems {

// Entrada: sistema, tolerância, limite e aproximação inicial opcional.
// Saída: aproximação final e diagnóstico. A atualização usa dois vetores para
// garantir que nenhuma componente recém-calculada seja usada na mesma iteração.
JacobiResult gaussJacobi(const Matrix& a, const Vector& b, double tolerance,
                        int maxIterations, const Vector& initial) {
    validateSystem(a, b);
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || maxIterations <= 0) {
        throw std::invalid_argument("Tolerância e limite de iterações devem ser positivos e finitos.");
    }
    const int n = static_cast<int>(a.size());
    Vector previous(n, 0.0);
    if (!initial.empty()) {
        if (initial.size() != b.size()) throw std::invalid_argument("Aproximação inicial incompatível.");
        infinityNorm(initial); // Também verifica se as componentes são finitas.
        previous = initial;
    }
    for (int i = 0; i < n; ++i) {
        if (a[i][i] == 0.0) throw std::invalid_argument("Diagonal nula: Jacobi não pode dividir por a_ii.");
        if (initial.empty()) previous[i] = b[i] / a[i][i];
        if (!std::isfinite(previous[i])) throw std::runtime_error("Aproximação inicial não representável.");
    }

    const double infinity = std::numeric_limits<double>::infinity();
    JacobiResult result{previous, JacobiStatus::IterationLimit, 0, infinity, infinity,
                        hasRowDominance(a), hasColumnDominance(a), {}};
    Vector current(n, 0.0);

    for (int step = 0; step < maxIterations; ++step) {
        const int k = step + 1;
        result.iterations = k;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) sum += a[i][j] * previous[j];
            }
            current[i] = (b[i] - sum) / a[i][i];
            if (!std::isfinite(current[i])) {
                result.status = JacobiStatus::NumericalFailure;
                return result; // result.x ainda contém a última aproximação finita.
            }
        }

        const double absolute = absoluteDistance(current, previous);
        if (!std::isfinite(absolute)) {
            result.status = JacobiStatus::NumericalFailure;
            return result;
        }
        const double relative = relativeDistance(absolute, current);
        result.x = current;
        result.absoluteDistance = absolute;
        result.relativeDistance = relative;
        result.history.push_back({k, absolute, relative});

        // O enunciado exige as distâncias absoluta E relativa: usa-se &&.
        if (absolute < tolerance && relative < tolerance) {
            result.status = JacobiStatus::Converged;
            return result;
        }
        previous.swap(current);
    }
    return result;
}

} // namespace linear_systems
