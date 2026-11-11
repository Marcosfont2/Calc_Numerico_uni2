#ifndef SISTEMAS_LINEARES_HPP
#define SISTEMAS_LINEARES_HPP

#include <string>
#include <vector>

namespace linear_systems {

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct GaussResult {
    Vector x;
    double determinant;
};

struct LUFactors {
    Matrix lower;
    Matrix upper;
};

enum class JacobiStatus { Converged, IterationLimit, NumericalFailure };

struct IterationInfo {
    int iteration;
    double absoluteDistance;
    double relativeDistance;
};

struct JacobiResult {
    Vector x;
    JacobiStatus status;
    int iterations;
    double absoluteDistance;
    double relativeDistance;
    bool rowCriterion;
    bool columnCriterion;
    std::vector<IterationInfo> history;
};

// Entrada: caminho do TXT e dimensão positiva esperada.
// Saída: vetor finito com n componentes; a primeira linha do TXT deve conter n.
// Lança exceção se o cabeçalho, os valores ou a quantidade forem inválidos.
Vector loadVector(const std::string& path, int expectedSize);

// Entrada: caminho e vetor finito. Saída: TXT com um valor por linha, sem
// cabeçalho, para soluções e resíduos. Lança exceção em caso de falha na escrita.
void saveVector(const std::string& path, const Vector& values);

// Entrada: A e b. Saída: nenhuma; verifica matriz não vazia, quadrada,
// dimensões compatíveis e valores finitos. Lança exceção se houver erro.
void validateSystem(const Matrix& a, const Vector& b);

// Entrada: matriz quadrada A e vetor x compatível. Saída: produto A*x.
// Dimensões inválidas ou resultados não finitos causam exceção.
Vector multiply(const Matrix& a, const Vector& x);

// Entrada: sistema original e aproximação x. Saída: r = b - A*x.
// Lança exceção para dados incompatíveis ou resultados não finitos.
Vector calculateResidual(const Matrix& a, const Vector& b, const Vector& x);

// Entrada: vetor finito. Saída: máximo dos módulos; para vetor vazio, zero.
double infinityNorm(const Vector& values);

// Entrada: sistema original, aproximação e seu resíduo.
// Saída: ||r||inf / (||A||inf * ||x||inf + ||b||inf).
// Denominador zero: retorna zero se o numerador for zero, infinito caso contrário.
double normalizedResidual(const Matrix& a, const Vector& b,
                          const Vector& x, const Vector& residual);

// Entrada: vetores finitos de mesmo tamanho. Saída: max |current[i]-previous[i]|.
double absoluteDistance(const Vector& current, const Vector& previous);

// Entrada: distância absoluta e aproximação atual. Saída: distância absoluta
// dividida por max |current[i]|. Trata denominador zero sem divisão por zero.
double relativeDistance(double absolute, const Vector& current);

// Entrada: matriz quadrada finita. Saída: true se todas as linhas (ou colunas)
// satisfizerem |a_ii| > soma dos módulos dos elementos fora da diagonal.
bool hasRowDominance(const Matrix& a);
bool hasColumnDominance(const Matrix& a);

// Entrada: matriz triangular inferior L e vetor b. Saída: y tal que L*y=b.
// Não modifica as entradas. Diagonal nula ou operação não finita causa exceção.
Vector forwardSubstitution(const Matrix& lower, const Vector& b);

// Entrada: matriz triangular superior U e vetor b. Saída: x tal que U*x=b.
// Não modifica as entradas. Diagonal nula ou operação não finita causa exceção.
Vector backSubstitution(const Matrix& upper, const Vector& b);

// Entrada: sistema quadrado A*x=b. Saída: solução por Gauss sem pivoteamento
// e determinante. Preserva A e b; pivô nulo ou falha numérica causa exceção.
// O determinante pode ser não representável mesmo quando a solução é finita.
GaussResult gaussianElimination(const Matrix& a, const Vector& b);

// Entrada: matriz quadrada A. Saída: L de diagonal unitária e U tais que A=L*U.
// Não modifica A, não troca linhas e lança exceção se não conseguir fatorar.
LUFactors factorLU(const Matrix& a);

// Entrada: fatores L e U e termo independente b. Saída: solução de L*U*x=b.
// Os fatores são preservados e podem ser reutilizados para outros vetores b.
Vector solveLU(const LUFactors& factors, const Vector& b);

// Entrada: sistema A*x=b, tolerância positiva e limite positivo de iterações.
// initial vazio escolhe x0[i]=b[i]/A[i][i]; caso contrário, usa o vetor informado.
// Saída: aproximação, estado, distâncias, critérios e histórico. Só converge
// quando AMBAS as distâncias forem menores que tolerance. A e b são preservados.
// Diagonal nula ou entradas inválidas causam exceção; falha durante a iteração
// devolve NumericalFailure e mantém a última aproximação finita.
JacobiResult gaussJacobi(const Matrix& a, const Vector& b,
                        double tolerance = 1e-6, int maxIterations = 10000,
                        const Vector& initial = {});

} // namespace linear_systems

#endif
