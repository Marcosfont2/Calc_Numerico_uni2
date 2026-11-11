#ifndef MATRIZ_MARKET_HPP
#define MATRIZ_MARKET_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

class MatrixMarketReader {
public:
    // Lê uma matriz no formato .rsa e a converte para representação densa (2D).
    static std::vector<std::vector<double>> carregarDensO(const std::string& caminho) {
        std::ifstream arquivo(caminho);
        if (!arquivo.is_open()) {
            throw std::runtime_error("Erro ao abrir o arquivo .rsa (Harwell-Boeing): " + caminho);
        }

        std::string linha;
        
        // Linha de Título e Chaves (Linha 1)
        if (!std::getline(arquivo, linha)) throw std::runtime_error("Erro ao ler cabeçalho do arquivo .rsa");

        //  Linha de Totais de Linhas/Ponteiros/Valores (Linha 2)
        //  TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD
        if (!std::getline(arquivo, linha)) throw std::runtime_error("Erro na linha 2 do arquivo .rsa");
        int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
        {
            std::stringstream ss(linha);
            ss >> totcrd >> ptrcrd >> indcrd >> valcrd >> rhscrd;
        }

        //  Linha de Tipo e Dimensões (Linha 3)
        //  MXTYPE, NROW, NCOL, NNZERO, NELTVL
        if (!std::getline(arquivo, linha)) throw std::runtime_error("Erro na linha 3 do arquivo .rsa");
        std::string mxtype;
        int nrow, ncol, nnzero, neltvl;
        {
            std::stringstream ss(linha);
            ss >> mxtype >> nrow >> ncol >> nnzero >> neltvl;
        }

        // Validação básica do tipo (deve ser 'R'eal e 'S'ymmetric)
        if (mxtype.empty() || mxtype[0] != 'R') {
            std::cout << "[AVISO] O tipo da matriz no arquivo .rsa pode não ser puramente Real (Tipo: " << mxtype << ")\n";
        }

        //  Linha de Formatos Fortran (Linha 4) - Apenas ignorada para leitura em fluxo livre
        if (!std::getline(arquivo, linha)) throw std::runtime_error("Erro na linha 4 do arquivo .rsa");

        // Se houver dados de lado direito (RHS) no arquivo, a linha 5 existe
        if (rhscrd > 0) {
            std::getline(arquivo, linha); 
        }

        // Leitura dos Ponteiros (col_ptr) - Ocupa 'ptrcrd' linhas
        std::vector<int> col_ptr(ncol + 1);
        for (int i = 0; i <= ncol; ) {
            if (!std::getline(arquivo, linha)) break;
            std::stringstream ss(linha);
            int val;
            while (ss >> val && i <= ncol) {
                col_ptr[i++] = val;
            }
        }

        // Leitura dos Índices de Linhas (row_ind) - Ocupa 'indcrd' linhas
        std::vector<int> row_ind(nnzero);
        for (int k = 0; k < nnzero; ) {
            if (!std::getline(arquivo, linha)) break;
            std::stringstream ss(linha);
            int val;
            while (ss >> val && k < nnzero) {
                row_ind[k++] = val;
            }
        }

        // Leitura dos Valores Numéricos (values) - Ocupa 'valcrd' linhas
        std::vector<double> valores(nnzero);
        for (int k = 0; k < nnzero; ) {
            if (!std::getline(arquivo, linha)) break;
            // O formato Fortran às vezes usa 'D' para expoentes em vez de 'E'
            for (char& c : linha) {
                if (c == 'D' || c == 'd') c = 'E';
            }
            std::stringstream ss(linha);
            double val;
            while (ss >> val && k < nnzero) {
                valores[k++] = val;
            }
        }

        // Montagem da matriz densa (NROW x NCOL) a partir da estrutura compactada CSC (Compressed Sparse Column)
        std::vector<std::vector<double>> A(nrow, std::vector<double>(ncol, 0.0));

        int indice_valor = 0;
        for (int j = 0; j < ncol; ++j) {
            // No padrão Harwell-Boeing, os índices em col_ptr começam em 1
            int inicio_coluna = col_ptr[j] - 1;
            int fim_coluna = col_ptr[j + 1] - 1;

            for (int p = inicio_coluna; p < fim_coluna; ++p) {
                int i = row_ind[indice_valor] - 1; // Índice de linha baseado em 1 -> 0
                double val = valores[indice_valor];
                indice_valor++;

                A[i][j] = val;
                // Como a matriz é simétrica (S in .rsa), preenchemos a posição transposta também
                if (i != j) {
                    A[j][i] = val;
                }
            }
        }

        std::cout << " Matriz carregada: " << nrow << "x" << ncol << " com " << nnzero << " não-nulos.\n";
        return A;
    }
};

#endif