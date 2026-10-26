#include "LinearSolver.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

// ============================ Конструктор и проверки ============================

LinearSolver::LinearSolver(const Matrix& a, const Vector& rhs)
    : A(a), b(rhs), n(static_cast<int>(a.size()))
{
    checkSquare();
}

void LinearSolver::checkSquare() const {
    if (A.size() != b.size()) {
        throw std::invalid_argument("Размеры матрицы A и вектора b не совпадают");
    }
    for (const auto& row : A) {
        if (row.size() != A.size()) {
            throw std::invalid_argument("Матрица A должна быть квадратной");
        }
    }
}

// ============================ Определитель ============================

double LinearSolver::det(Matrix M) {
    int m = static_cast<int>(M.size());
    if (m == 1) return M[0][0];
    if (m == 2) return M[0][0] * M[1][1] - M[0][1] * M[1][0];

    // Для m >= 3 — приведение к треугольному виду (Гаусс)
    double d = 1.0;
    for (int i = 0; i < m; ++i) {
        // Поиск опорного элемента
        int piv = i;
        for (int r = i + 1; r < m; ++r) {
            if (std::abs(M[r][i]) > std::abs(M[piv][i])) piv = r;
        }
        if (std::abs(M[piv][i]) < 1e-14) return 0.0;
        if (piv != i) {
            std::swap(M[i], M[piv]);
            d = -d;
        }
        d *= M[i][i];
        for (int r = i + 1; r < m; ++r) {
            double factor = M[r][i] / M[i][i];
            for (int c = i; c < m; ++c) {
                M[r][c] -= factor * M[i][c];
            }
        }
    }
    return d;
}

// ============================ Метод Крамера ============================

Vector LinearSolver::cramer() const {
    double Dz = det(A);
    if (std::abs(Dz) < 1e-14) {
        throw std::runtime_error("det(A) близок к нулю: метод Крамера неприменим");
    }

    Vector x(n, 0.0);
    for (int i = 0; i < n; ++i) {
        Matrix Ai = A;
        for (int r = 0; r < n; ++r) Ai[r][i] = b[r];
        x[i] = det(Ai) / Dz;
    }
    return x;
}

// ============================ Метод Гаусса ============================

Vector LinearSolver::gauss() const {
    Matrix M = A;
    Vector rhs = b;

    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента в k-м столбце
        int piv = k;
        for (int r = k + 1; r < n; ++r) {
            if (std::abs(M[r][k]) > std::abs(M[piv][k])) piv = r;
        }
        if (std::abs(M[piv][k]) < 1e-14) {
            throw std::runtime_error("Матрица вырождена (нулевой ведущий элемент)");
        }
        if (piv != k) {
            std::swap(M[k], M[piv]);
            std::swap(rhs[k], rhs[piv]);
        }
        // Исключение
        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j < n; ++j) M[i][j] -= factor * M[k][j];
            rhs[i] -= factor * rhs[k];
        }
    }

    // Обратный ход
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = rhs[i];
        for (int j = i + 1; j < n; ++j) sum -= M[i][j] * x[j];
        if (std::abs(M[i][i]) < 1e-14) {
            throw std::runtime_error("Нулевой диагональный элемент при обратном ходе");
        }
        x[i] = sum / M[i][i];
    }
    return x;
}

// ============================ Метод главных элементов ============================

Vector LinearSolver::gaussMainElement() const {
    Matrix M = A;
    Vector rhs = b;
    std::vector<int> colOrder(n);
    for (int i = 0; i < n; ++i) colOrder[i] = i;

    for (int k = 0; k < n; ++k) {
        // Поиск главного элемента в подматрице [k..n-1] x [k..n-1]
        int prow = k, pcol = k;
        double maxAbs = -1.0;
        for (int r = k; r < n; ++r) {
            for (int c = k; c < n; ++c) {
                if (std::abs(M[r][c]) > maxAbs) {
                    maxAbs = std::abs(M[r][c]);
                    prow = r;
                    pcol = c;
                }
            }
        }
        if (maxAbs < 1e-14) {
            throw std::runtime_error("Матрица вырождена (главный элемент = 0)");
        }
        // Перестановка строк
        if (prow != k) {
            std::swap(M[k], M[prow]);
            std::swap(rhs[k], rhs[prow]);
        }
        // Перестановка столбцов (с запоминанием порядка неизвестных)
        if (pcol != k) {
            for (int r = 0; r < n; ++r) std::swap(M[r][k], M[r][pcol]);
            std::swap(colOrder[k], colOrder[pcol]);
        }
        // Исключение
        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j < n; ++j) M[i][j] -= factor * M[k][j];
            rhs[i] -= factor * rhs[k];
        }
    }

    // Обратный ход
    Vector xtmp(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double sum = rhs[i];
        for (int j = i + 1; j < n; ++j) sum -= M[i][j] * xtmp[j];
        if (std::abs(M[i][i]) < 1e-14) {
            throw std::runtime_error("Нулевой диагональный элемент");
        }
        xtmp[i] = sum / M[i][i];
    }

    // Восстановление порядка неизвестных
    Vector x(n, 0.0);
    for (int i = 0; i < n; ++i) x[colOrder[i]] = xtmp[i];
    return x;
}

// ============================ Диагональное преобладание ============================

bool LinearSolver::hasDiagonalDominance() const {
    for (int i = 0; i < n; ++i) {
        double diag = std::abs(A[i][i]);
        double off = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) off += std::abs(A[i][j]);
        }
        if (diag <= off) return false;
    }
    return true;
}

// ============================ Метод простой итерации (Якоби) ============================

Vector LinearSolver::simpleIteration(double eps, int maxIter) const {
    // Построение C, D
    Matrix C(n, Vector(n, 0.0));
    Vector D(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (std::abs(A[i][i]) < 1e-14) {
            throw std::runtime_error("Нулевой диагональный элемент a_ii: переставьте уравнения");
        }
        D[i] = b[i] / A[i][i];
        for (int j = 0; j < n; ++j) {
            C[i][j] = (j == i) ? 0.0 : -A[i][j] / A[i][i];
        }
    }

    // Достаточное условие сходимости: ||C|| < 1
    double normC = 0.0;
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += std::abs(C[i][j]);
        normC = std::max(normC, s);
    }
    if (normC >= 1.0) {
        std::cerr << " [предупреждение] ||C|| = " << normC
                  << " >= 1, сходимость не гарантирована" << std::endl;
    }

    // Начальное приближение: x^(0) = D
    Vector x = D;
    for (int k = 0; k < maxIter; ++k) {
        Vector xNew(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double sum = D[i];
            for (int j = 0; j < n; ++j) sum += C[i][j] * x[j];
            xNew[i] = sum;
        }
        // Критерий (47)
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            if (std::abs(xNew[i] - x[i]) > eps) { converged = false; break; }
        }
        x = xNew;
        if (converged) return x;
    }
    throw std::runtime_error("Метод простой итерации не сошёлся за заданное число шагов");
}

// ============================ Метод Зейделя ============================

Vector LinearSolver::seidel(double eps, int maxIter) const {
    Matrix C(n, Vector(n, 0.0));
    Vector D(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (std::abs(A[i][i]) < 1e-14) {
            throw std::runtime_error("Нулевой диагональный элемент a_ii");
        }
        D[i] = b[i] / A[i][i];
        for (int j = 0; j < n; ++j) {
            C[i][j] = (i == j) ? 0.0 : -A[i][j] / A[i][i];
        }
    }

    Vector x = D; // начальное приближение
    for (int k = 0; k < maxIter; ++k) {
        Vector xOld = x;
        for (int i = 0; i < n; ++i) {
            double sum = D[i];
            for (int j = 0; j < n; ++j) sum += C[i][j] * x[j];
            // x[j] уже обновлён для j < i
            x[i] = sum;
        }
        // Критерий (47)
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            if (std::abs(x[i] - xOld[i]) > eps) { converged = false; break; }
        }
        if (converged) return x;
    }
    throw std::runtime_error("Метод Зейделя не сошёлся за заданное число шагов");
}

// ============================ Обратная матрица ============================

Matrix LinearSolver::inverse() const {
    Matrix aug(n, Vector(2 * n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) aug[i][j] = A[i][j];
        aug[i][n + i] = 1.0;
    }

    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int r = k + 1; r < n; ++r) {
            if (std::abs(aug[r][k]) > std::abs(aug[piv][k])) piv = r;
        }
        if (std::abs(aug[piv][k]) < 1e-14) {
            throw std::runtime_error("Матрица вырождена: обратной не существует");
        }
        if (piv != k) std::swap(aug[k], aug[piv]);

        double pivVal = aug[k][k];
        for (int j = 0; j < 2 * n; ++j) aug[k][j] /= pivVal;

        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            double factor = aug[i][k];
            for (int j = 0; j < 2 * n; ++j) aug[i][j] -= factor * aug[k][j];
        }
    }

    Matrix inv(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) inv[i][j] = aug[i][n + j];
    }
    return inv;
}

// ============================ Нормы и число обусловленности ============================

double LinearSolver::vecNorm(const Vector& v) {
    double s = 0.0;
    for (double x : v) s += std::abs(x);
    return s;
}

double LinearSolver::matNorm(const Matrix& M) {
    int m = static_cast<int>(M.size());
    if (m == 0) return 0.0;
    double best = 0.0;
    for (int j = 0; j < m; ++j) {
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += std::abs(M[i][j]);
        best = std::max(best, s);
    }
    return best;
}

double LinearSolver::conditionNumber() const {
    Matrix inv = inverse();
    return matNorm(A) * matNorm(inv);
}

// ============================ Уточнение корней ============================

Vector LinearSolver::refine(const Vector& x0, int iterations) const {
    Vector x = x0;
    for (int it = 0; it < iterations; ++it) {
        // Невязка ε = b - A*x
        Vector eps(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double ax = 0.0;
            for (int j = 0; j < n; ++j) ax += A[i][j] * x[j];
            eps[i] = b[i] - ax;
        }
        // Решаем A*delta = eps методом Гаусса
        LinearSolver ls(A, eps);
        Vector delta = ls.gauss();
        for (int i = 0; i < n; ++i) x[i] += delta[i];
    }
    return x;
}

// ============================ Вывод решения ============================

void LinearSolver::printSolution(const std::string& title, const Vector& x) const {
    std::cout << "=== " << title << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < n; ++i) {
        std::cout << " x" << (i + 1) << " = " << x[i] << std::endl;
    }

    // Проверка: A*x vs b
    std::cout << " Проверка (A*x - b):" << std::endl;
    double maxErr = 0.0;
    for (int i = 0; i < n; ++i) {
        double ax = 0.0;
        for (int j = 0; j < n; ++j) ax += A[i][j] * x[j];
        double err = ax - b[i];
        maxErr = std::max(maxErr, std::abs(err));
        std::cout << " уравнение " << (i + 1) << ": " << err << std::endl;
    }
    std::cout << " max|A*x - b| = " << maxErr << std::endl;
    std::cout << std::endl;
}