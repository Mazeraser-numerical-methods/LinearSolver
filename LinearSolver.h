#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <vector>
#include <string>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

/**
 * @brief Класс для решения систем линейных уравнений A*x = b.
 *
 * Реализует:
 *  - метод Крамера;
 *  - метод Гаусса (с частичным выбором главного элемента по столбцу);
 *  - метод главных элементов (полный выбор);
 *  - метод простой итерации (Якоби);
 *  - метод Зейделя;
 *  - число обусловленности cond(A);
 *  - уточнение корней по невязке.
 */
class LinearSolver {
public:
    /**
     * @brief Конструктор.
     * @param a матрица коэффициентов A
     * @param rhs столбец свободных членов b
     * @throws std::invalid_argument при несовпадении размеров или не квадратной A
     */
    LinearSolver(const Matrix& a, const Vector& rhs);

    /**
     * @brief Метод Крамера (формулы (7)-(9)): x_i = Δ_i / Δ.
     * @throws std::runtime_error если det(A) близок к нулю.
     */
    Vector cramer() const;

    /**
     * @brief Метод Гаусса с частичным выбором главного элемента по столбцу.
     * Прямой ход по (17), обратный ход по (18)-(19).
     */
    Vector gauss() const;

    /**
     * @brief Метод главных элементов (полный выбор).
     * На каждом шаге ищется максимальный по модулю элемент во всей рабочей подматрице.
     */
    Vector gaussMainElement() const;

    /**
     * @brief Проверка диагонального преобладания (44):
     * |a_ii| > sum_{j != i} |a_ij|.
     */
    bool hasDiagonalDominance() const;

    /**
     * @brief Метод простой итерации (Якоби).
     * Приведение к виду (35): x = D + C*x; итерации (41)-(42); критерий (47).
     * @param eps требуемая абсолютная погрешность
     * @param maxIter максимальное число итераций
     */
    Vector simpleIteration(double eps = 1e-6, int maxIter = 1000) const;

    /**
     * @brief Метод Зейделя (формулы (48)).
     * При вычислении x_i^(k+1) используются уже найденные x_1^(k+1), ..., x_{i-1}^(k+1).
     */
    Vector seidel(double eps = 1e-6, int maxIter = 1000) const;

    /// Обратная матрица методом Гаусса–Жордана.
    Matrix inverse() const;

    /// Норма вектора (30): сумма модулей элементов.
    static double vecNorm(const Vector& v);

    /// Норма матрицы (31): max из сумм модулей по столбцам.
    static double matNorm(const Matrix& M);

    /// Число обусловленности (32): cond(A) = ||A|| * ||A^{-1}||.
    double conditionNumber() const;

    /**
     * @brief Уточнение корней по невязке (22)-(24): x_new = x_old + delta, A*delta = b - A*x_old.
     */
    Vector refine(const Vector& x0, int iterations = 1) const;

    /// Вывод решения и проверки подстановкой.
    void printSolution(const std::string& title, const Vector& x) const;

private:
    Matrix A;   // матрица коэффициентов
    Vector b;   // столбец свободных членов
    int n;      // размерность

    /// Проверка размеров: A квадратная, размеры A и b совпадают.
    void checkSquare() const;

    /// Определитель матрицы (для небольших n — разложением, для больших — приведением к треугольному виду).
    static double det(Matrix M);
};

#endif // LINEAR_SOLVER_H