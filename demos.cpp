#include "demos.h"
#include "LinearSolver.h"

#include <iostream>
#include <iomanip>

// ============================ Пример 11 ============================

void demoExample11() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 11: Метод Гаусса" << std::endl;
    std::cout << "==============================================" << std::endl;

    Matrix A = { {2, 2, 4}, {2, -1, 3}, {3, -1, 2} };
    Vector b = {18, 9, 7};

    LinearSolver ls(A, b);
    Vector x = ls.gauss();
    ls.printSolution("Решение методом Гаусса (ожидается 1, 2, 3)", x);
}

// ============================ Пример 12 ============================

void demoExample12() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 12: метод главных элементов" << std::endl;
    std::cout << "==============================================" << std::endl;

    Matrix A = { { 0.34,  0.71,  0.63},
                 { 0.71, -0.65, -0.18},
                 { 1.17, -2.35,  0.75} };
    Vector b = {2.08, 0.17, 1.28};

    LinearSolver ls(A, b);
    Vector x = ls.gaussMainElement();
    ls.printSolution("Решение методом главных элементов (ожидается ~1.331, 0.693, 1.802)", x);
}

// ============================ Пример 13 ============================

void demoExample13() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 13: система с диагональным преобладанием" << std::endl;
    std::cout << "==============================================" << std::endl;

    // Преобразованная система из лекции
    Matrix A = { {10,  2, -1,  2},
                 { 1,  5,  1,  0},
                 { 1, -2, -5,  1},
                 { 3,  0,  0, -9} };
    Vector b = {-4, 1, 2, 10};

    LinearSolver ls(A, b);
    std::cout << "Диагональное преобладание: "
              << (ls.hasDiagonalDominance() ? "ДА" : "НЕТ") << std::endl;

    Vector x = ls.simpleIteration(1e-8, 5000);
    ls.printSolution("Решение методом простой итерации (Пример 13)", x);
}

// ============================ Пример 14 ============================

void demoExample14() {
    std::cout << "==============================================" << std::endl;
    std::cout << " Пример 14: метод простой итерации (Якоби)" << std::endl;
    std::cout << "==============================================" << std::endl;

    Matrix A = { {4, 0.24, -0.08},
                 {0.09, 3, -0.15},
                 {0.04, -0.08, 4} };
    Vector b = {8, 9, 20};

    LinearSolver ls(A, b);
    std::cout << "Диагональное преобладание: "
              << (ls.hasDiagonalDominance() ? "ДА" : "НЕТ") << std::endl;

    Vector x = ls.simpleIteration(1e-5, 1000);
    ls.printSolution("Решение методом простой итерации (ожидается ~1.9092, 3.19496, 5.04481)", x);
}

// ============================ Пример 15 ============================

void demoExample15() {
    std::cout << "========================" << std::endl;
    std::cout << " Пример 15: метод Зейделя" << std::endl;
    std::cout << "========================" << std::endl;

    Matrix A = { {9, -2, 1},
                 {2, -7, 1},
                 {1,  3, 8} };
    Vector b = {8, -4, 12};

    LinearSolver ls(A, b);
    std::cout << "Диагональное преобладание: "
              << (ls.hasDiagonalDominance() ? "ДА" : "НЕТ") << std::endl;

    Vector x = ls.seidel(1e-8, 1000);
    ls.printSolution("Решение методом Зейделя (ожидается 1, 1, 1)", x);
}

// ============================ Крамер и число обусловленности ============================

void demoCramerAndCond() {
    std::cout << "========================" << std::endl;
    std::cout << " Метод Крамера и число обусловленности" << std::endl;
    std::cout << "========================" << std::endl;

    Matrix A = { {2, 2, 4},
                 {2, -1, 3},
                 {3, -1, 2} };
    Vector b = {18, 9, 7};

    LinearSolver ls(A, b);
    Vector x = ls.cramer();
    ls.printSolution("Метод Крамера", x);

    std::cout << "Число обусловленности cond(A) = "
              << ls.conditionNumber() << std::endl;
    std::cout << std::endl;

    // Уточнение корней
    Vector xApprox = {1.001, 1.999, 3.002}; // искусственно возмущённое
    Vector xRefined = ls.refine(xApprox, 2);
    ls.printSolution("Уточнение корней по невязке (2 итерации)", xRefined);
}