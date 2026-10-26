#include "demos.h"

#include <iostream>
#include <exception>

int main() {
    try {
        std::cout << "==============================================" << std::endl;
        std::cout << " РЕШЕНИЕ СИСТЕМ ЛИНЕЙНЫХ УРАВНЕНИЙ (Лекция 6)" << std::endl;
        std::cout << "==============================================" << std::endl << std::endl;

        demoExample11();
        demoExample12();
        demoExample13();
        demoExample14();
        demoExample15();
        demoCramerAndCond();

        std::cout << "==============================================" << std::endl;
        std::cout << " Все демонстрационные примеры выполнены" << std::endl;
        std::cout << "==============================================" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}