// A simple example

#include <iostream>
#include "tinyexpr.h"

int main()
{
    int error;
    float f;

    std::string expr;
    std::cout << "Type arithmetic expression: ";
    std::cin >> expr;

    std::cout << "Let me calculate " << expr << " ... " << std::endl;

    f = te_interp(expr.c_str(), &error);

    if (error)
        std::cout << "Something wrong! Error code: " << error << std::endl;
    else
        std::cout << expr << " = " << f << std::endl;
}
