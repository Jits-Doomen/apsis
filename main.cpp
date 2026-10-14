#include "apsis/interpreter.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Apsis — scientific and astrodynamics language\n"
                  << "Usage: apsis <script.jdx>\n";
        return 2;
    }
    return apsis::runFile(argv[1]);
}
