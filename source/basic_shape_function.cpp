#include "basic_shape_function.hpp"

void BasicShapeFunction::calc_N(real x, real L) {
    N(0,0) = 1 - (x/L);
    N(0,3) = x / L;
    N(1,1) = 1 - 3*std::pow(x/L,2) + 2*std::pow(x/L,3);
    N(1,2) = x - 2*std::pow(x,2)/L + std::pow(x/L, 2)*x;
    N(1,4) = 3*std::pow(x/L, 2) - 2*std::pow(x/L, 3);
    N(1,5) = -x*(x/L) + x * std::pow(x/L,2);
}

void BasicShapeFunction::calc_B(real x, real L) {
    B(0,0) = -1/L;
    B(0,3) = 1/L;
    B(1,1) = -6*std::pow(1/L,2) + 12*x*std::pow(1/L,3);
    B(1,2) = - 4/L + 6*x*std::pow(1/L, 2);
    B(1,4) = 6*std::pow(1/L, 2) - 12*x*std::pow(1/L, 3);
    B(1,5) = -2/L + 6 * x* std::pow(1/L,2);
}

void BasicShapeFunction::calc_k(real L, BasicSection& sec) {
    real A = sec.get_A();
    real E = sec.get_E();
    real I = sec.get_I();
    // Row 1
    k(0,0) = E*A/L;
    k(0,3) = -E*A/L;
    // Row 2
    k(1,1) = 12*E*I/std::pow(L,3);
    k(1,2) = 6*E*I/std::pow(L,2);
    k(1,4) = -12*E*I/std::pow(L,3);
    k(1,5) = 6*E*I/std::pow(L,2);
    // Row 3
    k(2,1) = 6*E*I/std::pow(L,2);
    k(2,2) = 4*E*I/L;
    k(2,4) = -6*E*I/std::pow(L,2);
    k(2,5) = 2*E*I/L;
    // Row 4
    k(3,0) = -E*A/L;
    k(3,3) = E*A/L;
    // Row 5
    k(4,1) = -12*E*I/std::pow(L,3);
    k(4,2) = -6*E*I/std::pow(L,2);
    k(4,4) = 12*E*I/std::pow(L,3);
    k(4,5) = -6*E*I/std::pow(L,2);
    // Row 6
    k(5,1) = 6*E*I/std::pow(L,2);
    k(5,2) = 2*E*I/L;
    k(5,4) = -6*E*I/std::pow(L,2);
    k(5,5) = 4*E*I/L;
}
