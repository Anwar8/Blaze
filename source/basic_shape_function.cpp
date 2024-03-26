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
void BasicShapeFunction::calc_elem_mat_stiffness(real& L, BasicSection& sec, mat& k_mat)
{
    real A = sec.get_A();
    real E = sec.get_E();
    real I = sec.get_I();
    real EA = E*A;
    real EI = E*I;
    // Row 1
    k_mat(0,0) = EA/L;
    k_mat(0,3) = -EA/L;
    // Row 2
    k_mat(1,1) = 12*EI/std::pow(L,3);
    k_mat(1,2) = 6*EI/std::pow(L,2);
    k_mat(1,4) = -12*EI/std::pow(L,3);
    k_mat(1,5) = 6*EI/std::pow(L,2);
    // Row 3
    k_mat(2,1) = 6*EI/std::pow(L,2);
    k_mat(2,2) = 4*EI/L;
    k_mat(2,4) = -6*EI/std::pow(L,2);
    k_mat(2,5) = 2*EI/L;
    // Row 4
    k_mat(3,0) = -EA/L;
    k_mat(3,3) = EA/L;
    // Row 5
    k_mat(4,1) = -12*EI/std::pow(L,3);
    k_mat(4,2) = -6*EI/std::pow(L,2);
    k_mat(4,4) = 12*EI/std::pow(L,3);
    k_mat(4,5) = -6*EI/std::pow(L,2);
    // Row 6
    k_mat(5,1) = 6*EI/std::pow(L,2);
    k_mat(5,2) = 2*EI/L;
    k_mat(5,4) = -6*EI/std::pow(L,2);
    k_mat(5,5) = 4*EI/L;
}
void BasicShapeFunction::calc_elem_geom_stiffness(real& L, real& P, mat& k_g)
{
    real N = P/(30*L); /**< common factor of the geometric stiffness.*/
    k_g(1,1) = N*36;
    k_g(1,2) = N*3*L;
    k_g(1,4) = N*-36;
    k_g(1,5) = N*3*L;
    
    k_g(2,1) = N*3*L;
    k_g(2,2) = N*4*L*L;
    k_g(2,4) = N*-3*L;
    k_g(2,5) = N*-L*L;

    k_g(4,1) = N*-36;
    k_g(4,2) = N*-3*L;
    k_g(4,4) = N*36;
    k_g(4,5) = N*-3*L;

    k_g(5,1) = N*3*L;
    k_g(5,2) = N*-L*L;
    k_g(5,4) = N*-3*L;
    k_g(5,5) = N*L*L;
}
