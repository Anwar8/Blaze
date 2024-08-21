#ifndef MAIN_HPP
#define MAIN_HPP


#define VERBOSE 0
#define VERBOSE_STIFFNESSES 0
#define VERBOSE_NLB 0
#define LF_VERBOSE 0
// 1. original basic_beam, 2. Izzuddin's beam, 3. the new Linear Beam 2D, 4 the new Nonlinear Beam 2D, or 5 the plasatic version of the nonlinear beam-column.
#define ELEM 5


#define PI 3.14159265358979323846

enum ElementType {
    LinearElastic = 0,
    NonlinearElastic = 1,
    NonlinearPlastic = 2
};
#endif