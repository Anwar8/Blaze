#ifndef ELEMENT_TYPES
#define ELEMENT_TYPES
/**
 * @brief an enum that defines the types of beam-column elements
 * 
 */
enum ElementType {
    LinearElastic = 0,
    NonlinearElastic = 1,
    NonlinearPlastic = 2
};
#endif