#ifndef SECTION_BASE_CLASS
#define SECTION_BASE_CLASS

#include "maths_defaults.hpp"
/**
 * @file SectionBaseClass.hpp
 *
 * This class provides a foundation for creating sections of different types (basic and otherwise).
 */
class SectionBaseClass {
public:
    SectionBaseClass() = default;
    virtual ~SectionBaseClass() = default;
    virtual real get_E() {return 0.0;}
    virtual real get_A() {return 0.0;}
    virtual real get_I() {return 0.0;}
    virtual void update_section_state(vec eps) {};
    virtual mat get_D_t() {return make_xd_mat(2,2);}
};
#endif