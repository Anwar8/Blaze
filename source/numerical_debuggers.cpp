#include "numerical_debuggers.hpp"

// TODO: Fix BROKEN analyser
bool has_zero_row(spmat A) {
    int n = A.outerSize();
    // int* nnz = A.innerNonZeroPtr();
    // for(int i = 0; i < n; ++i)
    // {
    //     if(nnz[i] == 0)
    //     std::cout << "Row " << i << " is zero\n";
    //     return true;
    // }
    return false;
}

bool check_matrix(spmat A) {
    return has_zero_row(A);
}
