#ifndef TPETRA_WRAPPERS_HPP
#define  TPETRA_WRAPPERS_HPP

#include <Teuchos_DefaultMpiComm.hpp> 
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

typedef Tpetra::Vector<>::scalar_type tpetra_scalar;
typedef Tpetra::Vector<>::global_ordinal_type tpetra_global_ordinal;
typedef Tpetra::Vector<>::local_ordinal_type tpetra_local_ordinal;

Teuchos::RCP<int> test_rcp(int& an_int)
{
    return Teuchos::rcpFromRef(an_int);
} 


#endif