#ifndef TPETRA_WRAPPERS_HPP
#define  TPETRA_WRAPPERS_HPP

#include <Teuchos_DefaultMpiComm.hpp> 
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>

typedef Tpetra::Vector<>::scalar_type tpetra_scalar;
typedef Tpetra::Vector<>::global_ordinal_type tpetra_global_ordinal;
typedef Tpetra::Vector<>::local_ordinal_type tpetra_local_ordinal;

Teuchos::RCP<Tpetra::Vector<>> make_rcp_to_vector(Teuchos::Map<> vector_map)
{
    return (new rcp)
}

#endif