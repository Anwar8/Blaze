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

const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();

/**
 * @brief populate Tpetra::Vector<> by using the less efficient replaceGlobalValue method and iterating over a std::vector of Eigen::triplet<real> 
 * 
 * @tparam scalar_type the type of the Tpetra scalar
 * @param V Tpetra::Vector 
 * @param triplets a std::vector of Eigen::triplet<real> objects that were originally meant for populating an Eigen sparse matrix
 */
void set_from_triplets(Tpetra::Vector<real>& V, std::vector<spnz>& triplets)
{
    for (spnz& triplet : triplets)
    {
        V.replaceGlobalValue(triplet.row(), triplet.value());
    }
}

/**
 * @brief populate Tpetra::Vector<> by using the more efficient replaceLocalValue method and iterating over a std::vector of Eigen::triplet<real> 
 * 
 * @param V Tpetra::Vector 
 * @param triplets a std::vector of Eigen::triplet<real> objects that were originally meant for populating an Eigen sparse matrix
 * @param nz_i is the starting value for numbering the DoFs on a rank; used to convert from global to local numbering.
 * @warning only works because the mapping from global to local is a simple offset due to the renumbering of the nodes and DoFs on each rank.
 */
void set_from_triplets(Tpetra::Vector<real>& V, std::vector<spnz>& triplets, int const nz_i)
{
    for (spnz& triplet : triplets)
    {
        V.replaceLocalValue(triplet.row() - nz_i, triplet.value());
    }
}

/**
 * @brief Get a local memory (HostSpace) 1D view from a Tpetra::Vector
 * 
 * @param V 
 * @return Kokkos::View<const scalar_type*, Kokkos::HostSpace> 
 */
Kokkos::View<const real*, Kokkos::HostSpace> get_1d_view(Tpetra::Vector<real>& V)
{
    auto V_2d = V.getLocalViewHost(Tpetra::Access::ReadOnly);
    return Kokkos::subview (V_2d, Kokkos::ALL (), 0);
}

#endif