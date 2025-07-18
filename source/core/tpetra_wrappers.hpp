#ifndef TPETRA_WRAPPERS_HPP
#define  TPETRA_WRAPPERS_HPP

#include <Teuchos_DefaultMpiComm.hpp> 
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_FancyOStream.hpp>

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

std::map<tpetra_global_ordinal, std::pair<std::vector<tpetra_global_ordinal>, std::vector<real>>> map_triplets_to_rows(std::vector<spnz>& triplets)
{
    std::map<tpetra_global_ordinal, std::pair<std::vector<tpetra_global_ordinal>, std::vector<real>>> row_keyed_value_map;
    for (spnz& triplet : triplets)
    {
        row_keyed_value_map[tpetra_global_ordinal(triplet.row())].first.push_back(tpetra_global_ordinal(triplet.col())); 
        row_keyed_value_map[tpetra_global_ordinal(triplet.row())].second.push_back(triplet.value()); 
    }
    return row_keyed_value_map;
}

std::map<tpetra_global_ordinal, std::vector<tpetra_global_ordinal>> map_triplets_to_row_column_positions(std::vector<spnz>& triplets)
{
    std::map<tpetra_global_ordinal, std::vector<tpetra_global_ordinal>> row_keyed_position_map;
    for (spnz& triplet : triplets)
    {
        row_keyed_position_map[tpetra_global_ordinal(triplet.row())].push_back(tpetra_global_ordinal(triplet.col())); 
    }
    return row_keyed_position_map;
}
/**
 * @brief Replaces the existing values in a Tpetra::CrsMatrix based on triplets. Expects the matrix has already been initialised beforehand so will call fillResume.
 * @note I needed help with the reference and const conversions so I used GitHub Copilot running on GPT 4.1 on 19 July 2025 to help me with the conversions.
 * @param A a Tpetra::CrsMatrix that will be updated with new values from triplets.
 * @param triplets a std::vector of triplets that will be used to update the matrix A.
 */
void set_from_triplets(Teuchos::RCP<Tpetra::CrsMatrix<real>> A, std::vector<spnz>& triplets)
{
    A->resumeFill();
    auto row_map = map_triplets_to_rows(triplets);
    for (const auto& row_entry : row_map)
    {
        tpetra_global_ordinal row = row_entry.first;
        const std::vector<tpetra_global_ordinal>& cols = row_entry.second.first;
        const std::vector<real>& vals = row_entry.second.second;

        // Convert std::vector to Teuchos::ArrayView<const T> 
        Teuchos::ArrayView<const tpetra_global_ordinal> col_view(cols);
        Teuchos::ArrayView<const real> val_view(vals);

        // Replace existing values in CrsMatrix
        A->replaceGlobalValues(row, col_view, val_view);
    }
    A->fillComplete();
}

/**
 * @brief initialises an existing graph with global indices coming from triplets calculated from the contribution of finite elements.
 * 
 * @param A_graph a Teuchos RCP to Tpetra::CrsGraph that will be used to initialise the stiffness matrix.
 * @param triplets a std::vector of triplets that will be used to update the graph indices - their values will be ignored we only care about their column indices.
 */
void initialise_from_triplets(Teuchos::RCP<Tpetra::CrsGraph<>> A_graph, std::vector<spnz>& triplets)
{
    auto row_col_map = map_triplets_to_row_column_positions(triplets);
    std::cout << "there are " << row_col_map.size() << " entires in glob_mesh.update_elements_states()." << std::endl;
    for (const auto& row_entry : row_col_map)
    {
        tpetra_global_ordinal row = row_entry.first;

        // Convert std::vector to Teuchos::ArrayView<const T> 
        Teuchos::ArrayView<const tpetra_global_ordinal> col_view(row_entry.second);

        // Replace existing values in CrsMatrix
        A_graph->insertGlobalIndices(row, col_view);
    }
}


#endif