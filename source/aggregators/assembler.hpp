/**
 * @file assembler.hpp
 * @brief assembler object and its related functionality
 * 
 */
#ifndef ASSEMBLER
#define ASSEMBLER

#include "global_mesh.hpp"
#include "tpetra_wrappers.hpp"
#ifdef KOKKOS
    #include <Kokkos_Core.hpp>
#endif
/**
 * @brief helps assemble the global matrices.
 * 
 * @details supposed to help with organising the code by separating 
 * assembly; right now it is not helpful and only contains one big
 * function and redundunt member variables that appeared in \ref GlobalMesh.
 * 
 * @todo improve teh naming scheme of the variables. Right now they are either very
 * vague such as "K_glob" (now became K_glob_triplets) when it is actually referring to triplet vector, or way too
 * long such as "K_global_elem_triplet_contribution"
 * 
 * @todo break down the assembler into load and stiffness assembly steps in stead of combined together in a big mess.
 * 
 */
class Assembler {
    private:
        #ifndef WITH_MPI
        spmat K; /**< Global stiffness matrix \f$\boldsymbol{K}\f$.*/
        spmat P; /**< Nodal force vector \f$\boldsymbol{P}\f$.*/
        spmat R; /**< Resistance force vector \f$\boldsymbol{R}\f$.*/
        spmat G; /**< Out of balance force vector  \f$\boldsymbol{G}\f$.*/
        spvec U; /**< Global Nodal displacement vector \f$\boldsymbol{U}\f$ - also known as system state vector.*/
        spvec dU; /**< Change in global Nodal displacement vector \f$\Delta\boldsymbol{U}\f$ - also known as system state vector increment.*/
        #else
        Teuchos::RCP<Tpetra::Map<>> vector_map; /**<a map that explains which rows of the \f$\boldsymbol{P}\f$, \f$\boldsymbol{U}\f$, \f$\boldsymbol{R}\f$, etc. vectors go on which cores.*/
        Teuchos::RCP<Tpetra::CrsGraph<>> matrix_graph; /**<a graph that explains which rows of the \f$\boldsymbol{K}\f$ matrix go on which cores.*/
        Teuchos::RCP<Tpetra::Map<>> interface_map; /**< a map that is used for communicating the interface DoFs that are globally not locally allocated. */
        Teuchos::RCP<Tpetra::Import<>> interface_importer; /**< the importer object that is meant for communicating the \f$\boldsymbol{U}\f$ values on different nodes and ranks. */

        Teuchos::RCP<Tpetra::CrsMatrix<real>> K;
        Tpetra::MultiVector<real> P;
        Tpetra::MultiVector<real> R;
        Tpetra::MultiVector<real> G;
        Tpetra::MultiVector<real> U;
        Tpetra::MultiVector<real> dU;

        Tpetra::MultiVector<real> interface_U; /**< the interface version of the U vector.*/
            
        #endif

        realx2 G_max; /**< Max out-of-balance force retrieved by taking the square-root of the l2 norm of \f$ \sqrt{\norm{\boldsymbol{G}}}\f$.*/
        std::vector<spnz> K_global_triplets; /** The container for the triplets that are used for assembling the global stiffness matrix \f$ \boldsymbol{K}\f$ from element contributions.*/
        std::vector<spnz> R_global_triplets; /** The container for the triplets that are used for assembling the global resistance vector \f$ \boldsymbol{R}\f$ from element contributions.*/
        std::vector<spnz> P_global_triplets;  /** The container for the triplets that are used for assembling the global load matrix \f$ \boldsymbol{P}\f$ from nodal contributions.*/
    public:
        friend class BasicSolver;
        /**
         * @brief initialises the K sparse matrix to the size that correspond to the mesh being used.
         * 
         * @param glob_mesh  the global_mesh object which contains information about the number of nodes and degrees of freedom.
         */
        void initialise_stiffness_matrix(GlobalMesh& glob_mesh) {
            K_global_triplets.reserve(glob_mesh.nelems*36);
            #ifndef WITH_MPI
            K = make_spd_mat(glob_mesh.ndofs, glob_mesh.ndofs);
            #else
            setup_tpetra_crs_graph(glob_mesh);
            K = Teuchos::RCP(new Tpetra::CrsMatrix<real>(matrix_graph));
            #endif 
        }
        /**
         * @brief initialises the P, U, G, and dU vectors to the sizes that correspond to the mesh being used.
         * 
         * @param glob_mesh  the global_mesh object which contains information about the number of nodes and degrees of freedom.
         */
        void initialise_global_vectors(GlobalMesh& glob_mesh) {
            R_global_triplets.reserve(glob_mesh.ndofs);
            P_global_triplets.reserve(glob_mesh.ndofs);

            #ifndef WITH_MPI
            R = make_spd_mat(glob_mesh.ndofs, 1);
            G = make_spd_mat(glob_mesh.ndofs, 1);
            P = make_spd_mat(glob_mesh.ndofs, 1);
            U = make_spd_vec(glob_mesh.ndofs);
            dU = make_spd_vec(glob_mesh.ndofs);            
            #else
            setup_tpetra_vector_map(glob_mesh.get_ndofs(), glob_mesh.get_rank_ndofs());
            // The constructor for the vectors that follows prefills the vectors with zeros.
            R = Tpetra::MultiVector<real>(vector_map, 1);
            G = Tpetra::MultiVector<real>(vector_map, 1);
            P = Tpetra::MultiVector<real>(vector_map, 1);
            U = Tpetra::MultiVector<real>(vector_map, 1);
            dU = Tpetra::MultiVector<real>(vector_map, 1);
            setup_interface_import(glob_mesh);
            #endif 
        }
        /**
         * @brief Set the up tpetra vector map object which is used for initialising the distributed vectors.
         * 
         * @param ndofs the total number of active DoFs of the entire problem.
         * @param rank_ndofs the number of DoFs owned by the current rank. 
         */
        void setup_tpetra_vector_map(int ndofs, int rank_ndofs)
        {
            #ifdef WITH_MPI
            Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
            const size_t numLocalEntries = rank_ndofs;
            const Tpetra::global_size_t numGlobalEntries = ndofs;
            const tpetra_global_ordinal indexBase = 0;
            vector_map = Teuchos::rcp(new Tpetra::Map<>(numGlobalEntries, numLocalEntries, indexBase, comm));
            #endif
        }

        /**
         * @brief Setup tpetra crs graph object which is used to initialise the CrsMatrix used for stiffness. This is okay because the structure of the matrix will not change in the entire runtime so we should use the graph as it is much more efficient and allows for access using the local indices.
         * 
         * @param max_num_row_contributions the maximum number of entries per row of the stiffness matrix - retrieved from the GlobalMesh object.
         */
        void setup_tpetra_crs_graph(GlobalMesh& glob_mesh)
        {
            #ifdef WITH_MPI
            const size_t entriesPerRow = glob_mesh.max_num_stiffness_contributions;
            matrix_graph = Teuchos::rcp(new Tpetra::CrsGraph<>(vector_map, entriesPerRow));
            collect_global_K_triplets(glob_mesh);
            initialise_from_triplets(matrix_graph, K_global_triplets);
            matrix_graph->fillComplete();
            #endif
        }

        /**
         * @brief prints Tpetra::MultiVector or Tpetra::CrsMatrix or Tpetra::CrsGraph
         * 
         */
        void print_distributed_maths_object(std::string what, Teuchos::EVerbosityLevel verbosity = Teuchos::VERB_EXTREME)
        {
            #ifdef WITH_MPI
            Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
            if (what == "P")
            {
                P.describe(*out, verbosity);                
            }
            else if (what == "K")
            {
                K->describe(*out, verbosity);
            }
            else if (what == "U")
            {
                U.describe(*out, verbosity);                
            }
            else if (what == "dU")
            {
                dU.describe(*out, verbosity);           
            }
            else if (what == "R")
            {
                R.describe(*out, verbosity);                
            }
            else if (what == "G")
            {
                G.describe(*out, verbosity);                
            }
            else
            {
                std::cout << "Assembler:print_distributed_maths_object can only print U, dU, P, G, R, or K, but was asked to print " << what << std::endl;
                exit(1);
            }

            #endif
        }

        /**
         * @brief prints \ref stiffness_map to the output stream
         * 
         */
        void print_stiffness_graph()
        {
            #ifdef WITH_MPI
            Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
            matrix_graph->describe(*out, Teuchos::VERB_HIGH);
            #endif
        }

        /**
         * @brief retrieves global load contributions from all nodes.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the nodes and then uses the triplets to create the global load vector \f$ \boldsymbol{P}\f$
         * Note that as per the documentation of Eigen: "The initial contents of *this is destroyed." 
         * when calling [setFromTriplets](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0).
         * So this means that this function does not add to the existing force vector but recreates it from the contributions of the nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_P(GlobalMesh& glob_mesh)
        {
            P_global_triplets.clear();

            for (auto& node: glob_mesh.node_vector)
            {
                node->insert_load_triplets(this->P_global_triplets);
            }
            if (VERBOSE)
            {
                std::cout << "Assembler: all triplets are: " << std::endl;
                for (auto& triplet: P_global_triplets)
                {
                    std::cout << "row, col, val: " << triplet.row() << "," << triplet.col() << "," << triplet.value() << std::endl;
                }
                std::cout << "There are " << std::size(P_global_triplets) << " P_global contributions to add up." << std::endl;
            }
            #ifdef WITH_MPI
            set_from_triplets(P, P_global_triplets, glob_mesh.rank_starting_nz_i);
            if (VERBOSE)
            {
                print_distributed_maths_object("P");
            }
            #else
            P.setFromTriplets(P_global_triplets.begin(), P_global_triplets.end());
            P.makeCompressed();
            if (VERBOSE_NLB)
            {
                std::cout << "The P vector is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
            }
            #endif
        }

        void collect_global_K_triplets(GlobalMesh& glob_mesh)
        {
            // clear the members but keep the size reservation unchanged.
            K_global_triplets.clear();
            for (auto& elem: glob_mesh.elem_vector)
            {   
                elem->insert_global_stiffness_triplets(K_global_triplets);
            }

            // std::cout << "Assembler::collect_global_K_triplets: there are " << K_global_triplets.size() << " triplets that were collected from the elements." << std::endl;
        }

        /**
         * @brief retrieves global contributions to stiffness and resistance from all elements. For resistance, this function maps local element nodal forces \f$\boldsymbol{f}\f$ to the resistance vector \f$\boldsymbol{R}\f$ \ref R.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to create the global sparse stiffness
         * matrix, and the resistance vector. Note that as per the 
         * documentation of Eigen: "The initial contents of *this is destroyed." when calling [setFromTriplets](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0).
         * So this means that this function does not add to the existing matrix but recreates it from the contributions of the elements.
         * @todo add a function to calculate constrained nodes reactions.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_K_R(GlobalMesh& glob_mesh)
        {
            // clear the members but keep the size reservation unchanged.
            K_global_triplets.clear();
            R_global_triplets.clear();
            for (auto& elem: glob_mesh.elem_vector)
            {   
                elem->insert_global_stiffness_triplets(K_global_triplets);
                elem->insert_global_resistance_force_triplets(R_global_triplets);
            }

            #ifdef WITH_MPI
            set_from_triplets(R, R_global_triplets, glob_mesh.rank_starting_nz_i);
            set_from_triplets(K, K_global_triplets);
            if (VERBOSE)
            {
                print_distributed_maths_object("K");
            }
            if (VERBOSE_NLB)
            {
                std::cout << "The R vector is:" << std::endl;
                print_distributed_maths_object("R");
                std::cout << std::endl;
                std::cout << "and P vector is:" << std::endl;
                print_distributed_maths_object("P");
                std::cout << std::endl;
            }
            #else
            K.setFromTriplets(K_global_triplets.begin(), K_global_triplets.end());
            K.makeCompressed();
            R.setFromTriplets(R_global_triplets.begin(), R_global_triplets.end());
            R.makeCompressed();

            if (VERBOSE)
            {
                std::cout << "There are " << std::size(K_global_triplets) << " global_stiffness_triplets contributions to add up." << std::endl;
                std::cout << "The K_global_triplets is of size " << glob_mesh.ndofs << "x" << glob_mesh.ndofs << std::endl;
            }
            if (VERBOSE_NLB)
            {
                std::cout << "The R vector is:" << std::endl << Eigen::MatrixXd(R) << std::endl;
                std::cout << "KU, however, is:" << std::endl << Eigen::MatrixXd(K*U) << std::endl;
                std::cout << "and P is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
            }
            #endif
        }

        /**
         * @brief initialises \ref interface_map and \ref interface_U, and sets up \ref interface_importer for future communication.
         * 
         * @param glob_mesh the \ref GlobalMesh object.
         */
        void setup_interface_import(GlobalMesh& glob_mesh)
        {
            #ifdef WITH_MPI
            std::vector<unsigned> interface_dofs;
            // the presizing is too big as it does not account for inactive DoFs, but is better than trying to resize the vector many times over. The big size is not an issue because I expect only a few interface nodes anyway!
            interface_dofs.reserve(glob_mesh.interface_node_vector.size()*6);

            for (auto& node: glob_mesh.interface_node_vector)
            {
                int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
                int num_node_dofs = node->get_ndof(); // how many there are to loop over.
                std::set<int> node_active_dofs = node->get_active_dofs();
                int i = 0;
                for (auto& dof: node_active_dofs) {
                    interface_dofs.push_back(i + nzi);
                    ++i;
                }
            }
            Teuchos::Array<tpetra_global_ordinal> interface_dofs_array(interface_dofs.size());
            for (int k = 0; k < interface_dofs.size(); ++k) {
                interface_dofs_array[k] = interface_dofs[k];
            }
            // cyclicMap = rcp (new map_type (numGlobalEntries, elementList, indexBase, comm));
            interface_map = Teuchos::rcp(new Tpetra::Map<>(INVALID, interface_dofs_array, 0, vector_map->getComm()));
            interface_U = Tpetra::MultiVector<real>(interface_map, 1);
            interface_importer = Teuchos::rcp(new Tpetra::Import<>(vector_map, interface_map));
            // since it is initialisation, the Tpetra::CombineMode is INSERT.
            interface_U.doImport(U, *interface_importer, Tpetra::INSERT);
            #endif
        }

        /**
         * @brief maps state vector U back to nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void map_U_to_nodes(GlobalMesh& glob_mesh)
        {
            std::vector<std::shared_ptr<Node>>* nodes = &glob_mesh.node_vector;
            #ifdef WITH_MPI
            auto U_local_view = get_1d_view(U);
            int nzi = 0;
            for (auto& node: glob_mesh.node_vector)
                {
                    // int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = node->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = node->get_active_dofs();
                    int i = 0;
                    for (auto& dof: node_active_dofs) {
                        node->set_nodal_displacement(dof, U_local_view(i + nzi));
                        ++i;
                    }
                    nzi += i;
                }
                // since this is a call happening often, the Tpetra::CombineMode is REPLACE since the elements should already exist inplace.
                interface_U.doImport(U, *interface_importer, Tpetra::REPLACE);

                auto interface_U_local_view = get_1d_view(interface_U);
                nzi = 0;
                for (auto& node: glob_mesh.interface_node_vector)
                {
                    std::set<int> node_active_dofs = node->get_active_dofs();
                    int i = 0;
                    for (auto& dof: node_active_dofs) {
                        node->set_nodal_displacement(dof, interface_U_local_view(i + nzi));
                        ++i;
                    }
                    nzi += i;
                }
                if (VERBOSE_NLB)
                {
                    glob_mesh.read_nodal_U();
                }
            #else
            #ifdef KOKKOS
                Kokkos::parallel_for( "Assembler::map_U_to_nodes", glob_mesh.node_vector.size(), KOKKOS_LAMBDA (int i) {
                    int nzi = (*nodes)[i]->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = (*nodes)[i]->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = (*nodes)[i]->get_active_dofs();
                    int j = 0;
                    for (auto& dof: node_active_dofs) {
                        (*nodes)[i]->set_nodal_displacement(dof, U.coeff(j + nzi,0));
                        ++j;
                    }
                });
            #else
                for (auto& node: glob_mesh.node_vector)
                {
                    int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = node->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = node->get_active_dofs();
                    int i = 0;
                    
                    for (auto& dof: node_active_dofs) {
                        node->set_nodal_displacement(dof, U.coeff(i + nzi,0));
                        ++i;
                    }
                }
            #endif
                // There are only a few interface nodes and these are not worth the overhead of parallelising.
                for (auto& node: glob_mesh.interface_node_vector)
                {
                    int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = node->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = node->get_active_dofs();
                    int i = 0;
                    
                    for (auto& dof: node_active_dofs) {
                        node->set_nodal_displacement(dof, U.coeff(i + nzi,0));
                        ++i;
                    }
                }
            #endif
        }


        /**
         * @brief calculates out of balance forces from \f$\boldsymbol{G} =  \boldsymbol{R} - \boldsymbol{P}\f$.
         * 
         */
        void calculate_out_of_balance() {
            #ifdef WITH_MPI
            G.update(1.0, R, -1.0, P, 0.0);
            if (VERBOSE_NLB)
            {
                std::cout << "The G (out of balance) vector is:" << std::endl;
                print_distributed_maths_object("G");
                std::cout << std::endl;
            }
            #else
            G = R - P;
            if (VERBOSE_NLB)
            {
                std::cout << "The G (out of balance) vector is:" << std::endl << Eigen::MatrixXd(G) << std::endl;
            }
            #endif
        }
        /**
         * @brief updates \f$\boldsymbol{U}\f$ by incrementing with \f$ \Delta \boldsymbol{U}\f$.
         * 
         */
        void increment_U() {
            #ifdef WITH_MPI
            if (VERBOSE_NLB)
            {
                std::cout << "U before update is " << std::endl;
                print_distributed_maths_object("U", Teuchos::VERB_EXTREME);
                std::cout << std::endl;
            }
            U.update(-1.0, dU, 1.0);
            if (VERBOSE_NLB)
            {
                std::cout << "U after update is " << std::endl;
                print_distributed_maths_object("U", Teuchos::VERB_EXTREME);
                std::cout << std::endl;
            }
            #else
            if (VERBOSE_NLB)
                std::cout << "U before update is " << std::endl << U << std::endl;
            U += dU;
            if (VERBOSE_NLB)
                std::cout << "U after update is " << std::endl << U << std::endl;
            #endif
        }

        /**
         * @brief checks if the maximum square-root of the norm of out-of-balance is smaller than a tolerance.
         * 
         * @param tolerance maximum allowed out of balance force.
         * @return true if converged (i.e. tolerance more than max out-of-balance)
         * @return false if not converged (i.e. tolerance more than out-of-balance)
         */
        bool check_convergence(real tolerance)
        {
            #ifdef WITH_MPI
            // had some trouble with the types for the norms, so I used copilot for help with typing here.
            Teuchos::Array<typename Tpetra::MultiVector<>::mag_type> norms(1);
            G.norm2(norms);
            G_max = norms[0];
            #else
            G_max = std::sqrt(calc_l2_norm(G));
            #endif
            return G_max < tolerance;
        }

        realx2 get_G_max() {return G_max;}
};
#endif