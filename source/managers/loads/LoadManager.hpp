/**
 * @file LoadManager.hpp
 * @brief defines the LoadManager class which contains information about loads, and can apply them to the \ref GlobalMesh.
 */

#ifndef LOAD_MANAGER_HPP
#define LOAD_MANAGER_HPP
#include <vector>
#include <map>
#include "global_mesh.hpp"
#include "NodalLoad.hpp"
#include "node.hpp"
#include "maths_defaults.hpp"
/**
 * @brief 
 * @details this class should not only act as a manager, but also as a factory. 
 * It creates \ref NodalLoad objects, assigns them a nodes from the GlobalMesh, and then is used by the \ref SolutionProcedure to increment the load using these \ref NodalLoad objects.
 */
class LoadManager
{   
    protected:
        std::vector<NodalLoad> nodal_loads; /**< a vector of NodalLoad objects that are used to store information about the load conditions for a set of loaded nodes.*/
    public:
        /**
         * @brief creates a \ref NodalLoad object and adds it to the nodal loads controlled by this \ref LoadManager. Uses node shared_ptrs to assign nodes in stead of using ids.
         * @tparam NodePtrContainer STL container that can be dereferenced with the [] operator - used for the nodes to be loaded.
         * @tparam DofContainer STL container that can be dereferenced with the [] operator - used for the DoFs to be loaded. Can also be a std::set.
         * @tparam LoadContainer STL container that can be dereferenced with the [] operator - used for the loads corresponding to the DoFs.
         * @param loaded_nodes the nodes to be loaded.
         * @param loaded_dofs the DoFs to be loaded.
         * @param loads the loads corresponding to the DoFs.
         */
        template <typename NodePtrContainer, typename DofContainer, typename LoadContainer>
        void create_a_nodal_load_by_ptr(NodePtrContainer loaded_nodes, DofContainer loaded_dofs, LoadContainer loads)
        {
            NodalLoad nodal_load;
            nodal_load.assign_nodes_by_ptr(loaded_nodes);
            nodal_load.assign_dofs_loads(loaded_dofs, loads);
            nodal_loads.push_back(nodal_load);
        }

        /**
         * @brief creates a \ref NodalLoad object and adds it to the nodal loads controlled by this \ref LoadManager. Uses node IDs to assign nodes in stead of using ptrs.
         * @tparam DofContainer STL container that can be dereferenced with the [] operator - used for the DoFs to be loaded. Can also be a std::set.
         * @tparam LoadContainer STL container that can be dereferenced with the [] operator - used for the loads corresponding to the DoFs.
         * @param loaded_nodes the node IDs to be loaded.
         * @param loaded_dofs the DoFs to be loaded.
         * @param loads the loads corresponding to the DoFs.
         */
        template <typename DofContainer, typename LoadContainer>
        void create_a_nodal_load_by_id(std::vector<size_t> loaded_node_ids, DofContainer loaded_dofs, LoadContainer loads, GlobalMesh& glob_mesh)
        {
            NodalLoad nodal_load;
            nodal_load.assign_nodes_by_id(loaded_node_ids, glob_mesh);
            nodal_load.assign_dofs_loads(loaded_dofs, loads);
            nodal_loads.push_back(nodal_load);
        }

        /**
         * @brief creates a \ref NodalLoad object and adds it to the nodal loads controlled by this \ref LoadManager. Uses node IDs to assign nodes in stead of using ptrs. Only assigns loads to nodes owned by the curret rank.
         * @tparam DofContainer STL container that can be dereferenced with the [] operator - used for the DoFs to be loaded. Can also be a std::set.
         * @tparam LoadContainer STL container that can be dereferenced with the [] operator - used for the loads corresponding to the DoFs.
         * @param loaded_nodes the node IDs to be loaded.
         * @param loaded_dofs the DoFs to be loaded.
         * @param loads the loads corresponding to the DoFs.
         */
        template <typename DofContainer, typename LoadContainer>
        void create_a_distributed_nodal_load_by_id(std::vector<size_t> loaded_node_ids, DofContainer loaded_dofs, LoadContainer loads, GlobalMesh& glob_mesh)
        {
            std::set<size_t> loaded_nodes_on_rank = glob_mesh.filter_node_ids(loaded_node_ids, "rank_owned");
            NodalLoad nodal_load;
            nodal_load.assign_nodes_by_id(loaded_nodes_on_rank, glob_mesh);
            nodal_load.assign_dofs_loads(loaded_dofs, loads);
            nodal_loads.push_back(nodal_load);
        }
        
        /**
         * @brief initialise all load objects managed by this manager.
         * 
         */
        void initialise_loads()
        {
            for (auto& nodal_load : nodal_loads)
            {
                nodal_load.initialise_loads();
            }
        }

        /**
         * @brief increment the loads of all the nodal loads managed by this manager by the increment load_factor_increment.
         * 
         * @param load_factor_increment the multiplier by which to increment the loads: \f$ \Delta LF\f$.
         */
        void  increment_loads(real load_factor_increment)
        {
            for (auto&& nodal_load : nodal_loads)
            {
                nodal_load.increment_loads(load_factor_increment);
            }
        }

        /**
         * @brief remove all loads and clear loaded_dofs from the loaded nodes.
         * 
         */
        void remove_loads()
        {
            for (auto&& nodal_load : nodal_loads)
            {
                nodal_load.unload_loaded_nodes();
            }
        }

};

#endif 