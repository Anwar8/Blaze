/**
 * @file NodalRestraint.hpp
 * @brief defines the NodalRestraint class which contains information about the restraint conditions for a set of nodes.
 */

#ifndef NODAL_RESTRAINT_HPP
#define NODAL_RESTRAINT_HPP

#include "maths_defaults.hpp"
#include "node.hpp"
#include "global_mesh.hpp"

#include <vector>
#include <set>

/**
 * @brief an object that is used to store information about the restraint conditions for a std::set of similarly-restrained nodes.
 */
class NodalRestraint
{
    protected:
        std::set<std::shared_ptr<Node>> restrained_nodes; /**< a std::set of shared pointers to nodes that are loaded.*/
        std::set<int> restrained_dofs; /**< a std set of loaded DoFs; none at first, then those loaded are added.*/
    public:
        /**
         * @brief assigns dofs to be restrained by this instant of the \ref NodalRestraint object.
         * 
         * @tparam Container STL container that is compatible with STL iterators and is used for storing the DoFs to be restrained.
         * @param dofs the DoFs to be restrained by this restraint object.
         */
        template <typename Container>
        void assign_dofs_loads(Container dofs)
        {
            for (auto dof : dofs)
            {
                restrained_dofs.insert(dof);
            }
        }
        

        /**
         * @brief assigns nodes by ID to the \ref NodalRestraint object. That is, the nodes that this object will restrain.
         * 
         * @tparam Container STL container that is compatible with standard STL iterators and contains node IDs.
         * @param node_ids the IDs of the nodes to be loaded.
         * @param glob_mesh the global mesh object that contains all the nodes of the model.
         */
        template <typename Container>
        void assign_nodes_by_id(Container node_ids, GlobalMesh& glob_mesh)
        {
            for (auto node_id : node_ids)
            {
                restrained_nodes.insert(glob_mesh.get_node_by_id(node_id));
            }
        }

        /**
         * @brief assigns the shared pointer to the nodes directly to the \ref restrained_nodes container.
         * 
         * @param nodes a shared_ptr to a node object that will be loaded.
         */
        void assign_nodes_by_ptr(std::vector<std::shared_ptr<Node>> nodes)
        {
            for (auto node : nodes)
            {
                restrained_nodes.insert(node);
            }
        }

        /**
         * @brief applies the restraint conditions to the nodes restrained by this object.
         * 
         */
        void apply_restraints()
        {
            for (auto node : restrained_nodes)
            {
                node->fix_dofs(restrained_dofs);
            }
        }
        /**
         * @brief removes the restraints this \ref NodalRestraint object imposes on the nodes from the nodes it imposes them on.
         */
        void free_restraints()
        {
            for (auto node : restrained_nodes)
            {
                node->free_dofs(restrained_dofs);
            }
        }

        /**
         * @brief clears the nodes set. Used for unit testing.
         */
        void clear_restrained_nodes()
        {
            restrained_nodes.clear();
        }

        /**
         * @brief clears the DoFs set. Used for unit testing.
         */
        void clear_dofs()
        {
            restrained_dofs.clear();
        }

        /**
         * @brief clear both the loads and the loaded nodes. Basically resets state of the nodal load. Useful for unit testing.
         * 
         */
        void reset()
        {
            clear_dofs();
            clear_restrained_nodes();
        }
       
        /**
         * @name getters
         * @brief used for testing, mostly.
         * 
         */
        //@{
        std::set<std::shared_ptr<Node>> get_restrained_nodes() const {return restrained_nodes;}
        std::set<int> get_restrained_dofs() const {return restrained_dofs;};
        //@}
};
#endif