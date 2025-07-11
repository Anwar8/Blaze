/**
 * @file LoadManager.hpp
 * @brief defines the NodalLoad class which contains information about the load conditions for a set of nodes.
 */

#ifndef NODAL_LOAD_HPP
#define NODAL_LOAD_HPP

#include "maths_defaults.hpp"
#include "node.hpp"
#include "global_mesh.hpp"

#include <vector>
#include <set>

/**
 * @brief an object that is used to store information about the load conditions for a vector of loaded nodes.
 * @todo consider whether to change the way the loaded nodes are stored from a vector to a std::set so that the same node cannot be loaded twice by the same load object.
 */
class NodalLoad
{
    protected:
        std::vector<std::shared_ptr<Node>> loaded_nodes; /**< a vector of shared pointers to nodes that are loaded.*/
        std::set<int> loaded_dofs; /**< a std set of loaded DoFs; none at first, then those loaded are added.*/
        std::array<real, 6> nodal_loads = {0., 0., 0., 0., 0., 0.}; /**< a std array containing 6 slots to be filled with nodal loads corresponding to dofs; initialised to zero.*/
    public:
        /**
         * @brief adds a load to a particular DOF of the \ref NodalLoad object.
         * 
         * @tparam Container STL container that can be dereferenced with the [] operator - used for the DOFs.
         * @tparam Container2 STL container that can be dereferenced with the [] operator - used for the loads corresponding to the DOFs.
         * @param dofs the DoFs to be loaded and must correspond to the loads in the loads parameter.
         * @param loads the loads corresponding to the DoFs in the dofs parameter.
         * @warning dofs and loads must be the same size, and must correspond to each other one-to-one. That is, first load corresponds to first DoF, second load corresponds to second DoF, etc.
         */
        template <typename Container, typename Container2>
        void assign_dofs_loads(Container dofs, Container2 loads)
        {
            if (dofs.size() != loads.size()) {
                std::cout << "NodalLoad::load_dofs expects same number of DoFs as number of loads. Got " << dofs.size() << " DoFs, but " << loads.size() << " loads." << std::endl; 
                exit(1);
            }
            for (int i = 0; i < dofs.size(); ++i)
            {
                loaded_dofs.insert(dofs[i]);
                nodal_loads[dofs[i]] = loads[i];
            }
        }
        /**
         * @brief adds a load to a particular DOF of the \ref NodalLoad object. This is a particular override in case the dofs parameter is a std::set already.
         * 
         * @tparam Container STL container that can be dereferenced with the [] operator - used for the loads corresponding to the DOFs.
         * @param dofs a std::set of DoFs to be loaded.
         * @param loads the loads corresponding to the DoFs in the dofs parameter; must be the same size as the dofs parameter.
         */
        template<typename Container>
        void assign_dofs_loads(std::set<int> dofs, Container loads)
        {
            if (dofs.size() != loads.size()) {
                std::cout << "NodalLoad::load_dofs expects same number of DoFs as number of loads. Got " << dofs.size() << " DoFs, but " << loads.size() << " loads." << std::endl; 
                exit(1);
            }
            int i = 0;
            for (auto& dof : dofs)
            {
                loaded_dofs.insert(dof);
                nodal_loads[dof] = loads[i];
                ++i;
            }
        }

        /**
         * @brief assigns nodes by ID to the \ref NodalLoad object. That is, the nodes that this object will load.
         * 
         * @tparam Container STL container that is compatible with standard STL iterators.
         * @param node_ids the IDs of the nodes to be loaded.
         * @param glob_mesh the global mesh object that contains all the nodes of the model.
         */
        template <typename Container>
        void assign_nodes_by_id(Container node_ids, GlobalMesh& glob_mesh)
        {
            for (auto& node_id : node_ids)
            {
                loaded_nodes.push_back(glob_mesh.get_node_by_id(node_id, "rank_owned"));
            }
        }

        /**
         * @brief assigns the shared pointer to the nodes directly to the \ref loaded_nodes container.
         * 
         * @param nodes a shared_ptr to a node object that will be loaded.
         */
        void assign_nodes_by_ptr(std::vector<std::shared_ptr<Node>> nodes)
        {
            for (auto& node : nodes)
            {
                loaded_nodes.push_back(node);
            }
        }
        /**
         * @brief initialises the nodal load on all loaded nodes by calling the \ref Node::add_nodal_load function with a 0.0 load. 
         * This prevents having to have an if-statement during load incrementation.
         * 
         */
        void initialise_loads()
        {
            for (auto& node : loaded_nodes)
            {
                for (auto& dof : loaded_dofs)
                {
                    node->add_nodal_load(0.0, dof);
                }
            }
        }
        /**
         * @brief increments the load for \ref loaded_nodes by an amount of the load equivalent to the \f$ \Delta LF\f$ given. 
         * 
         * @param load_factor_increment the \f$ \Delta LF\f$ increment in the applied load at all loaded DoFs.
         */
        void increment_loads(real load_factor_increment)
        {
            for (auto& loaded_node: loaded_nodes)
            {
                for (auto& dof : loaded_dofs)
                {
                    loaded_node->increment_nodal_load(nodal_loads[dof]*load_factor_increment, dof);
                }
            }
        }
        /**
         * @brief removes all loaded DoFs and resets the nodal_loads to zero. Useful for unit testing.
         */
        void clear_loads()
        {
            loaded_dofs.clear();
            nodal_loads = {0., 0., 0., 0., 0., 0.};
        }
        /**
         * @brief clears the nodes vector. Used for unit testing.
         */
        void clear_loaded_nodes()
        {
            loaded_nodes.clear();
        }
        /**
         * @brief clear both the loads and the loaded nodes. Basically resets state of the nodal load. Useful for unit testing.
         * 
         */
        void reset()
        {
            clear_loads();
            clear_loaded_nodes();
        }
        /**
         * @brief removes ALL loads from the loaded_nodes. Calls the \ref Node::clear_nodal_loads function.
         * 
         */
        void unload_loaded_nodes()
        {
            for (auto& node: loaded_nodes)
            {
                node->clear_nodal_loads();
            }
        }
        /**
         * @name getters
         * @brief used for testing, mostly.
         * 
         */
        //@{
        std::vector<std::shared_ptr<Node>> get_loaded_nodes() const {return loaded_nodes;}
        int get_num_loaded_nodes() const {return loaded_nodes.size();}
        std::set<int> get_loaded_dofs() const {return loaded_dofs;};
        std::array<real, 6> get_nodal_loads() const {return nodal_loads;}
        //@}
};
#endif