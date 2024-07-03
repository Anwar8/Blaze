/**
 * @file LoadManager.hpp
 * @brief defines the LoadManager class which contains information about loads, and can apply them to the \ref GlobalMesh.
 */

#ifndef LOAD_MANAGER_HPP
#define LOAD_MANAGER_HPP
#include <vector>
#include <map>
#include "global_mesh.hpp"
#include "nodal_load.hpp"
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

        
    public:

        
        

};

#endif 