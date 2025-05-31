/**
 * @file Record.hpp
 * @brief defines the \ref Record class which is used by the scribe to store state data.
 */
#ifndef RECORD_HPP
#define RECORD_HPP
#include "maths_defaults.hpp"
#include "node.hpp"
#include "basic_utilities.hpp"
#include <set>
#include <utility>
#include <map>
#include <memory>


class Record {
    protected:
        std::shared_ptr<Node> tracked_node; /**< a shared pointer to the node that is being tracked by this record.*/
        unsigned tracked_node_id; /**< the ID of the node that is being tracked by this record.*/
        std::array<std::vector<real>, 6> recorded_data; /**< the data that is recorded in this record by the scribe.*/
        std::set<int> tracked_dofs; /**< a std set of tracked DoFs as decided by the scribe.*/
        
        bool full = false; /**< tells if the record is full and requires flushing. */
    
    public:
        Record() = default;

        // Record(std::shared_ptr<Node> node, std::set<int> dofs, int buffer_size)
        // {
        //     initialise_record(node, dofs, buffer_size);
        // } 

        /**
         * @brief initialises the record with the number of columns and rows that the recorded data will have, as well as the node being tracked and its ID.
         * 
         * @param node a shared_ptr to the node that is being tracked by this record.
         * @param dofs a std set of the DoFs that are being tracked by this record.
         */
        void initialise_record(std::shared_ptr<Node> node, std::set<int> dofs, int buffer_size)
        {
            
            tracked_node = node;
            tracked_node_id = node->get_id();
            tracked_dofs = dofs;
        }
        
        /**
         * @brief writes the current state of the tracked node to the record and checks if it is full.
         * 
         */
        void write_to_record(int row)
        {
            int i = 0;
            for (auto& dof : tracked_dofs)
            {
                real displacement = tracked_node->get_nodal_displacement(dof);
                (this->recorded_data[dof]).push_back(displacement);
                ++i;
            }
        }
        
        /**
         * @brief overloads the less than operator to compare records by their tracked node ID, allowing easy sorting of record libraries by \ref Scribe objects.
         * 
         */
        bool operator<(const Record& other_record) const
        { 
            return tracked_node_id < other_record.tracked_node_id; 
        } 

        /**
         * @brief reads the contents of the record to the output stream.
         * 
         */
        void read_record()
        {
            std::cout << "Record for node " << tracked_node_id << " tracking DoFs:";
            print_container(tracked_dofs);
            std::cout << "Record contents are: " << std::endl;
            for (auto& dof: tracked_dofs)
            {
                print_container(this->recorded_data[dof]);
            }
            
        }

        /**
         * @brief reads the contents of the record to the output stream for a particular row.
         * @param i the row number to read.
         */
        void read_record_at(int i)
        {
            std::cout << "Record for node " << tracked_node_id << " tracking DoFs:";
            print_container(tracked_dofs);
            std::cout << "Record contents at i = " << i << " are: " << std::endl;
            for (auto& dof: tracked_dofs)
            {
                std::cout << "dof: [" << dof << "] = " << this->recorded_data[dof][i] << ", ";
            }
            std::cout << std::endl;
        }
        /**
         * @brief Get the tracked node id.
         * 
         * @return unsigned 
         */
        unsigned get_tracked_node_id() const {return tracked_node_id;}
        
        /**
         * @brief Get the tracked node shared ptr.
         * 
         * @return std::shared_ptr<Node> tracked node.
         */
        std::shared_ptr<Node> get_tracked_node() const {return tracked_node;}

        /**
         * @brief Get the tracked dofs set.
         * 
         * @return std::set<int> tracked DoFs.
         */
        std::set<int> get_tracked_dofs() const {return tracked_dofs;}

        /**
         * @brief Get the recorded data.
         * 
         * @return mat recorded data.
         */
        std::array<std::vector<real>,6> get_recorded_data() const {return this->recorded_data;}

        /**
         * @brief overloads the less than operator to compare records by their tracked node ID, allowing easy sorting of node STL containers via \ref std::sort.
         * 
         */
        bool operator<(const Record& other_record) const
        { 
            return tracked_node_id < other_record.tracked_node_id; 
        } 

};

#endif