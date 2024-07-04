/**
 * @file Scribe.hpp
 * @brief defines the Scribe class which tracks state of the model and records output.
 */

#ifndef SCRIBE_HPP
#define SCRIBE_HPP

#include "maths_defaults.hpp"
#include "node.hpp"
#include "global_mesh.hpp"
#include "Record.hpp"

#include <vector>
#include <set>


/**
 * @brief  the size of the buffer used to store the data beyond which the data has to be flushed to file.
 */
#define BUFFER_SIZE 1024

/**
 * @brief A scribe manages the recording of the state of the model.
 */
class Scribe
{
    protected:
        // std::vector<std::shared_ptr<Node>> tracked_nodes; /**< a vector of shared pointers to nodes that are tracked.*/
        // std::set<int> tracked_dofs; /**< a std set of tracked DoFs; none at first, then those tracked are added.*/
        
        std::vector<Record> record_library; /**< a vector of records that are used to store the \ref Record objects for all tracked nodes.*/
        int current_row = 0; /**< the current row in the recorded data that is being filled. Used for deciding when the data needs flushing.*/
        int buffer_size = BUFFER_SIZE; /**< the size of the buffer used to store the data beyond which the data has to be flushed to file.*/

    public:


        /**
         * @brief assigns nodes by ID to the \ref Scribe object. That is, the nodes that this object will track.
         * 
         * @tparam Container STL container that is compatible with standard STL iterators.
         * @param node_ids the IDs of the nodes to be tracked.
         * @param glob_mesh the global mesh object that contains all the nodes of the model.
         */
        template <typename Container>
        void track_nodes_by_id(Container node_ids, std::set<int> dofs,  GlobalMesh& glob_mesh)
        {
            for (auto node_id : node_ids)
            {
                record_library.push_back(Record(glob_mesh.get_node_by_id(node_id)), dofs);
            }
            sort_record_library();
        }

        /**
         * @brief Creates a record directly with the shared pointer to the nodes.
         * 
         * @param nodes a shared_ptr to a node object that will be tracked by the created record.
         * @param dofs a std set of DoFs that are to be tracked by the record.
         */
        void track_nodes_by_ptr(std::vector<std::shared_ptr<Node>> nodes, std::set<int> dofs)
        {
            for (auto node : nodes)
            {
                record_library.push_back(Record(node, dofs, buffer_size));
            }
            sort_record_library();
        }

        /**
         * @brief sorts \ref record_library by the ID of the nodes being tracked.
         * 
         */
        void sort_record_library()
        {
            std::sort(record_library.begin(), record_library.end());
        }

        /**
         * @brief writes the current state of the tracked nodes to all the records.
         * 
         */
        void write_to_records()
        {
            for (auto record: record_library)
            {
                record.write_to_record(current_row);
            }
            ++current_row;

            // This check is only done once every time we write all the records.
            if (current_row >= buffer_size)
            {        
                flush_records();
                current_row = 0;
            }
        }

        /**
         * @brief flushes the records to file. Should use HDF5, but currently is not implemented.
         * 
         */
        void flush_records()
        {
            std::cout << "current_row = " << current_row << "Flushing records to file still not implemented. Exiting." << std::endl;
            exit(1);
        }
        /**
         * @brief reads the contents of a particular record corresponding to a particular node ID to the output stream.
         * 
         */
        void read_a_record(int node_id)
        {
            auto record_it =  get_id_iterator<std::vector<Record>::iterator, std::vector<Record>>(node_id, record_library);
            record_it->read_record();
        }
        /**
         * @brief reads the contents of all records in the \ref records_library to the output stream.
         * 
         */
        void read_all_records()
        {
            for (auto record: record_library)
            {
                record.read_record();
            }
        }

};
#endif