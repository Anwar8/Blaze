/**
 * @file FrameMesh.hpp
 * @brief Defines the FrameMesh class which creates meshes of portal frames for performance evaluation of Blaze.
 */
#ifndef FRAME_MESH
#define FRAME_MESH
#include "maths_defaults.hpp"
#include <set>
#include <utility>
#include <iostream>
#include <algorithm>

/**
 * @brief FrameMesh class for creating and managing frame meshes.
 * 
 * The FrameMesh class provides functionality to create frame meshes, calculate node coordinates, and map nodes and coordinates, and elements and nodes.
 */
class FrameMesh {
    protected:
        int nbays = -1; /**<Number of bays.*/
        int nfloors = -1; /**<Number of floors.*/
        real bay_length = 0.0;  /**<Length of each bay.*/
        real floor_height = 0.0; /**<Height of each floor.*/
        int column_divisions = 0.0; /**<Number of elements per column.*/
        int beam_divisions = 0.0; /**<Number of elements per beam.*/

        int nodes_per_column_line = 0; /**<Number of nodes per column line (all floors) including vertices.*/
        int nodes_per_full_bay = 0; /**<Number of nodes, excluding vertices, for all floor at a given bay. */
        int nodes_per_column = 0; /**<Number of nodes per column excluding vertices.*/
        int nodes_per_beam = 0; /**<Number of nodes per beam excluding vertices.*/
        int num_nodes = 0; /**<Total number of nodes in model.*/
        int num_elements = 0; /**<Total number of elements in model.*/
        real dx = 0.0; /**<Length of each element for the beams.*/
        real dy = 0.0; /**<Length of each element for the columns.*/

    public:

        FrameMesh() = default;
        /**
         * @brief Constructor for the FrameMesh class.
         * 
         * Initializes the FrameMesh with the given parameters.
         * 
         * @param nbays Number of bays in the frame.
         * @param nfloors Number of floors in the frame.
         * @param bay_length Length of each bay.
         * @param floor_height Height of each floor.
         * @param beam_divisions Number of divisions in each beam.
         * @param column_divisions Number of divisions in each column.
         */
        FrameMesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions)
        {
            this->nbays = nbays;
            this->nfloors = nfloors;
            this->bay_length = bay_length;
            this->floor_height = floor_height;
            this->beam_divisions = beam_divisions;
            this->column_divisions = column_divisions;
            basic_counts();
        }

        /**
         * @brief Calculates \ref nodes_per_column_line, \ref nodes_per_full_bay, \ref nodes_per_column, \ref nodes_per_beam, \ref dx, \ref dy, and \ref num_nodes
         * 
         */
        void basic_counts()
        {
            nodes_per_column_line = nfloors*column_divisions + 1;
            nodes_per_full_bay = nfloors*(beam_divisions - 1);

            nodes_per_column = column_divisions - 1;
            nodes_per_beam = beam_divisions - 1;

            dx = bay_length/beam_divisions;
            dy = floor_height/column_divisions;
            num_nodes = nodes_per_column_line * (nbays + 1) + nodes_per_full_bay*(nbays);
            num_elements = (nodes_per_column_line - 1) * (nbays + 1) + (beam_divisions * nfloors * nbays);
        }

        /**
         * @brief Gets the vertex ID for a given column line and floor.
         * 
         * @param column_line The column line index.
         * @param floor The floor index.
         * @return unsigned The vertex ID.
         */
        unsigned get_vertix_id(int column_line, int floor)
        {
            if (column_line >= 0 && column_line <= nbays && floor >= 0 && floor <= nfloors)
            {
                return 1 + floor*(column_divisions) + column_line*nodes_per_column_line + column_line*nodes_per_full_bay;
            } 
            else
            {
                std::cout << "FrameMesh::get_vertix_id received incorrect column_line or floor.  (column_line/nbays) = (" 
                << column_line << "/" << nbays << "), and (floor/nfloors) = (" << floor << "/" << nfloors << ")" << std::endl; 
                exit(1);
            }
        }

        /**
         * @brief Gets the vertex IDs for all vertices at a specific floor and all bays.
         * 
         * @param floor The floor at which to get the vertices - can be from 0 to \ref nfloors.
         * @return std::set<unsigned> A set of vertex IDs at the specified floor.
         */
        std::set<unsigned> get_vertices_ids_at_floor(int floor)
        {
            if (floor < 0 || floor > nfloors)
            {
                std::cout << "FrameMesh::get_vertices_ids_at_floor expects floor >= 0 and <= nfloors, but got (floor/nfloors) = (" << floor << "/" << nfloors << ")" << std::endl;
                exit(1); 
            }
            std::set<unsigned> vertices;
            for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
            {
                vertices.insert(get_vertix_id(column_line_i, floor));
            }
            return vertices;
        }

        /**
         * @brief Gets the vertices ids for the entire model.
         * 
         * @return std::set<unsigned> vertices_ids a std::set of vertices that are the ids of all vertices in the model.
         */
        std::set<unsigned> get_vertices_ids()
        {
            std::set<unsigned> vertices_ids;
            for (int floor_i = 0; floor_i <= nfloors; ++floor_i)
            {
                for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
                {
                    vertices_ids.insert(get_vertix_id(column_line_i, floor_i));
                }
            }
            return vertices_ids;
        }

        /**
         * @brief returns the column node ids at a given floor for a certain column_line including vertices. 
         * @param column_line the column line which the column we want is at. Ranges from 0 to \ref nbays.
         * @param floor floor at which the column starts. Ranges from 0 to < \ref nfloors.
         * @return std::set<unsigned> nodes a std::set of node ids for a column.
         */
        std::set<unsigned> get_column_node_ids(int column_line, int floor)
        {
            if (column_line < 0 || column_line > nbays || floor < 0 || floor >= nfloors)
            {
                std::cout << "FrameMesh::get_column_node_ids received incorrect column_line or floor. (column_line/nbays) = ("
                          << column_line << "/" << nbays << "), and (floor/nfloors) = (" << floor << "/" << nfloors << "). "
                          << "Valid column_line range: [0, " << nbays << "], valid floor range: [0, " << nfloors - 1 << "]" << std::endl;
                exit(1);
            }
            std::set<unsigned> nodes;

            unsigned starting_id = 1 + column_line*nodes_per_column_line + floor*column_divisions + column_line*nodes_per_full_bay;
            for (unsigned i = 1; i < nodes_per_column + 1; ++i)
            {
                nodes.insert(starting_id + i);
            }
            return nodes;
        }

        /**
         * @brief Returns the node ids for a given beam either including or excluding the vertices.
         * @param bay bay in which the beam is. Ranges from 1 to \ref nbays.
         * @param floor floor at which the beam is. Ranges from 1 to \ref nfloors.
         * @param include_vertices boolean to include or exclude vertices of the beam.
         * @return std::set<unsigned> A set of node IDs representing the nodes of the beam.
         */
        std::set<unsigned> get_beam_node_ids(int bay, int floor, bool include_vertices = false)
        {
            if ((bay < 1) || (floor < 1) || floor > nfloors || bay > nbays)
            {
                std::cout << "FrameMesh::get_beam_node_ids expects bay in range [1, " << nbays << "] and floor in range [1, " << nfloors << "], but got bay = " << bay << " and floor = " << floor << std::endl;
                exit(1);
            }
            std::set<unsigned> nodes;
            
            if (include_vertices)
            {
                nodes.insert(get_vertix_id(bay - 1, floor));
            }

            unsigned starting_id = bay*nodes_per_column_line + (bay - 1)*nodes_per_full_bay + (floor - 1) * nodes_per_beam;
            for (unsigned i = 1; i < nodes_per_beam + 1; ++i)
            {
                nodes.insert(starting_id + i);
            }
            if (include_vertices)
            {
                nodes.insert(get_vertix_id(bay, floor));
            }

            return nodes;
        }

        /**
         * @brief Gets the column line node IDs for a specific floor.
         * 
         * @param column_line The column_line index. Ranges from 0 to \ref nbays.
         * @return std::set<unsigned> A set of node IDs for the column line at the specified floor.
         */
        std::set<unsigned> get_column_line_node_ids(int column_line)
        {
            if (column_line < 0 || column_line > nbays)
            {
                std::cout << "FrameMesh::get_column_line_node_ids expects column_line >= 0 and <= nbays, but got (column_line/nbays) = (" << column_line << "/" << nbays << ")" << std::endl;
                exit(1);
            }
            std::set<unsigned> nodes;
            std::set<unsigned> column_node_ids;

            for (unsigned floor_i = 0; floor_i < nfloors; ++floor_i)
            {   
                nodes.insert(get_vertix_id(column_line, floor_i));
                column_node_ids = get_column_node_ids(column_line, floor_i);
                nodes.insert(column_node_ids.begin(), column_node_ids.end());
            }
            nodes.insert(get_vertix_id(column_line, nfloors));
            return nodes;
        }

        /**
         * @brief Gets the beam line node IDs for a specific floor. 
         * 
         * @param floor the floor at which the beam is located. Ranges from 1 to \ref nfloors.
         * @param include_vertices boolean to include or exclude vertices of the beam.
         * @return std::set<unsigned> A set of node IDs for the beam line at the specified beam line.
         */
        std::set<unsigned> get_beam_line_node_ids(int floor, bool include_vertices = false)
        {
            if (floor <= 0 || floor > nfloors)
            {
                std::cout << "FrameMesh::get_beam_line_node_ids expects floor > 0 and <= nfloors, but got (floor/nfloors) = (" << floor << "/" << nfloors << ")" << std::endl;
                exit(1); 
            }
            std::set<unsigned> nodes;
            std::set<unsigned> bay_node_ids;

            for (int bay_i = 1; bay_i <= nbays; ++bay_i)
            {
                bay_node_ids = get_beam_node_ids(bay_i, floor, include_vertices);
                nodes.insert(bay_node_ids.begin(), bay_node_ids.end());           
            }
            return nodes;
        }

        /**
         * @brief Gets the beam line node IDs for all floors. 
         * 
         * @param include_vertices boolean to include or exclude vertices of the beam.
         * @return std::set<unsigned> A set of node IDs for the beam lines of all floors.
         */
        std::set<unsigned> get_all_beam_line_node_ids(bool include_vertices = false)
        {
            std::set<unsigned> nodes;
            std::set<unsigned> floor_node_ids;
            for (int floor = 1; floor <= nfloors; ++floor)
            {
                floor_node_ids = get_beam_line_node_ids(floor, include_vertices);
                nodes.insert(floor_node_ids.begin(), floor_node_ids.end());           
            }
            return nodes;
        }

        /**
         * @brief Creates node coordinate pairs for column nodes.
         * 
         * This function generates a vector of pairs, where each pair consists of a node ID and its corresponding coordinates.
         * The nodes are created for each column line in the frame.
         * 
         * @return std::vector<std::pair<unsigned, coords>> A vector of pairs mapping node IDs to their coordinates for column nodes.
         */
        std::vector<std::pair<unsigned, coords>> create_column_node_coords_pairs()
        {
            std::vector<std::pair<unsigned, coords>> nodes_id_coords_pairs;
            for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
            {
                std::set<unsigned> node_ids = get_column_line_node_ids(column_line_i);
                real x_0 = column_line_i * bay_length;
                real y_0 = 0.0;
                real x = x_0;
                real y = y_0;
                for (int node_id : node_ids)
                {
                    nodes_id_coords_pairs.push_back(std::make_pair(node_id, coords(x, y, 0.0)));
                    y += dy;
                }
            }
            return nodes_id_coords_pairs;
        }
        /**
         * @brief Creates node coordinate pairs for beam nodes.
         * 
         * This function generates a vector of pairs, where each pair consists of a node ID and its corresponding coordinates.
         * The nodes are created for each beam line in the frame, *excluding* the vertices as these are created by \ref create_column_node_coords_pairs.
         * 
         * @return std::vector<std::pair<unsigned, coords>> A vector of pairs mapping node IDs to their coordinates for beam nodes.
         */
        std::vector<std::pair<unsigned, coords>> create_beam_node_coords_pairs()
        {
            std::vector<std::pair<unsigned, coords>> nodes_id_coords_pairs;
            std::set<unsigned> vertices = get_vertices_ids();
            for (int floor_i = 1; floor_i <= nfloors; ++floor_i)
            {
                std::set<unsigned> node_ids = get_beam_line_node_ids(floor_i, true);
                real x_0 = 0.0;
                real y_0 = floor_i * floor_height;
                real x = x_0;
                real y = y_0;
                for (unsigned node_id : node_ids)
                {
                    if (!vertices.count(node_id))
                    {
                        nodes_id_coords_pairs.push_back(std::make_pair(node_id, coords(x, y, 0.0)));
                    }
                    x += dx;
                }
            }
            return nodes_id_coords_pairs;
        }
        /**
         * @brief Gets the node coordinate pairs for the entire frame sorted by ID.
         * 
         * @details This function combines the node coordinate pairs for both column and beam nodes, sorts them by node ID, and returns the result.
         * 
         * @return std::vector<std::pair<unsigned, coords>> A sorted vector of pairs mapping node IDs to their coordinates for the entire frame.
         */
        std::vector<std::pair<unsigned, coords>> get_node_coords_pairs()
        {
           std::vector<std::pair<unsigned, coords>> nodes_coords_vector, beam_nodes_coords_vector;

           nodes_coords_vector = create_column_node_coords_pairs();
           beam_nodes_coords_vector = create_beam_node_coords_pairs();

           nodes_coords_vector.insert(nodes_coords_vector.end(), beam_nodes_coords_vector.begin(), beam_nodes_coords_vector.end());
           std::sort(nodes_coords_vector.begin(), nodes_coords_vector.end(), 
                                                      [](const auto& a, const auto& b) 
                                                      {
                                                        return a.first < b.first;
                                                      }
                                                    );
           return nodes_coords_vector;
        }

        /**
         * @brief Maps elements to nodes.
         * 
         * @details This function generates a vector of pairs, where each pair consists of an element ID and a vector of node IDs that the element connects to.
         * The elements are created for each bay and floor in the frame.
         * 
         * @return std::vector<std::pair<unsigned, std::vector<unsigned>>> A vector of pairs mapping element IDs to their corresponding node IDs.
         */
        std::vector<std::pair<unsigned, std::vector<unsigned>>> map_elements_to_nodes()
        {
            std::vector<std::pair<unsigned, std::vector<unsigned>>> elements_map;
            elements_map.reserve(nbays * nfloors * beam_divisions + (nbays + 1)*nfloors*column_divisions);
            unsigned element_id = 0;

            for (int bay_i = 1; bay_i <= nbays; ++bay_i) {
                // first column line of each bay
                std::set<unsigned> node_ids = get_column_line_node_ids(bay_i - 1);
                auto it = node_ids.begin();
                for (unsigned node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                    element_id++;
                    elements_map.push_back(std::make_pair(element_id, std::vector<unsigned>{*it, *(++it)}));
                }
                // beams on each floor
                for (int floor_i = 1; floor_i <= nfloors; ++floor_i) {
                    node_ids = get_beam_node_ids(bay_i, floor_i, true);
                    it = node_ids.begin();
                    for (unsigned node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                        element_id++;
                        elements_map.push_back(std::make_pair(element_id, std::vector<unsigned>{*it, *(++it)}));
                    }
                }
            }
            // last column line
            std::set<unsigned> node_ids = get_column_line_node_ids(nbays);
            auto it = node_ids.begin();
            for (unsigned node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                element_id++;
                elements_map.push_back(std::make_pair(element_id, std::vector<unsigned>{*it, *(++it)}));
            }
            return elements_map;
        }

        /**
         * @brief Gets the IDs of the frame that will need to be restrained out-of-plane, which is all except bases.
         * 
         * @return std::set<unsigned> A set of node IDs representing the nodes that will need out-of-plane restrain.
         */
        std::set<unsigned> get_out_of_plane_nodes()
        {
            std::set<unsigned> column_bases = get_column_bases();
            std::set<unsigned> all_node_ids;
            std::set<unsigned> temp_id_set;
            std::set<unsigned> out_of_plane_nodes;
            for (int floor_i = 1; floor_i <= nfloors; ++floor_i)
            {
                temp_id_set = get_beam_line_node_ids(floor_i);
                all_node_ids.insert(temp_id_set.begin(), temp_id_set.end());
            }
            for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
            {
                temp_id_set = get_column_line_node_ids(column_line_i);
                all_node_ids.insert(temp_id_set.begin(), temp_id_set.end());
            }
            
            std::set_difference(all_node_ids.begin(), all_node_ids.end(), 
                                column_bases.begin(), column_bases.end(),
                                std::inserter(out_of_plane_nodes, out_of_plane_nodes.begin()));
            return out_of_plane_nodes;
        }

        /**
         * @brief Get the node IDs that correspond to the bases of the column that will be fixed.
         * 
         * @return std::set<unsigned> nodes at the base of the column that will be fixed.
         */
        std::set<unsigned> get_column_bases()
        {
            return get_vertices_ids_at_floor(0);   
        }

        /**
         * @brief prints the \ref num_nodes and \ref num_elements to the standard output stream.
         * 
         */
        void read_frame_size()
        {
            std::cout << "Frame has " << num_nodes << " nodes and " << num_elements << " elements." << std::endl; 
        }
        
        /**
         * @brief returns the \ref num_nodes and \ref num_elements to the standard output stream.
         * 
         */
        std::pair<int, int> get_frame_size()
        {
            return std::make_pair(num_nodes, num_elements);
        }

        int get_num_elements() {return num_elements;}
        int get_num_nodes() {return num_nodes;}

};

/**
 * @brief prints out each element ID and its corresponding node IDs to the standard output stream.
 * 
 * @param element_map std::vector<std::pair<unsigned, std::vector<unsigned>>> that contains element ID followed by node IDs for this element.
 */
void inline read_element_map(std::vector<std::pair<unsigned, std::vector<unsigned>>> element_map)
{
    for (auto& pair : element_map)
    {
        std::cout << "Element " << pair.first << ", nodes: (" << pair.second[0] << ", " << pair.second[1] << ")" << std::endl;
    }
}

/**
 * @brief prints out the node IDs and their xyz coordinates to the standard output stream.
 * 
 * @param nodes_coords_vector std::vector<std::pair<unsigned, coords>> containing node IDs followed by coordinates of the node.
 */
void inline read_nodes_coords_vector(std::vector<std::pair<unsigned, coords>> nodes_coords_vector)
{
    for (auto& pair : nodes_coords_vector)
    {
        std::cout << "Node " << pair.first << ", xyz = (" << pair.second(0) << ", " << pair.second(1) << ", " << pair.second(2) << ")" << std::endl;
    }
}
#endif