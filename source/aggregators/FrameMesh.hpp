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
        }

        /**
         * @brief Gets the vertex ID for a given column line and floor.
         * 
         * @param column_line The column line index.
         * @param floor The floor index.
         * @return size_t The vertex ID.
         */
        size_t get_vertix_id(int column_line, int floor)
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
         * @return std::set<size_t> A set of vertex IDs at the specified floor.
         */
        std::set<size_t> get_vertices_ids_at_floor(int floor)
        {
            if (floor < 0 || floor > nfloors)
            {
                std::cout << "FrameMesh::get_vertices_ids_at_floor expects floor >= 0 and <= nfloors, but got (floor/nfloors) = (" << floor << "/" << nfloors << ")" << std::endl;
                exit(1); 
            }
            std::set<size_t> vertices;
            for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
            {
                vertices.insert(get_vertix_id(column_line_i, floor));
            }
            return vertices;
        }

        /**
         * @brief Gets the vertices ids for the entire model.
         * 
         * @return std::set<size_t> vertices_ids a std::set of vertices that are the ids of all vertices in the model.
         */
        std::set<size_t> get_vertices_ids()
        {
            std::set<size_t> vertices_ids;
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
         * @return std::set<size_t> nodes a std::set of node ids for a column.
         */
        std::set<size_t> get_column_node_ids(int column_line, int floor)
        {
            if (column_line < 0 || column_line > nbays || floor < 0 || floor >= nfloors)
            {
                std::cout << "FrameMesh::get_column_node_ids received incorrect column_line or floor. (column_line/nbays) = ("
                          << column_line << "/" << nbays << "), and (floor/nfloors) = (" << floor << "/" << nfloors << "). "
                          << "Valid column_line range: [0, " << nbays << "], valid floor range: [0, " << nfloors - 1 << "]" << std::endl;
                exit(1);
            }
            std::set<size_t> nodes;

            size_t starting_id = 1 + column_line*nodes_per_column_line + floor*column_divisions + column_line*nodes_per_full_bay;
            for (size_t i = 1; i < nodes_per_column + 1; ++i)
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
         * @return std::set<size_t> A set of node IDs representing the nodes of the beam.
         */
        std::set<size_t> get_beam_node_ids(int bay, int floor, bool include_vertices = false)
        {
            if ((bay < 1) || (floor < 1) || floor > nfloors || bay > nbays)
            {
                std::cout << "FrameMesh::get_beam_node_ids expects bay in range [1, " << nbays << "] and floor in range [1, " << nfloors << "], but got bay = " << bay << " and floor = " << floor << std::endl;
                exit(1);
            }
            std::set<size_t> nodes;
            
            if (include_vertices)
            {
                nodes.insert(get_vertix_id(bay - 1, floor));
            }

            size_t starting_id = bay*nodes_per_column_line + (bay - 1)*nodes_per_full_bay + (floor - 1) * nodes_per_beam;
            for (size_t i = 1; i < nodes_per_beam + 1; ++i)
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
         * @return std::set<size_t> A set of node IDs for the column line at the specified floor.
         */
        std::set<size_t> get_column_line_node_ids(int column_line)
        {
            if (column_line < 0 || column_line > nbays)
            {
                std::cout << "FrameMesh::get_column_line_node_ids expects column_line >= 0 and <= nbays, but got (column_line/nbays) = (" << column_line << "/" << nbays << ")" << std::endl;
                exit(1);
            }
            std::set<size_t> nodes;
            std::set<size_t> column_node_ids;

            for (size_t floor_i = 0; floor_i < nfloors; ++floor_i)
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
         * @return std::set<size_t> A set of node IDs for the beam line at the specified beam line.
         */
        std::set<size_t> get_beam_line_node_ids(int floor, bool include_vertices = false)
        {
            if (floor <= 0 || floor > nfloors)
            {
                std::cout << "FrameMesh::get_beam_line_node_ids expects floor > 0 and <= nfloors, but got (floor/nfloors) = (" << floor << "/" << nfloors << ")" << std::endl;
                exit(1); 
            }
            std::set<size_t> nodes;
            std::set<size_t> bay_node_ids;

            for (int bay_i = 1; bay_i <= nbays; ++bay_i)
            {
                bay_node_ids = get_beam_node_ids(bay_i, floor, include_vertices);
                nodes.insert(bay_node_ids.begin(), bay_node_ids.end());           
            }
            return nodes;
        }

        /**
         * @brief Creates node coordinate pairs for column nodes.
         * 
         * This function generates a vector of pairs, where each pair consists of a node ID and its corresponding coordinates.
         * The nodes are created for each column line in the frame.
         * 
         * @return std::vector<std::pair<size_t, coords>> A vector of pairs mapping node IDs to their coordinates for column nodes.
         */
        std::vector<std::pair<size_t, coords>> create_column_node_coords_pairs()
        {
            std::vector<std::pair<size_t, coords>> nodes_id_coords_pairs;
            for (int column_line_i = 0; column_line_i <= nbays; ++column_line_i)
            {
                std::set<size_t> node_ids = get_column_line_node_ids(column_line_i);
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
         * @return std::vector<std::pair<size_t, coords>> A vector of pairs mapping node IDs to their coordinates for beam nodes.
         */
        std::vector<std::pair<size_t, coords>> create_beam_node_coords_pairs()
        {
            std::vector<std::pair<size_t, coords>> nodes_id_coords_pairs;
            std::set<size_t> vertices = get_vertices_ids();
            for (int floor_i = 1; floor_i <= nfloors; ++floor_i)
            {
                std::set<size_t> node_ids = get_beam_line_node_ids(floor_i, true);
                real x_0 = 0.0;
                real y_0 = floor_i * floor_height;
                real x = x_0;
                real y = y_0;
                for (size_t node_id : node_ids)
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
         * @return std::vector<std::pair<size_t, coords>> A sorted vector of pairs mapping node IDs to their coordinates for the entire frame.
         */
        std::vector<std::pair<size_t, coords>> get_node_coords_pairs()
        {
           std::vector<std::pair<size_t, coords>> node_map, beam_node_map;

           node_map = create_column_node_coords_pairs();
           beam_node_map = create_beam_node_coords_pairs();

           node_map.insert(node_map.end(), beam_node_map.begin(), beam_node_map.end());
           std::sort(node_map.begin(), node_map.end(), 
                                                      [](const auto& a, const auto& b) 
                                                      {
                                                        return a.first < b.first;
                                                      }
                                                    );
           return node_map;
        }

        /**
         * @brief Maps elements to nodes.
         * 
         * @details This function generates a vector of pairs, where each pair consists of an element ID and a vector of node IDs that the element connects to.
         * The elements are created for each bay and floor in the frame.
         * 
         * @return std::vector<std::pair<size_t, std::vector<size_t>>> A vector of pairs mapping element IDs to their corresponding node IDs.
         */
        std::vector<std::pair<size_t, std::vector<size_t>>> map_elements_to_nodes()
        {
            std::vector<std::pair<size_t, std::vector<size_t>>> elements_map;
            elements_map.reserve(nbays * nfloors * beam_divisions + (nbays + 1)*nfloors*column_divisions);
            size_t element_id = 0;

            for (int bay_i = 1; bay_i <= nbays; ++bay_i) {
                // first column line of each bay
                std::set<size_t> node_ids = get_column_line_node_ids(bay_i - 1);
                auto it = node_ids.begin();
                for (size_t node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                    element_id++;
                    elements_map.push_back(std::make_pair(element_id, std::vector<size_t>{*it, *(++it)}));
                }
                // beams on each floor
                for (int floor_i = 1; floor_i <= nfloors; ++floor_i) {
                    node_ids = get_beam_node_ids(bay_i, floor_i, true);
                    it = node_ids.begin();
                    for (size_t node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                        element_id++;
                        elements_map.push_back(std::make_pair(element_id, std::vector<size_t>{*it, *(++it)}));
                    }
                }
            }
            // last column line
            std::set<size_t> node_ids = get_column_line_node_ids(nbays);
            auto it = node_ids.begin();
            for (size_t node_id_i = 0; node_id_i < node_ids.size() - 1; ++node_id_i) {
                element_id++;
                elements_map.push_back(std::make_pair(element_id, std::vector<size_t>{*it, *(++it)}));
            }
            return elements_map;
        }

        /**
         * @brief Gets the IDs of the frame that will need to be restrained out-of-plane, which is all except bases.
         * 
         * @return std::set<size_t> A set of node IDs representing the nodes that will need out-of-plane restrain.
         */
        std::set<size_t> get_out_of_plane_nodes()
        {
            std::set<size_t> column_bases = get_column_bases();
            std::set<size_t> all_node_ids;
            std::set<size_t> temp_id_set;
            std::set<size_t> out_of_plane_nodes;
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
         * @return std::set<size_t> nodes at the base of the column that will be fixed.
         */
        std::set<size_t> get_column_bases()
        {
            return get_vertices_ids_at_floor(0);   
        }

};

/**
 * @brief prints out each element ID and its corresponding node IDs to the standard output stream.
 * 
 * @param element_map std::vector<std::pair<size_t, std::vector<size_t>>> that contains element ID followed by node IDs for this element.
 */
void inline read_element_map(std::vector<std::pair<size_t, std::vector<size_t>>> element_map)
{
    for (auto& pair : element_map)
    {
        std::cout << "Element " << pair.first << ", nodes: (" << pair.second[0] << ", " << pair.second[1] << ")" << std::endl;
    }
}

/**
 * @brief prints out the node IDs and their xyz coordinates to the standard output stream.
 * 
 * @param node_map std::vector<std::pair<size_t, coords>> containing node IDs followed by coordinates of the node.
 */
void inline read_node_map(std::vector<std::pair<size_t, coords>> node_map)
{
    for (auto& pair : node_map)
    {
        std::cout << "Node " << pair.first << ", xyz = (" << pair.second(0) << ", " << pair.second(1) << ", " << pair.second(2) << ")" << std::endl;
    }
}
#endif