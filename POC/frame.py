import numpy as np
import matplotlib.pyplot as plt

NODE_COLOR = {"column":"g", "beam":"b", "vertix":"r"}

def find_node_by_id(nodes, id):
    for node in nodes:
        if node.id == id:
            return node
    print(f"Cannot find node with id = {id}")
    return None

class node:
    def __init__(self, id, x, y, node_type):
        self.id = id
        self.label = id
        self.x = x
        self.y = y
        self.node_type = node_type
    def print(self):
        print(f"node {self.label} at ({self.x},{self.y})")
    def set_rank(self, rank):
        self.rank = rank
        self.label = f"{self.id},{self.rank}"
    def plot(self, bay_length, floor_height):
        plt.scatter(self.x, self.y, c=NODE_COLOR[self.node_type])
        if (self.node_type == "column" or self.node_type == "vertix"):
            plt.text(self.x - bay_length*0.15, self.y, f"{self.label}", c=NODE_COLOR[self.node_type])
        else:
            plt.text(self.x, self.y + floor_height*0.1, f"{self.label}", c=NODE_COLOR[self.node_type])

class element:
    def __init__(self, id, node_i, node_f):
        self.id = id
        self.label = id
        self.nodes_ids = [node_i.id, node_f.id]
        self.xs = [node_i.x, node_f.x]
        self.ys = [node_i.y, node_f.y]
        self.length = np.sqrt((self.xs[1] - self.xs[0])**2 + (self.ys[1] - self.ys[0])**2)
        self.column = ((self.xs[1] - self.xs[0])**2) < 1e-6
    def set_rank(self, rank):
        self.rank = rank
        self.label = f"{self.id},{self.rank}"
    def plot(self):
        plt.plot(self.xs, self.ys, "-k")
        if self.column:
            plt.text(self.xs[0] + 0.15 * self.length, self.ys[0] + 0.5*self.length, f"{self.label}", c="k")
        else:
            plt.text(self.xs[0] + 0.5 * self.length, self.ys[0] - 0.35*self.length, f"{self.label}", c="k")


class frame:
    def __init__(self, num_bays, num_floors, bay_length, floor_height, elements_per_bay, elements_per_floor):
        self.num_bays = num_bays
        self.num_floors = num_floors
        self.bay_length = bay_length
        self.floor_height = floor_height
        self.frame_length = num_bays*bay_length
        self.frame_height = num_floors*floor_height
        self.elements_per_bay = elements_per_bay
        self.elements_per_floor = elements_per_floor
        self.basic_count()
        self.create_mesh()

    def basic_count(self):
        self.nodes_per_column_line = self.num_floors*self.elements_per_floor + 1
        self.nodes_per_full_bay = self.num_floors*(self.elements_per_bay - 1)

        self.nodes_per_column = self.elements_per_floor - 1
        self.nodes_per_bay = self.elements_per_bay - 1

        self.dx = self.bay_length/self.elements_per_bay
        self.dy = self.floor_height/self.elements_per_floor
        self.num_nodes = self.nodes_per_column_line * (self.num_bays + 1) + self.nodes_per_full_bay*(self.num_bays)
        
    # single-call initialisation
    def create_mesh(self):
        self.nodes = self.create_beam_nodes()
        self.nodes.extend(self.create_columns_nodes())
        self.elements = self.create_elements()

    # getters: they get which node number belongs where
    def get_bay_nodes(self, bay, floor):
        bay_nodes = []
        starting_node = 1 + floor*(self.nodes_per_column + 1) + (bay - 1)*self.nodes_per_column_line + (bay - 1)*self.nodes_per_full_bay
        ending_node = starting_node + self.nodes_per_column_line + self.nodes_per_full_bay
        bay_nodes.append(starting_node)
        
        starting_bay_node = bay*self.nodes_per_column_line + (bay - 1)*self.nodes_per_full_bay + (floor - 1)*self.nodes_per_bay 
        for node_i in np.arange(1, self.nodes_per_bay + 1):    
            bay_nodes.append(starting_bay_node + node_i)
        bay_nodes.append(ending_node)
        return bay_nodes
    
    def get_vertix_id(self, column_line, floor):
        return 1 + floor*(self.elements_per_floor) + column_line*self.nodes_per_column_line + column_line*self.nodes_per_full_bay
    
    def get_vertices_ids(self):
        vertices_ids = []
        for floor_i in np.arange(0, self.num_floors + 1):
            for column_line_i in np.arange(0, self.num_bays + 1):
                vertices_ids.append(self.get_vertix_id(column_line_i, floor_i))
        return vertices_ids

    def get_beam_node_ids(self, bay, floor, include_vertices = False):
        if ((floor < 1) or (bay < 1)):
            raise ValueError(f"Floor = {floor} and Bay = {bay}, but neither can be less than 1.")
        starting_id = bay*self.nodes_per_column_line + (bay - 1)*self.nodes_per_full_bay + (floor - 1) * self.nodes_per_bay
        nodes = []
        if include_vertices:
            nodes.append(self.get_vertix_id(bay - 1, floor))

        for node_i in np.arange(1,self.nodes_per_bay + 1):
            nodes.append(starting_id + node_i)

        if include_vertices:
            nodes.append(self.get_vertix_id(bay, floor))
        return nodes

    def get_column_node_ids(self, column_line, floor):
        starting_id = 1 + column_line*self.nodes_per_column_line + floor*self.elements_per_floor + column_line*self.nodes_per_full_bay
        nodes = []
        for node_i in np.arange(1,self.nodes_per_column + 1):
            nodes.append(starting_id + node_i)
        return nodes

    def get_column_line_node_ids(self, column_line):
        column_nodes = []
        for floor_i in np.arange(0, self.num_floors):
            column_nodes.append(self.get_vertix_id(column_line, floor_i))
            column_node_ids = self.get_column_node_ids(column_line, floor_i)
            column_nodes.extend(column_node_ids)
        column_nodes.append(self.get_vertix_id(column_line, self.num_floors))
        return column_nodes
    
    def get_beam_line_node_ids(self, floor, get_vertices = True):
        beam_nodes = []
        for bay_i in np.arange(1, self.num_bays + 1):
            if get_vertices:
                beam_nodes.append(self.get_vertix_id(bay_i - 1, floor))
            bay_node_ids = self.get_beam_node_ids(bay_i, floor)
            beam_nodes.extend(bay_node_ids)
            if get_vertices:
                beam_nodes.append(self.get_vertix_id(bay_i, floor))
        return sorted(list(set(beam_nodes)))
    
    # functions to create nodes and elements
    def create_columns_nodes(self):
        """
        @brief Creates nodes for columns in a frame structure.

        This function generates a list of nodes of type \ref node for each column line in a frame structure.
        It iterates over the range of column lines, calculates the x and y coordinates for each node,
        and appends the nodes to a list. The nodes here include the vertices.

        @return list of node objects for each column line in the frame structure.
        """
        nodes = []
        vertices_ids = self.get_vertices_ids()
        for column_line_i in np.arange(0, self.num_bays + 1):
            node_ids = self.get_column_line_node_ids(column_line_i)
            x_0 = column_line_i*self.bay_length
            y_0 = 0.0
            x = x_0
            y = y_0
            for node_id in node_ids:
                if not(node_id in vertices_ids):
                    nodes.append(node(node_id, x, y, "column"))
                else:
                    nodes.append(node(node_id, x, y, "vertix"))

                y = y + self.dy
        return nodes
    
    def create_beam_nodes(self):
        """
        @brief Creates beam nodes for each floor in the structure.
        This function generates beam nodes for each floor in the structure by iterating through 
        the floors and calculating the x and y coordinates for each node. It checks if the node 
        is not a vertex before appending it to the nodes list, so no vertices included.
        @return List of beam nodes with their respective IDs and coordinates.
        """
        nodes = []
        vertices = self.get_vertices_ids()
        for floor_i in np.arange(1, self.num_floors + 1):
            node_ids = self.get_beam_line_node_ids(floor_i, get_vertices=True)
            x_0 = 0.0
            y_0 = floor_i*self.floor_height
            x = x_0
            y = y_0
            for node_id in node_ids:
                if not(node_id in vertices):
                    nodes.append(node(node_id, x, y, "beam"))
                x = x + self.dx
        return nodes
    
    # element functions 
    def map_elements(self):
        self.elements_map = []
        element_id = 0

        for bay_i in np.arange(1, self.num_bays + 1):
            # first column line
            node_ids = self.get_column_line_node_ids(bay_i - 1)
            for node_id_i in np.arange(0, len(node_ids) - 1):
                element_id = element_id + 1
                self.elements_map.append((element_id, [node_ids[node_id_i], node_ids[node_id_i + 1]]))
            # beams on each floor
            for floor_i in np.arange(1, self.num_floors + 1):
                node_ids = self.get_beam_node_ids(bay_i, floor_i, include_vertices=True)
                for node_id_i in np.arange(0, len(node_ids) - 1):
                    element_id = element_id + 1    
                    self.elements_map.append((element_id, [node_ids[node_id_i], node_ids[node_id_i + 1]]))
        # last column line
        node_ids = self.get_column_line_node_ids(self.num_bays)
        for node_id_i in np.arange(0, len(node_ids) - 1):
            element_id = element_id + 1    
            self.elements_map.append((element_id, [node_ids[node_id_i], node_ids[node_id_i + 1]]))

    def create_elements(self):
        self.map_elements()
        elements = []
        for elem_node_pair in self.elements_map:
            node_i = find_node_by_id(self.nodes, elem_node_pair[1][0])
            node_f = find_node_by_id(self.nodes, elem_node_pair[1][1])
            elements.append(element(elem_node_pair[0], node_i, node_f))
        return elements  


    # plotting 
    def plot_mesh(self):
        if (not(self.nodes) or not(self.elements)):
            ValueError("self.nodes and/or self.elements is not initialised.")
        for node_ in self.nodes:
            node_.plot(self.bay_length, self.floor_height)
        for elem in self.elements:
            elem.plot()
        plt.gcf().set_size_inches(self.num_bays*2, self.num_floors*2)
        plt.tight_layout()