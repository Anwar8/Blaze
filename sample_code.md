# Using `gmsh` to assemble the global stiffness matrix
```
// Retrieve nodal coordinates
std::vector<double> coord;
std::vectorstd::size_t nodeTags;
gmsh::model::mesh::getNodes(nodeTags, coord);

// Retrieve element connectivity
std::vectorstd::size_t elTags;
std::vectorstd::size_t elNodes;
std::vectorstd::size_t elTypes;
gmsh::model::mesh::getElementsByType(3, elTypes, elTags, elNodes);

// Define global stiffness matrix
std::size_t numNodes = coord.size() / 3;
std::size_t numElems = elTags.size();
std::size_t numDofsPerNode = 3;
std::size_t numDofs = numNodes * numDofsPerNode;

Eigen::SparseMatrix<double> K(numDofs, numDofs);
K.reserve(numElems * 27); // Assuming 27 non-zero entries per element

// Compute element stiffness matrix
Eigen::Matrix<double, 9, 9> Ke;
// ... call your function to compute the element stiffness matrix ...

// Assemble element stiffness matrix into global stiffness matrix
for (std::size_t i = 0; i < elNodes[elTypes[el]]; ++i) {
    std::size_t node_i = elNodes[el * elNodes[elTypes[el]] + i] - 1;
    for (std::size_t j = 0; j < elNodes[elTypes[el]]; ++j) {
        std::size_t node_j = elNodes[el * elNodes[elTypes[el]] + j] - 1;
        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t l = 0; l < 3; ++l) {
                K.coeffRef(node_i * 3 + k, node_j * 3 + l) += Ke(i * 3 + k, j * 3 + l);
            }
        }
    }
}
```