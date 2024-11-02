# `Blaze` Libraries

`Blaze` currenlty consists of 6 libraries, each of which has its own directory inside `source`:
1. `core`: Contains the base functionality of `Blaze` including helper utility functions, and essential mathematics defaults. Also includes the file that implements the `Node` class. Depends on `Eigen3`. 
2. `materials`: Contains the material and different sections classes. Depends on `core`.
3. `elements`: Contains the interface for generic element classes, helper classes for elements such as those that deal with orientation and direction, and the classes for beam-column elements. Depends on `materials`.
4. `aggregators`: Definitions for the classes that *aggregate* others such `GlobalMesh` and `Assembler`. Depends on `elements`.
5. `managers`: Managers manage the various components of `Blaze` such as load and history, and this library contains the classes responsible for this. Depends on `aggregators`.
6. `solution`: Classes that perform the solution procedure are not `aggregators` nor `managers`, and so are placed here in their own library. Depends on `managers`.
Each of these libraries depends on the next library, and all are linked to an over-arching library called `model` which also implements the `Model` class. These libraries can all be treated as `INTERFACE` libraries as each only consists of header files. However, empty `.cpp` files were also added in order to allow for building `STATIC` versions of these libraries that are linked to a `main` `Blaze` function, or to test executable. See `CMakeLists.txt` in the parent directory of `Blaze` to see available build options.