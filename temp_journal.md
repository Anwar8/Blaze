# Temporary Journal - building Blaze on Precision 3680

## Updating `cmake`
Building with `Eigen3` version 5.0 and above [requires `cmake` version 3.5 or higher](https://libeigen.gitlab.io/eigen/docs-nightly/TopicCMakeGuide.html), while the default on this workstation is 3.28. A new version install was thus required. 
This is done via the command: 
```bash
sudo snap install cmake --classic
```
which has been retrieved from the [snap store](https://snapcraft.io/cmake), and installed version 4.2.1 of `cmake` in `/snap/bin/`.
This installation, however, does not override the original installation performed by `apt`, which can be located in `/usr/bin/`. To use the correct version of `cmake`, the path needs to be modified via the command:
```bash
export PATH=/snap/bin/:$PATH
```
which is placed in `~/.bashrc`. Additionally, to use the *correct* compilers for `C` and `C++`, another command needs to be added to `.bashrc`: 
```bash
export CC=icx
export CXX=icpx
```
noting that `icc` does not exist in the latest version of `oneAPI`, which I am assuming is installed on the system already. 

## Installing latest version of `eigen3`
A [new version of `eigen3` has been released](https://gitlab.com/libeigen/eigen/-/releases) recently going from 3.4 to 5.0. It is best to use this version when working with `Blaze` to ensure best performance possible, and to detect and correct any compatiblity errors early on. Unlike what the [documentation states](https://libeigen.gitlab.io/eigen/docs-5.0/GettingStarted.html), it *is* necessary to run `cmake` on `eigen` to get it properly installed and usable with `Blaze`. 
The following set of commands, executed in the directory where the comopressed release is located, will install `eigen3` in `/usr/local/include` using 24 cores:
```bash
tar -xvf eigen-5.0.1.tar
cd eigen-5.0.1
mkdir build
cd build
cmake .. 
sudo make install -j24
```
There will not be a need to tell `cmake` where to find `eigen3` when building `Blaze` as this installation process will add `eigen3` to the path, and will make itself discoverable by `cmake`'s `find_package`. 

## Installing `Trilinos`
For a machine with `oneAPI` and its `MKL`, the following command will install `Trilinos` with `Tpetra`, `Amesos2`, `Belos`, and all their dependencies. It will also enable the necessary links between `Belos` and `Tpetra`.
```bash
cmake -DTPL_ENABLE_MPI=ON -DTPL_ENABLE_MKL=ON -DMKL_LIBRARY_DIRS=$MKLROOT/lib/intel64 -DMKL_INCLUDE_DIRS=$MKLROOT/include -DBLAS_LIBRARY_DIRS=$MKLROOT/lib/intel64 -DLAPACK_LIBRARY_DIRS=$MKLROOT/lib/intel64 -DBLAS_LIBRARY_NAMES="mkl_intel_lp64;mkl_sequential;mkl_core" -DLAPACK_LIBRARY_NAMES="mkl_intel_lp64;mkl_sequential;mkl_core" -DTrilinos_ENABLE_Tpetra=ON -DTrilinos_ENABLE_Amesos2=ON -DTrilinos_ENABLE_Belos=ON -DBelos_ENABLE_Tpetra=ON -DBelos_ENABLE_KokkosKernels=ON -DCMAKE_CXX_COMPILER=/opt/intel/oneapi/mpi/2021.14/bin/mpicxx  -DCMAKE_INSTALL_PREFIX=../trilinos-install ..
```

I believe, however, that the flags `-DBelos_ENABLE_Tpetra=ON -DBelos_ENABLE_KokkosKernels=ON` are unnecessary as they should be turned on automatically as `Tpetra`, `KokkosKernels` and `Belos` are being built. 

There are other ways to direct `Trilinos` to `MKL`. Check journal entry for 23 July, which uses: 

```bash
-DTPL_ENABLE_MKL=ON -DTPL_MKL_INCLUDE_DIRS="${MKLROOT}/include" -DTPL_MKL_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" -DTPL_ENABLE_BLAS=ON -DTPL_ENABLE_LAPACK=ON -DTPL_BLAS_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" -DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" 
``` 
which avoids having to name the libraries explicitly using `-DBLAS_LIBRARY_NAMES="mkl_intel_lp64;mkl_sequential;mkl_core" -DLAPACK_LIBRARY_NAMES="mkl_intel_lp64;mkl_sequential;mkl_core"`.

Using `-DCMAKE_CXX_COMPILER=/opt/intel/oneapi/mpi/2021.14/bin/mpicxx` ensures we use the `MPI` wrapper of the `oneAPI` compiler. 

