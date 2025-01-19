if [ -z $1 ]; then
    cd bin
    ./Blaze "$@"
elif [ $1 = "mpi" ]; then
    shift
    cd bin
    if [ -z $2 ]; then
        mpirun -n 1 ./Blaze 
    else
        mpirun -n "$2" ./Blaze 
    fi
fi