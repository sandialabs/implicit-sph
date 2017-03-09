#!/bin/bash


rm -rf CMakeCache.txt CMakeFiles

cmake \
    -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
    -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
    -D Trilinos_ENABLE_Fortran:BOOL=OFF \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
    -D Trilinos_ENABLE_TESTS:BOOL=ON \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Epetra_ENABLE_TESTS:BOOL=ON \
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_NOX:BOOL=ON \
    -D NOX_ENABLE_TESTS:BOOL=ON \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Stratimikos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_Teko:BOOL=ON \
    -D Trilinos_ENABLE_Zoltan:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_BLAS:STRING=ON \
    -D TPL_BLAS_LIBRARIES="-L/usr/lib64 -lblas -llapack -lgfortran" \
    \
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    \
    -D CMAKE_INSTALL_PREFIX:PATH=/home/kyukim/Work/lib/trilinos/implicit-sph/sems \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
    -D CMAKE_CXX_FLAGS:STRING="-g -O3" \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    \
    /home/kyukim/Work/lib/trilinos/master
