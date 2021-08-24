#! /bin/bash
#
set -ex

FC=/usr/lib64/mpich/bin/mpif90

${FC} -c -Wall hello_mpi.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
${FC} hello_mpi.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hello_mpi.o
cp a.out $HOME/local/bin/hello_mpi
cp a.out hello_mpi
#
echo "Normal end of execution."

