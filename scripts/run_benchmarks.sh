#!/bin/bash

clean_up_dir () {
  echo "Directory ${1} exists!"
  echo -n "Remove contents of directory ${1} (y/n)? "
  read answer

  if [ "$answer" != "${answer#[Yy]}" ] ;then
    rm -r $1/*
  elif [ "$answer" != "${answer#[Nn]}" ] ;then
    echo "Directory will not be cleaned"
  else
    echo "Please enter (y/n)"
    exit 1
  fi
}


if [[ $# > 0 ]]; then
  while getopts ":i:p:b:r:" opt; do
    case $opt in
      i) instr_set="$OPTARG"
      ;;
      p) precision="$OPTARG"
      ;;
      b) build_dir="$OPTARG"
      ;;
      r) rm_build_dir="$OPTARG"
      ;;
    esac
  done
fi

if [ "$build_dir" = "" ]; then
  build_dir="../build"
  printf "Build dir %s\n" "$build_dir"
fi

# If no instruction set is specified, compilation 
# will use Vc architecture flags (in cmake)
if [ "$instr_set" = "sse" ]; then
  INSTRUCTION_SET="-DCOMPILE_SSE=ON"
  instr_set="SSE_"
  echo "Compiling for SSE"
elif [ "$instr_set" = "avx" ]; then
  INSTRUCTION_SET="-DCOMPILE_AVX=ON"
  instr_set="AVX_"
  echo "Compiling for AVX"
else
  echo "Compile with Vc Architecture flags"
  instr_set="Optimal_"
fi

# If precision is not set to double, it is already 
# automatically (in cmake)
if [ "$precision" = "double" ]; then
  COMPILE_OPTIONS="-DPRECISION_DOUBLE"
  echo "Precision double"
  precision="Double_"
else
  echo "Precision float"
  precision="Float_"
fi

printf "\n====================\nRun Configuration...\n====================\n\n"

# Clean build directory. You are on your own
if [ "$rm_build_dir" = "yes" ]; then
  rm -r $build_dir/*
# Be more polite
else
  clean_up_dir $build_dir
fi

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="g++" -DCMAKE_CXX_FLAGS=$COMPILE_OPTIONS $INSTRUCTION_SET -B $build_dir -S ..

echo "Compile command: cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=\"g++\" -DCMAKE_CXX_FLAGS=${COMPILE_OPTIONS} ${INSTRUCTION_SET} -B ${build_dir} -S .." >> ../visualization/${instr_set}${precision}ComplileCmd.out

printf "\n================\nBuild Project...\n================\n\n"

make -C $build_dir -j4

printf "\n=================\nRun Benchmarks...\n=================\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcHorizIntersection --benchmark_out_format=json --benchmark_out=../visualization/VcHoriz_bench.json 2> ../visualization/${instr_set}${precision}VcHoriz.out

printf "\n=========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcHorizIntersection_wres --benchmark_out_format=json --benchmark_out=../visualization/VcHoriz_wres_bench.json 2> ../visualization/${instr_set}${precision}VcHoriz_wres.out

printf "\n=========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcHybridIntersection --benchmark_out_format=json --benchmark_out=../visualization/VcHybrid_bench.json 2> ../visualization/${instr_set}${precision}VcHybrid.out

printf "\n========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcHybridIntersection_wres --benchmark_out_format=json --benchmark_out=../visualization/VcHybrid_wres_bench.json 2> ../visualization/${instr_set}${precision}VcHybrid_wres.out

printf "\n=========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcVertIntersection --benchmark_out_format=json --benchmark_out=../visualization/VcVert_bench.json 2> ../visualization/${instr_set}${precision}VcVert.out

printf "\n=========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkVcVertIntersection_wres --benchmark_out_format=json --benchmark_out=../visualization/VcVert_wres_bench.json 2> ../visualization/${instr_set}${precision}VcVert_wres.out

printf "\n=========\nnDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkEigenIntersection --benchmark_out_format=json --benchmark_out=../visualization/Eigen_bench.json 2> ../visualization/${instr_set}${precision}Eigen.out

printf "\n=========\nDone.\nSleeping to reduce system load for next benchmark...\n=========\n\n"
#sleep 5m
sleep 10

../build/src/BenchmarkEigenIntersection_wres --benchmark_out_format=json --benchmark_out=../visualization/Eigen_wres_bench.json 2> ../visualization/${instr_set}${precision}Eigen_wres.out

printf "\n=========\nDone.\n=========\n"