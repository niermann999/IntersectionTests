#!/bin/bash

clean_up_dir () {
  echo "Directory exists"
  echo -n "Remove directory ${1} (y/n)? "
  read answer

  if [ "$answer" != "${answer#[Yy]}" ] ;then
    rm -r $1
  elif [ "$answer" != "${answer#[Nn]}" ] ;then
    echo "Please clean directory manually!"
    exit 1
  else
    echo "Please enter (y/n)"
    exit 1
  fi
}


if [[ $# > 0 ]]; then
  while getopts ":r:" opt; do
    case $opt in
      r) rm_build_dir="$OPTARG"
      ;;
    esac
  done
fi

# Set output dir for this machine
name=$(hostname)
cpu=$(cat /proc/cpuinfo| grep 'model name' | head -n 1 | awk '{printf"%s_%s_%s", $4, $6, $9}')

name="${name//./_}"
cpu="${cpu//./_}"
cpu="${cpu//-/_}"
cpu="${cpu//(R)}"
out_dir="../visualization/${name}_${cpu}"

# Cleanup of previous runs
echo "Ouput dir: ${out_dir}"
[ -d $out_dir ] && clean_up_dir $out_dir
mkdir $out_dir 

# Remember Hardware details
cat /proc/cpuinfo >> ../visualization/$out_dir/cpuinfo


# Clean build directory. You are on your own
run_bench () {
  if [ "$rm_build_dir" = "yes" ]; then
    ./run_benchmarks.sh -r yes -i $1 -p $2
  # Be more polite
  else
    ./run_benchmarks.sh -i $1 -p $2
  fi

  # Plot results
  #python_init="../visualization/__init__.py"
  #[! -f $python_init ] && touch $python_init
  #python ../visualization/plot_benchmarks.py

  # Ouput dir has been cleaned before
  mkdir $3

  #mv ../visualization/aggregate_plots ../visualization/aggregate_plots $out_dir/Optimal_float
  mv ../visualization/*.out $3
  mv ../visualization/*.json $3
}

# Optimal float
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "" "" "${out_dir}/Optimal_float"

# Optimal double
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "" "double" "${out_dir}/Optimal_double"

# AVX float
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "avx" "" "${out_dir}/Avx_float"

# AVX double
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "avx" "double" "${out_dir}/Avx_double"

# SSE float
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "sse" "" "${out_dir}/Sse_float"

# SSE float
printf "\nRun benchmark for optimal compilation with single precision...\n"
run_bench "sse" "double" "${out_dir}/Sse_double"