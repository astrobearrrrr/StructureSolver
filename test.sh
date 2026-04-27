#!/usr/bin/env bash
set -euo pipefail

bash make.sh

compare_fieldstat() {
  local result_file=$1
  local baseline_file=$2
  local solver_name=$3

  if [ ! -f "${baseline_file}" ]; then
    echo "${solver_name} baseline missing: ${baseline_file}"
    return 2
  fi

  if diff -q "${result_file}" "${baseline_file}" >/dev/null 2>&1; then
    return 0
  fi

  echo "${solver_name} result differs from ${baseline_file}"
  return 1
}

cd examples
for testcase in */; do
  [ -d "${testcase}" ] || continue
  testcase="${testcase%/}"

  echo "Start checking ${testcase} case ..."
  start=$(date +%s)

  cd "${testcase}"
  rm -rf test
  mkdir test
  cp -r ./*.dat test/
  cd test

  ./../../../BeamStructure > /dev/null 2>&1
  ./../../../StructureSolver > /dev/null 2>&1

  beam_status=0
  structure_status=0

  compare_fieldstat fieldstat_BeamStructure.dat old_fieldstat_BeamStructure.dat BeamStructure || beam_status=$?
  compare_fieldstat fieldstat_StructureSolver.dat old_fieldstat_StructureSolver.dat StructureSolver || structure_status=$?

  end=$(date +%s)
  echo -n "Time used $(($end-$start)) seconds. "

  if [ "${beam_status}" -eq 0 ] && [ "${structure_status}" -eq 0 ]; then
    echo "Passed."
  else
    echo "Failed."
  fi

  cd ..
  rm -rf test
  cd ..
done
