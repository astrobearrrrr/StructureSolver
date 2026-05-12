#!/usr/bin/env bash
set -euo pipefail

bash make.sh

compare_fieldstat() {
  local result_file=$1
  local baseline_file=$2
  local solver_name=$3

  if [ ! -f "${baseline_file}" ]; then
    echo "${solver_name} no baseline file: ${baseline_file}"
    return 2
  fi

  if [ ! -f "${result_file}" ]; then
    echo "${solver_name} result missing: ${result_file}"
    return 1
  fi

  if awk '
    NR==FNR {
      n = n + 1
      base[n] = $NF
      next
    }
    {
      m = m + 1
      res[m] = $NF
    }
    END {
      if (n != m) exit 1
      for (i = 1; i <= n; i++) {
        if (substr(base[i], 1, 5) != substr(res[i], 1, 5)) exit 1
      }
    }
  ' "${baseline_file}" "${result_file}"; then
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

  beam_status=0
  structure_status=0

  if [ -f old_fieldstat_BeamStructure.dat ]; then
    ./../../../BeamStructure > /dev/null 2>&1
  fi
  compare_fieldstat fieldstat_BeamStructure.dat old_fieldstat_BeamStructure.dat BeamStructure || beam_status=$?

  if [ -f old_fieldstat_StructureSolver.dat ]; then
    ./../../../StructureSolver > /dev/null 2>&1
  fi
  compare_fieldstat fieldstat_StructureSolver.dat old_fieldstat_StructureSolver.dat StructureSolver || structure_status=$?

  end=$(date +%s)
  echo -n "Time used $(($end-$start)) seconds. "

  if [ "${beam_status}" -ne 1 ] && [ "${structure_status}" -ne 1 ]; then
    echo "Passed."
  else
    echo "Failed."
  fi

  cd ..
  rm -rf test
  cd ..
done
