gfortran -ffree-line-length-none -o StructureSolver StructureSolver.f90
gfortran -ffree-line-length-none -o BeamStructure BeamStrucutre.f90
cd examples
for frqcase in */; do
  if [ -d "${frqcase}" ]; then
    frqcase="${frqcase%/}"
    echo "Start checking ${frqcase} case ..."
    start=$(date +%s)
    cd "${frqcase}"
    rm -rf test
    mkdir test
    cp -r ./*.dat test/
    cd test

    ./../../../BeamStructure > /dev/null 2>&1
    ./../../../StructureSolver > /dev/null 2>&1

    beam_result=$(diff -q disp_ele_2_BeamStructure.dat old_disp_ele_2_BeamStructure.dat || true)
    structure_result=$(diff -q disp_ele_2_StructureSolver.dat old_disp_ele_2_StructureSolver.dat || true)

    end=$(date +%s)
    echo -n "Time used $(($end-$start)) seconds. "
    if [ -z "${beam_result}" ] && [ -z "${structure_result}" ]; then
      echo "Passed."
    else
      echo "Failed."
      if [ -n "${beam_result}" ]; then
        echo "BeamStructure result differs from old_disp_ele_2_BeamStructure.dat"
      fi
      if [ -n "${structure_result}" ]; then
        echo "StructureSolver result differs from old_disp_ele_2_StructureSolver.dat"
      fi
    fi

    cd ..
    rm -rf test
    cd ..
  fi
done
