#!/bin/bash
if [ -z "$2" ]
then
  echo "Invalid arguments"
  exit
fi

fidelity=$1
angle_id=$2
n_discretisation=$((2 ** ("${fidelity}" + 5)))
angle=$(../../system/angle "${angle_id}")
name="${1}_${2}"
path="../data/${name}"

if [ ! -d "template${n_discretisation}" ]
then
  echo "Generating mesh..."
  cp -r template "template${n_discretisation}"
  cd "template${n_discretisation}/system" || exit
  python3 generateblockMeshDict.py "$n_discretisation"
  cd ..
  surfaceFeatureExtract > log.surfaceFeatureExtract
  blockMesh > log.blockMesh
  snappyHexMesh -overwrite > log.snappyHexMesh
  extrudeMesh > log.extrudeMesh
  renumberMesh -overwrite > log.renumberMesh
  touch template_done
  echo "Mesh generated."
  cd ..
fi

# For if another process is generating the template in parallel
while [ ! -f "template${n_discretisation}/template_done" ]
do
  true
done

start=$(date +%s)

if [ -d "${path}" ]
then
    echo "Snapshot ${path} already exists. Aborting."
    exit
fi

cp -r "template${n_discretisation}" "$path"
cd "${path}" || exit
../../../system/mesh_displacement/mesh_displacement -a "${angle}" > log.meshDisplacement
cp 0/polyMesh/points constant/polyMesh/points
checkMesh > log.checkMesh
simpleFoam > log.simpleFoam

cd ../../system || exit

qoi=$(../../system/qoi "$fidelity" "$angle_id")
end=$(date +%s)

cost=$((end - start))

../../system/write_snapshot_data "$fidelity" "$angle_id" "$qoi" "$cost"
../../system/collect_data ../

cd "$path" || exit

echo "${qoi} $cost" >> "data_ascii"

