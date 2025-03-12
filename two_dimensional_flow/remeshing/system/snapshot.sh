#!/bin/bash
if [ -z "$2" ]
then
  echo "Invalid arguments"
  exit
fi

start=$(date +%s)

fidelity=$1
angle_id=$2
n_discretisation=$((2 ** ("${fidelity}" + 5)))
angle=$(../../system/angle "${angle_id}")
name="${1}_${2}"
path="../data/${name}"

if [ -d "../data/${path}" ]
then
  echo "Snapshot = ${path} already exists. Aborting."
  exit
fi

cp -r template "$path"
cd "${path}/system" || exit
python3 generateblockMeshDict.py "$n_discretisation"
cd ../constant/triSurface || exit
python3 rotate.py "$angle"
cd ../..
surfaceFeatureExtract > log.surfaceFeatureExtract
blockMesh > log.blockMesh
snappyHexMesh -overwrite > log.snappyHexMesh
extrudeMesh > log.extrudeMesh
renumberMesh -overwrite > log.renumberMesh
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

