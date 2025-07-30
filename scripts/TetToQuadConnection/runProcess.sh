#!/bin/zsh
#Run this file from the directory where it is located, the parent directory.
#The input tet mesh must be named: “inputTetMesh.obj” and placed in the “input” directory 
#   which resides in the TetToQuadConnection directory.
set -e
echo "Starting script"
echo "Creating quadrilateral mesh..."
cd build
if ! ./quadriflow ../input/inputTetMesh.obj ./output/quad.vtk setting.config; then
    echo "Quadriflow failed after mesh generation, continuing with mesh generated..."
fi
echo "Copying output to input directory..."
cp -i ./tempdir/quad.vtk ../input
cd ../src
echo "Converting vtk file to obj..."
echo "Connecting the tetrahedral mesh to the quadrilateral mesh..."
python3 main.py ../input/inputTetMesh.obj ../input/convertedVTK.obj
echo "output file saved to: TetToQuadConnection/output/quadVertices.txt"


#cleanup?
#delete quad.vtk from input and the tempdirectory? I could unless people want to see
# the midpoints. It isnt harming much they should simply be overwritten each time. As long as i am rewriting files never appending.
