mol new <filename>.psf
mol addfile <binaryfilename>.coor
set sel [atomselect top all]
$sel writepdb <finalframe>.pdb

