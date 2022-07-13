use CTypes;

use octree;

require "index_space.h";
require "grgeometry.h";
require "input_database.h";
extern type Point = 3*int(32);

extern record Box {
    var lo: Point;
    var up: Point;
};

extern record BoxArray {
  var boxes : c_ptr(Box);
  var boxlimits : c_ptr(c_int);
  var size : c_uint;
}

extern record GrGeomSolid {
    var data: c_ptr(GrGeomOctree);
    var patches: c_ptr(c_ptr(GrGeomOctree));
    var octree_bg_level, octree_ix, octree_iy, octree_iz: int;
    
    var interior_boxes: c_ptr(BoxArray);
    var surface_boxes : c_ptr(c_ptr(BoxArray));
    var patch_boxes : c_ptr(c_ptr(c_ptr(BoxArray)));
    
}


iter GrGeomInLoopBoxes_iter(ref grgeom: GrGeomSolid, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {

        var boxes = grgeom.interior_boxes[0];
        //writeln("boxes:", boxes);
        if boxes.size == 0 then return;
        for i in 0..<boxes.size {
            var box = boxes.boxes[i];
            var xlo = max(ix, box.lo[0]);
            var ylo = max(iy, box.lo[1]);
            var zlo = max(iz, box.lo[2]);
            var xhi = min(ix + nx -1, box.up[0]);
            var yhi = min(iy + ny -1, box.up[1]);
            var zhi = min(iz + nz -1, box.up[2]);
            //writeln(xlo, " ", xhi, " ", ylo, " ", yhi, " ", zlo," ", zhi);
            yield(xlo,xhi,ylo,yhi,zlo,zhi);
        }
    //writeln("done yielding");
}

iter GrGeomSurfLoopBoxes_iter(ref grgeom: GrGeomSolid, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    for f in 0..<GrGeomOctreeNumFaces {
        var fdir: [0..2] int = create_fdir(f);
        var boxes = grgeom.surface_boxes[f][0];
        //writeln("boxes:", boxes);
        if boxes.size == 0 then continue;
        for i in 0..<boxes.size {
            var box = boxes.boxes[i];
            var xlo = max(ix, box.lo[0]);
            var ylo = max(iy, box.lo[1]);
            var zlo = max(iz, box.lo[2]);
            var xhi = min(ix + nx -1, box.up[0]);
            var yhi = min(iy + ny -1, box.up[1]);
            var zhi = min(iz + nz -1, box.up[2]);
            //writeln(xlo, " ", xhi, " ", ylo, " ", yhi, " ", zlo," ", zhi);
            yield(xlo,xhi,ylo,yhi,zlo,zhi,fdir);
        }
    }
    //writeln("done yielding");
}


iter GrGeomInLoop_iter(ref grgeom: GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    if(r == 0 && grgeom.interior_boxes != nil) {

        for space in GrGeomInLoopBoxes_iter(grgeom,ix,iy,iz,nx,ny,nz) {
            yield space;
        }
    }
    else {
        
        var offset = 2 ** r;
        var i = grgeom.octree_ix * offset;
        var j = grgeom.octree_iy * offset;
        var k = grgeom.octree_iz * offset;
        for space in GrGeomOctreeInteriorLoop(i,j,k, grgeom.data[0], r + grgeom.octree_bg_level, ix, iy, iz, nx, ny, nz) {
            yield space;
        }
    }
}

iter GrGeomSurfLoop_iter(ref grgeom: GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    if(r == 0 && grgeom.surface_boxes[5] != nil) {
        //writeln("boxes.");
        for space in GrGeomSurfLoopBoxes_iter(grgeom,ix,iy,iz,nx,ny,nz) {
            yield space;
        }
    }
    else {
        //writeln("octree");
        var offset = 2 ** r;
        var i = grgeom.octree_ix * offset;
        var j = grgeom.octree_iy * offset;
        var k = grgeom.octree_iz * offset;
        for space in GrGeomOctreeFaceLoop_iter(i,j,k, grgeom.data[0], r + grgeom.octree_bg_level, ix, iy, iz, nx, ny, nz) {
            yield space;
        }
    }
} 

export proc GrGeomSurfLoop_chapel(ref grgeom: GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) { 
    //writeln("About to iterate.");
    for space in GrGeomSurfLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        writeln(space);
    }
}


iter GrGeomPatchLoopBoxesNoFdir(ref grgeom: GrGeomSolid, patch_num: int, ovrlnd: int,  r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    for f in 0..<GrGeomOctreeNumFaces {
        var boxes = grgeom.patch_boxes[f][patch_num][0];
        for box_index in 0..<boxes.size {
            var box = boxes.boxes[box_index];
            var xlo = max(ix, box.lo[0]);
            var ylo = max(iy, box.lo[1]);
            var zlo = max(iz, box.lo[2]);
            var xhi = min(ix + nx -1, box.up[0]);
            var yhi = min(iy + ny -1, box.up[1]);
            var zhi = min(iz + nz -1, box.up[2]);
            //writeln(xlo, " ", xhi, " ", ylo, " ", yhi, " ", zlo," ", zhi);
            yield(xlo,xhi,ylo,yhi,zlo,zhi,f);
        }
    }
}


iter GrGeomPatchLoopNoFdir_iter(ref grgeom: GrGeomSolid, patch_num: int, ovrlnd: int,  r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    if(r == 0 && grgeom.patch_boxes[5][patch_num] != nil) {
        for space in GrGeomPatchLoopBoxesNoFdir(grgeom, patch_num, ovrlnd, r, ix, iy, iz, nx, ny, nz) {
            yield space;
        }
    } else {
        var offset = 2 ** r;
        var i = grgeom.octree_ix * offset;
        var j = grgeom.octree_iy * offset;
        var k = grgeom.octree_iz * offset; 
        for space in GrGeomOctreeFaceLoopNoFdir_iter(i,j,k, grgeom.patches[patch_num][0], r + grgeom.octree_bg_level, ix, iy, iz, nx, ny, nz) {
            yield space;
        }   
    }
}