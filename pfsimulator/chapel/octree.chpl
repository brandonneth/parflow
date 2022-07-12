use CTypes;
require "grgeom_octree.h";

extern record GrGeomOctree {
  var flags : c_uchar;
  var faces : c_uchar;
  var parent : c_ptr(GrGeomOctree);
  var children : c_ptr(c_ptr(GrGeomOctree));
}

const GrGeomOctreeNumFaces = 6;
const GrGeomOctreeNumChildren = 8;

const GrGeomOctreeNodeEmpty: uint(8) = 1;
const GrGeomOctreeNodeOutside: uint(8) = 2;
const GrGeomOctreeNodeInside: uint(8) = 4;
const GrGeomOctreeNodeFull: uint(8) = 8;
const GrGeomOctreeNodeLeaf: uint(8) = 16;



const GrGeomOctreeFaceL = 0;
const GrGeomOctreeFaceR = 1;
const GrGeomOctreeFaceD = 2;
const GrGeomOctreeFaceU = 3;
const GrGeomOctreeFaceB = 4;
const GrGeomOctreeFaceF = 5;

//corresponds to GrGeomOctreeFaceValue
proc faceValue(face) {
    return 1 << face:uint;
}

enum GrGeomOctreeOctant {LDB, RDB, LUB, RUB, LDF, RDF,
                         LUF, RUF};

proc flagIs(node: GrGeomOctree, flag) {
    return (node.flags & flag) == flag;
}

proc isInside(node: GrGeomOctree) {
    return flagIs(node, GrGeomOctreeNodeInside);
}

proc isEmpty(node: GrGeomOctree) {
    return flagIs(node, GrGeomOctreeNodeEmpty);
}

proc isOutside(node: GrGeomOctree) {
    return flagIs(node, GrGeomOctreeNodeOutside);
}

proc isLeaf(node: GrGeomOctree) {
    return flagIs(node, GrGeomOctreeNodeLeaf);
}

proc isFull(node: GrGeomOctree) {
    return flagIs(node, GrGeomOctreeNodeFull);
}

proc GrGeomOctreeHasFace(node:GrGeomOctree, face_index) {
    return (node.faces & faceValue(face_index)) != 0 ;
}

proc create_fdir(f: int) {
    var fdir: [0..2] int = 0;
    select f {
        when GrGeomOctreeFaceL do fdir[0] = -1;
        when GrGeomOctreeFaceR do fdir[0] = 1;
        when GrGeomOctreeFaceD do fdir[1] = -1;
        when GrGeomOctreeFaceU do fdir[1] = 1;
        when GrGeomOctreeFaceB do fdir[2] = -1;
        when GrGeomOctreeFaceF do fdir[2] = 1;
        otherwise do fdir[0] = -999;
    }
    return fdir;
}

iter GrGeomOctreeInsideLoop(in i: int, in j: int, in k: int, root: GrGeomOctree, level_of_interest: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    var l = 0;
    var increment = 1 << level_of_interest;
    var node = root;

    var visit_child = false;
    var PV_visiting: [-1..level_of_interest + 1] int = 0;

    while (l >= 0) {
        if (l == level_of_interest) {
            if(isInside(node)) {
                if((i >= ix) && (i < (ix + nx)) && (j >= iy) && (j < (iy + ny)) && (k >= iz) && (k < (iz + nz))) {
                    for f in 0..<GrGeomOctreeNumFaces {
                        if (GrGeomOctreeHasFace(node,f)) {
                            yield(i,i,j,j,k,k,f);
                        }
                    }
                }
            }
            visit_child = false;
        } else if(isLeaf(node)){
            if (isInside(node)) {
                var xlo = max(ix,i);
                var ylo = max(iy,j);
                var zlo = max(iz,k);
                var xhi = min(ix+nx,i+increment);
                var yhi = min(iy+ny,j+increment);
                var zhi = min(iz+nz,k+increment);
                for f in 0..<GrGeomOctreeNumFaces {
                        if (GrGeomOctreeHasFace(node,f)) {
                            //var fdir: [0..2] int = create_fdir(f);
                            //the C code does < checks rather than <= checks.
                            yield(xlo,xhi-1,ylo,yhi-1,zlo,zhi-1,f);
                        }
                    }
            }
            visit_child = false;
        } else if (!isInside(node)) {
            visit_child = false;
        } else if (PV_visiting[l] < GrGeomOctreeNumChildren) {
            visit_child = true;
        } else {visit_child = false;}

        if (visit_child) {
            node = node.children[PV_visiting[l]][0];
            increment = increment >> 1;
            if (PV_visiting[l] & 1) {
                i += increment;
            }
            if (PV_visiting[l] & 2) {
                j += increment;
            }
            if (PV_visiting[l] & 4) {
                k += increment;
            }
            l += 1;
            PV_visiting[l] = 0;
        } else {
            l -= 1;
            if (PV_visiting[l] & 1) {
                i -= increment;
            }
            if (PV_visiting[l] & 2) {
                j -= increment;
            }
            if (PV_visiting[l] & 4) {
                k -= increment;
            }
            increment = increment << 1;
            if(node.parent != nil) {
                node = node.parent[0];
            }
            PV_visiting[l] += 1;
        }
    }
}

iter GrGeomOctreeInteriorLoop(in i: int, in j: int, in k: int, root: GrGeomOctree, level_of_interest: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    var l = 0;
    var increment = 1 << level_of_interest;
    var node = root;

    var visit_child = false;
    var PV_visiting: [-1..level_of_interest + 1] int = 0;

    while (l >= 0) {
        if (l == level_of_interest) {
            if(isInside(node) || isFull(node)) {
                if((i >= ix) && (i < (ix + nx)) && (j >= iy) && (j < (iy + ny)) && (k >= iz) && (k < (iz + nz))) {
                    yield(i,i,j,j,k,k);
                }
            }
            visit_child = false;
        } else if(isLeaf(node) || isFull(node)){
            if (isInside(node) || isFull(node)) {
                var xlo = max(ix,i);
                var ylo = max(iy,j);
                var zlo = max(iz,k);
                var xhi = min(ix+nx,i+increment);
                var yhi = min(iy+ny,j+increment);
                var zhi = min(iz+nz,k+increment);
                //the C code does < checks rather than <= checks.
                yield(xlo,xhi-1,ylo,yhi-1,zlo,zhi-1);
            }
            visit_child = false;
        } else if (!isInside(node)) {
            visit_child = false;
        } else if (PV_visiting[l] < GrGeomOctreeNumChildren) {
            visit_child = true;
        } else {visit_child = false;}

        if (visit_child) {
            node = node.children[PV_visiting[l]][0];
            increment = increment >> 1;
            if (PV_visiting[l] & 1) {
                i += increment;
            }
            if (PV_visiting[l] & 2) {
                j += increment;
            }
            if (PV_visiting[l] & 4) {
                k += increment;
            }
            l += 1;
            PV_visiting[l] = 0;
        } else {
            l -= 1;
            if (PV_visiting[l] & 1) {
                i -= increment;
            }
            if (PV_visiting[l] & 2) {
                j -= increment;
            }
            if (PV_visiting[l] & 4) {
                k -= increment;
            }
            increment = increment << 1;
            if(node.parent != nil) {
                node = node.parent[0];
            }
            PV_visiting[l] += 1;
        }
    }
}

iter GrGeomOctreeFaceLoop_iter(i: int, j: int, k: int, in root: GrGeomOctree, level_of_interest: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int) {
    for (il,iu,jl,ju,kl,ku,f) in GrGeomOctreeInsideLoop(i,j,k,root,level_of_interest,ix,iy,iz,nx,ny,nz) {
        var fdir: [0..2] int = create_fdir(f);
        yield (il,iu,jl,ju,kl,ku,fdir);
   }
}

