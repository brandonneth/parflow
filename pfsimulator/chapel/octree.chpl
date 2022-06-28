extern record GrGeomOctree {
  var flags : c_uchar;
  var faces : c_uchar;
  var parent : c_ptr(grgeom_octree);
  var children : c_ptr(c_ptr(grgeom_octree));
}

const GrGeomOctreeNumFaces = 6;
const GrGeomOctreeNumChildren = 8;

const GrGeomOctreeNodeEmpty: uint(8) = 1;
const GrGeomOctreeNodeOutside: uint(8) = 2;
const GrGeomOctreeNodeInside: uint(8) = 4;
const GrGeomOctreeNodeFull: uint(8) = 8;
const GrGeomOctreeNodeLeaf: uint(8) = 16;


enum GrGeomOctreeFace {L,R,D,U,B,F};

//corresponds to GrGeomOctreeFaceValue
proc faceValue(face: GrGeomOctreeFace) {
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
