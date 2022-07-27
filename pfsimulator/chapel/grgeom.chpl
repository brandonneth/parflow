use CTypes;

use octree;

require "index_space.h";
require "grgeometry.h";
require "input_database.h";
extern type Point = 3*int(32);

extern record Box {
    var lo: Point;
    var up: Point;

    proc dom() { return {lo[0]..up[0],lo[1]..up[1],lo[2]..up[2]}; }

    iter points(param tag: iterKind) where tag == iterKind.standalone {
        forall point in dom() do yield point;
    }
};

extern record BoxArray {
  var boxes : c_ptr(Box);
  var boxlimits : c_ptr(c_int);
  var size : c_uint;

  iter these(param tag: iterKind) where tag == iterKind.standalone
  {
    if boxes == nil then return;
    forall i in 0..<size do
        yield boxes[i];
  }
  iter these()
  {
    if boxes == nil then return;
    for i in 0..<size do
        yield boxes[i];
  }


}

extern record GrGeomSolid {
    var data: c_ptr(GrGeomOctree);
    var patches: c_ptr(c_ptr(GrGeomOctree));
    var octree_bg_level, octree_ix, octree_iy, octree_iz: int;
    var interior_boxes: c_ptr(BoxArray);
    var surface_boxes : c_ptr(c_ptr(BoxArray));
    var patch_boxes : c_ptr(c_ptr(c_ptr(BoxArray)));

    proc interiorBoxes() { return interior_boxes[0]; }

    proc surfaceBoxes(face: int) { 
        return surface_boxes[face][0]; 
    }

    proc patchBoxes(face: int, patchNum: int) { return patch_boxes[face][patchNum][0]; }
}


iter groundGeometryPatchBoxes(param tag: iterKind, ref groundGeometry: GrGeomSolid, patchNum: int, minPoint: Point, maxPoint: Point)
    where tag == iterKind.standalone
{
    forall f in 0..<GrGeomOctreeNumFaces {
        const boxArray = groundGeometry.patchBoxes(f,patchNum);
        forall boxIndex in 0..<boxArray.size {
            const box = boxArray.boxes[boxIndex];
            var low: Point = max(minPoint, box.lo);
            var high: Point = min(maxPoint, box.up);
            var points: domain(3) = {low[0]..high[0],low[1]..high[1],low[2]..high[2]};
            forall (i,j,k) in points do
                yield (i,j,k,f);
        }
    }
}

iter groundGeometryInteriorBoxes(ref groundGeometry: GrGeomSolid, outerDom: domain(3, int(32))) {
    for box in groundGeometry.interiorBoxes() do
        for point in box.dom()[outerDom] do
            yield point;
}
iter groundGeometryInteriorBoxes(param tag: iterKind, ref groundGeometry: GrGeomSolid, outerDom: domain(3, int(32))) 
  where tag == iterKind.standalone {
    forall box in groundGeometry.interiorBoxes() do
        forall point in box.dom()[outerDom] do
            yield point;
}

iter groundGeometrySurfaceBoxes(ref groundGeometry: GrGeomSolid, outerDom: domain(3, int(32))) {
    for face in 0..<GrGeomOctreeNumFaces {
        for box in groundGeometry.surfaceBoxes(face) {
            for (i,j,k) in box.dom()[outerDom] do
                yield (i,j,k,create_fdir(face));
        }
    }
}
