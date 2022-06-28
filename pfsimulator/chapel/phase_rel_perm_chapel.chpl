use grgeom;
use CTypes;

export proc compute_vang_curve(ref grgeom:GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, gravity: real) {
    writeln("computing vang curve.");
    for (xl,xu,yl,yu,zl,zu,fdir) in GrGeomSurfLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..<xu, yl..<yu, zl..<zu};
        writeln("computing for:", (xl,xu,yl,yu,zl,zu));
        for (i,j,k) in dom {
            //do thing
        }
    }
}
