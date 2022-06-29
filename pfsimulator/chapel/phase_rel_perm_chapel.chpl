use grgeom;
use CTypes;

proc subvector_elt_index(x,y,z,ix,iy,iz,nx,ny) {
    return (((x) - ix) + (((y) - iy) + (((z) - iz)) *  ny) * nx);
}

export proc compute_vang_curve(ref grgeom:GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, gravity: real, region_idx: int, ixv: int, iyv: int, izv: int, nxv: int, nyv: int) {
    //writeln("computing vang curve.");
    var sum1 = 0.0;
    for i in 0..<pr_sub.size {
        sum1 += pr_sub[i];
    }
    //writeln("pr_sub sum in chapel, start: ", sum1);
    for (xl,xu,yl,yu,zl,zu,fdir) in GrGeomSurfLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        //writeln("computing for:", (xl,xu,yl,yu,zl,zu));
        for (i,j,k) in dom {

            //do thing
            var idx = subvector_elt_index(i + fdir[0], j + fdir[1], k + fdir[2], ixv, iyv, izv, nxv, nyv);

            if (pp_sub[idx] >= 0.0) {
                pr_sub[idx] = 1.0;
            } else {
                var alpha = alphas[region_idx];
                var n = ns[region_idx];
                var m = 1.0 - (1.0 / n);

                const head = abs(pp_sub[idx]) / (pd_sub[idx] * gravity);
                const opahn = 1.0 + (alpha * head) ** n;
                const ahnm1 = (alpha * head) ** (n-1);
                const top = (1.0 - ahnm1 / (opahn ** m)) ** 2;
                const bottom = opahn ** (m / 2);
                //writeln("head, opahn, ahnm1, prsub val:", head, " ", opahn, " ",ahnm1, " ",top / bottom);
                pr_sub[idx] = top / bottom;
            }
        }
    }
    var sum2 = 0.0;
    for i in 0..<pr_sub.size {
        sum2 += pr_sub[i];
    }
    //writeln("pr_sub sum in chapel, end: ", sum2);
}
