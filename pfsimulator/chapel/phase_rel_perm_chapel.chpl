use grgeom;
use CTypes;
config param call_only = 0;
proc subvector_elt_index(x,y,z,ix,iy,iz,nx,ny) {
    return ((x - ix) + ((y - iy) + ((z - iz) *  ny)) * nx);
}

export proc calcfcn_compute_vang_curve_surface(ref grgeom:GrGeomSolid, r: int, 
    ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, 
    pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, 
    gravity: real, region_idx: int, ixv: int, iyv: int, izv: int, nxv: int, nyv: int) {
    
    writeln("calcfcn_compute_vang_curve_surface");
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
}


export proc calcder_compute_vang_curve_surface(ref grgeom:GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, gravity: real, region_idx: int, ixv: int, iyv: int, izv: int, nxv: int, nyv: int) {
    if(call_only) {
        return;
    }
    for (xl,xu,yl,yu,zl,zu,fdir) in GrGeomSurfLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        //writeln("domain size:", dom.size);
        for (i,j,k) in dom {

            //do thing
            var idx = subvector_elt_index(i + fdir[0], j + fdir[1], k + fdir[2], ixv, iyv, izv, nxv, nyv);

            if (pp_sub[idx] >= 0.0) {
                pr_sub[idx] = 0.0;
            } else {
                var alpha = alphas[region_idx];
                var n = ns[region_idx];
                var m = 1.0 - (1.0 / n);

                const head = abs(pp_sub[idx]) / (pd_sub[idx] * gravity);
                const opahn = 1.0 + (alpha * head) ** n;
                const ahnm1 = (alpha * head) ** (n-1);
                const coeff = 1.0 - ahnm1 * (opahn ** -m);

                //writeln("head, opahn, ahnm1, prsub val:", head, " ", opahn, " ",ahnm1, " ",top / bottom);
                pr_sub[idx] =  2.0 * (coeff / ((opahn ** (m / 2))))
                                 * ((n - 1) * ((alpha * head) ** (n - 2)) * alpha
                                    * (opahn ** -m)
                                    - ahnm1 * m * (opahn ** -(m + 1)) * n * alpha * ahnm1)
                                 + (coeff ** 2) * (m / 2) * (opahn ** (-(m + 2) / 2))
                                 * n * alpha * ahnm1;
            }
        }
    }
}
export proc calcder_compute_vang_curve_interior(ref grgeom:GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, gravity: real, region_idx: int, ixv: int, iyv: int, izv: int, nxv: int, nyv: int) {
    if(call_only) {
        return;
    }
    
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        //writeln("computing for:", (xl,xu,yl,yu,zl,zu));
        for (i,j,k) in dom {

            //do thing
            var idx = subvector_elt_index(i , j , k, ixv, iyv, izv, nxv, nyv);

            if (pp_sub[idx] >= 0.0) {
                pr_sub[idx] = 0.0;
            } else {
                var alpha = alphas[region_idx];
                var n = ns[region_idx];
                var m = 1.0 - (1.0 / n);

                const head = abs(pp_sub[idx]) / (pd_sub[idx] * gravity);
                const opahn = 1.0 + (alpha * head) ** n;
                const ahnm1 = (alpha * head) ** (n-1);
                const coeff = 1.0 - ahnm1 * (opahn ** -m);

                //writeln("head, opahn, ahnm1, prsub val:", head, " ", opahn, " ",ahnm1, " ",top / bottom);
                pr_sub[idx] =  2.0 * (coeff / ((opahn ** (m / 2))))
                                 * ((n - 1) * ((alpha * head) ** (n - 2)) * alpha
                                    * (opahn ** -m)
                                    - ahnm1 * m * (opahn ** -(m + 1)) * n * alpha * ahnm1)
                                 + (coeff ** 2) * (m / 2) * (opahn ** (-(m + 2) / 2))
                                 * n * alpha * ahnm1;
            }
        }
    }
}
export proc calcfcn_compute_vang_curve_interior(ref grgeom:GrGeomSolid, r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, pr_sub: [] real, pp_sub: [] real, pd_sub: [] real, alphas: []real, ns: [] real, gravity: real, region_idx: int, ixv: int, iyv: int, izv: int, nxv: int, nyv: int) {
    if(call_only) {
        return;
    }

    //writeln("pr_sub sum in chapel, start: ", sum1);
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(grgeom, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        //writeln("computing for:", (xl,xu,yl,yu,zl,zu));
        for (i,j,k) in dom {

            //do thing
            var idx = subvector_elt_index(i, j, k, ixv, iyv, izv, nxv, nyv);

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
}
