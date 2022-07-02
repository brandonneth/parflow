use grgeom;
use CTypes;
use phase_rel_perm_chapel;


//Corresponds to the GrGeomInLoop call on line 599 in richards_jacobian_eval.c
export proc richards_gravity_and_second_order_derivative_interior(ref gr_domain: GrGeomSolid,
r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, //args for the iteration
pp: [] real, dp: [] real, rpp: [] real, ddp: [] real, rpdp: [] real, //presure, density, relperm, derivatives therein
permxp: [] real, permyp: [] real, permzp: [] real, //permeabilities
fb_x : [] real, fb_y: [] real, fb_z: [] real, //FB data (not exactly sure what FB means)
x_ssl: [] real, y_ssl: [] real, z_mult: [] real, //x and y slopes and variable dz
cp: [] real, wp: [] real, ep: [] real, sop: [] real, np: [] real, lp: [] real, up: [] real, // jacobian submatrix stencils
grid2d_iz: int, dx: int, dy: int, dz: int, dt: int, sy_v: int, sz_v: int, //scalar parameters, sizes, etc
tfgupwind: int, gravity: real, viscosity: real // other scalars
)
{
    writeln("chapel richards gravity interior.");
    var ffx = dy * dz;
    var ffy = dx * dz;
    var ffz = dx * dy;
    
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(gr_domain, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        for (i,j,k) in dom {

            
        } // for ijk in dom
    } //for in GrGoemInLoop_iter
}//richards_gravity_and_second_order_derivative_interior
