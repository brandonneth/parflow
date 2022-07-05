use grgeom;
use CTypes;
use phase_rel_perm_chapel;
//config param call_only = 0;
proc mean(a,b) { return (a+b) / 2;}

proc harmonic_mean(a,b) {
    if (a + b == 0) {
        return 0;
    } else {
        return (2 * a * b) / (a + b);
    }
}
proc harmonic_mean_dz(a,b,c,d) {
    if((c*b) + (a*d) != 0) {
        return (((c+d) * a*b) / ((b*c) + (a*d)));
    } else {
        return 0;
    }
}
proc upstream_mean(a,b,c,d) {if (a - b) >= 0 then return c; else return d;}

proc Mean(a,b) {
    return mean(a,b);
}
proc PMean(a,b,c,d) {
    return harmonic_mean(c,d);
}
proc PMeanDZ(a,b,c,d) {
    return harmonic_mean_dz(a,b,c,d);
}

proc RPMean(a,b,c,d) {
    return upstream_mean(a,b,c,d);
}



//Corresponds to the GrGeomInLoop call on line 599 in richards_jacobian_eval.c
export proc richards_gravity_and_second_order_derivative_interior(ref gr_domain: GrGeomSolid,
r: int, ix: int, iy: int, iz: int, nx: int, ny: int, nz: int, //args for the iteration
pp: [] real, dp: [] real, rpp: [] real, ddp: [] real, rpdp: [] real, //presure, density, relperm, derivatives therein
permxp: [] real, permyp: [] real, permzp: [] real, //permeabilities
fb_x : [] real, fb_y: [] real, fb_z: [] real, //FB data (not exactly sure what FB means)
x_ssl: [] real, y_ssl: [] real, z_mult: [] real, //x and y slopes and variable dz
J: [] real, cp: [] real, wp: [] real, ep: [] real, sop: [] real, np: [] real, lp: [] real, up: [] real, // jacobian submatrix stencils
grid2d_iz: int, dx: int, dy: int, dz: int, dt: int, sy_v: int, sz_v: int, sy_m: int, sz_m: int,//scalar parameters, sizes, etc
tfgupwind: int, gravity: real, viscosity: real, // other scalars
pix: int, piy: int, piz: int, pnx: int, pny: int,
jix: int, jiy: int, jiz: int, jnx: int, jny: int,
six: int, siy: int, siz: int, snx: int, sny: int,
symm_part: int
)
{
    if(call_only) {
        return;
    }
    writeln("chapel richards gravity interior.");
    var ffx = dy * dz;
    var ffy = dx * dz;
    var ffz = dx * dy;
    
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(gr_domain, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        forall (i,j,k) in dom {
            var ip = subvector_elt_index(i,j,k,pix,piy,piz,pnx,pny);
            var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
            var ioo = subvector_elt_index(i,j,grid2d_iz,six,siy,siz,snx,sny);

            var prod = rpp[ip] * dp[ip];
            var prod_der = rpdp[ip] * dp[ip] + rpp[ip] * ddp[ip];

            var prod_rt = rpp[ip + 1] * dp[ip + 1];
            var prod_rt_der = rpdp[ip + 1] * dp[ip + 1] + rpp[ip + 1] * ddp[ip + 1];

            var prod_no = rpp[ip + sy_v] * dp[ip + sy_v];
            var prod_no_der = rpdp[ip + sy_v] * dp[ip + sy_v]
                    + rpp[ip + sy_v] * ddp[ip + sy_v];

            var prod_up = rpp[ip + sz_v] * dp[ip + sz_v];
            var prod_up_der = rpdp[ip + sz_v] * dp[ip + sz_v]
                    + rpp[ip + sz_v] * ddp[ip + sz_v];
            
            var x_dir_g: real = NAN;
            var x_dir_g_c: real = NAN;
            var y_dir_g: real = NAN;
            var y_dir_g_c: real = NAN;

            select tfgupwind {
                when 0 {
                    x_dir_g = mean(gravity * sin(atan(x_ssl[ioo])), gravity * sin(atan(x_ssl[ioo + 1])));
                    x_dir_g_c = mean(gravity * cos(atan(x_ssl[ioo])), gravity * cos(atan(x_ssl[ioo + 1])));
                    y_dir_g = mean(gravity * sin(atan(y_ssl[ioo])), gravity * sin(atan(y_ssl[ioo + sy_v])));
                    y_dir_g_c = mean(gravity * cos(atan(y_ssl[ioo])), gravity * cos(atan(y_ssl[ioo + sy_v])));
                }
                when 1 {
                    x_dir_g = gravity * sin(atan(x_ssl[ioo]));
                    x_dir_g_c = gravity * cos(atan(x_ssl[ioo]));
                    y_dir_g = gravity * sin(atan(y_ssl[ioo]));
                    y_dir_g_c = gravity * cos(atan(y_ssl[ioo]));
                }
                when 2 {
                    x_dir_g = x_ssl[ioo];
                    x_dir_g_c = 1.0;
                    y_dir_g = y_ssl[ioo];
                    y_dir_g_c = 1.0;
                }
            }

            var diff = pp[ip] - pp[ip + 1];
            var updir = (diff / dx) * x_dir_g_c - x_dir_g;

            var x_coeff = fb_x[ip] * dt * ffx * (1.0 / dx) * z_mult[ip]
                * PMean(pp[ip], pp[ip+1], permxp[ip], permxp[ip + 1])
                / viscosity;
            var sym_west_temp = (-x_coeff
                       * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM TFG contributions, sym


            var west_temp = (-x_coeff * diff
                        * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g_c
                        + sym_west_temp;

            west_temp += (x_coeff * dx * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g; //@RMM TFG contributions, non sym

            var sym_east_temp = (-x_coeff
                            * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM added sym TFG contributions

            var east_temp = (x_coeff * diff
                        * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g_c
                        + sym_east_temp;

            east_temp += -(x_coeff * dx * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g; //@RMM  TFG contributions non sym

            /* diff >= 0 implies flow goes south to north */
            diff = pp[ip] - pp[ip + sy_v];
            updir = (diff / dy) * y_dir_g_c - y_dir_g;


            /* multiply y_coeff by FB in y */
            var y_coeff = fb_y[ip] * dt * ffy * (1.0 / dy) * z_mult[ip]
                        * PMean(pp[ip], pp[ip + sy_v], permyp[ip], permyp[ip + sy_v])
                        / viscosity;

            var sym_south_temp = -y_coeff
                            * RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM TFG contributions, SYMM

            var south_temp = -y_coeff * diff
                        * RPMean(updir, 0.0, prod_der, 0.0) * y_dir_g_c
                        + sym_south_temp;

            south_temp += (y_coeff * dy * RPMean(updir, 0.0, prod_der, 0.0)) * y_dir_g; //@RMM TFG contributions, non sym


            var sym_north_temp = y_coeff
                            * -RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM  TFG contributions non SYMM

            var north_temp = y_coeff * diff
                        * RPMean(updir, 0.0, 0.0,
                                    prod_no_der) * y_dir_g_c
                        + sym_north_temp;

            north_temp += -(y_coeff * dy * RPMean(updir, 0.0, 0.0, prod_no_der)) * y_dir_g; //@RMM  TFG contributions non sym

            var sep = (dz * Mean(z_mult[ip], z_mult[ip + sz_v]));
            /* diff >= 0 implies flow goes lower to upper */


            var lower_cond = pp[ip] / sep - (z_mult[ip] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip] * gravity;

            var upper_cond = pp[ip + sz_v] / sep + (z_mult[ip + sz_v] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip + sz_v] * gravity;


            diff = lower_cond - upper_cond;

            /* multiply z_coeff by FB in z */
            var z_coeff = fb_z[ip] * dt * ffz
                        * PMeanDZ(permzp[ip], permzp[ip + sz_v], z_mult[ip], z_mult[ip + sz_v])
                        / viscosity;

            var sym_lower_temp = -z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                            * RPMean(lower_cond, upper_cond, prod,
                                        prod_up);

            var lower_temp = -z_coeff
                        * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                            + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip]
                                * RPMean(lower_cond, upper_cond, prod,
                                        prod_up)))
                        + sym_lower_temp;

            var sym_upper_temp = z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                            * -RPMean(lower_cond, upper_cond, prod,
                                        prod_up);

            var upper_temp = z_coeff
                        * (diff * RPMean(lower_cond, upper_cond, 0.0,
                                            prod_up_der)
                            + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip + sz_v]
                                * RPMean(lower_cond, upper_cond, prod,
                                        prod_up)))
                        + sym_upper_temp;

            cp[im] += -(west_temp + south_temp + lower_temp);
            cp[im + 1] += -east_temp;
            cp[im + sy_m] += -north_temp;
            cp[im + sz_m] += -upper_temp;   

            if(!symm_part) {
                ep[im] += east_temp;
                np[im] += north_temp;
                up[im] += upper_temp;

                wp[im + 1] += west_temp;
                sop[im + sy_m] += south_temp;
                lp[im + sz_m] += lower_temp;
            } else {
                ep[im] += sym_east_temp;
                np[im] += sym_north_temp;
                up[im] += sym_upper_temp;
            }         
        } // for ijk in dom
    } //for in GrGoemInLoop_iter
}//richards_gravity_and_second_order_derivative_interior
