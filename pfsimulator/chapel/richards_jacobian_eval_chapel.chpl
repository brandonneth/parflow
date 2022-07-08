use grgeom;
use CTypes;
use phase_rel_perm_chapel;
//config param call_only = 0;
proc mean(a,b) { return (a+b) / 2;}

//config const dataParMinGranularity = 1000;
//config const dataParTasksPerLocale = 32;
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
    writeln("sy_m and sz_m:", sy_m, " ", sz_m);
    writeln("im params:", jix, " ", jiy, " ", jiz, " ", jnx, " ", jny);
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(gr_domain, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        var mod: domain(3) = {zl..zu,yl..yu,xl..xu};


        forall (k,j,i) in mod {
            
            
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
            /*
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
            }    */     
        }
    } //for in GrGoemInLoop_iter
}//richards_gravity_and_second_order_derivative_interior

export proc richards_gravity_and_second_order_derivative_interior_split(ref gr_domain: GrGeomSolid,
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
    writeln("chapel richards gravity interior split.");
    var ffx = dy * dz;
    var ffy = dx * dz;
    var ffz = dx * dy;
    
    for (xl,xu,yl,yu,zl,zu) in GrGeomInLoop_iter(gr_domain, r, ix, iy, iz, nx, ny, nz) {
        var dom: domain(3) = {xl..xu, yl..yu, zl..zu};
        var north_temp: [dom] real;
        var east_temp: [dom] real;
        var south_temp: [dom] real;
        var west_temp: [dom] real;
        var upper_temp: [dom] real;
        var lower_temp: [dom] real;
        var sym_upper_temp: [dom] real;
        var sym_lower_temp: [dom] real;
        var sym_east_temp: [dom] real;
        var sym_north_temp: [dom] real;
        var sym_south_temp: [dom] real;
        var sym_west_temp: [dom] real;
        //writeln("domain size:", dom.size);
        //for (i,j,k) in {0..1,0..1,0..1} {
        forall (i,j,k) in dom {

            //writeln("iteration: ", i, j, k);
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
            sym_west_temp[i,j,k] = (-x_coeff
                       * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM TFG contributions, sym


            west_temp[i,j,k] = (-x_coeff * diff
                        * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g_c
                        + sym_west_temp[i,j,k];

            west_temp[i,j,k] += (x_coeff * dx * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g; //@RMM TFG contributions, non sym

            sym_east_temp[i,j,k] = (-x_coeff
                            * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM added sym TFG contributions

            east_temp[i,j,k] = (x_coeff * diff
                        * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g_c
                        + sym_east_temp[i,j,k];

            east_temp[i,j,k] += -(x_coeff * dx * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g; //@RMM  TFG contributions non sym

            /* diff >= 0 implies flow goes south to north */
            diff = pp[ip] - pp[ip + sy_v];
            updir = (diff / dy) * y_dir_g_c - y_dir_g;


            /* multiply y_coeff by FB in y */
            var y_coeff = fb_y[ip] * dt * ffy * (1.0 / dy) * z_mult[ip]
                        * PMean(pp[ip], pp[ip + sy_v], permyp[ip], permyp[ip + sy_v])
                        / viscosity;

            sym_south_temp[i,j,k] = -y_coeff
                            * RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM TFG contributions, SYMM

            south_temp[i,j,k] = -y_coeff * diff
                        * RPMean(updir, 0.0, prod_der, 0.0) * y_dir_g_c
                        + sym_south_temp[i,j,k];

            south_temp[i,j,k] += (y_coeff * dy * RPMean(updir, 0.0, prod_der, 0.0)) * y_dir_g; //@RMM TFG contributions, non sym


            sym_north_temp[i,j,k] = y_coeff
                            * -RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM  TFG contributions non SYMM

            north_temp[i,j,k] = y_coeff * diff
                        * RPMean(updir, 0.0, 0.0,
                                    prod_no_der) * y_dir_g_c
                        + sym_north_temp[i,j,k];

            north_temp[i,j,k] += -(y_coeff * dy * RPMean(updir, 0.0, 0.0, prod_no_der)) * y_dir_g; //@RMM  TFG contributions non sym

            var sep = (dz * Mean(z_mult[ip], z_mult[ip + sz_v]));
            /* diff >= 0 implies flow goes lower to upper */


            var lower_cond = pp[ip] / sep - (z_mult[ip] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip] * gravity;

            var upper_cond = pp[ip + sz_v] / sep + (z_mult[ip + sz_v] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip + sz_v] * gravity;


            diff = lower_cond - upper_cond;

            /* multiply z_coeff by FB in z */
            var z_coeff = fb_z[ip] * dt * ffz
                        * PMeanDZ(permzp[ip], permzp[ip + sz_v], z_mult[ip], z_mult[ip + sz_v])
                        / viscosity;

            sym_lower_temp[i,j,k] = -z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                            * RPMean(lower_cond, upper_cond, prod,
                                        prod_up);

            lower_temp[i,j,k] = -z_coeff
                        * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                            + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip]
                                * RPMean(lower_cond, upper_cond, prod,
                                        prod_up)))
                        + sym_lower_temp[i,j,k];

            sym_upper_temp[i,j,k] = z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                            * -RPMean(lower_cond, upper_cond, prod,
                                        prod_up);

            upper_temp[i,j,k] = z_coeff
                        * (diff * RPMean(lower_cond, upper_cond, 0.0,
                                            prod_up_der)
                            + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip + sz_v]
                                * RPMean(lower_cond, upper_cond, prod,
                                        prod_up)))
                        + sym_upper_temp[i,j,k];

               
        } // for ijk in dom

        forall (i,j,k) in dom {
            var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
            cp[im] += -(west_temp[i,j,k] + south_temp[i,j,k] + lower_temp[i,j,k]);
             

            if(!symm_part) {
                ep[im] += east_temp[i,j,k];
                np[im] += north_temp[i,j,k];
                up[im] += upper_temp[i,j,k];

                wp[im + 1] += west_temp[i,j,k];
                sop[im + sy_m] += south_temp[i,j,k];
                lp[im + sz_m] += lower_temp[i,j,k];
            } else {
                ep[im] += sym_east_temp[i,j,k];
                np[im] += sym_north_temp[i,j,k];
                up[im] += sym_upper_temp[i,j,k];
            }      
        }
        forall (i,j,k) in dom {
            var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
            cp[im + 1] += -east_temp[i,j,k];
        }
        forall (i,j,k) in dom {
            var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
            cp[im + sy_m] += -north_temp[i,j,k];
        }
        forall (i,j,k) in dom {
            var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
            cp[im + sz_m] += -upper_temp[i,j,k];  
        }
        
            
    } //for in GrGoemInLoop_iter
}//richards_gravity_and_second_order_derivative_interior


export proc richards_gravity_and_second_order_derivative_interior_overlapped_tiled(ref gr_domain: GrGeomSolid,
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
        
        
        var x_dist = xu - xl + 1;
        var y_dist = yu - yl + 1;
        var z_dist = zu - zl + 1;
        const nti = 4;
        const ntj = 4;
        const ntk = 4;
        writeln("entire bounds:", xl, " ", xu, " ", yl, " ", yu, " ", zl, " ", zu, " ");
        var tile_size_x = x_dist / nti;
        var tile_size_y = y_dist / ntj;
        var tile_size_z = z_dist / ntk;
        writeln("tile sizes:", tile_size_x, " ", tile_size_y, " ", tile_size_z);
        var tile_domain: domain(3) = {0..nti,0..ntj,0..ntk};
        forall (ti,tj,tk) in tile_domain {
            //tile
            var x_start = max(xl, xl + ti * (x_dist / nti));
            var x_stop = min(xu+1, xl + (ti+1) * (x_dist / nti) );

            var y_start = max(yl, yl + tj * (y_dist / ntj));
            var y_stop = min(yu+1, yl + (tj+1) * (y_dist / ntj) );

            var z_start = max(zl, zl + tk * (z_dist / ntk));
            var z_stop = min(zu+1, zl + (tk+1) * (z_dist / ntk) );

            
            var base_domain: domain(3) = {x_start..<x_stop, y_start..<y_stop, z_start..<z_stop};

            var overlap_domain :domain(3) = {x_start-1..<x_stop, y_start-1..<y_stop, z_start-1..<z_stop};
            
            var north_temp: [overlap_domain] real;
            var east_temp: [overlap_domain] real;
            var south_temp: [overlap_domain] real;
            var west_temp: [overlap_domain] real;
            var upper_temp: [overlap_domain] real;
            var lower_temp: [overlap_domain] real;
            var sym_upper_temp: [overlap_domain] real;
            var sym_lower_temp: [overlap_domain] real;
            var sym_east_temp: [overlap_domain] real;
            var sym_north_temp: [overlap_domain] real;
            var sym_south_temp: [overlap_domain] real;
            var sym_west_temp: [overlap_domain] real;
            //writeln("domain size:", dom.size);
            //for (i,j,k) in {0..1,0..1,0..1} {
            for (i,j,k) in overlap_domain {
                //writeln("iteration: ", i, j, k);
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
                sym_west_temp[i,j,k] = (-x_coeff
                        * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM TFG contributions, sym


                west_temp[i,j,k] = (-x_coeff * diff
                            * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g_c
                            + sym_west_temp[i,j,k];

                west_temp[i,j,k] += (x_coeff * dx * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g; //@RMM TFG contributions, non sym

                sym_east_temp[i,j,k] = (-x_coeff
                                * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM added sym TFG contributions

                east_temp[i,j,k] = (x_coeff * diff
                            * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g_c
                            + sym_east_temp[i,j,k];

                east_temp[i,j,k] += -(x_coeff * dx * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g; //@RMM  TFG contributions non sym

                /* diff >= 0 implies flow goes south to north */
                diff = pp[ip] - pp[ip + sy_v];
                updir = (diff / dy) * y_dir_g_c - y_dir_g;


                /* multiply y_coeff by FB in y */
                var y_coeff = fb_y[ip] * dt * ffy * (1.0 / dy) * z_mult[ip]
                            * PMean(pp[ip], pp[ip + sy_v], permyp[ip], permyp[ip + sy_v])
                            / viscosity;

                sym_south_temp[i,j,k] = -y_coeff
                                * RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM TFG contributions, SYMM

                south_temp[i,j,k] = -y_coeff * diff
                            * RPMean(updir, 0.0, prod_der, 0.0) * y_dir_g_c
                            + sym_south_temp[i,j,k];

                south_temp[i,j,k] += (y_coeff * dy * RPMean(updir, 0.0, prod_der, 0.0)) * y_dir_g; //@RMM TFG contributions, non sym


                sym_north_temp[i,j,k] = y_coeff
                                * -RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM  TFG contributions non SYMM

                north_temp[i,j,k] = y_coeff * diff
                            * RPMean(updir, 0.0, 0.0,
                                        prod_no_der) * y_dir_g_c
                            + sym_north_temp[i,j,k];

                north_temp[i,j,k] += -(y_coeff * dy * RPMean(updir, 0.0, 0.0, prod_no_der)) * y_dir_g; //@RMM  TFG contributions non sym

                var sep = (dz * Mean(z_mult[ip], z_mult[ip + sz_v]));
                /* diff >= 0 implies flow goes lower to upper */


                var lower_cond = pp[ip] / sep - (z_mult[ip] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip] * gravity;

                var upper_cond = pp[ip + sz_v] / sep + (z_mult[ip + sz_v] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip + sz_v] * gravity;


                diff = lower_cond - upper_cond;

                /* multiply z_coeff by FB in z */
                var z_coeff = fb_z[ip] * dt * ffz
                            * PMeanDZ(permzp[ip], permzp[ip + sz_v], z_mult[ip], z_mult[ip + sz_v])
                            / viscosity;

                sym_lower_temp[i,j,k] = -z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                                * RPMean(lower_cond, upper_cond, prod,
                                            prod_up);

                lower_temp[i,j,k] = -z_coeff
                            * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                                + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip]
                                    * RPMean(lower_cond, upper_cond, prod,
                                            prod_up)))
                            + sym_lower_temp[i,j,k];

                sym_upper_temp[i,j,k] = z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                                * -RPMean(lower_cond, upper_cond, prod,
                                            prod_up);

                upper_temp[i,j,k] = z_coeff
                            * (diff * RPMean(lower_cond, upper_cond, 0.0,
                                                prod_up_der)
                                + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip + sz_v]
                                    * RPMean(lower_cond, upper_cond, prod,
                                            prod_up)))
                            + sym_upper_temp[i,j,k];
            }
            forall (i,j,k) in base_domain {
                var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
                cp[im] += -(west_temp[i,j,k] + south_temp[i,j,k] + lower_temp[i,j,k]);
                cp[im+1-1] += -east_temp[i-1,j,k];
                cp[im + sy_m - sy_m] += -north_temp[i,j-1,k];
                cp[im + sz_m - sz_m] += -upper_temp[i,j,k-1];  
                if(!symm_part) {
                    ep[im] += east_temp[i,j,k];
                    np[im] += north_temp[i,j,k];
                    up[im] += upper_temp[i,j,k];

                    wp[im + 1] += west_temp[i,j,k];
                    sop[im + sy_m] += south_temp[i,j,k];
                    lp[im + sz_m] += lower_temp[i,j,k];
                } else {
                    ep[im] += sym_east_temp[i,j,k];
                    np[im] += sym_north_temp[i,j,k];
                    up[im] += sym_upper_temp[i,j,k];
                }      
            }
        }
            
    }
}//richards_gravity_and_second_order_derivative_interior_overlapped_tiled

export proc richards_gravity_and_second_order_derivative_interior_overlapped_tiled_vectorized(ref gr_domain: GrGeomSolid,
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
        
        
        var x_dist = xu - xl + 1;
        var y_dist = yu - yl + 1;
        var z_dist = zu - zl + 1;
        const nti = 4;
        const ntj = 4;
        const ntk = 4;
        writeln("entire bounds:", xl, " ", xu, " ", yl, " ", yu, " ", zl, " ", zu, " ");
        var tile_size_x = x_dist / nti;
        var tile_size_y = y_dist / ntj;
        var tile_size_z = z_dist / ntk;
        writeln("tile sizes:", tile_size_x, " ", tile_size_y, " ", tile_size_z);
        var tile_domain: domain(3) = {0..nti,0..ntj,0..ntk};
        forall (ti,tj,tk) in tile_domain {
            //tile
            var x_start = max(xl, xl + ti * (x_dist / nti));
            var x_stop = min(xu+1, xl + (ti+1) * (x_dist / nti) );

            var y_start = max(yl, yl + tj * (y_dist / ntj));
            var y_stop = min(yu+1, yl + (tj+1) * (y_dist / ntj) );

            var z_start = max(zl, zl + tk * (z_dist / ntk));
            var z_stop = min(zu+1, zl + (tk+1) * (z_dist / ntk) );

            
            var base_domain: domain(3) = {x_start..<x_stop, y_start..<y_stop, z_start..<z_stop};

            var overlap_domain :domain(3) = {x_start-1..<x_stop, y_start-1..<y_stop, z_start-1..<z_stop};
            
            var north_temp: [overlap_domain] real;
            var east_temp: [overlap_domain] real;
            var south_temp: [overlap_domain] real;
            var west_temp: [overlap_domain] real;
            var upper_temp: [overlap_domain] real;
            var lower_temp: [overlap_domain] real;
            var sym_upper_temp: [overlap_domain] real;
            var sym_lower_temp: [overlap_domain] real;
            var sym_east_temp: [overlap_domain] real;
            var sym_north_temp: [overlap_domain] real;
            var sym_south_temp: [overlap_domain] real;
            var sym_west_temp: [overlap_domain] real;
            //writeln("domain size:", dom.size);
            //for (i,j,k) in {0..1,0..1,0..1} {
            for k in z_start-1..<z_stop {
                for j in y_start-1..<y_stop {
                    var ip_start = subvector_elt_index(x_start-1,j,k,pix,piy,piz,pnx,pny);
                    var ip_stop = subvector_elt_index(x_stop,j,k,pix,piy,piz,pnx,pny);

                    var im_start = subvector_elt_index(x_start-1,j,k,jix,jiy,jiz,jnx,jny);
                    var im_stop = subvector_elt_index(x_stop,j,k,jix,jiy,jiz,jnx,jny);
                    var i_dom = x_start-1..<x_stop;
                    var ip_dom = ip_start..<ip_stop;
                    var im_dom = im_start..<im_stop;
                    foreach (i,ip,im) in zip(i_dom, ip_dom, im_dom) {
                        
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
                sym_west_temp[i,j,k] = (-x_coeff
                        * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM TFG contributions, sym


                west_temp[i,j,k] = (-x_coeff * diff
                            * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g_c
                            + sym_west_temp[i,j,k];

                west_temp[i,j,k] += (x_coeff * dx * RPMean(updir, 0.0, prod_der, 0.0)) * x_dir_g; //@RMM TFG contributions, non sym

                sym_east_temp[i,j,k] = (-x_coeff
                                * RPMean(updir, 0.0, prod, prod_rt)) * x_dir_g_c; //@RMM added sym TFG contributions

                east_temp[i,j,k] = (x_coeff * diff
                            * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g_c
                            + sym_east_temp[i,j,k];

                east_temp[i,j,k] += -(x_coeff * dx * RPMean(updir, 0.0, 0.0, prod_rt_der)) * x_dir_g; //@RMM  TFG contributions non sym

                /* diff >= 0 implies flow goes south to north */
                diff = pp[ip] - pp[ip + sy_v];
                updir = (diff / dy) * y_dir_g_c - y_dir_g;


                /* multiply y_coeff by FB in y */
                var y_coeff = fb_y[ip] * dt * ffy * (1.0 / dy) * z_mult[ip]
                            * PMean(pp[ip], pp[ip + sy_v], permyp[ip], permyp[ip + sy_v])
                            / viscosity;

                sym_south_temp[i,j,k] = -y_coeff
                                * RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM TFG contributions, SYMM

                south_temp[i,j,k] = -y_coeff * diff
                            * RPMean(updir, 0.0, prod_der, 0.0) * y_dir_g_c
                            + sym_south_temp[i,j,k];

                south_temp[i,j,k] += (y_coeff * dy * RPMean(updir, 0.0, prod_der, 0.0)) * y_dir_g; //@RMM TFG contributions, non sym


                sym_north_temp[i,j,k] = y_coeff
                                * -RPMean(updir, 0.0, prod, prod_no) * y_dir_g_c; //@RMM  TFG contributions non SYMM

                north_temp[i,j,k] = y_coeff * diff
                            * RPMean(updir, 0.0, 0.0,
                                        prod_no_der) * y_dir_g_c
                            + sym_north_temp[i,j,k];

                north_temp[i,j,k] += -(y_coeff * dy * RPMean(updir, 0.0, 0.0, prod_no_der)) * y_dir_g; //@RMM  TFG contributions non sym

                var sep = (dz * Mean(z_mult[ip], z_mult[ip + sz_v]));
                /* diff >= 0 implies flow goes lower to upper */


                var lower_cond = pp[ip] / sep - (z_mult[ip] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip] * gravity;

                var upper_cond = pp[ip + sz_v] / sep + (z_mult[ip + sz_v] / (z_mult[ip] + z_mult[ip + sz_v])) * dp[ip + sz_v] * gravity;


                diff = lower_cond - upper_cond;

                /* multiply z_coeff by FB in z */
                var z_coeff = fb_z[ip] * dt * ffz
                            * PMeanDZ(permzp[ip], permzp[ip + sz_v], z_mult[ip], z_mult[ip + sz_v])
                            / viscosity;

                sym_lower_temp[i,j,k] = -z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                                * RPMean(lower_cond, upper_cond, prod,
                                            prod_up);

                lower_temp[i,j,k] = -z_coeff
                            * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)
                                + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip]
                                    * RPMean(lower_cond, upper_cond, prod,
                                            prod_up)))
                            + sym_lower_temp[i,j,k];

                sym_upper_temp[i,j,k] = z_coeff * (1.0 / (dz * Mean(z_mult[ip], z_mult[ip + sz_v])))
                                * -RPMean(lower_cond, upper_cond, prod,
                                            prod_up);

                upper_temp[i,j,k] = z_coeff
                            * (diff * RPMean(lower_cond, upper_cond, 0.0,
                                                prod_up_der)
                                + (-gravity * 0.5 * dz * (Mean(z_mult[ip], z_mult[ip + sz_v])) * ddp[ip + sz_v]
                                    * RPMean(lower_cond, upper_cond, prod,
                                            prod_up)))
                            + sym_upper_temp[i,j,k];
                    }
                }
            }
            
            forall (i,j,k) in base_domain {
                var im = subvector_elt_index(i,j,k,jix,jiy,jiz,jnx,jny);
                cp[im] += -(west_temp[i,j,k] + south_temp[i,j,k] + lower_temp[i,j,k]);
                cp[im+1-1] += -east_temp[i-1,j,k];
                cp[im + sy_m - sy_m] += -north_temp[i,j-1,k];
                cp[im + sz_m - sz_m] += -upper_temp[i,j,k-1];  
                if(!symm_part) {
                    ep[im] += east_temp[i,j,k];
                    np[im] += north_temp[i,j,k];
                    up[im] += upper_temp[i,j,k];

                    wp[im + 1] += west_temp[i,j,k];
                    sop[im + sy_m] += south_temp[i,j,k];
                    lp[im + sz_m] += lower_temp[i,j,k];
                } else {
                    ep[im] += sym_east_temp[i,j,k];
                    np[im] += sym_north_temp[i,j,k];
                    up[im] += sym_upper_temp[i,j,k];
                }      
            }
        }
            
    }
}//richards_gravity_and_second_order_derivative_interior_overlapped_tiled