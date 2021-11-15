#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265

// LETS PROTOTYPE THIS FOR 2 BODIES

// Hacking quad and tidal forces in for the 2-body system
void additional_forces(struct reb_simulation* sim){

  struct reb_particle* p1 = &(sim->particles[0]);
  struct reb_particle* p2 = &(sim->particles[1]);

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.2;// rebx_get_param(sim->extras, p1->ap, "k");
  // const double* const q1 = rebx_get_param(sim->extras, p1->ap, "Q");
  // const double* const sx1 = rebx_get_param(sim->extras, p1->ap, "sx");
  // const double* const sy1 = rebx_get_param(sim->extras, p1->ap, "sy");
  // const double* const sz1 = rebx_get_param(sim->extras, p1->ap, "sz");
  const double moi1 = 0.1;// rebx_get_param(sim->extras, p1->ap, "MOI"); // Are there 3 MOIs????
  const double s1 = p1->r;
  const double m1 = p1->m;

  // for p2
  const double k2 = 0.3; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q2 = 10000; //rebx_get_param(sim->extras, p2->ap, "Q");
  //const double* const sx2 = rebx_get_param(sim->extras, p2->ap, "sx");
  //const double* const sy2 = rebx_get_param(sim->extras, p2->ap, "sy");
  //const double* const sz2 = rebx_get_param(sim->extras, p2->ap, "sz");
  const double moi2 = 0.3;//rebx_get_param(sim->extras, p2->ap, "MOI");
  const double s2 = p2->r;
  const double m2 = p2->m;

  // Dummy particles to track spin motions.  vx, vy, vz correspond to sx, sy, sz
  struct reb_particle* dummy1 = &(sim->particles[2]); // star
  struct reb_particle* dummy2 = &(sim->particles[3]); // planet

  // extract spin info
  double sx1 = dummy1->vx; // set using initial conditions
  double sy1 = dummy1->vy;
  double sz1 = dummy1->vz;

  double sx2 = dummy2->vx; // set using initial conditions
  double sy2 = dummy2->vy;
  double sz2 = dummy2->vz;

  // Orbital elements for the mutual orbit
  struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *p1, *p2);
  double n = o.n;
  double a = o.a;

  // distance vector
  double dx = p2->x - p1->x;
  double dy = p2->y - p1->y;
  double dz = p2->z - p1->z;
  double dist = sqrt(dx * dx + dy * dy + dz * dz);
  // Unit distance vector
  double dx_hat = dx / dist;
  double dy_hat = dy / dist;
  double dz_hat = dz / dist;

  // quadrupole moment of body 1
  double quad_prefactor_1 = (pow(s1, 5) * (1 + (m1 / m2)) * k1) / pow(dist, 4);
  double omega_dot_rhat_1 = sx1 * dx_hat + sy1 * dy_hat + sz1 * dz_hat;
  double omega_squared_1 = sx1 * sx1 + sy1 * sy1 + sz1 * sz1;

  double qx1 = quad_prefactor_1 * ((5 * pow(omega_dot_rhat_1, 2) - omega_squared_1 - (12. * sim->G * m2 / pow(dist, 3))) * dx_hat - (2 * omega_dot_rhat_1 * (sx1)));
  double qy1 = quad_prefactor_1 * ((5 * pow(omega_dot_rhat_1, 2) - omega_squared_1 - (12. * sim->G * m2 / pow(dist, 3))) * dy_hat - (2 * omega_dot_rhat_1 * (sy1)));
  double qz1 = quad_prefactor_1 * ((5 * pow(omega_dot_rhat_1, 2) - omega_squared_1 - (12. * sim->G * m2 / pow(dist, 3))) * dz_hat - (2 * omega_dot_rhat_1 * (sz1)));

  // quadrupole moment of body 2
  double quad_prefactor_2 = (pow(s2, 5) * (1 + (m2 / m1) * k2)) / pow(dist, 4);
  double omega_dot_rhat_2 = sx2 * dx_hat + sy2 * dy_hat + sz2 * dz_hat;
  double omega_squared_2 = sx2 * sx2 + sy2 * sy2 + sz2 * sz2;

  double qx2 = quad_prefactor_2 * ((5 * pow(omega_dot_rhat_2, 2) - omega_squared_2 - (12. * sim->G * m1 / pow(dist, 3))) * dx_hat - (2 * omega_dot_rhat_2 * sx2));
  double qy2 = quad_prefactor_2 * ((5 * pow(omega_dot_rhat_2, 2) - omega_squared_2 - (12. * sim->G * m1 / pow(dist, 3))) * dy_hat - (2 * omega_dot_rhat_2 * sy2));
  double qz2 = quad_prefactor_2 * ((5 * pow(omega_dot_rhat_2, 2) - omega_squared_2 - (12. * sim->G * m1 / pow(dist, 3))) * dz_hat - (2 * omega_dot_rhat_2 * sz2));

  // TIDAL DAMPING
  // Relative velocities of two bodies
  double dvx = p2->vx - p1->vx;
  double dvy = p2->vy - p1->vy;
  double dvz = p2->vz - p1->vz;

  // relevant dot and cross products
  double rhat_dot_vel = dx_hat * dvx + dy_hat * dvy + dz_hat * dvz;
  double rhat_cross_v_x = dy_hat * dvz - dz_hat * dvy;
  double rhat_cross_v_y = dz_hat * dvx - dx_hat * dvz;
  double rhat_cross_v_z = dx_hat * dvy - dy_hat * dvx;

  // first bracketed vector - constant across both bodies
  double vec1_x = 3 * rhat_dot_vel * dx_hat;
  double vec1_y = 3 * rhat_dot_vel * dy_hat;
  double vec1_z = 3 * rhat_dot_vel * dz_hat;

  // FOR BODY 1 IGNORE TIDES ON STAR FOR NOW
  // double prefactor_1 = -(6 * n * (*k1) / (*q1)) * (m2 / m1) * pow((s1 / a), 5) * pow((a / dist));

  // Build vector 2 - this is the parenthesis term in (4)
  // double temp_x1 = rhat_cross_v_x - dist * sx1;
  // double temp_y1 = rhat_cross_v_y - dist * sy1;
  // double temp_z1 = rhat_cross_v_z - dist * sz1;

  // double vec2_x1 = temp_y1 * dz_hat - temp_z1 * dy_hat;
  // double vec2_y1 = temp_z1 * dx_hat - temp_x1 * dz_hat;
  // double vec2_z1 = temp_x1 * dy_hat - temp_y1 * dx_hat;

  double tx1 = 0; //vec1_x + vec2_x1;
  double ty1 = 0; //vec1_y + vec2_y1;
  double tz1 = 0; //vec1_z + vec2_z1;

  // FOR BODY 2
  double prefactor_2 = -((6 * n * k2) / q2) * (m1 / m2) * pow((s2 / a), 5) * pow((a / dist), 8);

  // Build vector 2 - this is the parenthesis term in (4)
  double temp_x2 = rhat_cross_v_x - dist * sx2;
  double temp_y2 = rhat_cross_v_y - dist * sy2;
  double temp_z2 = rhat_cross_v_z - dist * sz2;

  double vec2_x2 = temp_y2 * dz_hat - temp_z2 * dy_hat;
  double vec2_y2 = temp_z2 * dx_hat - temp_x2 * dz_hat;
  double vec2_z2 = temp_x2 * dy_hat - temp_y2 * dx_hat;

  double tx2 = prefactor_2 * (vec1_x + vec2_x2);
  double ty2 = prefactor_2 * (vec1_y + vec2_y2);
  double tz2 = prefactor_2 * (vec1_z + vec2_z2);

  // total forces for convenience
  double totx_1 = qx1 + tx1;
  double toty_1 = qy1 + ty1;
  double totz_1 = qz1 + tz1;

  double totx_2 = qx2 + tx2;
  double toty_2 = qy2 + ty2;
  double totz_2 = qz2 + tz2;

  // Force on planet 1
  p1->ax += (m2 / (m1 + m2)) * (totx_1 + totx_2);
  p1->ay += (m2 / (m1 + m2)) * (toty_1 + toty_2);
  p1->az += (m2 / (m1 + m2)) * (totz_1 + totz_2);

  // Force on planet 2
  p2->ax += (m1 / (m1 + m2)) * (totx_1 + totx_2);
  p2->ay += (m1 / (m1 + m2)) * (toty_1 + toty_2);
  p2->az += (m1 / (m1 + m2)) * (totz_1 + totz_2);

  // printf("Quad Force 1: %f\n", sqrt(qx1 * qx1 + qy1 * qy1 + qz1 * qz1));
  // printf("Quad Force 2: %f\n", sqrt(qx2 * qx2 + qy2 * qy2 + qz2 * qz2));

  // printf("Tidal Force 1: %f\n", sqrt(tx1 * tx1 + ty1 * ty1 + tz1 * tz1));
  // printf("Tidal Force 2: %f\n", sqrt(tx2 * tx2 + ty2 * ty2 + tz2 * tz2));

  // Rewuired for convergence criteria only
  dummy1->x =  1.;
  dummy1->y =  1.;
  dummy1->z =  1.;
  dummy2->x =  1.;
  dummy2->y =  1.;
  dummy2->z =  1.;

  // Overwrite gravity acceleration
  dummy1->ax = 0.;
  dummy1->ay = 0.;
  dummy1->az = 0.;
  dummy2->ax = 0.;
  dummy2->ay = 0.;
  dummy2->az = 0.;

  // EOMs for spin vector
  double mu12 = -(p1->m * p2->m) / ((p1->m + p2->m));
  double rvec_x = mu12 * dx;
  double rvec_y = mu12 * dy;
  double rvec_z = mu12 * dz;

  // Now take the cross product
  dummy1->ax = (rvec_y * totz_1 - rvec_z * toty_1) / moi1;
  dummy1->ay = (rvec_z * totx_1 - rvec_x * totz_1) / moi1;
  dummy1->az = (rvec_x * toty_1 - rvec_y * totx_1) / moi1;

  dummy2->ax = (rvec_y * totz_2 - rvec_z * toty_2) / moi2;
  dummy2->ay = (rvec_z * totx_2 - rvec_x * totz_2) / moi2;
  dummy2->az = (rvec_x * toty_2 - rvec_y * totx_2) / moi2;

  // double dsx = rvec_y2 * totz_2 - rvec_z2 * toty_2;
  // double dsy = rvec_z2 * totx_2 - rvec_x2 * totz_2;
  // double dsz = rvec_x1 * toty_1 - rvec_y1 * totx_1;

  // if (sim->t < 0.01)
  // printf("%f, %.15f, %.15f, %.15f, %15f\n", sim->t, rvec_z2, totx_2, rvec_x2, totz_2);
  // printf("%.15f, %.15f, %.15f\n", dsx, dsy, dsz);
  //printf("%f, %.15f, %.15f, %.15f, %.6f, %.6f, %.6f\n", sim->t, rvec_x2, rvec_y2, rvec_z2, dx, dy, dz);

  // if (sim->t < 100)
  //printf("%f, %f, %f, %f, %f, %f, %f\n", sim->t, dsx, dsy, dsz, dummy2->vx, dummy2->vy, dummy2->vz);
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m r", 1., 0.00465);                // Central object


    reb_add_fmt(r, "m a e r inc", 1e-3, 0.04072, 0.01, 0.000476, 0.087); // Hot Jupiter
    reb_move_to_com(r);

    r->N_active = 2;
    // Tot rot angular momentum is 0.8 * pi * m * r^2 / trot
    // take solar obliquity = 6 degrees
    double solar_spin_period = 27 * 2 * PI / 365;
    double solar_spin = (0.8 * PI * 0.00465 * 0.00465) / solar_spin_period;
    reb_add_fmt(r, "m vx vy vz", 0., 0., solar_spin * sin(6. * (PI / 180.)), solar_spin * cos(6. * (PI / 180.)));         // Dummy particle
                                     // Mass is not relevant

    double jup_spin_period = 5 * 2 * PI / 365;
    double jup_spin = (0.8 * PI * 1e-3 * 0.000476 * 0.000476) / jup_spin_period;
    double obliquity = 10 * (PI / 180);
    double res_angle = 15 * (PI / 180);
    double spin_x = jup_spin * cos(res_angle) * sin(obliquity);
    double spin_y = jup_spin * sin(res_angle) * sin(5. * (PI / 180.));
    double spin_z = jup_spin * cos(5. * (PI / 180.));
    // printf("%f\n", jup_spin);
    reb_add_fmt(r, "m vx vy vz", 0., spin_x, spin_y, spin_z);
    // printf("%.15f, %.15f, %.15f\n", spin_y, spin_z, sqrt(spin_y * spin_y + spin_z * spin_z));
    // printf("%f\n", r->particles[3].vx);


    r->additional_forces = additional_forces;
    r->integrator = REB_INTEGRATOR_IAS15;
    r->force_is_velocity_dependent = 1;


   FILE* f = fopen("out.txt","w");

    for (int i=0; i<300; i++){
        reb_integrate(r,r->t+100);

        struct reb_particle sun = r->particles[0];
        struct reb_particle p = r->particles[1];
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, p, sun);
        struct reb_particle dummy1 = r->particles[2];
        struct reb_particle dummy2 = r->particles[3];

        struct reb_vec3d energy = reb_tools_angular_momentum(r);
        double mag = sqrt(dummy2.vx * dummy2.vx + dummy2.vy * dummy2.vy + dummy2.vz * dummy2.vz);
        double tot_angular_momentum = sqrt(energy.x * energy.x + energy.y * energy.y + energy.z * energy.z);// + mag + sqrt(dummy1.vx * dummy1.vx + dummy1.vy * dummy1.vy + dummy1.vz * dummy1.vz);

        double torb = ((0.8 * PI * 1e-3 * 0.000476 * 0.000476) / (mag * 2 * PI)) * 365;
        double obliquity = acos(dummy2.vz / mag) * (180 / PI);
        printf("torb=%.10f \t t=%.4f\t planet = %6.3f \t spin axis = %.10f %.10f %.10f \t ob = %.3f \t ang = %.10f\n", torb, r->t, o.a, dummy2.vx, dummy2.vy, dummy2.vz, obliquity, tot_angular_momentum);
        fprintf(f, "%e %e %e\n", r->t, p.x, dummy1.vx);
    }
    fclose(f);

    reb_free_simulation(r);
}
