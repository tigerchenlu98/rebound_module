#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265

// LETS PROTOTYPE THIS FOR 3 BODIES
// Hacking quad and tidal forces in for the 3-body system
void additional_forces(struct reb_simulation* sim){

  struct reb_particle* p1 = &(sim->particles[0]);
  struct reb_particle* p2 = &(sim->particles[1]);
  struct reb_particle* p3 = &(sim->particles[2]);

  // Migration forces
  if (sim->t <= 2e6 * 2 * PI){
    p2->vx *= (1 - 1e-7 * sim->dt);
    p2->vy *= (1 - 1e-7 * sim->dt);
    p2->vz *= (1 - 1e-7 * sim->dt);

    p3->vx *= (1 - 1.1 * 1e-7 * sim->dt);
    p3->vy *= (1 - 1.1 * 1e-7 * sim->dt);
    p3->vz *= (1 - 1.1 * 1e-7 * sim->dt);
  }

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.035;// rebx_get_param(sim->extras, p1->ap, "k");
  const double q1 = 100000.; //rebx_get_param(sim->extras, p1->ap, "Q");
  // const double* const sx1 = rebx_get_param(sim->extras, p1->ap, "sx");
  // const double* const sy1 = rebx_get_param(sim->extras, p1->ap, "sy");
  // const double* const sz1 = rebx_get_param(sim->extras, p1->ap, "sz");
  const double moi1 = 0.4 * p1->m * p1->r * p1->r;// rebx_get_param(sim->extras, p1->ap, "MOI"); // Are there 3 MOIs????

  // for p2
  const double k2 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q2 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");
  const double moi2 = 0.4 * p2->m * p2->r * p2->r;//rebx_get_param(sim->extras, p2->ap, "MOI");

  // for p3
  const double k3 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q3 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");
  const double moi3 = 0.4 * p3->m * p3->r * p3->r;//rebx_get_param(sim->extras, p2->ap, "MOI");

  // Dummy particles to track spin motions.  vx, vy, vz correspond to sx, sy, sz
  struct reb_particle* dummy1 = &(sim->particles[3]); // star
  struct reb_particle* dummy2 = &(sim->particles[4]); // planet 1
  struct reb_particle* dummy3 = &(sim->particles[5]); // planet 2

  // Pack these into arrays
  const double apsidals[3] = {k1, k2, k3};
  const double tidals[3] = {q1, q2, q3};
  const double mois[3] = {moi1, moi2, moi3};

  // Required for convergence criteria only
  dummy1->x =  1.;
  dummy1->y =  1.;
  dummy1->z =  1.;
  dummy2->x =  1.;
  dummy2->y =  1.;
  dummy2->z =  1.;
  dummy3->x =  1.;
  dummy3->y =  1.;
  dummy3->z =  1.;

  // Overwrite gravity acceleration
  dummy1->ax = 0.;
  dummy1->ay = 0.;
  dummy1->az = 0.;
  dummy2->ax = 0.;
  dummy2->ay = 0.;
  dummy2->az = 0.;
  dummy3->ax = 0.;
  dummy3->ay = 0.;
  dummy3->az = 0.;

  // Calculates the effect of perturbing body j on body i
  int nb = 3;
  for (int i = 0; i < nb; i++){
    for (int j = 0; j < nb; j++){
      if (i != j){
        // Extract relevant physical parameters
        struct reb_particle* pi = &(sim->particles[i]);
        struct reb_particle* pj = &(sim->particles[j]);

        double ri = pi->r;
        double mi = pi->m;
        double mj = pj->m;

        // Orbital elements for the mutual orbit
        // This is a bit tricky b/c sometimes the star acts as the perturber: need to check against this
        double n;
        double a;
        if (j == 0) {
          struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pi);
          n = o.n;
          a = o.a;
        }
        else {
          struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pj);
          n = o.n;
          a = o.a;
        }

        // Extract spin info from dummy
        struct reb_particle* dummy = &(sim->particles[i + 3]);
        double sx = dummy->vx;
        double sy = dummy->vy;
        double sz = dummy->vz;

        // distance vector FROM j TO i
        double dx = pj->x - pi->x;
        double dy = pj->y - pi->y;
        double dz = pj->z - pi->z;
        double dist = sqrt(dx * dx + dy * dy + dz * dz);
        // Unit distance vector
        double dx_hat = dx / dist;
        double dy_hat = dy / dist;
        double dz_hat = dz / dist;

        // acceleration on i due to quadrupole of j
        double quad_prefactor = (pow(ri, 5) * (1 + (mi / mj)) * apsidals[i]) / pow(dist, 4);
        double omega_dot_rhat = sx * dx_hat + sy * dy_hat + sz * dz_hat;
        double omega_squared = sx * sx + sy * sy + sz * sz;

        double qx = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dx_hat - (2 * omega_dot_rhat * sx));
        double qy = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dy_hat - (2 * omega_dot_rhat * sy));
        double qz = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dz_hat - (2 * omega_dot_rhat * sz));

        // Velocity vecto: j to i
        double dvx = pj->vx - pi->vx;
        double dvy = pj->vy - pi->vy;
        double dvz = pj->vz - pi->vz;

        // relevant dot and cross products
        double rhat_dot_vel = dx_hat * dvx + dy_hat * dvy + dz_hat * dvz;
        double rhat_cross_v_x = dy_hat * dvz - dz_hat * dvy;
        double rhat_cross_v_y = dz_hat * dvx - dx_hat * dvz;
        double rhat_cross_v_z = dx_hat * dvy - dy_hat * dvx;

        // first bracketed vector
        double vec1_x = 3 * rhat_dot_vel * dx_hat;
        double vec1_y = 3 * rhat_dot_vel * dy_hat;
        double vec1_z = 3 * rhat_dot_vel * dz_hat;

        // Tidal damping on i due to j
        double prefactor = -((6 * n * apsidals[i]) / tidals[i]) * (mj / mi) * pow((ri / a), 5) * pow((a / dist), 8);

        // Build vector 2 - this is the parenthesis term in (4)
        double temp_x = rhat_cross_v_x - dist * sx;
        double temp_y = rhat_cross_v_y - dist * sy;
        double temp_z = rhat_cross_v_z - dist * sz;

        double vec2_x = temp_y * dz_hat - temp_z * dy_hat;
        double vec2_y = temp_z * dx_hat - temp_x * dz_hat;
        double vec2_z = temp_x * dy_hat - temp_y * dx_hat;

        double tx = prefactor * (vec1_x + vec2_x);
        double ty = prefactor * (vec1_y + vec2_y);
        double tz = prefactor * (vec1_z + vec2_z);

        // total forces for convenience
        double totx = qx + tx;
        double toty = qy + ty;
        double totz = qz + tz;

      /*
        double uqx1 = qx1 / sqrt(qx1 * qx1 + qy1 * qy1 + qz1 * qz1);
        double uqy1 = qy1 / sqrt(qx1 * qx1 + qy1 * qy1 + qz1 * qz1);
        double uqz1 = qy1 / sqrt(qx1 * qx1 + qy1 * qy1 + qz1 * qz1);

        double utx1 = tx1 / sqrt(tx1 * tx1 + ty1 * ty1 + tz1 * tz1);
        double uty1 = ty1 / sqrt(tx1 * tx1 + ty1 * ty1 + tz1 * tz1);
        double utz1 = ty1 / sqrt(tx1 * tx1 + ty1 * ty1 + tz1 * tz1);

        double uqx2 = qx2 / sqrt(qx2 * qx2 + qy2 * qy2 + qz2 * qz2);
        double uqy2 = qy2 / sqrt(qx2 * qx2 + qy2 * qy2 + qz2 * qz2);
        double uqz2 = qz2 / sqrt(qx2 * qx2 + qy2 * qy2 + qz2 * qz2);

        double utx2 = tx2 / sqrt(tx2 * tx2 + ty2 * ty2 + tz2 * tz2);
        double uty2 = ty2 / sqrt(tx2 * tx2 + ty2 * ty2 + tz2 * tz2);
        double utz2 = ty2 / sqrt(tx2 * tx2 + ty2 * ty2 + tz2 * tz2);
      */

      //  if (sim->t < 0.001){
          //printf("%.20f\n", prefactor_2);
          //printf("Unit quad moment 1: %f %f %f\n", uqx1, uqy1, uqz1);
          //printf("Unit quad moment 2: %f %f %f\n", uqx2, uqy2, uqz2);
          //printf("Unit tidal moment 1: %f %f %f\n", utx1, uty1, utz1);
          //printf("Unit tidal moment 2: %f %f %f\n", utx2, uty2, utz2);
          /*printf("time: %f\n", sim->t);
          printf("vel: %f %f %f\n", p1->vx, p1->vy, p1->vz);
          printf("a_1: %f %f %f\n", (m2 / (m1 + m2)) * (totx_1 + totx_2), (m2 / (m1 + m2)) * (toty_1 + toty_2), (m2 / (m1 + m2)) * (totz_1 + totz_2));
          printf("p_1: %.10f %.10f %.10f\n", p1->x, p1->y, p1->z);
          printf("p_2: %.10f %.10f %.10f\n", p2->x, p2->y, p2->z);
          printf("r: %f %f %f %.10f\n", dx, dy, dz, sqrt(dx*dx+dy*dy+dz*dz));*/
      //  }

        // Force on planet 1
        pi->ax -= ((mj / (mi + mj)) * totx);
        pi->ay -= ((mj / (mi + mj)) * toty);
        pi->az -= ((mj / (mi + mj)) * totz);

        // Force on planet 2
        pj->ax += ((mi / (mi + mj)) * totx);
        pj->ay += ((mi / (mi + mj)) * toty);
        pj->az += ((mi / (mi + mj)) * totz);

        // printf("Quad Force 1: %f\n", sqrt(qx1 * qx1 + qy1 * qy1 + qz1 * qz1));
        // printf("Quad Force 2: %f\n", sqrt(qx2 * qx2 + qy2 * qy2 + qz2 * qz2));

        // printf("Tidal Force 1: %f\n", sqrt(tx1 * tx1 + ty1 * ty1 + tz1 * tz1));
        // printf("Tidal Force 2: %f\n", sqrt(tx2 * tx2 + ty2 * ty2 + tz2 * tz2));

        // EOMs for spin vector
        double mu_ij = -(pi->m * pj->m) / ((pi->m + pj->m));

        // Now take the cross product
        dummy->ax += ((dy * totz - dz * toty) * (mu_ij / mois[i]));
        dummy->ay += ((dz * totx - dx * totz) * (mu_ij / mois[i]));
        dummy->az += ((dx * toty - dy * totx) * (mu_ij / mois[i]));

        //printf("t = %f, x = %.10f, y = %.10f, z = %.10f\n", sim->t, dummy->vx, dummy->vy, dummy->vz);

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
    }
  }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m r", 1., 0.00465);                // Central object

    reb_add_fmt(r, "m a e r inc", 5. * 3.0e-6, 0.175, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.)); // Planet 1
    reb_add_fmt(r, "m a e r inc", 5.5 * 3.0e-6, 0.235, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.));
    reb_move_to_com(r);

    //printf("%f, %f, %f\n", r->particles[1].vx, r->particles[1].vy, r->particles[1].vz);

    r->N_active = 3;
    // Tot rot angular momentum is 0.8 * pi * m * r^2 / trot
    // take solar obliquity = 6 degrees
    double solar_spin_period = 20 / (2 * PI * 365);
    double solar_spin = (2 * PI) / solar_spin_period;
    double obliquity_solar = 0. * (PI / 180.);
    double res_angle_solar = 60. * (PI / 180.);
    double spin_x_solar = solar_spin * sin(obliquity_solar) * cos(res_angle_solar);
    double spin_y_solar = solar_spin * sin(obliquity_solar) * sin(res_angle_solar);
    double spin_z_solar = solar_spin * cos(obliquity_solar);
    reb_add_fmt(r, "m vx vy vz", 0., 0., spin_x_solar, spin_y_solar, spin_z_solar);         // Dummy particle
                                     // Mass is not relevant

    double spin_period_1 = 5 / (365 * 2 * PI); // 5 days in reb years
    double spin_1 = (2 * PI) / spin_period_1; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_1 = 1. * (PI / 180.);
    double res_angle_1 = 120. * (PI / 180.);
    double spin_x1 = spin_1 * sin(obliquity_1) * cos(res_angle_1);
    double spin_y1 = spin_1 * sin(obliquity_1) * sin(res_angle_1);
    double spin_z1 = spin_1 * cos(obliquity_1);
    reb_add_fmt(r, "m vx vy vz", 0., spin_x1, spin_y1, spin_z1);
    // printf("%.15f, %.15f, %.15f\n", spin_y, spin_z, sqrt(spin_y * spin_y + spin_z * spin_z));
    // printf("%f\n", r->particles[3].vx);

    double spin_period_2 = 3 / (365 * 2 * PI); // 5 days in reb years
    double spin_2 = (2 * PI) / spin_period_2; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_2 = 1. * (PI / 180.);
    double res_angle_2 = 180. * (PI / 180.);
    double spin_x2 = spin_2 * sin(obliquity_2) * cos(res_angle_2);
    double spin_y2 = spin_2 * sin(obliquity_2) * sin(res_angle_2);
    double spin_z2 = spin_2 * cos(obliquity_2);
    reb_add_fmt(r, "m vx vy vz", 0., spin_x2, spin_y2, spin_z2);

    r->additional_forces = additional_forces;
    r->integrator = REB_INTEGRATOR_IAS15;
    r->force_is_velocity_dependent = 1;


   FILE* f = fopen("out.txt","w");
    for (int i=0; i<500; i++){

        struct reb_particle sun = r->particles[0];
        struct reb_particle p1 = r->particles[1];
        struct reb_particle p2 = r->particles[2];
        struct reb_particle dummy1 = r->particles[3];
        struct reb_particle dummy2 = r->particles[4];
        struct reb_particle dummy3 = r->particles[5];

        struct reb_orbit o1 = reb_tools_particle_to_orbit(r->G, p1, sun);
        double a1 = o1.a;
        struct reb_orbit o2 = reb_tools_particle_to_orbit(r->G, p2, sun);
        double a2 = o2.a;

        double mag1 = sqrt(dummy2.vx * dummy2.vx + dummy2.vy * dummy2.vy + dummy2.vz * dummy2.vz);
        double ob1 = acos(dummy2.vz / mag1) * (180 / PI);
        double mag2 = sqrt(dummy3.vx * dummy3.vx + dummy3.vy * dummy3.vy + dummy3.vz * dummy3.vz);
        double ob2 = acos(dummy3.vz / mag2) * (180 / PI);

        //printf("torb=%.10f \t t=%.4f\t planet = %6.3f \t spin axis = %.10f %.10f %.10f \t ob = %.3f \t ang = %.10f\n", torb, r->t, o.a, dummy2.vx, dummy2.vy, dummy2.vz, obliquity, tot_angular_momentum);
        printf("t=%f\t a1 = %.6f\t a2 = %.6f\t ob1 = %.10f\t ob2 = %.10f\n", r->t / (2 * PI), a1, a2, ob1, ob2);
        //printf("%f, %f, %f, %f, %f, %f\n", sun.vx, sun.vy, sun.vz, p.vx, p.vy, p.vz);
        // fprintf(f, "%.4f %.10f %.10f %.10f %.10f\n", r->t, o.a, torb, o.e, obliquity);

        reb_integrate(r,r->t+(10000 * 2 * PI));
    }
    fclose(f);

    reb_free_simulation(r);
}
