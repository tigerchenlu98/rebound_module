// EXTRA TLu File 1/11/22

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static struct reb_vec3d rebx_calculate_spin_magnitude(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* p, struct reb_particle* source, struct reb_particle* dummy){
  // This calculates ALL forces on the ith particle
  const int _N_real = sim->N - sim->N_var; // Still not sure what a variational particle is

  // Extract physical constants
  const double ri = p->r;
  const double mi = p->m;
  const double moi = 0.4 * mi * ri * ri;

  double ki = 0.0;
  double qi = 0.0;
  double sx = dummy->vx;
  double sy = dummy->vy;
  double sz = dummy->vz;

  // Extract rebx constants
  const double* const ki_ptr = rebx_get_param(sim->extras, p->ap, "k");
  const double* const qi_ptr = rebx_get_param(sim->extras, p->ap, "q");

  if(ki_ptr != NULL){
      ki = *ki_ptr;
  }

  if(qi_ptr != NULL){
      qi = *qi_ptr;
  }

  struct reb_vec3d tot_force = {0};

  // This calculates the effect of the jth particle
  for (int j = 0; j < 3; j++){ // MAGIC NUMBER
    struct reb_particle* pj = &(sim->particles[j]);
    double rj = pj->r;
    double mj = pj->m;

    if (p != pj){
      // Need to calculate mutual orbit first
      // This is a bit tricky b/c sometimes the star acts as the perturber: need to check against this
      double n;
      double a;
      if (j == 0) { // If the star is the perturber
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *source, *p);
        n = o.n;
        a = o.a;
      }
      else {
        struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *source, *pj);
        n = o.n;
        a = o.a;
      }

      // Forces which i is the PERTURBED by j
      // distance vector FROM j TO i
      double dx = pj->x - p->x;
      double dy = pj->y - p->y;
      double dz = pj->z - p->z;
      double dist = sqrt(dx * dx + dy * dy + dz * dz);
      // Unit distance vector
      double dx_hat = dx / dist;
      double dy_hat = dy / dist;
      double dz_hat = dz / dist;

      // acceleration on i due to quadrupole of j
      double quad_prefactor = (pow(ri, 5) * (1 + (mi / mj)) * ki) / pow(dist, 4);
      double omega_dot_rhat = sx * dx_hat + sy * dy_hat + sz * dz_hat;
      double omega_squared = sx * sx + sy * sy + sz * sz;

      double qx = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dx_hat - (2 * omega_dot_rhat * sx));
      double qy = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dy_hat - (2 * omega_dot_rhat * sy));
      double qz = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dz_hat - (2 * omega_dot_rhat * sz));

      // Velocity vecto: j to i
      double dvx = pj->vx - p->vx;
      double dvy = pj->vy - p->vy;
      double dvz = pj->vz - p->vz;

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
      double prefactor = -((6 * n * ki) / qi) * (mj / mi) * pow((ri / a), 5) * pow((a / dist), 8);

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

      double mu_ij = -(mi * mj) / ((mi + mj));

      // Now take the cross product
      tot_force.x += ((dy * totz - dz * toty) * (mu_ij / moi));
      tot_force.y += ((dz * totx - dx * totz) * (mu_ij / moi));
      tot_force.z += ((dx * toty - dy * totx) * (mu_ij / moi));
    }
  }
  return tot_force;
}

// Spin vector using dummy particles
// This is probably a really stupid implementation but lets see if it works
void rebx_spin_eom(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    for (int i = 0; i < 6; i++){
      struct reb_particle* const p = &sim->particles[i];
      const int* const dummy_ptr = rebx_get_param(sim->extras, p->ap, "dummy");
      if (*dummy_ptr == 1 && dummy_ptr != NULL){
        // Required for convergence criteria only
        struct reb_particle* dummy = &sim->particles[i];
        dummy->x =  1.;
        dummy->y =  1.;
        dummy->z =  1.;

        // Overwrite gravity acceleration
        dummy->ax = 0.;
        dummy->ay = 0.;
        dummy->az = 0.;

        // This IS a dummy particle, now let's go back and find the associated real particle
        struct reb_particle* real = &sim->particles[i - 3]; // MAGIC NUMBER
        struct reb_vec3d tot_force = rebx_calculate_spin_magnitude(sim, force, real, &sim->particles[0], dummy);

        // Set relevant update_accelerations
        dummy->ax = tot_force.x;
        dummy->ay = tot_force.y;
        dummy->az = tot_force.z;
      }
    }
}
