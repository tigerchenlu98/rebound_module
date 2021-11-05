#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_tidal_damping(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
  const int _N_real = sim->N - sim->N_var; // Still not sure what a variational particle is

  for (int i = 0; i < _N_real; i++){ // This loop goes through each particle ACTED ON
    struct reb_particle* p1 = &sim->particles[i];
    const double* const k1 = rebx_get_param(sim->extras, p1->ap, "k");
    const double* const q1 = rebx_get_param(sim->extras, p1->ap, "Q");
    const double* const sx = rebx_get_param(sim->extras, p1->ap, "sx");
    const double* const sy = rebx_get_param(sim->extras, p1->ap, "sy");
    const double* const sz = rebx_get_param(sim->extras, p1->ap, "sz");

    // Particle only feels quadrupole acceleration if apsidal constant and love number are defined
    if (k1 == NULL || sx == NULL || sy == NULL || sz == NULL) continue;

    const double s1 = p1->r;
    const double m1 = p1->m;

    for (int j = 0; j < _N_real; j++){ // This loop calculates the contribution of particle j on particle i
      if (i == j) continue;

      struct reb_particle* const p2 = &sim->particles[j];
      const double m2 = p2->m;

      // Unclear on generalization to N-bodies - for now I consider effect of all bodies
      double dx = p1->x - p2->x;
      double dy = p1->y - p2->y;
      double dz = p1->z - p2->z;
      double dist = sqrt(dx * dx + dy * dy + dz * dz);

      // Unit distance vector
      double dx_hat = dx / dist;
      double dy_hat = dy / dist;
      double dz_hat = dz / dist;

      // mean motion and semimajor axis
      struct reb_orbit o =  reb_tools_particle_to_orbit(sim->G, sim->particles[j], sim->particles[i]);
      double n = o.n;
      double a = o.a;

      // Relative velocities of two bodies
      double dvx = p1->vx - p2->vx;
      double dvy = p1->vy - p2->vy;
      double dvz = p1->vz - p2->vz;

      // relevant dot and cross products
      double rhat_dot_vel = dx_hat * dvx + dy_hat * dvy + dz_hat * dvz;
      double rhat_cross_v_x = dy_hat * dvz - dz_hat * dvy;
      double rhat_cross_v_y = dz_hat * dvx - dx_hat * dvz;
      double rhat_cross_v_z = dx_hat * dvy - dy_hat * dvx;

      double prefactor = -(6 * n * k1 / q1) * (m2 / m1) * pow((s1 / a), 5) * pow((a / dist));

      double vec1_x = 3 * rhat_dot_vel * dx_hat;
      double vec1_y = 3 * rhat_dot_vel * dy_hat;
      double vec1_z = 3 * rhat_dot_vel * dz_hat;

      // Build vector 2
      double temp_x = rhat_cross_v_x - dist * sx;
      double temp_y = rhat_cross_v_y - dist * sy;
      double temp_z = rhat_cross_v_z - dist * sz;

      double vec2_x = temp_y * dz_hat - temp_z * dy_hat;
      double vec2_y = temp_z * dx_hat - temp_x * dz_hat;
      double vec2_z = temp_x * dy_hat - temp_y * dx_hat;

      // Force on planet 1
      p1->ax += (m2 / (m1 + m2)) * prefactor * (vec1_x + vec2_x);
      p1->ay += (m2 / (m1 + m2)) * prefactor * (vec1_y + vec2_y);
      p1->az += (m2 / (m1 + m2)) * prefactor * (vec1_z + vec2_z);

      // Force on planet 2
      p1->ax += (m1 / (m1 + m2)) * prefactor * (vec1_x + vec2_x);
      p1->ay += (m1 / (m1 + m2)) * prefactor * (vec1_y + vec2_y);
      p1->az += (m1 / (m1 + m2)) * prefactor * (vec1_z + vec2_z);
    }
  }
}
