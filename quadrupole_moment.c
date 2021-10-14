#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_quadrupole_moment(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
  const int _N_real = sim->N - sim->N_var; // Still not sure what a variational particle is

  for (int i = 0; i < _N_real; i++){ // This loop goes through each particle ACTED ON
    struct reb_particle* p1 = &sim->particles[i];
    const double* const k1 = rebx_get_param(sim->extras, p1->ap, "k");
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
      const double mu = sim->G * p2->m;

      // Unclear on generalization to N-bodies - for now I consider effect of all bodies
      double dx = p1->x - p2->x;
      double dy = p1->y - p2->y;
      double dz = p1->z - p2->z;
      double norm = sqrt(dx * dx + dy * dy + dz * dz);

      double dx_hat = dx / norm;
      double dy_hat = dy / norm;
      double dz_hat = dz / norm;

      double prefactor = pow(s1, 5) * (1 + m1 / m2) * (*k1) / pow(norm, 4);
      double o_dot_r = (*sx) * dx_hat + (*sy) * dy_hat + (*sz) * dz_hat;
      double o_squared = (*sx) * (*sx) + (*sy) * (*sy) + (*sz) * (*sz);

      p1->ax += prefactor * ((5 * pow(o_dot_r, 2) - o_squared - (12. * mu / pow(norm, 3))) * dx_hat - (2 * o_dot_r * (*sx)));
      p1->ay += prefactor * ((5 * pow(o_dot_r, 2) - o_squared - (12. * mu / pow(norm, 3))) * dy_hat - (2 * o_dot_r * (*sy)));
      p1->az += prefactor * ((5 * pow(o_dot_r, 2) - o_squared - (12. * mu / pow(norm, 3))) * dz_hat - (2 * o_dot_r * (*sz)));
    }
  }
}
