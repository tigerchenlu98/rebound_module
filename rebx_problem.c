#include "rebound.h"
#include "reboundx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m r", 1., 0.00465);
    r->particles[0].hash = reb_hash("star");                // Central object
    reb_add_fmt(r, "m a e r inc", 5. * 3.0e-6, 0.175, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.));
    r->particles[1].hash = reb_hash("p1");   // Planet 1
    reb_add_fmt(r, "m a e r inc", 5.5 * 3.0e-6, 0.235, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.));
    r->particles[2].hash = reb_hash("p2");   // Planet 2
    reb_move_to_com(r);

    r->N_active = 3;
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->dt = 1;
    // Add spin status of planets and stars by adding dummy particles
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


    double spin_period_2 = 3 / (365 * 2 * PI); // 5 days in reb years
    double spin_2 = (2 * PI) / spin_period_2; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_2 = 1. * (PI / 180.);
    double res_angle_2 = 180. * (PI / 180.);
    double spin_x2 = spin_2 * sin(obliquity_2) * cos(res_angle_2);
    double spin_y2 = spin_2 * sin(obliquity_2) * sin(res_angle_2);
    double spin_z2 = spin_2 * cos(obliquity_2);
    reb_add_fmt(r, "m vx vy vz", 0., spin_x2, spin_y2, spin_z2);

    // Add migration via reboundx
    struct rebx_extras* rebx = rebx_attach(r);
    struct rebx_operator* mo = rebx_load_operator(rebx, "modify_orbits_direct");
    rebx_add_operator(rebx, mo);

    double tau_a1 = 5e6 * 2 * PI;
    double tau_a2 = tau_a1 / 1.1;
    rebx_set_param_double(rebx, &r->particles[1].ap, "tau_a", -tau_a1);
    rebx_set_param_double(rebx, &r->particles[2].ap, "tau_a", -tau_a2);

    // Add quad forces and spin EOMs
    struct rebx_force* qt = rebx_load_force(rebx, "quad_and_tidal_forces");
    //struct rebx_force* spin = rebx_load_force(rebx, "spin_eom");
    rebx_add_force(rebx, qt);
    //rebx_add_operator(rebx, spin);

    rebx_set_param_double(rebx, &r->particles[0].ap, "k", 0.035);
    rebx_set_param_double(rebx, &r->particles[0].ap, "Q", 100000.);

    //rebx_set_param_double(rebx, &r->particles[0].ap, "dummy", 0);

    rebx_set_param_double(rebx, &r->particles[1].ap, "k", 0.2);
    rebx_set_param_double(rebx, &r->particles[1].ap, "Q", 10000.);
    //rebx_set_param_double(rebx, &r->particles[1].ap, "dummy", 0);

    rebx_set_param_double(rebx, &r->particles[2].ap, "k", 0.2);
    rebx_set_param_double(rebx, &r->particles[2].ap, "Q", 10000.);
    //rebx_set_param_double(rebx, &r->particles[2].ap, "dummy", 0);

    //rebx_set_param_double(rebx, &r->particles[3].ap, "dummy", 1);
    //rebx_set_param_double(rebx, &r->particles[4].ap, "dummy", 1);
    //rebx_set_param_double(rebx, &r->particles[5].ap, "dummy", 1);

   FILE* f = fopen("out.txt","w");

   double mig_timescale = 2e6 * 2 * PI;
    for (int i=0; i<500; i++){

        if (r->t >= mig_timescale){
          rebx_set_param_double(rebx, &r->particles[1].ap, "tau_a", INFINITY);
          rebx_set_param_double(rebx, &r->particles[2].ap, "tau_a", INFINITY);
        }

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
        printf("t=%f\t a1 = %.3f\t a2 = %.3f\t ob1 = %.3f\t ob2 = %.3f\n", r->t / (2 * PI), a1, a2, ob1, ob2);
        //printf("%f, %f, %f, %f, %f, %f\n", sun.vx, sun.vy, sun.vz, p.vx, p.vy, p.vz);
        // fprintf(f, "%.4f %.10f %.10f %.10f %.10f\n", r->t, o.a, torb, o.e, obliquity);

        reb_integrate(r,r->t+(1 * 2 * PI));
    }
    fclose(f);

    rebx_free(rebx);
    reb_free_simulation(r);
}
