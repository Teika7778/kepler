#pragma once

struct kepler_orbit_denorm {
    double a;     // Semi-major axis (arcsec)
    double e;     // Eccentricity
    double w;     // Argument of periapsis (deg)
    double omega; // Longitude of ascending node (deg)
    double i;     // Inclination (deg)
    double T0;    // Time of periapsis t_0
    double t0;
} typedef kepler_orbit_denorm;


struct kepler_orbit {
    double a;     // Semi-major axis (m)
    double e;     // Eccentricity (-)
    double w;     // Argument of periapsis (rad)
    double omega; // Longitude of ascending node (rad)
    double i;     // Inclination (rad)
    double M0;    // Mean anomaly at t_0 (rad)
    double t0;    // Epoch (JD)
} typedef kepler_orbit;
