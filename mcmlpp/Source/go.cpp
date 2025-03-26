/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright M.H.J. Lam, 2025.
 *	Launch, move, and record photon weight.
 ****/


#include "mcml.hpp"

#include <limits>
#include <numbers>


/*******************************************************************************
 *  Choose a new direction for photon propagation by sampling:
 *	1. The polar (deflection) angle θ, measuring from downwards z-axis
 *	2. The azimuthal angle Ψ, measuring rotation around the z-axis in the xy-plane
 * 
 *  θ range:    0 - π  (0 to 180 degrees)    [sin(θ) range:  0 to 1]
 *  Ψ range:    0 - 2π (0 to 360 degrees)    [cos(Ψ) range: -1 to 1]
 * 
 *  Since cos²(θ) + sin²(θ) = 1, and if sin(θ) is known; cos(θ) = sqrt(1-sin(θ)²)
 * 
 *  In the Henyey-Greenstein phase function, g is the asymmetry parameter and 
 *  typically takes values in the range [-1, 1]:
 *      g = 0: isotropic scattering (uniform in all directions).
 *      g > 0: forward scattering (continue in the same direction).
 *      g < 0: backward scattering (scatter in the opposite direction).
 * 
 ****/
void Spin(Photon& photon, double g) // g: Henyey-Greenstein asymmetry parameter
{
    // Cosine and sine of Ψ
    double ux = photon.direction.x;
    double uy = photon.direction.y;
    double uz = photon.direction.z;

    // Sample θ
    // Generate a random variable based on the Henyey-Greenstein distribution
    double theta = (1 - g * g) / (1 - g + 2 * g * g_rand.next());
    double cos_theta = (g == 0.0) ? 2 * g_rand.next() - 1 : (1 + g * g - theta * theta) / (2 * g);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    // Sample Ψ
    double psi = 2.0 * std::numbers::pi * g_rand.next();
    double cos_psi = std::cos(psi);
    double sin_psi = (psi < std::numbers::pi) ? std::sqrt(1.0 - cos_psi * cos_psi) : -std::sqrt(1.0 - cos_psi * cos_psi);

    // Update directional cosines

    // Close to perpendicular
    if (1.0 - std::abs(uz) <= COS_0_TOLERANCE) {
        photon.direction.x = sin_theta * cos_psi;
        photon.direction.y = sin_theta * sin_psi;
        photon.direction.z = cos_theta * std::copysign(1.0, uz);
    }
    else {
        double temp = std::sqrt(1.0 - uz * uz);
        photon.direction.x = sin_theta * (ux * uz * cos_psi - uy * sin_psi) / temp + ux * cos_theta;
        photon.direction.y = sin_theta * (uy * uz * cos_psi + ux * sin_psi) / temp + uy * cos_theta;
        photon.direction.z = -sin_theta * cos_psi * temp + uz * cos_theta;
    }
}

/*******************************************************************************
 *	Initialize a photon packet.
 *
 *	If an isotropic source is launched inside a glass layer, check whether the 
 *  photon will be total-internally reflected. If so, kill the photon to avoid 
 *  infinite travelling inside the glass layer.
 ****/
void LaunchPhoton(RunParams& params, Tracer& tracer, Photon& photon)
{
    photon.weight = 1.0 - tracer.R_spec;
    photon.alive = 1;
    photon.current_layer = (params.source.layer_index != 0) ? params.source.layer_index : 1;
    photon.step_size = 0;
    photon.step_size_left = 0;
    photon.num_scatters = 0;
    photon.flight_time = 0;

    photon.position.x = 0.0;
    photon.position.y = 0.0;
    photon.position.z = params.source.z;
    photon.direction.x = 0.0;
    photon.direction.y = 0.0;
    photon.direction.z = 1.0;

    tracer.A_i = 0.0;
    tracer.Tb_i = 0.0;
    tracer.T_i = 0.0;
    tracer.Rb_i = 0.0;
    tracer.R_i = 0.0;

    if (params.source.beam  == BeamType::Isotropic) {
        Layer layer = params.layers[photon.current_layer];

        // To avoid scoring into reflectance or transmittance
        photon.num_scatters++;

        // Isotropically scatter the photon
        Spin(photon, 0.0);

        // Glass layer.
        if (layer.mu_a == 0.0 && layer.mu_s == 0.0) {
            // Total internal reflection
            if (std::abs(photon.direction.z) <= layer.cos_theta_c0 && 
                std::abs(photon.direction.z) <= layer.cos_theta_c1) {
                photon.alive = 0;
            }
        }
    }
}

/*******************************************************************************
 *	Move the photon away in the current layer.
 ****/
void Hop(Photon& photon, double dist, double n)
{
    photon.position.x += dist * photon.direction.x;
    photon.position.y += dist * photon.direction.y;
    photon.position.z += dist * photon.direction.z;
    photon.flight_time += dist * n * SPEED_OF_LIGHT_INV;
}

/*******************************************************************************
 *	Pick a step size in dimensionless unit for a photon packet. If step_size is 
 *  zero, make a new step size: -log(rnd). Otherwise, finish the leftover.
 ****/
void SetStepSize(Photon& photon)
{
    // Make a new step
    if (photon.step_size == 0.0) {
        double rnd; // Avoid zero
        while ((rnd = g_rand.next()) <= 0.0);

        photon.step_size = -std::log(rnd);
    }
}

/*******************************************************************************
 *	Return the distance between the photon position to the boundary along the 
 *  photon direction.
 ****/
double PathToBoundary(Photon& photon, RunParams& params)
{
    // Length to boundary
    short layer = photon.current_layer;
    double uz = photon.direction.z;

    // Distance to the boundary
    double path = std::numeric_limits<double>::max(); // infinity

    // Path > 0
    if (uz > 0.0) {
        path = (params.layers[layer].z1 - photon.position.z) / uz;
    }
    // Path > 0
    else if (uz < 0.0) {
        path = (params.layers[layer].z0 - photon.position.z) / uz;
    }

    return path;
}

/*******************************************************************************
 *	Drop photon weight inside the tissue. The photon is assumed alive.
 *	The weight drop is dw = w * μa / (μa + μs).
 *	The dropped weight is assigned to the absorption array elements.
 ****/
void Drop(RunParams& params, Photon& photon, Tracer& tracer)
{
    // Absorbed weight
    double x = photon.position.x;
    double y = photon.position.y;

    short layer = photon.current_layer;

    // Update photon weight
    double mua = params.layers[layer].mu_a;
    double mus = params.layers[layer].mu_s;
    double dwa = photon.weight * mua / (mua + mus);
    photon.weight -= dwa;

    // Compute array indices
    std::size_t iz, ir, it;

    if (params.A_rzt || params.A_zt || params.A_z || params.A_rz) {
        if (photon.position.z >= params.max_z) {
            iz = params.num_z - 1;
        }
        else {
            iz = static_cast<short>(photon.position.z / params.grid_z);
        }
    }

    if (params.A_rzt || params.A_zt || params.A_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_t - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_t));
        }
    }

    // Assign dwa to an absorption array element

    // Scattered
    if (photon.num_scatters) {
        if (params.A_rzt || params.A_rz) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = static_cast<short>(temp / params.grid_r);
            }
        }
        if (params.A_rzt) {
            tracer.A_rzt[ir][iz][it] += dwa;
        }
        if (params.A_rz) {
            tracer.A_rz[ir][iz] += dwa;
        }
    }

    // Ballistic
    else {
        if (params.A_rzt) {
            tracer.Ab_zt[iz][it] += dwa;
        }
        if (params.A_rz) {
            tracer.Ab_z[iz] += dwa;
        }
    }

    if (params.A_zt) {
        tracer.A_zt[iz][it] += dwa;
    }
    if (params.A_z) {
        tracer.A_z[iz] += dwa;
    }

    if (params.A_t) {
        tracer.A_t[it] += dwa;
    }
    tracer.A_i += dwa;
}

/*******************************************************************************
 *  Determine photon survival via roulette when photon weight becomes too small.
 ****/
void Roulette(Photon& photon)
{
    // Already dead
    if (photon.weight == 0.0) {
        photon.alive = 0;
    }

    // Survived the roulette
    else if (g_rand.next() < ROULETTE_SURVIVAL) {
        photon.weight /= ROULETTE_SURVIVAL;
    }
    else {
        photon.alive = 0;
    }
}

/*******************************************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle ai is positive, and the case 
 *  when the angle is greater than the critical angle is ruled out.
 *
 *	cos_ai: cosine of the incident angle ai.
 *	cos_at: pointer to the cosine of the transmission angle at.
 *
 * 	Avoid trigonometric function operations as much as possible, because they 
 *  are computationally intensive.
 *
 *  eta_i:  incident refractive index
 *  eta_t:  transmit refractive index
 *  cos_ai: cosine of angle ai
 ****/
double FresnelReflectance(double eta_i, double eta_t, double cos_ai, double& cos_at)
{
    // Matched boundary
    if (eta_i == eta_t) {
        cos_at = cos_ai;
        return 0.0;
    }
    // Normal incidence
    else if (1.0 - cos_ai <= COS_0_TOLERANCE) {
        cos_at = cos_ai;
        double r = (eta_t - eta_i) / (eta_t + eta_i);
        return r * r;
    }
    // Very slant incidence
    else if (cos_ai < COS_90_TOLERANCE) {
        cos_at = 0.0;
        return 1.0;
    }
    // General incidence
    else {
        // Sine of the angles ai & at
        double sin_ai = std::sqrt(1 - cos_ai * cos_ai);
        double sin_at = eta_i * sin_ai / eta_t;

        // Total internal reflection
        if (sin_at >= 1.0) {
            cos_at = 0.0;
            return 1.0;
        }

        // Cosine of at
        cos_at = std::sqrt(1 - sin_at * sin_at);

        // Cosines of (ai + at) and (ai - at)
        double cos_ap = cos_ai * cos_at - sin_ai * sin_at;  // c+ = cc - ss.
        double cos_am = cos_ai * cos_at + sin_ai * sin_at;  // c- = cc + ss.

        // Sines of (ai + at) and (ai - at)
        double sin_ap = sin_ai * cos_at + cos_ai * sin_at;	// s+ = sc + cs.
        double sin_am = sin_ai * cos_at - cos_ai * sin_at;	// s- = sc - cs.

        // Arranged for speed
        return (0.5 * sin_am * sin_am) * 
               (cos_am * cos_am + cos_ap * cos_ap) / 
               (sin_ap * sin_ap * cos_am * cos_am);
    }
    return 0.0;
}

/*******************************************************************************
 *	Record photon weight exiting the first layer (uz < 0) to the reflectance 
 *  array and update its weight.
 ****/
void RecordReflectance(RunParams& params, Photon& photon, Tracer& tracer, double reflectance)
{
    double x = photon.position.x;
    double y = photon.position.y;

    // Index to r & angle
    std::size_t ir, ia, it;

    if (params.R_rat || params.R_at || params.R_rt || params.R_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_t - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_t));
        }
    }

    // Scattered
    if (photon.num_scatters) {
        if (params.R_rat || params.R_rt || params.R_ra || params.R_r) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = static_cast<short>(temp / params.grid_r);
            }
        }
        if (params.R_rat || params.R_at || params.R_ra || params.R_a) {
            if ((ia = static_cast<short>(std::acos(-photon.direction.z) / params.grid_a) > params.num_a - 1)) {
                ia = params.num_a - 1;
            }
        }

        // Assign photon weight to the reflection array element
        if (params.R_rat) {
            tracer.R_rat[ir][ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.R_ra) {
            tracer.R_ra[ir][ia] += photon.weight * (1.0 - reflectance);
        }

        if (params.R_rt) {
            tracer.R_rt[ir][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.R_r) {
            tracer.R_r[ir] += photon.weight * (1.0 - reflectance);
        }

        if (params.R_at) {
            tracer.R_at[ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.R_a) {
            tracer.R_a[ia] += photon.weight * (1.0 - reflectance);
        }
        if (params.R_t) {
            tracer.R_t[it] += photon.weight * (1.0 - reflectance);
        }
        tracer.R_i += photon.weight * (1.0 - reflectance);
    }

    // Ballistic reflection
    else {
        tracer.Rb_i += photon.weight * (1.0 - reflectance);
    }

    photon.weight *= reflectance;
}

/*******************************************************************************
 *	Record the photon weight exiting the last layer (uz > 0), no matter whether 
 *  the layer is glass or not, to the transmittance array.
 ****/
void RecordTransmittance(RunParams& params, Photon& photon, Tracer& tracer, double reflectance)
{
    double x = photon.position.x;
    double y = photon.position.y;

    // Index to r & angle
    std::size_t ir, ia, it;
    if (params.T_rat || params.T_at || params.T_rt || params.T_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_t - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_t));
        }
    }

    // Scattered
    if (photon.num_scatters) {
        if (params.T_rat || params.T_rt || params.T_ra || params.T_r) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = (short)(temp / params.grid_r);
            }
        }
        if (params.T_rat || params.T_at || params.T_ra || params.T_a) {
            if ((ia = static_cast<short>(std::acos(photon.direction.z) / params.grid_a) > params.num_a - 1)) {
                ia = params.num_a - 1;
            }
        }

        // Assign photon weight to the transmittance array element
        if (params.T_rat) {
            tracer.T_rat[ir][ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.T_ra) {
            tracer.T_ra[ir][ia] += photon.weight * (1.0 - reflectance);
        }

        if (params.T_rt) {
            tracer.T_rt[ir][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.T_r) {
            tracer.T_r[ir] += photon.weight * (1.0 - reflectance);
        }

        if (params.T_at) {
            tracer.T_at[ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (params.T_a) {
            tracer.T_a[ia] += photon.weight * (1.0 - reflectance);
        }

        if (params.T_t) {
            tracer.T_t[it] += photon.weight * (1.0 - reflectance);
        }
        tracer.T_i += photon.weight * (1.0 - reflectance);
    }

    // Collimated
    else {
        tracer.Tb_i += photon.weight * (1.0 - reflectance);
    }

    photon.weight *= reflectance;
}

/*******************************************************************************
 *	Decide whether the photon will be transmitted or reflected on the upper
 *  boundary (uz < 0) of the current layer.
 *
 *	If current_layer is the first layer, the photon packet will be partially 
 *  transmitted and partially reflected if PARTIAL_REFLECTION active, or 
 *  the photon packet will be either transmitted or reflected determined 
 *  statistically if PARTIAL_REFLECTION is inactive.
 *
 *	Record the transmitted photon weight as reflection.
 *
 *	If the current_layer is not the first layer and the photon packet is 
 *  transmitted, move the photon to the previous layer.
 *
 *	Update photon parameters.
 ****/
void CrossUp(RunParams& params, Photon& photon, Tracer& tracer)
{
    // Z directional cosine
    double uz = photon.direction.z;

    // Cosines of transmission alpha.
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double eta_i = params.layers[layer].eta;
    double eta_t = params.layers[layer - 1].eta;

    // Get reflectance; if 1.0, then total internal reflection
    double refl = (-uz <= params.layers[layer].cos_theta_c0) ? 1.0 : FresnelReflectance(eta_i, eta_t, -uz, uz1);

    if (PARTIAL_REFLECTION) {
        // Partially transmitted
        if (layer == 1 && refl < 1.0) {
            photon.direction.z = -uz1; // Escaped photon
            RecordReflectance(params, photon, tracer, refl);
            photon.direction.z = -uz; // Reflected photon
        }
        // Transmitted to current_layer - 1
        else if (g_rand.next() > refl) {
            photon.current_layer--;
            photon.direction.x *= eta_i / eta_t;
            photon.direction.y *= eta_i / eta_t;
            photon.direction.z = -uz1;
        }
        // Reflected
        else {
            photon.direction.z = -uz;
        }
    }
    else {
        // Transmitted to current_layer - 1
        if (g_rand.next() > refl) {
            // Escaped
            if (layer == 1) {
                photon.direction.z = -uz1;
                RecordReflectance(params, photon, tracer, 0.0);
                photon.alive = 0;
            }
            else {
                photon.current_layer--;
                photon.direction.x *= eta_i / eta_t;
                photon.direction.y *= eta_i / eta_t;
                photon.direction.z = -uz1;
            }
        }
        // Reflected
        else {
            photon.direction.z = -uz;
        }
    }
}

/*******************************************************************************
 *	Decide whether the photon will be transmitted or reflected on the bottom 
 *  boundary (uz > 0) of the current layer.
 *
 *	If the photon is transmitted, move the photon to current_layer + 1. 
 *  If current_layer is the last layer, record the weight as transmittance. 
 *
 *	Update the photon parmameters.
 ****/
void CrossDown(RunParams& params, Photon& photon, Tracer& tracer)
{
    // Z directional cosine
    double uz = photon.direction.z;

    // Cosines of transmission alpha
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double eta_i = params.layers[layer].eta;
    double eta_t = params.layers[layer + 1].eta;

    // Get reflectance; if 1.0, then total internal reflection
    double refl = (uz <= params.layers[layer].cos_theta_c1) ? 1.0 : FresnelReflectance(eta_i, eta_t, uz, uz1);

    if (PARTIAL_REFLECTION)
    {
        if (layer == params.num_layers && refl < 1.0) {
            photon.direction.z = uz1;
            RecordTransmittance(params, photon, tracer, refl);
            photon.direction.z = -uz;
        }
        // Transmitted to current_layer + 1
        else if (g_rand.next() > refl) {
            photon.current_layer++;
            photon.direction.x *= eta_i / eta_t;
            photon.direction.y *= eta_i / eta_t;
            photon.direction.z = uz1;
        }
        // Reflected
        else
        {
            photon.direction.z = -uz;
        }
    }
    else 
    {
        // Transmitted to current_layer + 1
        if (g_rand.next() > refl) {
            if (layer == params.num_layers) {
                photon.direction.z = uz1;
                RecordTransmittance(params, photon, tracer, 0.0);

                // Escaped
                photon.alive = 0;
            }
            else {
                photon.current_layer++;
                photon.direction.x *= eta_i / eta_t;
                photon.direction.y *= eta_i / eta_t;
                photon.direction.z = uz1;
            }
        }
        // Reflected
        else {
            photon.direction.z = -uz;
        }
    }
}

/*******************************************************************************
 *	Set a step size if the previous step has finished.
 *
 *	If the step size fits in the current layer, move the photon, drop weight, 
 *  and choose a new photon direction for propagation.
 *
 *	If the step size is long enough for the photon to hit an interface, this step 
 *  is divided into three steps:
 * 
 *	First, move the photon to the boundary free of absorption or scattering.
 *	Second, update the step size to the unfinished step size.
 *	Third, decide whether the photon is reflected or transmitted.
 ****/
void HopDropSpin(RunParams& params, Photon& photon, Tracer& tracer)
{
    Layer& layer = params.layers[photon.current_layer];
    double mu_t = layer.mu_a + layer.mu_s;
    
    SetStepSize(photon);

    // Distance between photon and boundary cm
    double path = PathToBoundary(photon, params);

    // Hit boundary
    if (path * mu_t <= photon.step_size) {
        // Move to boundary plane
        Hop(photon, path, layer.eta);

        // Update s
        photon.step_size -= path * mu_t;

        if (photon.direction.z < 0.0) {
            CrossUp(params, photon, tracer);
        }
        else {
            CrossDown(params, photon, tracer);
        }
    }
    // Fit in current_layer
    else {
        Hop(photon, photon.step_size / mu_t, layer.eta);

        // Update s
        photon.step_size = 0;
        Drop(params, photon, tracer);
        Spin(photon, params.layers[photon.current_layer].g);
        photon.num_scatters++;
    }
}

/*******************************************************************************
 *	Trace a photon, then compute the absorption, transmittance, and reflection 
 *  constants, including their standard errors.
 ****/
void TracePhoton(RunParams& params, Photon& photon, Tracer& tracer)
{
    do {
        HopDropSpin(params, photon, tracer);
        if (photon.alive && photon.weight < params.weight_treshold) {
            Roulette(photon);
        }
    } while (photon.alive);

    tracer.A_total += tracer.A_i;
    tracer.A_error += tracer.A_i * tracer.A_i;

    tracer.Tb_total += tracer.Tb_i;
    tracer.Tb_error += tracer.Tb_i * tracer.Tb_i;
    tracer.T_total += tracer.T_i;
    tracer.T_error += tracer.T_i * tracer.T_i;

    tracer.Rb_total += tracer.Rb_i;
    tracer.Rb_error += tracer.Rb_i * tracer.Rb_i;
    tracer.R_total += tracer.R_i;
    tracer.R_error += tracer.R_i * tracer.R_i;
}
