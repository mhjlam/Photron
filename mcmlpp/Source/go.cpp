/*******************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *  Copyright M.H.J. Lam, 2025.
 *	Launch, move, and record photon weight.
 ****/


#include "mcml.hpp"

#include <limits>
#include <numbers>


/*******************************************************************************
 *	Compute the specular reflectance.
 *
 *	If the first layer is a turbid name, use the Fresnel reflection from the 
 *  boundary between the top abmient name and the first layer as the specular 
 *  reflectance.
 ****/
double Rspecular(std::vector<Layer>& layers)
{
    double r = (layers[0].eta - layers[1].eta) / (layers[0].eta + layers[1].eta);
    return (r * r);
}

/*******************************************************************************
 *  Choose a new direction for photon propagation by sampling:
 *	1. The polar (deflection) angle θ, measuring from downwards z-axis
 *	2. The azimuthal angle Ψ, measuring rotation around the z-axis in the xy-plane
 * 
 *  θ range:    0 - π  (0 to 180 degrees)    [sin(θ) range:  0 to 1]
 *  Ψ range:    0 - 2π (0 to 360 degrees)    [cos(Ψ) range: -1 to 1]
 * 
 *  Since cos²(θ) + sin²(θ) = 1, and if sin(θ) is known, then: cos(θ) = sqrt(1-sin(θ)²)
 ****/
void Spin(double g, Photon& photon)
{
    // Cosine and sine of Ψ
    double ux = photon.ux;
    double uy = photon.uy;
    double uz = photon.uz;

    // Sample θ
    double theta = (1 - g * g) / (1 - g + 2 * g * Rand.next());
    double cos_theta = (g == 0.0) ? 2 * Rand.next() - 1 : (1 + g * g - theta * theta) / (2 * g);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    // Sample Ψ
    double psi = 2.0 * std::numbers::pi * Rand.next();
    double cos_psi = std::cos(psi);
    double sin_psi = (psi < std::numbers::pi) ? std::sqrt(1.0 - cos_psi * cos_psi) : -std::sqrt(1.0 - cos_psi * cos_psi);

    // Update directional cosines

    // Close to perpendicular
    if (1 - std::abs(uz) <= COS_0_TOLERANCE) {
        photon.ux = sin_theta * cos_psi;
        photon.uy = sin_theta * sin_psi;
        photon.uz = cos_theta * std::copysign(1.0, uz);
    }
    else {
        double temp = std::sqrt(1.0 - uz * uz);
        photon.ux = sin_theta * (ux * uz * cos_psi - uy * sin_psi) / temp + ux * cos_theta;
        photon.uy = sin_theta * (uy * uz * cos_psi + ux * sin_psi) / temp + uy * cos_theta;
        photon.uz = -sin_theta * cos_psi * temp + uz * cos_theta;
    }
}

/*******************************************************************************
 *	Initialize a photon packet.
 *
 *	If an isotropic source is launched inside a glass layer, check whether the 
 *  photon will be total-internally reflected. If so, kill the photon to avoid 
 *  infinite travelling inside the glass layer.
 ****/
void LaunchPhoton(double Rsp, RunParams& params, Tracer& tracer, Photon& photon)
{
    photon.weight = 1.0 - Rsp;
    photon.alive = 1;
    photon.current_layer = (params.source_layer != 0) ? params.source_layer : 1;
    photon.step_size = 0;
    photon.step_size_left = 0;
    photon.num_scatters = 0;
    photon.flight_time = 0;

    photon.x = 0.0;
    photon.y = 0.0;
    photon.z = params.source_z;
    photon.ux = 0.0;
    photon.uy = 0.0;
    photon.uz = 1.0;

    tracer.A.ai = 0.0;
    tracer.T.bi = 0.0;
    tracer.T.di = 0.0;
    tracer.R.bi = 0.0;
    tracer.R.di = 0.0;

    if (params.source == BeamType::Isotropic) {
        Layer layer = params.layers[photon.current_layer];

        // To avoid scoring into reflectance or transmittance
        photon.num_scatters++;

        // Isotropically scatter the photon
        Spin(0.0, photon);

        // Glass layer.
        if (layer.mua == 0.0 && layer.mus == 0.0) {
            // Total internal reflection
            if (std::abs(photon.uz) <= layer.cos_crit0 && std::abs(photon.uz) <= layer.cos_crit1) {
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
    photon.x += dist * photon.ux;
    photon.y += dist * photon.uy;
    photon.z += dist * photon.uz;
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
        while ((rnd = Rand.next()) <= 0.0);

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
    double uz = photon.uz;

    // Distance to the boundary
    double path = std::numeric_limits<double>::max(); // infinity

    // Path > 0
    if (uz > 0.0) {
        path = (params.layers[layer].bot_z - photon.z) / uz;
    }
    // Path > 0
    else if (uz < 0.0) {
        path = (params.layers[layer].top_z - photon.z) / uz;
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
    double x = photon.x;
    double y = photon.y;

    short layer = photon.current_layer;
    Record record = params.record;

    // Update photon weight
    double mua = params.layers[layer].mua;
    double mus = params.layers[layer].mus;
    double dwa = photon.weight * mua / (mua + mus);
    photon.weight -= dwa;

    // Compute array indices
    std::size_t iz, ir, it;

    if (record.A_rzt || record.A_zt || record.A_z || record.A_rz) {
        if (photon.z >= params.max_z) {
            iz = params.num_z - 1;
        }
        else {
            iz = static_cast<short>(photon.z / params.grid_z);
        }
    }

    if (record.A_rzt || record.A_zt || record.A_t) {
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
        if (record.A_rzt || record.A_rz) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = static_cast<short>(temp / params.grid_r);
            }
        }
        if (record.A_rzt) {
            tracer.A.rzt[ir][iz][it] += dwa;
        }
        if (record.A_rz) {
            tracer.A.rz[ir][iz] += dwa;
        }
    }

    // Ballistic
    else {
        if (record.A_rzt) {
            tracer.A.bzt[iz][it] += dwa;
        }
        if (record.A_rz) {
            tracer.A.bz[iz] += dwa;
        }
    }

    if (record.A_zt) {
        tracer.A.zt[iz][it] += dwa;
    }
    if (record.A_z) {
        tracer.A.z[iz] += dwa;
    }

    if (record.A_t) {
        tracer.A.t[it] += dwa;
    }
    tracer.A.ai += dwa;
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
    else if (Rand.next() < ROULETTE_SURVIVAL) {
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
double RFresnel(double eta_i, double eta_t, double cos_ai, double& cos_at)
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
void RecordR(double reflectance, RunParams& params, Photon& photon, Tracer& tracer)
{
    double x = photon.x;
    double y = photon.y;
    Record record = params.record;

    // Index to r & angle
    std::size_t ir, ia, it;

    if (record.Rd_rat || record.Rd_at || record.Rd_rt || record.Rd_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_t - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_t));
        }
    }

    // Scattered
    if (photon.num_scatters) {
        if (record.Rd_rat || record.Rd_rt || record.Rd_ra || record.Rd_r) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = static_cast<short>(temp / params.grid_r);
            }
        }
        if (record.Rd_rat || record.Rd_at || record.Rd_ra || record.Rd_a) {
            if ((ia = static_cast<short>(std::acos(-photon.uz) / params.grid_a) > params.num_a - 1)) {
                ia = params.num_a - 1;
            }
        }

        // Assign photon weight to the reflection array element
        if (record.Rd_rat) {
            tracer.R.rat[ir][ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Rd_ra) {
            tracer.R.ra[ir][ia] += photon.weight * (1.0 - reflectance);
        }

        if (record.Rd_rt) {
            tracer.R.rt[ir][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Rd_r) {
            tracer.R.r[ir] += photon.weight * (1.0 - reflectance);
        }

        if (record.Rd_at) {
            tracer.R.at[ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Rd_a) {
            tracer.R.a[ia] += photon.weight * (1.0 - reflectance);
        }
        if (record.Rd_t) {
            tracer.R.t[it] += photon.weight * (1.0 - reflectance);
        }
        tracer.R.di += photon.weight * (1.0 - reflectance);

    }

    // Ballistic reflection
    else {
        tracer.R.bi += photon.weight * (1.0 - reflectance);
    }

    photon.weight *= reflectance;
}

/*******************************************************************************
 *	Record the photon weight exiting the last layer (uz > 0), no matter whether 
 *  the layer is glass or not, to the transmittance array.
 ****/
void RecordT(double reflectance, RunParams& params, Photon& photon, Tracer& tracer)
{
    double x = photon.x;
    double y = photon.y;
    Record record = params.record;

    // Index to r & angle
    std::size_t ir, ia, it;
    if (record.Td_rat || record.Td_at || record.Td_rt || record.Td_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_t - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_t));
        }
    }

    // Scattered
    if (photon.num_scatters) {
        if (record.Td_rat || record.Td_rt || record.Td_ra || record.Td_r) {
            double temp = std::sqrt(x * x + y * y);
            if (temp >= params.max_r) {
                ir = params.num_r - 1;
            }
            else {
                ir = (short)(temp / params.grid_r);
            }
        }
        if (record.Td_rat || record.Td_at || record.Td_ra || record.Td_a) {
            if ((ia = static_cast<short>(std::acos(photon.uz) / params.grid_a) > params.num_a - 1)) {
                ia = params.num_a - 1;
            }
        }

        // Assign photon weight to the transmittance array element
        if (record.Td_rat) {
            tracer.T.rat[ir][ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Td_ra) {
            tracer.T.ra[ir][ia] += photon.weight * (1.0 - reflectance);
        }

        if (record.Td_rt) {
            tracer.T.rt[ir][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Td_r) {
            tracer.T.r[ir] += photon.weight * (1.0 - reflectance);
        }

        if (record.Td_at) {
            tracer.T.at[ia][it] += photon.weight * (1.0 - reflectance);
        }
        if (record.Td_a) {
            tracer.T.a[ia] += photon.weight * (1.0 - reflectance);
        }

        if (record.Td_t) {
            tracer.T.t[it] += photon.weight * (1.0 - reflectance);
        }
        tracer.T.di += photon.weight * (1.0 - reflectance);
    }

    // Collimated
    else {
        tracer.T.bi += photon.weight * (1.0 - reflectance);
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
    double uz = photon.uz;

    // Cosines of transmission alpha.
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double ni = params.layers[layer].eta;
    double nt = params.layers[layer - 1].eta;

    // Get reflectance; if 1.0, then total internal reflection
    double r = (-uz <= params.layers[layer].cos_crit0) ? 1.0 : RFresnel(ni, nt, -uz, uz1);

    if (PARTIAL_REFLECTION) {
        // Partially transmitted
        if (layer == 1 && r < 1.0) {
            // Escaped photon
            photon.uz = -uz1;

            RecordR(r, params, photon, tracer);

            // Reflected photon
            photon.uz = -uz;
        }
        // Transmitted to current_layer - 1
        else if (Rand.next() > r) {
            photon.current_layer--;
            photon.ux *= ni / nt;
            photon.uy *= ni / nt;
            photon.uz = -uz1;
        }
        // Reflected
        else {
            photon.uz = -uz;
        }
    }
    else {
        // Transmitted to current_layer - 1
        if (Rand.next() > r) {
            // Escaped
            if (layer == 1) {
                photon.uz = -uz1;
                RecordR(0.0, params, photon, tracer);
                photon.alive = 0;
            }
            else {
                photon.current_layer--;
                photon.ux *= ni / nt;
                photon.uy *= ni / nt;
                photon.uz = -uz1;
            }
        }
        // Reflected
        else {
            photon.uz = -uz;
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
    double uz = photon.uz;

    // Cosines of transmission alpha
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double eta_i = params.layers[layer].eta;
    double eta_t = params.layers[layer + 1].eta;

    // Get reflectance; if 1.0, then total internal reflection
    double r = (uz <= params.layers[layer].cos_crit1) ? 1.0 : RFresnel(eta_i, eta_t, uz, uz1);

    if (PARTIAL_REFLECTION)
    {
        if (layer == params.num_layers && r < 1.0) {
            photon.uz = uz1;
            RecordT(r, params, photon, tracer);
            photon.uz = -uz;
        }
        // Transmitted to current_layer + 1
        else if (Rand.next() > r) {
            photon.current_layer++;
            photon.ux *= eta_i / eta_t;
            photon.uy *= eta_i / eta_t;
            photon.uz = uz1;
        }
        // Reflected
        else
        {
            photon.uz = -uz;
        }
    }
    else 
    {
        // Transmitted to current_layer + 1
        if (Rand.next() > r) {
            if (layer == params.num_layers) {
                photon.uz = uz1;
                RecordT(0.0, params, photon, tracer);

                // Escaped
                photon.alive = 0;
            }
            else {
                photon.current_layer++;
                photon.ux *= eta_i / eta_t;
                photon.uy *= eta_i / eta_t;
                photon.uz = uz1;
            }
        }
        // Reflected
        else {
            photon.uz = -uz;
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
    double mut = layer.mua + layer.mus;
    
    SetStepSize(photon);

    // Distance between photon and boundary cm
    double path = PathToBoundary(photon, params);

    // Hit boundary
    if (path * mut <= photon.step_size) {
        // Move to boundary plane
        Hop(photon, path, layer.eta);

        // Update s
        photon.step_size -= path * mut;

        if (photon.uz < 0.0) {
            CrossUp(params, photon, tracer);
        }
        else {
            CrossDown(params, photon, tracer);
        }
    }
    // Fit in current_layer
    else {
        Hop(photon, photon.step_size / mut, layer.eta);

        // Update s
        photon.step_size = 0;
        Drop(params, photon, tracer);
        Spin(params.layers[photon.current_layer].ani, photon);
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

    tracer.A.ab += tracer.A.ai;
    tracer.A.ae += tracer.A.ai * tracer.A.ai;

    tracer.T.br += tracer.T.bi;
    tracer.T.be += tracer.T.bi * tracer.T.bi;
    tracer.T.dr += tracer.T.di;
    tracer.T.de += tracer.T.di * tracer.T.di;

    tracer.R.br += tracer.R.bi;
    tracer.R.be += tracer.R.bi * tracer.R.bi;
    tracer.R.dr += tracer.R.di;
    tracer.R.de += tracer.R.di * tracer.R.di;
}
