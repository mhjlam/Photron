/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Launch, move, and record photon min_weight.
 ****/


#include "mcml.hpp"

#include <numbers>


/**************************************************************************
 *	Compute the specular reflectance.
 *
 *	If the first current_layer is a turbid name, use the Fresnel reflection from
 *  the boundary between the top abmient name and the first current_layer as the
 *  specular reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly initialized.
 ****/
double Rspecular(std::vector<Layer>& layers)
{
    double r = (layers[0].eta - layers[1].eta) / (layers[0].eta + layers[1].eta);
    return (r * r);
}

/**************************************************************************
 *  Choose a new direction for photon propagation by sampling
 *	1. the polar deflection angle theta
 *	2. the azimuthal angle psi.
 *
 *  Note:
 *  theta: 0 - pi so sin(theta) is always positive
 *  Feel free to use sqrt() for cos(theta).
 *
 *  psi: 0 - 2pi
 *  for 0-pi:   sin(psi) is +
 *  for pi-2pi: sin(psi) is -
 ****/
void Spin(double g, Photon& photon)
{
    // cosine and sine of psi.
    double ux = photon.ux;
    double uy = photon.uy;
    double uz = photon.uz;

    // sample theta.
    double temp = (1 - g * g) / (1 - g + 2 * g * Rand.next());
    double cost = (g == 0.0) ? 2 * Rand.next() - 1 : (1 + g * g - temp * temp) / (2 * g);

    // sqrt is faster than sin.
    double sint = std::sqrt(1.0 - cost * cost);

    // sample psi.
    double psi = 2.0 * std::numbers::pi * Rand.next();
    double cosp = std::cos(psi);

    // sqrt is faster than sin.
    double sinp = (psi < std::numbers::pi) ? std::sqrt(1.0 - cosp * cosp) : -std::sqrt(1.0 - cosp * cosp);

    // update directional cosines.

    // close to perpendicular.
    if (1 - std::abs(uz) <= ONE_MINUS_COS_ZERO) {
        photon.ux = sint * cosp;
        photon.uy = sint * sinp;
        photon.uz = cost * std::copysign(1.0, uz);
    }
    else {
        double temp = std::sqrt(1.0 - uz * uz);
        photon.ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
        photon.uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
        photon.uz = -sint * cosp * temp + uz * cost;
    }
}

/**************************************************************************
 *	Initialize a photon packet.
 *
 *	If an Isotropic source is launched inside a glass current_layer, we check
 *	whether the photon will be total-internally reflected.  If it
 *	does, the photon is killed to avoid a infinite travelling inside
 *	the glass layer.
 ****/
void LaunchPhoton(double Rsp, RunParams& params, Tracer& tracer, Photon& photon)
{
    photon.min_weight = 1.0 - Rsp;
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

        // to avoid scoring into Rb or Tb.
        photon.num_scatters++;

        // isotropically scatter the photon.
        Spin(0.0, photon);

        // glass current_layer.
        if (layer.mua == 0.0 && layer.mus == 0.0) {
            // total internal reflection
            if (std::abs(photon.uz) <= layer.cos_crit0 && std::abs(photon.uz) <= layer.cos_crit1) {
                photon.alive = 0;
            }
        }
    }
}

/**************************************************************************
 *	Move the photon S away in the current current_layer of name.
 ****/
void Hop(Photon& photon, double S, double n)
{
    photon.x += S * photon.ux;
    photon.y += S * photon.uy;
    photon.z += S * photon.uz;
    photon.flight_time += S * n * SPEED_OF_LIGHT_INV;
}

/**************************************************************************
 *	Pick a step size in dimensionless unit for a photon packet.
 *	If the member s is zero, make a new step size
 *	with: -log(rnd).
 *	Otherwise, finish the leftover in s.
 ****/
void SetStepSize(Photon& photon)
{
    // make a new step.
    if (photon.step_size == 0.0) {
        // avoid zero.
        double rnd;
        while ((rnd = Rand.next()) <= 0.0);
        photon.step_size = -std::log(rnd);
    }
}

/**************************************************************************
 *	Return the distance between the photon position to the boundary along
 *  the photon direction.
 ****/
double PathToBoundary(Photon& photon, RunParams& params)
{
    // length to boundary.
    short layer = photon.current_layer;
    double uz = photon.uz;

    // Distance to the boundary.
    double path = DBL_MAX; // infinity

    // path > 0.
    if (uz > 0.0) {
        path = (params.layers[layer].bot_z - photon.z) / uz;
    }
    // path > 0.
    else if (uz < 0.0) {
        path = (params.layers[layer].top_z - photon.z) / uz;
    }

    return (path);
}

/**************************************************************************
 *	Drop photon min_weight inside the tissue (not glass).
 *
 *      The photon is assumed alive.
 *
 *	The min_weight drop is dw = w*mua/(mua+mus).
 *
 *	The dropped min_weight is assigned to the absorption array
 *	elements.
 ****/
void Drop(RunParams& params, Photon& photon, Tracer& tracer)
{
    // Absorbed min_weight.
    double x = photon.x;
    double y = photon.y;

    short layer = photon.current_layer;
    Record record = params.record;

    // Update photon min_weight.
    double mua = params.layers[layer].mua;
    double mus = params.layers[layer].mus;
    double dwa = photon.min_weight * mua / (mua + mus);
    photon.min_weight -= dwa;

    // Compute array indices.
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
            it = params.num_time - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_time));
        }
    }

    // Assign dwa to an absorption array element.

    // Scattered.
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

    // Ballistic.
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

/**************************************************************************
 *	The photon min_weight is small, and the photon packet tries
 *	to survive a roulette.
 ****/
void Roulette(Photon& photon)
{
    // already dead.
    if (photon.min_weight == 0.0) {
        photon.alive = 0;
    }

    // survived the roulette.
    else if (Rand.next() < ROULETTE_SURVIVAL) {
        photon.min_weight /= ROULETTE_SURVIVAL;
    }
    else {
        photon.alive = 0;
    }
}

/**************************************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle ai
 *	is positive, and the case when the angle is greater
 *	than the critical angle is ruled out.
 *
 *	cai: cosine of the incident angle ai.
 *	cos_at: pointer to the cosine of the transmission angle at.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 *
 * ni:  incident refractive index
 * num_time:  transmit refractive index
 * cai: cosine of angle ai. +
 ****/
double RFresnel(double ni, double nt, double cai, double& cos_at)
{
    // Matched boundary.
    if (ni == nt) {
        cos_at = cai;
        return 0.0;
    }
    // Normal incidence.
    else if (1 - cai <= ONE_MINUS_COS_ZERO) {
        cos_at = cai;
        double r = (nt - ni) / (nt + ni);
        return r * r;
    }
    // Very slant incidence.
    else if (cai < COS_90_DEGREES) {
        cos_at = 0.0;
        return 1.0;
    }
    // General incidence.
    else {
        // Sine of the angles ai & at.
        double sai = std::sqrt(1 - cai * cai);
        double sat = ni * sai / nt;

        // Total internal reflection.
        if (sat >= 1.0) {
            cos_at = 0.0;
            return 1.0;
        }

        // cosine of at.
        double cat = std::sqrt(1 - sat * sat);
        cos_at = cat;

        // cosines of ai+at & ai-at.
        double cap = cai * cat - sai * sat;	// c+ = cc - ss.
        double cam = cai * cat + sai * sat;	// c- = cc + ss.

        // sines of ai+at & ai-at.
        double sap = sai * cat + cai * sat;	// s+ = sc + cs.
        double sam = sai * cat - cai * sat;	// s- = sc - cs.

        // arranged for speed.
        return 0.5 * sam * sam * (cam * cam + cap * cap) / (sap * sap * cam * cam);
    }
    return 0.0;
}

/**************************************************************************
 *	Record the photon min_weight exiting the first current_layer(uz<0),
 *	to the reflection array.
 *
 *	Update the photon min_weight as well.
 ****/
void RecordR(double Refl, RunParams& params, Photon& photon, Tracer& tracer) // reflectance.
{
    double x = photon.x;
    double y = photon.y;
    Record record = params.record;

    std::size_t ir, ia, it;	// index to r & angle.
    if (record.Rd_rat || record.Rd_at || record.Rd_rt || record.Rd_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_time - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_time));
        }
    }

    // Scattered.
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
            if ((ia = static_cast<short>(std::acos(-photon.uz) / params.grid_alpha) > params.num_alpha - 1)) {
                ia = params.num_alpha - 1;
            }
        }

        // Assign photon min_weight to the reflection array element.
        if (record.Rd_rat) {
            tracer.R.rat[ir][ia][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Rd_ra) {
            tracer.R.ra[ir][ia] += photon.min_weight * (1.0 - Refl);
        }

        if (record.Rd_rt) {
            tracer.R.rt[ir][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Rd_r) {
            tracer.R.r[ir] += photon.min_weight * (1.0 - Refl);
        }

        if (record.Rd_at) {
            tracer.R.at[ia][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Rd_a) {
            tracer.R.a[ia] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Rd_t) {
            tracer.R.t[it] += photon.min_weight * (1.0 - Refl);
        }
        tracer.R.di += photon.min_weight * (1.0 - Refl);

    }

    // Ballistic.
    else {
        tracer.R.bi += photon.min_weight * (1.0 - Refl);
    }

    photon.min_weight *= Refl;
}

/**************************************************************************
 *	Record the photon min_weight exiting the last current_layer (uz > 0), no matter
 *  whether the current_layer is glass or not, to the transmittance array.
 *
 *	Update the photon min_weight as well.
 ****/
void RecordT(double Refl, RunParams& params, Photon& photon, Tracer& tracer)
{
    double x = photon.x;
    double y = photon.y;
    Record record = params.record;

    // Index to r & angle.
    std::size_t ir, ia, it;
    if (record.Td_rat || record.Td_at || record.Td_rt || record.Td_t) {
        if (photon.flight_time >= params.max_time) {
            it = params.num_time - 1;
        }
        else {
            it = static_cast<short>(std::floor(photon.flight_time / params.grid_time));
        }
    }

    // Scattered.
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
            if ((ia = static_cast<short>(std::acos(photon.uz) / params.grid_alpha) > params.num_alpha - 1)) {
                ia = params.num_alpha - 1;
            }
        }

        // Assign photon min_weight to the transmittance array element.
        if (record.Td_rat) {
            tracer.T.rat[ir][ia][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Td_ra) {
            tracer.T.ra[ir][ia] += photon.min_weight * (1.0 - Refl);
        }

        if (record.Td_rt) {
            tracer.T.rt[ir][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Td_r) {
            tracer.T.r[ir] += photon.min_weight * (1.0 - Refl);
        }

        if (record.Td_at) {
            tracer.T.at[ia][it] += photon.min_weight * (1.0 - Refl);
        }
        if (record.Td_a) {
            tracer.T.a[ia] += photon.min_weight * (1.0 - Refl);
        }

        if (record.Td_t) {
            tracer.T.t[it] += photon.min_weight * (1.0 - Refl);
        }
        tracer.T.di += photon.min_weight * (1.0 - Refl);
    }

    // Collimated.
    else {
        tracer.T.bi += photon.min_weight * (1.0 - Refl);
    }

    photon.min_weight *= Refl;
}

/**************************************************************************
 *	Decide whether the photon will be transmitted or reflected on the upper
 *  boundary (uz<0) of the current current_layer.
 *
 *	If "current_layer" is the first current_layer, the photon packet will
 *	be partially transmitted and partially reflected if
 *	PARTIAL_REFLECTION is set to 1, or the photon packet will be either 
 *  transmitted or reflected determined statistically if PARTIAL_REFLECTION
 *	is set to 0.
 *
 *	Record the transmitted photon min_weight as reflection.
 *
 *	If the "current_layer" is not the first current_layer and the photon
 *	packet is transmitted, move the photon to "current_layer-1".
 *
 *	Update the photon parmameters.
 ****/
void CrossUpOrNot(RunParams& params, Photon& photon, Tracer& tracer)
{
    // z directional cosine.
    double uz = photon.uz;

    // cosines of transmission alpha. always * +.
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double ni = params.layers[layer].eta;
    double nt = params.layers[layer - 1].eta;

    // Get reflectance. If 1.0, then total internal reflection.
    double r = (-uz <= params.layers[layer].cos_crit0) ? 1.0 : RFresnel(ni, nt, -uz, uz1);

    if (PARTIAL_REFLECTION) {
        // partially transmitted.
        if (layer == 1 && r < 1.0) {
            photon.uz = -uz1; // escaped photon.
            RecordR(r, params, photon, tracer);
            photon.uz = -uz; // reflected photon.
        }
        // transmitted to current_layer-1.
        else if (Rand.next() > r) {
            photon.current_layer--;
            photon.ux *= ni / nt;
            photon.uy *= ni / nt;
            photon.uz = -uz1;
        }
        // reflected.
        else {
            photon.uz = -uz;
        }
    }
    else {
        // transmitted to current_layer-1.
        if (Rand.next() > r) {
            if (layer == 1) {
                photon.uz = -uz1;
                RecordR(0.0, params, photon, tracer);
                photon.alive = 0;	// escaped.
            }
            else {
                photon.current_layer--;
                photon.ux *= ni / nt;
                photon.uy *= ni / nt;
                photon.uz = -uz1;
            }
        }
        // reflected.
        else {
            photon.uz = -uz;
        }
    }
}

/**************************************************************************
 *	Decide whether the photon will be transmitted  or be
 *	reflected on the bottom boundary (uz>0) of the current
 *	current_layer.
 *
 *	If the photon is transmitted, move the photon to
 *	"current_layer+1". If "current_layer" is the last current_layer, record the
 *	transmitted min_weight as transmittance. See comments for
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void CrossDnOrNot(RunParams& params, Photon& photon, Tracer& tracer)
{
    // z directional cosine.
    double uz = photon.uz;

    // cosines of transmission alpha.
    double uz1 = 0.0;

    short layer = photon.current_layer;
    double ni = params.layers[layer].eta;
    double nt = params.layers[layer + 1].eta;

    // Get reflectance. If 1.0, then total internal reflection.
    double r = (uz <= params.layers[layer].cos_crit1) ? 1.0 : RFresnel(ni, nt, uz, uz1);

    if (PARTIAL_REFLECTION)
    {
        if (layer == params.num_layers && r < 1.0) {
            photon.uz = uz1;
            RecordT(r, params, photon, tracer);
            photon.uz = -uz;
        }
        // transmitted to current_layer+1.
        else if (Rand.next() > r) {
            photon.current_layer++;
            photon.ux *= ni / nt;
            photon.uy *= ni / nt;
            photon.uz = uz1;
        }
        // reflected.
        else
        {
            photon.uz = -uz;
        }
    }
    else 
    {
        // transmitted to current_layer+1.
        if (Rand.next() > r) {
            if (layer == params.num_layers) {
                photon.uz = uz1;
                RecordT(0.0, params, photon, tracer);

                // escaped.
                photon.alive = 0;
            }
            else {
                photon.current_layer++;
                photon.ux *= ni / nt;
                photon.uy *= ni / nt;
                photon.uz = uz1;
            }
        }
        // reflected.
        else {
            photon.uz = -uz;
        }
    }
}

/**************************************************************************
 *	Set a step size if the previous step has finished.
 *
 *	If the step size fits in the current current_layer, move the photon,
 *	drop some min_weight, choose a new photon direction for propagation.
 *
 *	If the step size is long enough for the photon to
 *	hit an interface, this step is divided into three steps.
 *	First, move the photon to the boundary free of
 *	absorption or scattering.
 *	Second, update the step size to the unfinished step size.
 *	Third, decide whether the photon is reflected or transmitted.
 ****/
void HopDropSpin(RunParams& params, Photon& photon, Tracer& tracer)
{
    Layer& layer = params.layers[photon.current_layer];
    double mut = layer.mua + layer.mus;
    
    SetStepSize(photon);

    // distance between photon and boundary cm.
    double path = PathToBoundary(photon, params);

    // hit boundary.
    if (path * mut <= photon.step_size) {
        // move to boundary plane.
        Hop(photon, path, layer.eta);

        // update s.
        photon.step_size -= path * mut;

        if (photon.uz < 0.0) {
            CrossUpOrNot(params, photon, tracer);
        }
        else {
            CrossDnOrNot(params, photon, tracer);
        }
    }
    // fit in current_layer.
    else {
        Hop(photon, photon.step_size / mut, layer.eta);

        // update s.
        photon.step_size = 0;
        Drop(params, photon, tracer);
        Spin(params.layers[photon.current_layer].aniso, photon);
        photon.num_scatters++;
    }
}

/**************************************************************************
 *	Trace a photon, then compute the 0D constants including A, R, T and
 *	their standard errors.
 ****/
void TracePhoton(RunParams& params, Photon& photon, Tracer& tracer)
{
    do {
        HopDropSpin(params, photon, tracer);
        if (photon.alive && photon.min_weight < params.min_weight) {
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
