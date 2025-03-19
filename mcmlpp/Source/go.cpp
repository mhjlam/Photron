/**************************************************************************
 *	Copyright Univ. of Texas M.D. Anderson Cancer Center, 1992-1996.
 *	Program name: MCML.
 *	Launch, move, and record photon min_weight.
 ****/


#include "mcml.hpp"

#include <random>


 // 1=split photon, 0=statistical reflection.
constexpr int PARTIAL_REFLECTION = 0;

// If 1-cos(theta) <= ONE_MINUS_COS_ZERO, fabs(theta) <= 1e-6 rad.
// If 1+cos(theta) <= ONE_MINUS_COS_ZERO, fabs(PI-theta) <= 1e-6 rad.
constexpr double ONE_MINUS_COS_ZERO = 1.0E-12;

// If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad.
constexpr double COS_90_DEGREES = 1.0E-6;

/**************************************************************************
 *	A random number generator that generates uniformly
 *	distributed random numbers between 0 and 1 inclusive.
 *	The algorithm is based on:
 *	W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *	Flannery, "Numerical Recipes in C," Cambridge University
 *	Press, 2nd edition, (1992).
 *	and
 *	D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *	of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *	When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *	When Type is 1, returns a random number.
 *	When Type is 2, gets the status of the generator.
 *	When Type is 3, restores the status of the generator.
 *
 *	The status of the generator is represented by Status[0..56].
 *
 *	Make sure you initialize the seed before you get random
 *	numbers.
 ****/
double RandomGen()
{
    std::random_device random_device;                               // Non-deterministic seed
    std::mt19937 generator(random_device());                        // Mersenne Twister generator
    std::uniform_real_distribution<double> distribution(0.0, 1.0);  // Uniform distribution [0,1]
    return distribution(generator);
}

/**************************************************************************
 *	Compute the specular reflectance.
 *
 *	If the first layer is a turbid medium, use the Fresnel reflection from
 *  the boundary between the top abmient medium and the first layer as the
 *  specular reflectance.
 *
 *	The subroutine assumes the Layerspecs array is correctly initialized.
 ****/
double Rspecular(std::vector<Layer>& layers)
{
    double r = (layers[0].n - layers[1].n) / (layers[0].n + layers[1].n);
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
    double temp = (1 - g * g) / (1 - g + 2 * g * RandomGen());
    double cost = (g == 0.0) ? 2 * RandomGen() - 1 : (1 + g * g - temp * temp) / (2 * g);

    // sqrt is faster than sin.
    double sint = sqrt(1.0 - cost * cost);

    // sample psi.
    double psi = 2.0 * PI * RandomGen();
    double cosp = cos(psi);

    // sqrt is faster than sin.
    double sinp = (psi < PI) ? sqrt(1.0 - cosp * cosp) : -sqrt(1.0 - cosp * cosp);

    // update directional cosines.

    // close to perpendicular.
    if (1 - fabs(uz) <= ONE_MINUS_COS_ZERO) {
        photon.ux = sint * cosp;
        photon.uy = sint * sinp;

        // SIGN() is faster than division.
        photon.uz = cost * sign(uz);
    }
    else {
        double temp = sqrt(1.0 - uz * uz);
        photon.ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
        photon.uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
        photon.uz = -sint * cosp * temp + uz * cost;
    }
}

/**************************************************************************
 *	Initialize a photon packet.
 *
 *	If an Isotropic source is launched inside a glass layer, we check
 *	whether the photon will be total-internally reflected.  If it
 *	does, the photon is killed to avoid a infinite travelling inside
 *	the glass layer.
 ****/
void LaunchPhoton(double Rsp, RunParams& run_params, Tracer& tracer, Photon& photon)
{
    photon.w = 1.0 - Rsp;
    photon.alive = 1;
    photon.layer = (run_params.slayer != 0) ? run_params.slayer : 1;
    photon.s = 0;
    photon.sleft = 0;
    photon.scatters = 0;
    photon.time = 0;

    photon.x = 0.0;
    photon.y = 0.0;
    photon.z = run_params.source_z;
    photon.ux = 0.0;
    photon.uy = 0.0;
    photon.uz = 1.0;

    tracer.A.i = 0.0;
    tracer.T.bi = 0.0;
    tracer.T.di = 0.0;
    tracer.R.bi = 0.0;
    tracer.R.di = 0.0;

    if (run_params.source == BeamType::Isotropic) {
        Layer layer = run_params.layers[photon.layer];

        // to avoid scoring into Rb or Tb.
        photon.scatters++;

        // isotropically scatter the photon.
        Spin(0.0, photon);

        // glass layer.
        if (layer.mua == 0.0 && layer.mus == 0.0) {
            // total internal reflection
            if (fabs(photon.uz) <= layer.cos_crit0 && fabs(photon.uz) <= layer.cos_crit1) {
                photon.alive = 0;
            }
        }
    }
}

/**************************************************************************
 *	Move the photon S away in the current layer of medium.
 ****/
void Hop(Photon& photon, double S, double n)
{
    photon.x += S * photon.ux;
    photon.y += S * photon.uy;
    photon.z += S * photon.uz;
    photon.time += S * n * ONE_OVER_C;
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
    if (photon.s == 0.0) {
        // avoid zero.
        double rnd;
        while ((rnd = RandomGen()) <= 0.0);
        photon.s = -log(rnd);
    }
}

/**************************************************************************
 *	Return the distance between the photon position to the boundary along
 *  the photon direction.
 ****/
double PathToBoundary(Photon& photon, RunParams& run_params)
{
    // length to boundary.
    short layer = photon.layer;
    double uz = photon.uz;

    // Distance to the boundary.
    double path = DBL_MAX; // infinity

    // path > 0.
    if (uz > 0.0) {
        path = (run_params.layers[layer].z1 - photon.z) / uz;
    }
    // path > 0.
    else if (uz < 0.0) {
        path = (run_params.layers[layer].z0 - photon.z) / uz;
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
void Drop(RunParams& run_params, Photon& photon, Tracer& tracer)
{
    // Absorbed min_weight.
    double x = photon.x;
    double y = photon.y;

    short layer = photon.layer;
    Record record = run_params.record;

    // Update photon min_weight.
    double mua = run_params.layers[layer].mua;
    double mus = run_params.layers[layer].mus;
    double dwa = photon.w * mua / (mua + mus);
    photon.w -= dwa;

    // Compute array indices.
    short iz, ir, it;
    if (record.A_rzt || record.A_zt || record.A_z || record.A_rz) {
        if (photon.z >= run_params.zm) {
            iz = run_params.nz - 1;
        }
        else {
            iz = (short)(photon.z / run_params.dz);
        }
    }

    if (record.A_rzt || record.A_zt || record.A_t) {
        if (photon.time >= run_params.tm) {
            it = run_params.nt - 1;
        }
        else {
            it = (short)floor(photon.time / run_params.dt);
        }
    }

    // Assign dwa to an absorption array element.

    // Scattered.
    if (photon.scatters) {
        if (record.A_rzt || record.A_rz) {
            double temp = sqrt(x * x + y * y);
            if (temp >= run_params.rm) {
                ir = run_params.nr - 1;
            }
            else {
                ir = (short)(temp / run_params.dr);
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
    tracer.A.i += dwa;
}

/**************************************************************************
 *	The photon min_weight is small, and the photon packet tries
 *	to survive a roulette.
 ****/
void Roulette(Photon& photon)
{
    // already dead.
    if (photon.w == 0.0) {
        photon.alive = 0;
    }

    // survived the roulette.
    else if (RandomGen() < CHANCE) {
        photon.w /= CHANCE;
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
 *	cat_Ptr: pointer to the cosine of the transmission
 *			 angle at.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 *
 * ni:  incident refractive index
 * nt:  transmit refractive index
 * cai: cosine of angle ai. +
 ****/
double RFresnel(double ni, double nt, double cai, double* cat_Ptr)
{
    // Matched boundary.
    if (ni == nt) {
        *cat_Ptr = cai;
        return 0.0;
    }
    // Normal incidence.
    else if (1 - cai <= ONE_MINUS_COS_ZERO) {
        *cat_Ptr = cai;
        double r = (nt - ni) / (nt + ni);
        return r * r;
    }
    // Very slant incidence.
    else if (cai < COS_90_DEGREES) {
        *cat_Ptr = 0.0;
        return 1.0;
    }
    // General incidence.
    else {
        // Sine of the angles ai & at.
        double sai = sqrt(1 - cai * cai);
        double sat = ni * sai / nt;

        // Total internal reflection.
        if (sat >= 1.0) {
            *cat_Ptr = 0.0;
            return 1.0;
        }

        // cosine of at.
        double cat = sqrt(1 - sat * sat);
        *cat_Ptr = cat;

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
 *	Record the photon min_weight exiting the first layer(uz<0),
 *	to the reflection array.
 *
 *	Update the photon min_weight as well.
 ****/
void RecordR(double Refl, RunParams& run_params, Photon& photon, Tracer& tracer) // reflectance.
{
    double x = photon.x;
    double y = photon.y;
    Record record = run_params.record;

    short ir, ia, it;	// index to r & angle.
    if (record.Rd_rat || record.Rd_at || record.Rd_rt || record.Rd_t) {
        if (photon.time >= run_params.tm) {
            it = run_params.nt - 1;
        }
        else {
            it = (short)floor(photon.time / run_params.dt);
        }
    }

    // Scattered.
    if (photon.scatters) {
        if (record.Rd_rat || record.Rd_rt || record.Rd_ra || record.Rd_r) {
            double temp = sqrt(x * x + y * y);
            if (temp >= run_params.rm) {
                ir = run_params.nr - 1;
            }
            else {
                ir = (short)(temp / run_params.dr);
            }
        }
        if (record.Rd_rat || record.Rd_at || record.Rd_ra || record.Rd_a) {
            if ((ia = (short)(acos(-photon.uz) / run_params.da) > run_params.na - 1)) {
                ia = run_params.na - 1;
            }
        }

        // Assign photon min_weight to the reflection array element.
        if (record.Rd_rat) {
            tracer.R.rat[ir][ia][it] += photon.w * (1.0 - Refl);
        }
        if (record.Rd_ra) {
            tracer.R.ra[ir][ia] += photon.w * (1.0 - Refl);
        }

        if (record.Rd_rt) {
            tracer.R.rt[ir][it] += photon.w * (1.0 - Refl);
        }
        if (record.Rd_r) {
            tracer.R.r[ir] += photon.w * (1.0 - Refl);
        }

        if (record.Rd_at) {
            tracer.R.at[ia][it] += photon.w * (1.0 - Refl);
        }
        if (record.Rd_a) {
            tracer.R.a[ia] += photon.w * (1.0 - Refl);
        }
        if (record.Rd_t) {
            tracer.R.t[it] += photon.w * (1.0 - Refl);
        }
        tracer.R.di += photon.w * (1.0 - Refl);

    }

    // Ballistic.
    else {
        tracer.R.bi += photon.w * (1.0 - Refl);
    }

    photon.w *= Refl;
}

/**************************************************************************
 *	Record the photon min_weight exiting the last layer (uz > 0), no matter
 *  whether the layer is glass or not, to the transmittance array.
 *
 *	Update the photon min_weight as well.
 ****/
void RecordT(double Refl, RunParams& run_params, Photon& photon, Tracer& tracer)
{
    double x = photon.x;
    double y = photon.y;
    Record record = run_params.record;

    // Index to r & angle.
    short ir, ia, it;
    if (record.Td_rat || record.Td_at || record.Td_rt || record.Td_t) {
        if (photon.time >= run_params.tm) {
            it = run_params.nt - 1;
        }
        else {
            it = (short)floor(photon.time / run_params.dt);
        }
    }

    // Scattered.
    if (photon.scatters) {
        if (record.Td_rat || record.Td_rt || record.Td_ra || record.Td_r) {
            double temp = sqrt(x * x + y * y);
            if (temp >= run_params.rm) {
                ir = run_params.nr - 1;
            }
            else {
                ir = (short)(temp / run_params.dr);
            }
        }
        if (record.Td_rat || record.Td_at || record.Td_ra || record.Td_a) {
            if ((ia = (short)(acos(photon.uz) / run_params.da) > run_params.na - 1)) {
                ia = run_params.na - 1;
            }
        }

        // Assign photon min_weight to the transmittance array element.
        if (record.Td_rat) {
            tracer.T.rat[ir][ia][it] += photon.w * (1.0 - Refl);
        }
        if (record.Td_ra) {
            tracer.T.ra[ir][ia] += photon.w * (1.0 - Refl);
        }

        if (record.Td_rt) {
            tracer.T.rt[ir][it] += photon.w * (1.0 - Refl);
        }
        if (record.Td_r) {
            tracer.T.r[ir] += photon.w * (1.0 - Refl);
        }

        if (record.Td_at) {
            tracer.T.at[ia][it] += photon.w * (1.0 - Refl);
        }
        if (record.Td_a) {
            tracer.T.a[ia] += photon.w * (1.0 - Refl);
        }

        if (record.Td_t) {
            tracer.T.t[it] += photon.w * (1.0 - Refl);
        }
        tracer.T.di += photon.w * (1.0 - Refl);
    }

    // Collimated.
    else {
        tracer.T.bi += photon.w * (1.0 - Refl);
    }

    photon.w *= Refl;
}

/**************************************************************************
 *	Decide whether the photon will be transmitted or
 *	reflected on the upper boundary (uz<0) of the current
 *	layer.
 *
 *	If "layer" is the first layer, the photon packet will
 *	be partially transmitted and partially reflected if
 *	PARTIAL_REFLECTION is set to 1,
 *	or the photon packet will be either transmitted or
 *	reflected determined statistically if PARTIAL_REFLECTION
 *	is set to 0.
 *
 *	Record the transmitted photon min_weight as reflection.
 *
 *	If the "layer" is not the first layer and the photon
 *	packet is transmitted, move the photon to "layer-1".
 *
 *	Update the photon parmameters.
 ****/
void CrossUpOrNot(RunParams& run_params, Photon& photon, Tracer& tracer)
{
    // z directional cosine.
    double uz = photon.uz;

    // cosines of transmission alpha. always * +.
    double uz1 = 0.0;

    short layer = photon.layer;
    double ni = run_params.layers[layer].n;
    double nt = run_params.layers[layer - 1].n;

    // Get reflectance. If 1.0, then total internal reflection.
    double r = (-uz <= run_params.layers[layer].cos_crit0) ? 1.0 : RFresnel(ni, nt, -uz, &uz1);

#if PARTIAL_REFLECTION
    if (layer == 1 && r < 1.0) {	// partially transmitted.
        photon.uz = -uz1;	// escaped photon.
        RecordR(r, run_params, Photon_Ptr, Out_Ptr);
        photon.uz = -uz;	// reflected photon.
    }
    else if (RandomNum > r) {	// transmitted to layer-1.
        photon.layer--;
        photon.ux *= ni / nt;
        photon.uy *= ni / nt;
        photon.uz = -uz1;
    }
    else			// reflected.
        photon.uz = -uz;
#else
    // transmitted to layer-1.
    if (RandomGen() > r) {
        if (layer == 1) {
            photon.uz = -uz1;
            RecordR(0.0, run_params, photon, tracer);
            photon.alive = 0;	// escaped.
        }
        else {
            photon.layer--;
            photon.ux *= ni / nt;
            photon.uy *= ni / nt;
            photon.uz = -uz1;
        }
    }
    // reflected.
    else {
        photon.uz = -uz;
    }
#endif
}

/**************************************************************************
 *	Decide whether the photon will be transmitted  or be
 *	reflected on the bottom boundary (uz>0) of the current
 *	layer.
 *
 *	If the photon is transmitted, move the photon to
 *	"layer+1". If "layer" is the last layer, record the
 *	transmitted min_weight as transmittance. See comments for
 *	CrossUpOrNot.
 *
 *	Update the photon parmameters.
 ****/
void CrossDnOrNot(RunParams& run_params, Photon& photon, Tracer& tracer)
{
    // z directional cosine.
    double uz = photon.uz;

    // cosines of transmission alpha.
    double uz1 = 0.0;

    short layer = photon.layer;
    double ni = run_params.layers[layer].n;
    double nt = run_params.layers[layer + 1].n;

    // Get reflectance. If 1.0, then total internal reflection.
    double r = (uz <= run_params.layers[layer].cos_crit1) ? 1.0 : RFresnel(ni, nt, uz, &uz1);

#if PARTIAL_REFLECTION
    if (layer == run_params.layers.size() && r < 1.0) {
        photon.uz = uz1;
        RecordT(r, run_params, Photon_Ptr, Out_Ptr);
        photon.uz = -uz;
    }
    else if (RandomNum > r) {	// transmitted to layer+1.
        photon.layer++;
        photon.ux *= ni / nt;
        photon.uy *= ni / nt;
        photon.uz = uz1;
    }
    else			// reflected.
        photon.uz = -uz;
#else
    // transmitted to layer+1.
    if (RandomGen() > r) {
        if (layer == run_params.layers.size()) {
            photon.uz = uz1;
            RecordT(0.0, run_params, photon, tracer);

            // escaped.
            photon.alive = 0;
        }
        else {
            photon.layer++;
            photon.ux *= ni / nt;
            photon.uy *= ni / nt;
            photon.uz = uz1;
        }
    }
    // reflected.
    else {
        photon.uz = -uz;
    }
#endif
}

/**************************************************************************
 *	Set a step size if the previous step has finished.
 *
 *	If the step size fits in the current layer, move the photon,
 *	drop some min_weight, choose a new photon direction for propagation.
 *
 *	If the step size is long enough for the photon to
 *	hit an interface, this step is divided into three steps.
 *	First, move the photon to the boundary free of
 *	absorption or scattering.
 *	Second, update the step size to the unfinished step size.
 *	Third, decide whether the photon is reflected or transmitted.
 ****/
void HopDropSpin(RunParams& run_params, Photon& photon, Tracer& tracer)
{
    Layer layer_struct = run_params.layers[photon.layer];
    double mut = layer_struct.mua + layer_struct.mus;
    SetStepSize(photon);

    // distance between photon and boundary cm.
    double path = PathToBoundary(photon, run_params);

    // hit boundary.
    if (path * mut <= photon.s) {
        // move to boundary plane.
        Hop(photon, path, layer_struct.n);

        // update s.
        photon.s -= path * mut;

        if (photon.uz < 0.0) {
            CrossUpOrNot(run_params, photon, tracer);
        }
        else {
            CrossDnOrNot(run_params, photon, tracer);
        }
    }
    // fit in layer.
    else {
        Hop(photon, photon.s / mut, layer_struct.n);

        // update s.
        photon.s = 0;
        Drop(run_params, photon, tracer);
        Spin(run_params.layers[photon.layer].g, photon);
        photon.scatters++;
    }
}

/**************************************************************************
 *	Trace a photon, then compute the 0D constants including A, R, T and
 *	their standard errors.
 ****/
void TracePhoton(RunParams& run_params, Photon& photon, Tracer& tracer)
{
    do {
        HopDropSpin(run_params, photon, tracer);
        if (photon.alive && photon.w < run_params.min_weight) {
            Roulette(photon);
        }
    } while (photon.alive);

    tracer.A.a += tracer.A.i;
    tracer.A.e += tracer.A.i * tracer.A.i;

    tracer.T.b += tracer.T.bi;
    tracer.T.be += tracer.T.bi * tracer.T.bi;
    tracer.T.d += tracer.T.di;
    tracer.T.de += tracer.T.di * tracer.T.di;

    tracer.R.b += tracer.R.bi;
    tracer.R.be += tracer.R.bi * tracer.R.bi;
    tracer.R.d += tracer.R.di;
    tracer.R.de += tracer.R.di * tracer.R.di;
}
