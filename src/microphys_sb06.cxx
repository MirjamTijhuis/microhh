/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>
#include <vector>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "stats.h"
#include "cross.h"
#include "column.h"
#include "thermo.h"
#include "thermo_moist_functions.h"
#include "constants.h"

#include "microphys.h"

#include "microphys_sb06.h"         // Class definition
#include "microphys_sb_common.h"    // Warm/cold common kernels
#include "microphys_sb_cold.h"      // Kernels from ICON

namespace
{
    using namespace Constants;
    namespace fm = Fast_math;
    namespace tmf = Thermo_moist_functions;
}

namespace warm
{
    /* These are the "old" kernels used in the old `2mom_warm` scheme,
       ported (mainly) from UCLA-LES and DALES. */

    template<typename TF> constexpr TF pi       = 3.14159265359;
    template<typename TF> constexpr TF K_t      = 2.5e-2;              // Conductivity of heat [J/(sKm)]
    template<typename TF> constexpr TF D_v      = 3.e-5;               // Diffusivity of water vapor [m2/s]
    template<typename TF> constexpr TF rho_w    = 1.e3;                // Density water
    template<typename TF> constexpr TF rho_0    = 1.225;               // SB06, p48
    template<typename TF> constexpr TF pirhow   = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF mc_min   = 4.2e-15;             // Min mean mass of cloud droplet
    template<typename TF> constexpr TF mc_max   = 2.6e-10;             // Max mean mass of cloud droplet
    template<typename TF> constexpr TF mr_min   = mc_max<TF>;          // Min mean mass of precipitation drop
    template<typename TF> constexpr TF mr_max   = 3e-6;                // Max mean mass of precipitation drop // as in UCLA-LES
    template<typename TF> constexpr TF ql_min   = 1.e-6;               // Min cloud liquid water for which calculations are performed
    template<typename TF> constexpr TF qr_min   = 1.e-15;              // Min rain liquid water for which calculations are performed
    template<typename TF> constexpr TF cfl_min  = 1.e-5;               // Small non-zero limit at the CFL number


    template<typename TF>
    inline TF tanh2(const TF x)
    {
        // Rational `tanh` approximation
        return x * (TF(27) + x * x) / (TF(27) + TF(9) * x * x);
    }

    template<typename TF>
    inline TF calc_rain_mass(const TF qr, const TF nr)
    {
        // Calculate mean mass of rain droplets (kg)
        TF mr = qr / std::max(nr, TF(1.));
        return std::min(std::max(mr, mr_min<TF>), mr_max<TF>);
    }

    template<typename TF>
    inline TF calc_rain_diameter(const TF mr)
    {
        // Given mean mass rain drop, calculate mean diameter (m)
        return pow(mr/pirhow<TF>, TF(1.)/TF(3.));
    }

    template<typename TF>
    inline TF calc_mu_r(const TF dr)
    {
        // Calculate shape parameter mu_r
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return TF(10) * (TF(1) + tanh2(TF(1200) * (dr - TF(0.0015)))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    template<typename TF> CUDA_MACRO
    inline TF calc_lambda_r(const TF mur, const TF dr)
    {
        // Calculate Slope parameter lambda_r
        return pow((mur+3)*(mur+2)*(mur+1), TF(1.)/TF(3.)) / dr;
    }

    template<typename TF>
    void prepare_microphysics_slice(
            TF* const restrict rain_mass,
            TF* const restrict rain_diameter,
            TF* const restrict mu_r,
            TF* const restrict lambda_r,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict rho,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jstride + k*kstride;
                const int ij  = i + j*jstride;

                if (qr[ij] > qr_min<TF> * rho[k])  // TODO: remove *rho...
                {
                    rain_mass[ij]     = calc_rain_mass(qr[ij], nr[ij]);
                    rain_diameter[ij] = calc_rain_diameter(rain_mass[ij]);
                    mu_r[ij]          = calc_mu_r(rain_diameter[ij]);
                    lambda_r[ij]      = calc_lambda_r(mu_r[ij], rain_diameter[ij]);
                }
                else
                {
                    rain_mass[ij]     = TF(0);
                    rain_diameter[ij] = TF(0);
                    mu_r[ij]          = TF(0);
                    lambda_r[ij]      = TF(0);
                }
            }
    }


    template<typename TF>
    void autoconversion(
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF Nc0, const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /* Formation of rain by coagulating cloud droplets */

        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES (kg)
        const TF k_cc = 9.44e9;          // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48, same as ICON (m3 kg-2 s-1)
        const TF nu_c = 1;               // SB06, Table 1., same as UCLA-LES (-)
        const TF kccxs = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / fm::pow2(nu_c+1);

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (ql[ijk] > ql_min<TF> * rho[k])  // TODO: remove *rho[k]...
                {
                    // Mean mass of cloud drops (kg):
                    const TF xc = ql[ijk] / Nc0;
                    // Dimensionless internal time scale (SB06, Eq 5):
                    const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ij] + dsmall);
                    // Collection kernel (SB06, Eq 6). Constants are from UCLA-LES and differ from SB06:
                    const TF phi_au = TF(600.) * pow(tau, TF(0.68)) * fm::pow3(TF(1.) - pow(tau, TF(0.68)));
                    // Autoconversion tendency (SB06, Eq 4, kg m-3 s-1):
                    const TF au_tend = rho_0<TF>/rho[k] * kccxs * fm::pow2(ql[ijk]) * fm::pow2(xc) *
                                       (TF(1.) + phi_au / fm::pow2(TF(1.)-tau)); // SB06, eq 4

                    // Update 2D slices:
                    qrt[ij] += au_tend;
                    nrt[ij] += au_tend / x_star;
                }
            }
    }


    template<typename TF>
    void accretion(
            TF* const restrict qrt,
            const TF* const restrict qr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /* Accretion: growth of raindrops collecting cloud droplets */

        const TF k_cr = 5.25; // SB06, p49 (m3 kg-1 s-1)

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (ql[ijk] > ql_min<TF> * rho[k] && qr[ij] > qr_min<TF> * rho[k])  // TODO: remove *rho[k]...
                {
                    // Dimensionless internal time scale (SB06, Eq 5):
                    const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ij]);
                    // Collection kernel (SB06, Eq 8):
                    const TF phi_ac = fm::pow4(tau / (tau + TF(5e-5)));
                    // Accreation tendency (SB06, Eq 7, kg m-3 s-1):
                    const TF ac_tend = k_cr * ql[ijk] *  qr[ij] * phi_ac * sqrt(rho_0<TF> / rho[k]);

                    // Update 2D slice:
                    qrt[ij] += ac_tend;
                }
            }
    }


    template<typename TF>
    void evaporation(
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict qi,
            const TF* const restrict qt,
            const TF* const restrict T,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF* const restrict p,
            const TF* const restrict rain_mass,
            const TF* const restrict rain_diameter,
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Evaporation: evaporation of rain drops in unsaturated environment
        const TF lambda_evap = TF(1.); // 1.0 in UCLA, 0.7 in DALES

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (qr[ij] > qr_min<TF> * rho[k]) // TODO: remove *rho..
                {
                    // Supersaturation over water (-, KP97 Eq 4-37).
                    // TMP BvS, for comparison with old scheme:
                    //const TF qv = qt[ijk] - ql[ijk] - qi[ijk];
                    const TF qv = qt[ijk] - ql[ijk];

                    const TF qs = tmf::qsat_liq(p[k], T[ijk]) * rho[k];
                    const TF S  = qv / qs - TF(1.);

                    if (S < TF(0))
                    {
                        const TF mr  = rain_mass[ij];
                        const TF dr  = rain_diameter[ij];

                        // ...Condensation/evaporation rate...?
                        const TF es = tmf::esat_liq(T[ijk]);
                        const TF Glv = TF(1.) / (Rv<TF> * T[ijk] / (es * D_v<TF>) +
                                                 (Lv<TF> / (K_t<TF> * T[ijk])) * (Lv<TF> / (Rv<TF> * T[ijk]) - TF(1.)));

                        // Ventilation factor. UCLA-LES=1, calculated in SB06 = TODO..
                        const TF F   = TF(1.);

                        // Evaporation tendency (kg m-3 s-1).
                        const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ij];

                        // Update 2D slices:
                        qrt[ij] += ev_tend;
                        nrt[ij] += lambda_evap * ev_tend / mr;
                    }
                }
            }
    }


    template<typename TF>
    void selfcollection_breakup(
            TF* const restrict nrt,
            const TF* const restrict nr,
            const TF* const restrict qr,
            const TF* const restrict rho,
            const TF* const restrict rain_mass,
            const TF* const restrict rain_diameter,
            const TF* const restrict lambda_r,
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (qr[ij] > qr_min<TF> * rho[k]) // TODO: remove *rho..
                {
                    // Short-cuts...
                    const TF dr = rain_diameter[ij];

                    // Selfcollection tendency:
                    // NOTE: this is quite different in ICON, UCLA-LES had 4 different versions, ...
                    const TF sc_tend = -k_rr * nr[ij] * qr[ij] /
                                       pow(TF(1.) + kappa_rr /
                                                    lambda_r[ij] * pow(pirhow<TF>, TF(1.)/TF(3.)), -9) * sqrt(rho_0<TF> / rho[k]);

                    // Update 2D slice:
                    nrt[ij] += sc_tend;

                    // Breakup by collisions:
                    const TF dDr = dr - D_eq;
                    if (dr > TF(0.35e-3))
                    {
                        TF phi_br;
                        if (dr <= D_eq)
                            phi_br = k_br1 * dDr;
                        else
                            phi_br = TF(2.) * exp(k_br2 * dDr) - TF(1.);

                        const TF br_tend = -(phi_br + TF(1.)) * sc_tend;

                        // Update 2D slice:
                        nrt[ij] += br_tend;
                    }
                }
            }
    }
}



template<typename TF>
Microphys_sb06<TF>::Microphys_sb06(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swmicrophys = Microphys_type::SB06;

    // Read microphysics switches and settings
    cfl_max = inputin.get_item<TF>("micro", "cflmax", "", 1.2);
    sw_warm = inputin.get_item<bool>("micro", "swwarm", "", true);
    Nc0 = inputin.get_item<TF>("micro", "Nc0", "");

    auto add_type = [&](
            const std::string& symbol,
            const std::string& short_name,
            const std::string& long_name,
            const std::string& units,
            const bool is_mass)
    {
        Hydro_type<TF> tmp = {short_name, long_name, units, is_mass};
        hydro_types.emplace(symbol, tmp);
    };

    // Switch between mass and density for precipitation types.
    const bool is_mass = true;

    if (sw_warm)
    {
        add_type("qr", "rain", "rain specific humidity", "kg kg-1", is_mass);
        add_type("nr", "rain", "number density rain", "kg-1", !is_mass);
    }
    else
    {
        //add_type("qi", "ice", "ice specific humidity", "kg kg-1", is_mass);
        //add_type("ni", "ice", "number density ice", "kg-1", !is_mass);

        add_type("qr", "rain", "rain specific humidity", "kg kg-1", is_mass);
        add_type("nr", "rain", "number density rain", "kg-1", !is_mass);

        //add_type("qs", "snow", "snow specific humidity", "kg kg-1", is_mass);
        //add_type("ns", "snow", "number density snow", "kg-1", !is_mass);

        //add_type("qg", "graupel", "graupel specific humidity", "kg kg-1", is_mass);
        //add_type("ng", "graupel", "number density graupel", "kg-1", !is_mass);

        //add_type("qh", "hail", "hail specific humidity", "kg kg-1", is_mass);
        //add_type("nh", "hail", "number density hail", "kg-1", !is_mass);
    }

    const std::string group_name = "thermo";
    for (auto& it : hydro_types)
    {
        fields.init_prognostic_field(it.first, it.second.long_name, it.second.units, group_name, gd.sloc);
        fields.sp.at(it.first)->visc = inputin.get_item<TF>("fields", "svisc", it.first);
    }

    // Set the rain and rain_coeff types to the default provided values
    cloud = cloud_nue1mue1;
    rain = rainSBB;
    rain_coeffs = rainSBBcoeffs;

    // bulk ventilation coefficient, Eq. (88) of SB2006
    auto vent_coeff_a = [](const Particle<TF>& parti, const int n)
    {
        const TF vent_coeff_a = parti.a_ven
                * std::tgamma((parti.nu+n+parti.b_geo)/parti.mu)
                / std::tgamma((parti.nu+1.0)/parti.mu)
                * std::pow( std::tgamma((parti.nu+1.0)/parti.mu)
                / std::tgamma((parti.nu+2.0)/parti.mu), (parti.b_geo+n-1.0) );

        return vent_coeff_a;
    };

    // bulk ventilation coefficient, Eq. (89) of SB2006
    auto vent_coeff_b = [](const Particle<TF>& parti, const int n)
    {
        constexpr TF m_f = 0.500; // see PK, S.541. Do not change.

        const TF vent_coeff_b = parti.b_ven
             * std::tgamma((parti.nu+n+(m_f+1.0)*parti.b_geo+m_f*parti.b_vel)/parti.mu)
                         / std::tgamma((parti.nu+1.0)/parti.mu)
                       * std::pow( std::tgamma((parti.nu+1.0)/parti.mu)
                         / std::tgamma((parti.nu+2.0)/parti.mu),
                         ((m_f+1.0)*parti.b_geo+m_f*parti.b_vel+n-1.0) );

        return vent_coeff_b;
    };

    auto moment_gamma = [](const Particle<TF>& p, const int n)
    {
        const TF moment_gamma = std::tgamma((n+p.nu+1.0)/p.mu) / std::tgamma((p.nu+1.0)/p.mu)
                * std::pow( std::tgamma((  p.nu+1.0)/p.mu) / std::tgamma((p.nu+2.0)/p.mu), n);

        return moment_gamma;
    };

    // Setup the subset of the coefficients that does not follow from the copy.
    auto setup_particle_coeffs = [&vent_coeff_a, &vent_coeff_b, &moment_gamma](const Particle<TF>& ptype, Particle_coeffs<TF>& pcoeffs)
    {
        constexpr TF N_sc = 0.710;  // Schmidt-Zahl (PK, S.541)
        constexpr TF n_f = 0.333;   // Exponent von N_sc im Vent-koeff. (PK, S.541)
        constexpr TF nu_l = 1.5e-5; // Kinematic viscosity of air (added by CvH).

        pcoeffs.c_i = 1. / ptype.cap;
        pcoeffs.a_f = vent_coeff_a(ptype, 1);
        pcoeffs.b_f = vent_coeff_b(ptype, 1) * std::pow(N_sc, n_f) / std::sqrt(nu_l);
        pcoeffs.c_z = moment_gamma(ptype, 2);
    };

    // Moved function to after assignment to prevent overwriting.
    setup_particle_coeffs(cloud, cloud_coeffs);
    setup_particle_coeffs(rain, rain_coeffs);

    // Setup autoconversion and selfcollection constants.
    Sb_cold::setup_cloud_autoconversion(cloud, cloud_coeffs);

    // CvH: ICON overrides using the following code, but they are as far as I see the same values.
    // rain_coeffs%cmu0 = cfg_params%rain_cmu0
    // rain_coeffs%cmu1 = cfg_params%rain_cmu1
    // rain_coeffs%cmu3 = cfg_params%rain_cmu3kkk
}

template<typename TF>
Microphys_sb06<TF>::~Microphys_sb06()
{
}


template<typename TF>
void Microphys_sb06<TF>::init()
{
    auto& gd = grid.get_grid_data();

    for (auto& it : hydro_types)
        if (it.second.is_mass)
            it.second.precip_rate.resize(gd.ijcells);
}

template<typename TF>
void Microphys_sb06<TF>::create(
        Input& inputin, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump, Column<TF>& column)
{
    const std::string group_name = "thermo";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        for (auto& it : hydro_types)
        {
            // Time series
            if (it.second.is_mass)
                stats.add_time_series("r" + it.first.substr(1,1), "Mean surface " + it.second.name + "rate", "kg m-2 s-1", group_name);

            // Profiles
            stats.add_prof("v" + it.first, "Fall velocity " + it.second.name + "mass density", "m s-1", "z" , group_name);

            // Tendencies
            stats.add_tendency(*fields.st.at(it.first), "z", tend_name, tend_longname);
        }

        // Thermo tendencies
        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
    }

    if (column.get_switch())
    {
        for (auto& it : hydro_types)
            if (it.second.is_mass)
                column.add_time_series("r" + it.first.substr(1,1), "Surface " + it.second.name + " rate", "kg m-2 s-1");
    }

    // Create cross-sections
    // Variables that this class can calculate/provide:
    std::vector<std::string> allowed_crossvars;
    for (auto& it : hydro_types)
    {
        allowed_crossvars.push_back("r" + it.first.substr(1,1) + "_bot");
    }

    // Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Get thermodynamic variables
    bool cyclic = false;
    bool is_stat = false;

    auto ql = fields.get_tmp();
    auto qi = fields.get_tmp();
    auto T = fields.get_tmp();

    thermo.get_thermo_field(*ql, "ql", cyclic, is_stat);
    thermo.get_thermo_field(*qi, "qi", cyclic, is_stat);
    thermo.get_thermo_field(*T, "T", cyclic, is_stat);

    const std::vector<TF>& p = thermo.get_basestate_vector("p");
    const std::vector<TF>& exner = thermo.get_basestate_vector("exner");

    // TMP/HACK BvS
    //const std::vector<TF>& rho = thermo.get_basestate_vector("rho");
    const std::vector<TF>& rho = fields.rhoref;

    // 2D slices for quantities shared between different kernels.
    auto rain_mass = fields.get_tmp_xy();
    auto rain_diameter = fields.get_tmp_xy();
    auto mu_r = fields.get_tmp_xy();
    auto lambda_r = fields.get_tmp_xy();

    // Setup 2D slices for implicit solver
    const int n_slices = hydro_types.size() * 8;
    if (n_slices > gd.kcells)
        throw std::runtime_error("TODO.... :-)");

    auto tmp_slices = fields.get_tmp();
    std::fill(tmp_slices->fld.begin(), tmp_slices->fld.end(), TF(0));
    int n = 0;
    for (auto& it : hydro_types)
    {
        it.second.v_sed_now = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.v_sed_new = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.flux_now = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.flux_new = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.sum = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.impl = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.slice = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.conversion_tend = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
    }

    // Help functions to zero XY fields.
    auto zero_tmp_xy = [&](std::shared_ptr<std::vector<TF>>& fld_xy)
    {
        std::fill((*fld_xy).begin(), (*fld_xy).end(), TF(0));
    };

    auto convert_units_short = [&](TF* data_ptr, const bool is_to_kgm3)
    {
        Sb_common::convert_unit(
                data_ptr,
                rho.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells,
                is_to_kgm3);
    };

    // Convert all units from `kg kg-1` to `kg m-3` (mass) and `kg-1` to `m-3` (density).
    const bool to_kgm3 = true;
    convert_units_short(ql->fld.data(), to_kgm3);
    convert_units_short(qi->fld.data(), to_kgm3);
    convert_units_short(fields.ap.at("qt")->fld.data(), to_kgm3);
    for (auto& it : hydro_types)
        convert_units_short(fields.ap.at(it.first)->fld.data(), to_kgm3);

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        for (auto& it : hydro_types)
        {
            // Copy 3D fields to 2D slices, and do partial
            // integration of dynamics tendencies.
            Sb_common::copy_slice_and_integrate(
                    it.second.slice,
                    fields.sp.at(it.first)->fld.data(),
                    fields.st.at(it.first)->fld.data(),
                    rho.data(),
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Zero slice which gathers the conversion tendencies
            for (int n=0; n<gd.ijcells; ++n)
                it.second.conversion_tend[n] = TF(0);
        }

        // Density correction fall speeds
        // In ICON, `rhocorr` is written into the cloud/rain/etc particle types as `rho_v`.
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
        const TF rho_corr = std::exp(-Sb_cold::rho_vel<TF>*hlp);

        // Sedimentation velocity rain species.
        Sb_cold::sedi_vel_rain<TF>(
                hydro_types.at("qr").v_sed_now,
                hydro_types.at("nr").v_sed_now,
                hydro_types.at("qr").slice,
                hydro_types.at("nr").slice,
                ql->fld.data(),
                rho.data(),
                rain, rain_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k, use_ql_sedi_rain);

        const TF rdzdt = TF(0.5) * gd.dzi[k] * dt;

        for (auto& it : hydro_types)
            Sb_common::implicit_core(
                    it.second.slice,
                    it.second.sum,
                    it.second.impl,
                    it.second.v_sed_new,
                    it.second.v_sed_now,
                    it.second.flux_new,
                    it.second.flux_now,
                    rdzdt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);

        if (sw_warm)
        {
            // Calculate microphysics processes.
            // These are the "old" `2mom_warm` kernels ported from UCLA-LES/DALES.

            // Autoconversion; formation of rain drop by coagulating cloud droplets.
            warm::autoconversion(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    ql->fld.data(),
                    rho.data(),
                    exner.data(),
                    Nc0, dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Accretion; growth of rain droplets by collecting cloud droplets
            warm::accretion(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("qr").slice,
                    ql->fld.data(),
                    rho.data(),
                    exner.data(),
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Calculate quantities used by multiple kernels.
            warm::prepare_microphysics_slice(
                    (*rain_mass).data(),
                    (*rain_diameter).data(),
                    (*mu_r).data(),
                    (*lambda_r).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    rho.data(),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Evaporation of rain droplets.
            warm::evaporation(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    ql->fld.data(),
                    // CvH, this is temporarily switched back to the sat-adjust qi until ice is enabled.
                    // fields.sp.at("qi")->fld.data(),
                    qi->fld.data(),
                    fields.sp.at("qt")->fld.data(),
                    T->fld.data(),
                    rho.data(),
                    exner.data(),
                    p.data(),
                    (*rain_mass).data(),
                    (*rain_diameter).data(),
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Selfcollection & breakup: growth of raindrops by mutual (rain-rain)
            // coagulation, and breakup by collisions.
            warm::selfcollection_breakup(
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("nr").slice,
                    hydro_types.at("qr").slice,
                    rho.data(),
                    (*rain_mass).data(),
                    (*rain_diameter).data(),
                    (*lambda_r).data(),
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);
        }
        else
        {
            // Calculate microphysics processes.
            // These are the new kernels ported from ICON.

            // Autoconversion; formation of rain drop by coagulating cloud droplets.
            Sb_cold::autoconversionSB(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    ql->fld.data(),
                    cloud_coeffs,
                    cloud, rain,
                    rho_corr,
                    Nc0,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);

            Sb_cold::accretionSB(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("qr").slice,
                    ql->fld.data(),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);

            Sb_cold::rain_selfcollectionSB(
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    rain,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
        }

        for (auto& it : hydro_types)
        {
            // Integrate conversion tendencies into qr/Nr slices before implicit step.
            Sb_common::integrate_process(
                    it.second.slice,
                    it.second.conversion_tend,
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);

            // Implicit sedimentation step
            Sb_common::implicit_time(
                    it.second.slice,
                    it.second.sum,
                    it.second.impl,
                    it.second.v_sed_new,
                    it.second.v_sed_now,
                    it.second.flux_new,
                    rdzdt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);

            // Evaluate total tendencies back from implicit solver.
            Sb_common::diagnose_tendency(
                    fields.st.at(it.first)->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    it.second.slice,
                    rho.data(),
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
        }

        // Calculate thermodynamic tendencies `thl` and `qt`,
        // from microphysics tendencies excluding sedimentation.
        Sb_common::calc_thermo_tendencies(
                fields.st.at("thl")->fld.data(),
                fields.st.at("qt")->fld.data(),
                hydro_types.at("qr").conversion_tend,
                rho.data(),
                exner.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k);
    }

    for (auto& it : hydro_types)
    {
        // Store surface precipitation rates.
        if (it.second.is_mass)
            for (int n=0; n<gd.ijcells; ++n)
                it.second.precip_rate[n] = it.second.flux_new[n];

        // Convert all units back from `kg m-3` to `kg kg-1` (mass) and `m-3` to `kg-1` (density)
        convert_units_short(fields.ap.at(it.first)->fld.data(), !to_kgm3);
    }

    // Convert specific humidity from `kg m-3` to `kg kg-`
    convert_units_short(fields.ap.at("qt")->fld.data(), !to_kgm3);

    // Calculate tendencies.
    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt" ), tend_name);
    for (auto& it : hydro_types)
        stats.calc_tend(*fields.st.at(it.first), tend_name);

    // Release temporary fields.
    fields.release_tmp(ql);
    fields.release_tmp(qi);
    fields.release_tmp(T);
    fields.release_tmp(tmp_slices);

    fields.release_tmp_xy(rain_mass);
    fields.release_tmp_xy(rain_diameter);
    fields.release_tmp_xy(mu_r);
    fields.release_tmp_xy(lambda_r);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    const bool is_stat = true;
    const bool cyclic = false;

    // Time series
    for (auto& it : hydro_types)
        if (it.second.is_mass)
            stats.calc_stats_2d("r" + it.first.substr(1,1), it.second.precip_rate, no_offset);

    // Profiles
    auto vq = fields.get_tmp();
    auto vn = fields.get_tmp();

    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Sedimentation
        // Density correction fall speeds
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
        const TF rho_corr = std::exp(-Sb_cold::rho_vel<TF> * hlp);

        Sb_cold::sedi_vel_rain<TF>(
                &vq->fld.data()[k * gd.ijcells],
                &vn->fld.data()[k * gd.ijcells],
                &fields.sp.at("qr")->fld.data()[k*gd.ijcells],
                &fields.sp.at("nr")->fld.data()[k*gd.ijcells],
                nullptr,
                rho.data(),
                rain, rain_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k, use_ql_sedi_rain);
    }

    stats.calc_stats("vqr", *vq, no_offset, no_threshold);
    stats.calc_stats("vnr", *vn, no_offset, no_threshold);

    fields.release_tmp(vq);
    fields.release_tmp(vn);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;

    for (auto& it : hydro_types)
        if (it.second.is_mass)
            column.calc_time_series("r" + it.first.substr(1,1), it.second.precip_rate.data(), no_offset);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (cross.get_switch())
    {
        for (auto& it : crosslist)
        {
            // Yikes (BvS), this is easy to break.... Any better ideas? Should we stick to
            // using the full name, like e.g. `rain_rate` or `qr_flux` or ...?
            const std::string letter = it.substr(1,1);
            cross.cross_plane(hydro_types.at("q" + letter).precip_rate.data(), it, iotime);
        }
    }
}

#ifndef USECUDA
template<typename TF>
unsigned long Microphys_sb06<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();

    double cfl = 0.;

    // TO-DO

    // Get maximum CFL across all MPI tasks
    master.max(&cfl, 1);
    fields.release_tmp(tmp);

    // Prevent zero division.
    cfl = std::max(cfl, 1.e-5);

    return idt * this->cfl_max / cfl;
}
#endif

template<typename TF>
bool Microphys_sb06<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_sb06<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    std::string message = "SB06 microphysics scheme can not provide mask: \"" + mask_name +"\"";
    throw std::runtime_error(message);
}

template<typename TF>
void Microphys_sb06<TF>::get_surface_rain_rate(std::vector<TF>& field)
{
    auto& gd = grid.get_grid_data();

    for (int n=0; n<gd.ijcells; ++n)
        field[n] = TF(0);

    for (auto& it : hydro_types)
    {
        if (it.second.is_mass)
        {
            for (int n=0; n<gd.ijcells; ++n)
                field[n] += it.second.precip_rate[n];
        }
    }
}

template class Microphys_sb06<double>;
template class Microphys_sb06<float>;
