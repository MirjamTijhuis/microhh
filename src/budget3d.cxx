/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include <algorithm>
#include <stdexcept>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "budget3d.h"

namespace
{
    template<typename TF>
    void accumulate_tendency(
            TF* const restrict acc,
            const TF* const restrict tend,
            const TF sub_dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    acc[ijk] += sub_dt * tend[ijk];
                }
    }


    template<typename TF>
    void scale_field(TF* const restrict fld, const TF fac, const int ncells)
    {
        #pragma ivdep
        for (int n=0; n<ncells; ++n)
            fld[n] *= fac;
    }
}

template<typename TF>
Budget3d<TF>::Budget3d(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin),
    elapsed_time(0.)
{
    swbudget3d = inputin.get_item<bool>("budget3d", "swbudget3d", "", false);

    if (swbudget3d)
    {
        tendencylist = inputin.get_list<std::string>("budget3d", "tendencylist", "", std::vector<std::string>());

        if (tendencylist.empty())
            throw std::runtime_error("[budget3d] tendencylist is empty");

        for (auto& varname : tendencylist)
        {
            if (fields.at.find(varname) == fields.at.end())
                throw std::runtime_error("budget3d variable \"" + varname + "\" in tendencylist is not a prognostic variable");

            auto& tend = *fields.at.at(varname);
            fields.init_diagnostic_field(
                    varname + "_tend", "Accumulated mean tendency of " + varname,
                    tend.unit, "budget3d", tend.loc);
        }
    }
}

template<typename TF>
Budget3d<TF>::~Budget3d()
{
}

#ifndef USECUDA
template<typename TF>
void Budget3d<TF>::exec(const double sub_dt)
{
    if (!swbudget3d)
        return;

    auto& gd = grid.get_grid_data();
    const TF sub_dt_tf = static_cast<TF>(sub_dt);

    for (auto& varname : tendencylist)
    {
        auto& acc  = *fields.sd.at(varname + "_tend");
        auto& tend = *fields.at.at(varname);

        accumulate_tendency(
                acc.fld.data(),
                tend.fld.data(),
                sub_dt_tf,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    elapsed_time += sub_dt;
}
#endif

template<typename TF>
void Budget3d<TF>::prepare_dump()
{
    if (!swbudget3d || elapsed_time <= 0.)
        return;

    // At this point `fld` on the host is up to date: the generic dump/cross/stats
    // block in Model::exec() always calls fields->backward_device() first.
    auto& gd = grid.get_grid_data();
    const TF dti = TF(1.) / static_cast<TF>(elapsed_time);

    for (auto& varname : tendencylist)
    {
        auto& acc = *fields.sd.at(varname + "_tend");
        scale_field(acc.fld.data(), dti, gd.ncells);
    }
}

template<typename TF>
void Budget3d<TF>::reset()
{
    if (!swbudget3d)
        return;

    for (auto& varname : tendencylist)
    {
        auto& acc = *fields.sd.at(varname + "_tend");
        std::fill(acc.fld.begin(), acc.fld.end(), TF(0));

        #ifdef USECUDA
        fields.forward_field_device_3d(acc.fld_g, acc.fld.data());
        #endif
    }

    elapsed_time = 0.;
}


#ifdef FLOAT_SINGLE
template class Budget3d<float>;
#else
template class Budget3d<double>;
#endif
