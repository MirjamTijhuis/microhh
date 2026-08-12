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

#ifndef BUDGET3D_H
#define BUDGET3D_H

#include <string>
#include <vector>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

/**
 * Accumulates the total tendency of a set of prognostic variables into
 * dedicated 3D diagnostic fields ("<var>_tend"), time-weighted over the
 * RK substeps so the result is the true mean tendency over the interval
 * between two dumps. The accumulated fields are ordinary diagnostic
 * fields, so they are written out through the existing [dump] mechanism
 * (dumplist / dumplist_coarse) with no additional plumbing.
 */
template<typename TF>
class Budget3d
{
    public:
        Budget3d(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Budget3d();

        bool get_switch() const { return swbudget3d; }

        void exec(const double sub_dt);   ///< Accumulate one RK substep's tendency.
        void prepare_dump();              ///< Normalize accumulators into a mean tendency before writing.
        void reset();                     ///< Zero the accumulators after writing.

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool swbudget3d;
        std::vector<std::string> tendencylist;
        double elapsed_time;
};
#endif
