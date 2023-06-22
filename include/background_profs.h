//
// Created by Mirjam Tijhuis on 27/03/2023.
//

#ifndef BACKGROUND_PROFS_H
#define BACKGROUND_PROFS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <memory>
#include "Types.h"
#include "Source_functions.h"
#include "Gas_concs.h"

using Aerosol_concs = Gas_concs;

class Master;
class Input;
class Netcdf_handle;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timedep;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Thermo;
template<typename> class Field3d;

template<typename TF>
class Background
{
public:
    Background(Master&, Grid<TF>&, Fields<TF>&, Input&);
    ~Background();

    void init(Netcdf_handle&, Timeloop<TF>&);
    void create(Input&, Netcdf_handle&, Stats<TF>&);
    void exec(Thermo<TF>&);
    void exec_stats(Stats<TF>&);
    void update_time_dependent(Timeloop<TF>&);

    // GPU functions and variables
    void prepare_device();
    void clear_device();
    std::map<std::string, TF*> gasprofs_g;    ///< Map of profiles with gasses stored by its name.

    void get_tpm(Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Gas_concs&);
    void get_gasses(Gas_concs&);
    void get_aerosols(Aerosol_concs&);

private:
    Master& master;
    Grid<TF>& grid;
    Fields<TF>& fields;

    // Case switches
    bool sw_update_background;
    bool sw_aerosol;
    bool sw_aerosol_timedep;
    double dt_rad;
    unsigned long idt_rad;

//    const TF n_era_layers = 136;
//    const TF n_era_levels = 137;
    TF n_era_layers;
    TF n_era_levels;

    std::map<std::string, Timedep<TF>*> tdep_gases;
    std::vector<std::string> gaslist;        ///< List of gases that have timedependent background profiles.
    std::map<std::string, std::vector<TF>> gasprofs; ///< Map of profiles with gases stored by its name.

    // Arrays
    // to fill with input values
    // temperature, pressure and moisture
    std::vector<TF> t_lay;
    std::vector<TF> t_lev;
    std::vector<TF> p_lay;
    std::vector<TF> p_lev;
    std::vector<TF> h2o;
    //aerosols
    std::vector<TF> aermr01;
    std::vector<TF> aermr02;
    std::vector<TF> aermr03;
    std::vector<TF> aermr04;
    std::vector<TF> aermr05;
    std::vector<TF> aermr06;
    std::vector<TF> aermr07;
    std::vector<TF> aermr08;
    std::vector<TF> aermr09;
    std::vector<TF> aermr10;
    std::vector<TF> aermr11;

    std::unique_ptr<Timedep<TF>> tdep_t_lay;
    std::unique_ptr<Timedep<TF>> tdep_t_lev;
    std::unique_ptr<Timedep<TF>> tdep_p_lay;
    std::unique_ptr<Timedep<TF>> tdep_p_lev;
    std::unique_ptr<Timedep<TF>> tdep_h2o;
    std::unique_ptr<Timedep<TF>> tdep_o3;
    std::unique_ptr<Timedep<TF>> tdep_aermr01;
    std::unique_ptr<Timedep<TF>> tdep_aermr02;
    std::unique_ptr<Timedep<TF>> tdep_aermr03;
    std::unique_ptr<Timedep<TF>> tdep_aermr04;
    std::unique_ptr<Timedep<TF>> tdep_aermr05;
    std::unique_ptr<Timedep<TF>> tdep_aermr06;
    std::unique_ptr<Timedep<TF>> tdep_aermr07;
    std::unique_ptr<Timedep<TF>> tdep_aermr08;
    std::unique_ptr<Timedep<TF>> tdep_aermr09;
    std::unique_ptr<Timedep<TF>> tdep_aermr10;
    std::unique_ptr<Timedep<TF>> tdep_aermr11;

    // GPU functions and variables
    TF* t_lay_g;    ///< Pointer to GPU array.
    TF* t_lev_g;
    TF* p_lay_g;
    TF* p_lev_g;
    TF* h2o_g;
    TF* aermr01_g;
    TF* aermr02_g;
    TF* aermr03_g;
    TF* aermr04_g;
    TF* aermr05_g;
    TF* aermr06_g;
    TF* aermr07_g;
    TF* aermr08_g;
    TF* aermr09_g;
    TF* aermr10_g;
    TF* aermr11_g;
};

#endif //BACKGROUND_PROFS_H