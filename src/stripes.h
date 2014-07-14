//
//  stripes.h
//  ealife
//
//  Created by Heather Goldsby on 11/22/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//


#ifndef _EALIFE_STRIPES_H_
#define _EALIFE_STRIPES_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>


#include <ea/digital_evolution.h>
#include <ea/digital_evolution/hardware.h>
#include <ea/digital_evolution/instruction_set.h>
#include <ea/digital_evolution/discrete_spatial_environment.h>
#include <ea/digital_evolution/utils/resource_consumption.h>
#include <ea/digital_evolution/utils/task_switching.h>
#include <ea/digital_evolution/population_founder.h>
#include <ea/datafiles/reactions.h>
#include <ea/cmdline_interface.h>
#include <ea/metapopulation.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/selection/tournament.h>
#include <ea/mutation.h>
#include <ea/recombination.h>

using namespace ealib;


LIBEA_MD_DECL(STRIPE_FIT, "ea.stripes.fit", int);
LIBEA_MD_DECL(STRIPE_FIT_PT, "ea.strip    es.fit.pt", int);

LIBEA_MD_DECL(ANCESTOR, "ea.stripes.ancestor", int);
LIBEA_MD_DECL(NUM_PROPAGULE_CELL, "ea.stripes.num_propagule_cell", int);
LIBEA_MD_DECL(STRIPE_FIT_FUNC, "ea.stripes.fit_func", int); // 0 = sum; 1 = prod
LIBEA_MD_DECL(MULTICELL_REP_TIME, "ea.stripes.mc_rep_time", int);
LIBEA_MD_DECL(MC_RESOURCE_UNITS, "ea.stripes.mc_res_units", int);
LIBEA_MD_DECL(FIT_MAX, "ea.stripes.fit_max", int);
LIBEA_MD_DECL(FIT_MIN, "ea.stripes.fit_min", int);
LIBEA_MD_DECL(FIT_GAMMA, "ea.stripes.fit_gamma", double);
LIBEA_MD_DECL(RES_UPDATE, "ea.stripes.res_update", double);



template <typename EA>
void eval_permute_stripes(EA& ea) {
    // vert stripes
    double one_fit_not = 0;
    double one_fit_nand = 0;
    double two_fit_not = 0;
    double two_fit_nand = 0;
    // horizontal stripes
    double three_fit_not = 0;
    double three_fit_nand = 0;
    double four_fit_not = 0;
    double four_fit_nand = 0;
    // diagonal stripes
    double five_fit_not = 0;
    double five_fit_nand = 0;
    double six_fit_not = 0;
    double six_fit_nand = 0;
    
    
    for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
        for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
            typename EA::environment_type::location_ptr_type l = ea.env().location(x,y);
            if (!l->occupied()) {
                continue;
            }
            
            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");
            // Vertical stripes!
            if((y % 2) == 0) {
                if (lt == "nand") { ++one_fit_nand; }
                if (lt == "not") { ++two_fit_not; }
            } else {
                if(lt == "not") { ++one_fit_not; }
                if (lt == "nand") { ++two_fit_nand; }
            }
            
            // Horizontal stripes
            if ((x % 2) == 0) {
                if (lt == "nand") { ++three_fit_nand; }
                if (lt == "not") { ++four_fit_not; }
            } else {
                if(lt == "not") { ++three_fit_not; }
                if (lt == "nand") { ++four_fit_nand; }
            }
            
            // Diagonal stripes
            if (((x % 2) == 0) && ((y % 2) == 0 )) {
                if(lt == "not") { ++five_fit_not; }
                if (lt == "nand") { ++six_fit_nand; }
            } else {
                if(lt == "nand") { ++five_fit_nand; }
                if (lt == "not") { ++six_fit_not; }
            }
        }
    }
    
    // fix to use fit func for prod or sum
    double tmp_one_fit;
    double tmp_two_fit;
    double tmp_three_fit;
    double tmp_four_fit;
    double tmp_five_fit;
    double tmp_six_fit;
    double min_fit;
    double max_fit;
    
    switch (get<STRIPE_FIT_FUNC>(ea,1)) {
        case 0: {
            tmp_one_fit = one_fit_not +  one_fit_nand;
            tmp_two_fit = two_fit_not + two_fit_nand;
            tmp_three_fit = three_fit_not + three_fit_nand;
            tmp_four_fit = four_fit_not + four_fit_nand;
            tmp_five_fit = five_fit_not + five_fit_nand;
            tmp_six_fit = six_fit_not + six_fit_nand;
            min_fit = get<POPULATION_SIZE>(ea) / 2;
            max_fit = get<POPULATION_SIZE>(ea);
            break;
        }
        case 1:
        default:
            tmp_one_fit = (one_fit_not + 1)  * (one_fit_nand + 1);
            tmp_two_fit = (two_fit_not + 1)  * (two_fit_nand + 1);
            tmp_three_fit = (three_fit_not + 1)  * (three_fit_nand + 1);
            tmp_four_fit = (four_fit_not + 1)  * (four_fit_nand + 1);
            tmp_five_fit = (five_fit_not + 1)  * (five_fit_nand + 1);
            tmp_six_fit = (six_fit_not + 1)  * (six_fit_nand + 1);
            min_fit = 0;
            max_fit = pow(((get<POPULATION_SIZE>(ea) / 2) + 1), 2);
            break;
    }
    
    
    double tmp_fit = std::max(tmp_one_fit, tmp_two_fit);
    tmp_fit = std::max(tmp_fit, tmp_three_fit);
    tmp_fit = std::max(tmp_fit, tmp_four_fit);
    tmp_fit = std::max(tmp_fit, tmp_five_fit);
    tmp_fit = std::max(tmp_fit, tmp_six_fit);
    
    if (tmp_fit < min_fit) {
        tmp_fit = min_fit;
    }
    
    
    put<STRIPE_FIT_PT>(tmp_fit, ea);
    // rescale fitness!
    double rescaled_fit = (get<FIT_MAX>(ea) - get<FIT_MIN>(ea)) * pow (((tmp_fit - min_fit) / (max_fit - min_fit)), (get<FIT_GAMMA>(ea))) + get<FIT_MIN>(ea);
    
    
    put<STRIPE_FIT>(rescaled_fit,ea);
}





/*! Compete to evolve stripes -- even number rows nand; odd number rows not
 */
template <typename EA>
struct permute_stripes : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    permute_stripes(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("permute_stripes.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness");
    }
    
    virtual ~permute_stripes() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        
        
        // calculate "fitness":
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            
            eval_permute_stripes(i->ea());
            
            // copy the stripe fit to the accumulator and also the subpop
            double sf =get<STRIPE_FIT>(i->ea());
            fit(sf);
            put<STRIPE_FIT>(sf, *i);
            
        }
        
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .endl();
        
        std::size_t n=get<META_POPULATION_SIZE>(ea);
        typename EA::population_type offspring; // container of (pointers to) subpopulations
        recombine_n(ea.population(), offspring,
                    selection::tournament < access::meta_data<STRIPE_FIT> > (n, ea.population(), ea),
                    recombination::propagule_without_replacement(),
                    n, ea);
        
        
        
        configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
        int s = get<POPULATION_SIZE>(ea);
        
        // Mutate and fill each offspring multicell.
        for(typename EA::population_type::iterator i=offspring.begin(); i!=offspring.end(); ++i) {
            assert((*i)->ea().population().size() == 1);
            
            // clear founders...
            (*i)->ea().founder().clear();
            
            // mutate it:
            mutate(**((*i)->ea().population().begin()),m,(*i)->ea());
            typename EA::individual_type::ea_type::individual_type g = (**((*i)->ea().population().begin()));
            
            // add first org as founder
            (*i)->ea().founder().insert((*i)->ea().founder().end(), (*i)->ea().copy_individual(g));
            
            // and fill up the offspring population with copies of the germ:
            for (int k=1; k<get<NUM_PROPAGULE_CELL>(ea); ++k) {
                typename EA::individual_type::ea_type::individual_ptr_type o = (*i)->ea().copy_individual(g);
                (*i)->insert((*i)->end(), o);
                
                // move to random location
                std::size_t pos = (*i)->ea().rng()(s);
                (*i)->ea().env().swap_locations(k, pos);
                
                // add org as founders
                (*i)->ea().founder().insert((*i)->ea().founder().end(), (*i)->ea().copy_individual(*o));
            }
        }
        
        // swap populations
        std::swap(ea.population(), offspring);
        
        
    }
    
    datafile _df;
};

//! Performs multicell replication using germ lines.
template <typename EA>
struct stripes_replication : end_of_update_event<EA> {
    //! Constructor.
    stripes_replication(EA& ea) : end_of_update_event<EA>(ea), _df("stripes_rep.dat") {
        _df.add_field("update")
        .add_field("mean_rep_time")
        .add_field("last_fitness")
        .add_field("last_fitness_pt")
        .add_field("multicell_size")
        .add_field("replication_count");
        num_rep = 0;
    }
    
    
    //! Destructor.
    virtual ~stripes_replication() {
    }
    
    //! Perform germline replication among populations.
    virtual void operator()(EA& ea) {
        
            // See if any subpops have exceeded the resource threshold
            typename EA::population_type offspring;
            for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
                
                // Calc fitness for each subpop
                eval_permute_stripes(i->ea());
                
                // copy the stripe fit to the accumulator and also the subpop
                double sf =get<STRIPE_FIT>(i->ea());
                put<STRIPE_FIT>(sf, *i);
                get<MC_RESOURCE_UNITS>(*i,0) += sf;
                
                // track time since group rep
                get<MULTICELL_REP_TIME>(*i,0) +=1;
                
                if (get<MC_RESOURCE_UNITS>(*i) > get<GROUP_REP_THRESHOLD>(*i)){
                    
                    // find a germ -- we are picking the first cell, since they are genetically identical, it is ok.
                    typename EA::individual_type::ea_type::individual_type& germ= **(i->ea().population().begin());
                    germ.repr().resize(germ.hw().original_size());
                    germ.hw().initialize();
                    
                    // setup the population (really, an ea):
                    typename EA::individual_ptr_type p = ea.make_individual();
                    
                    // mutate it:
                    configurable_per_site m(get<GERM_MUTATION_PER_SITE_P>(ea));
                    mutate(germ,m,p->ea());
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o=p->ea().copy_individual(germ.repr());
                    inherits_from(germ, *o, p->ea());
                    
                    // clear founders...
                    p->ea().founder().clear(); // offspring
                    
                    
                    int s = get<POPULATION_SIZE>(ea);
                    // and fill up the offspring population with copies of the germ:
                    for (int k=0; k<get<NUM_PROPAGULE_CELL>(ea); ++k) {
                        typename EA::individual_type::ea_type::individual_ptr_type o2 = p->ea().copy_individual(*o);
                        p->insert(p->end(), o2);
                        
                        // move to random location
                        std::size_t pos = p->ea().rng()(s);
                        p->ea().env().swap_locations(k, pos);
                        
                        // add org as founders
                        //p->ea().founder().insert(p->ea().founder().end(), p->ea().copy_individual(*o2));
                    }
                    

                    
                    offspring.push_back(p);
                    
                    // track last fitness AND time to rep
                    multicell_rep.push_back(get<MULTICELL_REP_TIME>(*i));
                    multicell_last_fitness.push_back(get<STRIPE_FIT>(*i));
                    multicell_size.push_back(i->ea().size());
                    multicell_last_fitness_pt.push_back(get<STRIPE_FIT_PT>(i->ea()));
                    ++num_rep;
                    
                    // reset parent multicell
                    i->ea().env().reset_resources();
                    put<MC_RESOURCE_UNITS>(0,*i);
                    put<MULTICELL_REP_TIME>(0,*i);
                    
                    i->ea().clear(); // kills existing population
                    i->ea().env().clear_env(i->ea());
                    int count =0;
                    
                    for(typename EA::individual_type::ea_type::population_type::iterator j=i->ea().founder().begin(); j!=i->ea().founder().end(); ++j) {
                        int s = get<POPULATION_SIZE>(i->ea());
                        std::size_t pos = i->ea().rng()(s);
                        
                        typename EA::individual_type::ea_type::individual_ptr_type o1 = i->ea().copy_individual(**j);
                        o1->hw().initialize();
                        i->ea().insert(i->ea().end(), o1);
                        i->ea().env().swap_locations(count, pos);
                        ++count;
                    }
                    

                    // i == parent individual;
                    typename EA::population_type parent_pop, offspring_pop;
                    parent_pop.push_back(*i.base());
                    offspring_pop.push_back(p);
                    inherits(parent_pop, offspring_pop, ea);
                    
                    if (multicell_rep.size() > 100) {
                        multicell_rep.pop_front();
                        multicell_last_fitness.pop_front();
                        multicell_last_fitness_pt.pop_front();
                        multicell_size.pop_front();
                    }
                    
                }
            }
            
            
            // select surviving parent groups
            if (offspring.size() > 0) {
                int n = get<META_POPULATION_SIZE>(ea) - offspring.size();
                
                typename EA::population_type survivors;
                select_n<selection::random< > >(ea.population(), survivors, n, ea);
                
                // add the offspring to the list of survivors:
                survivors.insert(survivors.end(), offspring.begin(), offspring.end());
                
                // and swap 'em in for the current population:
                std::swap(ea.population(), survivors);
            }
            
        
        if ((ea.current_update() % 100) == 0) {
            if (multicell_rep.size() > 0) {
                _df.write(ea.current_update())
                .write(std::accumulate(multicell_rep.begin(), multicell_rep.end(), 0.0)/multicell_rep.size())
                .write(std::accumulate(multicell_last_fitness.begin(), multicell_last_fitness.end(), 0.0)/multicell_last_fitness.size())
                .write(std::accumulate(multicell_last_fitness_pt.begin(), multicell_last_fitness_pt.end(), 0.0)/multicell_last_fitness_pt.size())
                .write(std::accumulate(multicell_size.begin(), multicell_size.end(), 0.0)/multicell_size.size())
                .write(num_rep)
                .endl();
                num_rep = 0;
            } else {
                _df.write(ea.current_update())
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(0.0)
                .write(num_rep)
                .endl();
            }
        }
        
        }
    
    datafile _df;
    std::deque<double> multicell_rep;
    std::deque<double> multicell_last_fitness;
    std::deque<double> multicell_last_fitness_pt;
    std::deque<double> multicell_size;
    
    int num_rep;
    
    
};

#endif
