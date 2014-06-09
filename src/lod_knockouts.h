
#ifndef _EALIFE_LOD_KNOCKOUTS_H_
#define _EALIFE_LOD_KNOCKOUTS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>

#include "stripes.h"



namespace ealib {
    namespace analysis {
        
                
        /*! lod_knockouts reruns each subpopulation along a line of descent and records how the subpopulation
         fares with key coordination instructions removed.
         
         */
        LIBEA_ANALYSIS_TOOL(lod_knockouts) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("lod_knockouts.dat");
            df.add_field("lod_depth")
            .add_field("no_knockouts")
            .add_field("rx_knockedout");
//            .add_field("location_knockedout");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control = ea.make_individual();
                control->ea().rng().reset(get<RNG_SEED>(i->ea()));
                
                
                typename EA::individual_ptr_type knockout_rx = ea.make_individual();
                knockout_rx->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<instructions::rx_msg,instructions::nop_x>(knockout_rx->ea());
                
                /*
                typename EA::individual_ptr_type knockout_location = ea.make_individual();
                knockout_location->ea().rng().reset(get<RNG_SEED>(i->ea()));
                knockout<instructions::get_xy,instructions::nop_x>(knockout_location->ea());
                */
                
                // Setup founders!
                int count =0;
                
                for(typename EA::individual_type::ea_type::population_type::iterator j=i->ea().founder().begin(); j!=i->ea().founder().end(); ++j) {
                    
                    
                    int s = get<POPULATION_SIZE>(i->ea());
                    std::size_t pos = i->ea().rng()(s);
                    
                    
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o1 = i->ea().copy_individual(**j);
                    o1->hw().initialize();
                    control->ea().insert(control->ea().end(), o1);
                    control->ea().env().swap_locations(count, pos);

                    
                    typename EA::individual_type::ea_type::individual_ptr_type o2 = i->ea().copy_individual(**j);
                    o2->hw().initialize();
                    knockout_rx->ea().insert(knockout_rx->ea().end(), o2);

                    /*
                    typename EA::individual_type::ea_type::individual_ptr_type o3 = i->ea().copy_individual(**j);
                    o3->hw().initialize();
                    knockout_location->ea().insert(knockout_location->ea().end(), o3);
                     */
                    ++count;
                    
                }
                
                int update_max = get<METAPOP_COMPETITION_PERIOD>(control->ea());
                
                
                for (int i=0; i<=update_max; ++i) {
                    control->ea().update();
                    knockout_rx->ea().update();
                    //knockout_location->ea().update();
                }
                
                // re-eval fitness here... write fitness output
                eval_permute_stripes(control->ea());
                eval_permute_stripes(knockout_rx->ea());
                //eval_permute_stripes(knockout_location->ea());
                
                df.write(get<STRIPE_FIT>(control->ea(),0));
                df.write(get<STRIPE_FIT>(knockout_rx->ea(),0));
                //df.write(get<STRIPE_FIT>(knockout_location->ea(),0));
      
                
                df.endl();
                
                
                ++lod_depth;
            }
            
            
        }
        
        
    }
}
#endif
