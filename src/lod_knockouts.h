
#ifndef _EALIFE_LOD_KNOCKOUTS_H_
#define _EALIFE_LOD_KNOCKOUTS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>


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
            .add_field("rx_knockedout")
            .add_field("location_knockedout");
            
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control_ea = ea.make_individual();
                control_ea->rng().reset(get<RNG_SEED>(**i));
                
                typename EA::individual_ptr_type knockout_rx_ea = ea.make_individual();
                knockout_rx_ea->rng().reset(get<RNG_SEED>(**i));
                knockout<instructions::rx_msg,instructions::nop_x>(*knockout_rx_ea);
                
                typename EA::individual_ptr_type knockout_location_ea = ea.make_individual();
                knockout_location_ea->rng().reset(get<RNG_SEED>(**i));
                knockout<instructions::get_xy,instructions::nop_x>(*knockout_location_ea);
                
                // setup the founder
                typename EA::individual_type::individual_ptr_type o= (*i)->make_individual((*i)->founder().repr());
                o->hw().initialize();
                control_ea->append(o);
                
                // setup send knockout founder
                typename EA::individual_type::individual_ptr_type ko_s = (*i)->make_individual((*i)->founder().repr());
                ko_s->hw().initialize();
                knockout_rx_ea->append(ko_s);
                
                // setup location knockout founder
                typename EA::individual_type::individual_ptr_type ko_l = (*i)->make_individual((*i)->founder().repr());
                ko_l->hw().initialize();
                knockout_location_ea->append(ko_l);
                
                int update_max = get<METAPOP_COMPETITION_PERIOD>(*control_ea);
                
                for (int i=0; i<=update_max; ++i) {
                    control_ea->update();
                    knockout_rx_ea->update();
                    knockout_location_ea->update();
                }
                
                // reeval fitness here... write fitness output
      
                
                df.endl();
                
                ++lod_depth;
            }
        }
    }
}
#endif
