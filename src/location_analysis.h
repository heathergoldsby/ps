#ifndef _PS_LOC_ANALYSIS_H_
#define _PS_LOC_ANALYSIS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ea/datafile.h>
#include <ea/line_of_descent.h>
#include <ea/analysis.h>
#include <ea/digital_evolution/utils/task_switching.h>


namespace ealib {
    namespace analysis {
        
        
        /*! lod_movie
         */
        LIBEA_ANALYSIS_TOOL(location_analysis) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::iterator i=lod.begin(); ++i;
            
            datafile df("location_analysis");
            
            int lod_depth = 0;
            // skip def ancestor (that's what the +1 does)
            for( ; i!=lod.end(); ++i) {
                
                df.write(lod_depth);
                
                // **i is the EA, AS OF THE TIME THAT IT DIED!
                
                // To replay, need to create new eas for each knockout exper.
                // setup the population (really, an ea):
                typename EA::individual_ptr_type control = ea.make_individual();
                control->ea().rng().reset(get<RNG_SEED>(i->ea()));
                
                
                // Setup founders!
                
                for(typename EA::individual_type::ea_type::population_type::iterator j=i->ea().founder().begin(); j!=i->ea().founder().end(); ++j) {
                    
                    typename EA::individual_type::ea_type::individual_ptr_type o1 = i->ea().copy_individual(**j);
                    o1->hw().initialize();
                    control->ea().insert(control->ea().end(), o1);
                }
                
                int update_max = get<METAPOP_COMPETITION_PERIOD>(control->ea());
                
                
                for (int i=0; i<=update_max; ++i) {
                    control->ea().update();
                }
                
                
                df.write(control->ea().size());
                if (control->ea().size() >= 36) {
                for(typename EA::individual_type::ea_type::population_type::iterator j=control->ea().population().begin(); j!=control->ea().population().end(); ++j) {
                        int x =(control->ea().env().location((**j).position())->x);
                        int y =(control->ea().env().location((**j).position())->y);
                    
                        df.write(x);
                        df.write(y);
                        df.endl();
                    }
                }
                df.endl();
                df.endl();
                
            }
          
           

            
            df.endl();
            
        }
    }
}
#endif
