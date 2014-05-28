#ifndef _PS_MOVIE_H_
#define _PS_MOVIE_H_

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
        LIBEA_ANALYSIS_TOOL(movie) {
            
            line_of_descent<EA> lod = lod_load(get<ANALYSIS_INPUT>(ea), ea);
            
            typename line_of_descent<EA>::reverse_iterator i=lod.rbegin(); ++i;
            datafile df("movie.dat");
            typename EA::individual_ptr_type control_sp = ea.make_individual();
            control_sp->ea().rng().reset(get<RNG_SEED>((i->ea())));
            
            
            
            typename EA::individual_type::ea_type::individual_type g (i->ea().founder());
            typename EA::individual_type::ea_type::individual_ptr_type o = i->ea().copy_individual(g);

//            typename EA::individual_type::ea_type::individual_ptr_type o= i->ea().make_individual(i->ea().founder().repr());
            o->hw().initialize();
//            inherits_from(i->ea().founder(), *q, p->ea());

            
            control_sp->ea().insert(control_sp->ea().end(), o);
            

            
            int update_max = get<METAPOP_COMPETITION_PERIOD>(ea);

            df.write(get<SPATIAL_X>(ea));
            df.write(get<SPATIAL_Y>(ea));
            df.endl();
            for (int i=0; i<=update_max; ++i) {
                control_sp->ea().update();
                df.write(i);
                // grab info based on location...
                for (int x=0; x < get<SPATIAL_X>(ea); ++x) {
                    for (int y=0; y<get<SPATIAL_Y>(ea); ++y){
                        typename EA::individual_type::ea_type::environment_type::location_ptr_type l = control_sp->ea().env().location(x,y);
                        if (l->occupied()) {
                            std::string lt = get<LAST_TASK>(*(l->inhabitant()),"");

                            if(lt == "not") {
                                df.write("1");
                            }
                            if (lt == "nand") {
                                df.write("2");
                            }
                            if (lt == "") {
                                df.write("0");
                            }
                            
            
                        } else {
                            df.write("-1");
                        }
                    
                    }
                }
                df.endl();

            }
            
            df.endl();
            
             }
    }
}
#endif
