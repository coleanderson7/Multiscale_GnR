//Models the G&R of a TEVG
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include "vessel.h"
#include "functions.h"

using std::string;
using std::vector;
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//All units meters, kg, days
//Conversions
double um_to_m = pow(10, -6);
double mm_to_m = pow(10, -3);
double kPa_to_Pa = pow(10, 3);

int main(int ac, char* av[]) {
    //Parse program options
    try {

        //Input argument variables to store to
        //Set name, number steps, and new pressure and wss
        int restart_arg;
        int step_arg;
        double step_size;
        double P_arg;
        double tauw_arg;
        int gnr_arg;
        string name_arg;
        int opt;
        int portnum;
        int num_days;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("restart,r", po::value<int>(&restart_arg)->default_value(0),
                "restarting simulation")
            ("step,s", po::value<int>(&step_arg), "number of timesteps to run simulation")
            ("time step size,d", po::value<double>(&step_size)->default_value(1.0), "size of each time step in days")
            ("name,n", po::value<string>(&name_arg)->default_value(""),
                "suffix of simulation files.")
            ("pressure,p", po::value<double>(&P_arg), "input pressure")
            ("wss,w", po::value<double>(&tauw_arg), "input wall shear stress")
            ("simulate,u", po::value<int>(&gnr_arg)->default_value(1), "execute simulation")
            ("max_days,m", po::value<int>(&num_days)->default_value(361), "maximum days to simulate")
        ;

        po::positional_options_description p;
        p.add("name", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << "Usage: options_description [options]\n";
            cout << desc;
            return 0;
        }

        std::cout << "Restarting simulation: " << restart_arg << "\n";
        std::cout << "Filename suffix: " << name_arg << "\n";

        string immune_file = "Immune_in_" + name_arg;
        string scaffold_file = "Scaffold_in_" + name_arg;

        //Initialize the reference vessel for the simulation
        string native_file = "Native_in_" + name_arg;
        vessel native_vessel;
        native_vessel.initializeNative(native_file,1,1.0);
        
        //Initialize TEVG
        vessel TEVG;
        TEVG.initializeTEVG(scaffold_file,immune_file,native_vessel,num_days,step_size);
        
        if (!vm.count("step")) {
            step_arg = int( num_days / step_size );
        }
        std::cout << "Steps to simulate: " << step_arg << "\n";
        if (vm.count("pressure")){
            TEVG.P = P_arg;
            std::cout << "Updated pressure: " << P_arg << std::endl;
        }
        //Apparent viscosity vars
        double mu = 0.0;
        if (vm.count("wss")){
            TEVG.bar_tauw = tauw_arg;
            mu = get_app_visc(&TEVG, 0);
            TEVG.Q = native_vessel.bar_tauw / (4 * mu / (3.14159265 * pow(TEVG.a_h*100, 3)));
            std::cout << "Updated WSS: " << tauw_arg << std::endl;
        }

        //Get all other input arguements and apply to TEVG
        TEVG.gnr_name = TEVG.gnr_name + "_" + name_arg;
        TEVG.exp_name = TEVG.exp_name + "_" + name_arg;
        TEVG.file_name = TEVG.file_name + "_" + name_arg;

        //------------------------------------------------------------------------

        //PD Tests
        double P_low = 0.0 * 133.33;
        double P_high = 10.0 * 133.33;
        double lambda_z_test = 1.0;
        //int pd_test_check = 0;
        //int pd_test_count = 0;
        //vector<int> pd_test_ind = { 1, 11, 91, 151 };
        //Pressure ramping
        TEVG.P_prev = TEVG.P;

        if(!restart_arg)
        {
            std::cout << "Initializing new simulation." << std::endl;
            //Setup file I/O for G&R output
            TEVG.GnR_out.open(TEVG.gnr_name);
            TEVG.Exp_out.open(TEVG.exp_name);

            //Write initial state to file
            int sn = 0;
            TEVG.printTEVGOutputs();
            //Simulate number of timesteps
            //TEVG.nts = step_arg;

            if (!(vm.count("wss"))){
                TEVG.wss_calc_flag = 1;
            }

            //Run the G&R time stepping
            for (int sn = 1; sn < std::min(step_arg,TEVG.nts); sn++) {

                if(gnr_arg){
                    TEVG.s = TEVG.dt * sn;
                    TEVG.sn = sn;
                    update_time_step(TEVG);
                    printf("%s \n", "---------------------------");
                    fflush(stdout);
                    //Store axial stretch history
                    TEVG.lambda_z_tau[sn] = TEVG.lambda_z_curr;
                    TEVG.P_prev = TEVG.P;
                }
                else {
                    //TEVG.s = 0;
                    //TEVG.sn = 0;
                    //Print current state
                    printf("%s %f %s %f %s %f %s %f %s %f %s %f\n", "Time:", TEVG.s, "a: ", TEVG.a[TEVG.sn], "h:", TEVG.h[TEVG.sn],
                   "Cbar_r:", TEVG.Cbar[0], "Cbar_t:", TEVG.Cbar[1], "Cbar_z:", TEVG.Cbar[2]);
                    printf("%s \n", "xxxxxxxxxxxxxxxxxxxxxxxxxxx");
                    fflush(stdout);    
                }
                
                //system("read");

                //To check to run a pressure-diameter test
                //if (pd_test_count < pd_test_ind.size() && pd_test_check > 0) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
                TEVG.printTEVGOutputs();
            }

            //Print vessel to file
            TEVG.save();

        }
        else
        {
            std::cout << "Continuing simulation from file." << std::endl;
            //Setup file I/O for G&R output
            //std::ofstream GnR_out;
            TEVG.GnR_out.open(TEVG.gnr_name, std::ofstream::out | std::ofstream::app);
            TEVG.Exp_out.open(TEVG.exp_name, std::ofstream::out | std::ofstream::app);

            //Read vessel from file
            TEVG.load();
            //Load pressure and WSS
            if (vm.count("pressure")){
                TEVG.P = P_arg;
            }
            if (vm.count("wss")){
                TEVG.bar_tauw = tauw_arg;
                mu = get_app_visc(&TEVG, TEVG.sn);
                TEVG.Q = TEVG.bar_tauw / (4 * mu / (3.14159265 * pow(TEVG.a[TEVG.sn]*100, 3)));
            }

            if (TEVG.sn == 0){
                //Write initial state to file
                int sn = 0;
                TEVG.printTEVGOutputs();
            }

            if (!(vm.count("wss"))){
                TEVG.wss_calc_flag = 1;
            }

            int csn = TEVG.sn+1;

            //Run the G&R time stepping
            for (int sn = csn; sn < std::min(csn+step_arg,TEVG.nts); sn++) {

                if(gnr_arg){
                    TEVG.s = TEVG.dt * sn;
                    TEVG.sn = sn;
                    update_time_step(TEVG);
                    printf("%s \n", "---------------------------");
                    fflush(stdout);
                    //Store axial stretch history
                    TEVG.lambda_z_tau[sn] = TEVG.lambda_z_curr;
                    TEVG.P_prev = TEVG.P;
                }
                else {
                    //Print current state
                    printf("%s %f %s %f %s %f %s %f %s %f %s %f\n", "Time:", TEVG.s, "a: ", TEVG.a[TEVG.sn], "h:", TEVG.h[TEVG.sn],
                   "Cbar_r:", TEVG.Cbar[0], "Cbar_t:", TEVG.Cbar[1], "Cbar_z:", TEVG.Cbar[2]);
                    printf("%s \n", "xxxxxxxxxxxxxxxxxxxxxxxxxxx");
                    fflush(stdout);    
                }
                
                //To check to run a pressure-diameter test
                //if (pd_test_count < pd_test_ind.size()) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
                TEVG.printTEVGOutputs();

            }

            //Print vessel to file
            TEVG.save();

        }

        TEVG.GnR_out.close();
        TEVG.Exp_out.close();
    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;

}
