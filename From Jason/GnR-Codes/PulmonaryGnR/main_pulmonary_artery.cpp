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

int main( int ac, char* av[] ) {

    try{

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
        double gamma_p;
        double gamma_q;
        double gamma_act;
        int gnr_equil_arg;
        int gnr_iter_flag;
        int gnr_out_flag;
        int mech_infl_flag;
        int mech_exp_flag;

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
            ("simulate_equil", po::value<int>(&gnr_equil_arg)->default_value(1), "execute equilibrated simulation")
            ("max_days,m", po::value<int>(&num_days)->default_value(361), "maximum days to simulate")
            ("gamma_p", po::value<double>(&gamma_p)->default_value(0.0), "fold change in pressure from baseline")
            ("gamma_q", po::value<double>(&gamma_q)->default_value(0.0), "fold change in flow from baseline")
            ("gamma_act", po::value<double>(&gamma_act)->default_value(0.0), "fold change in active stress from baseline")
            ("gnr_iter_flag", po::value<int>(&gnr_iter_flag)->default_value(0), "flag for iterating at fixed time with new inputs")
            ("gnr_out_flag", po::value<int>(&gnr_out_flag)->default_value(1), "flag for outputting to GnR_out")
            ("mech_infl_flag", po::value<int>(&mech_infl_flag)->default_value(0), "flag for controlling mech-mediated infl.")
            ("mech_exp_flag", po::value<int>(&mech_exp_flag)->default_value(0), "flag for controlling mech exp")
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

        //Initialize the reference vessel for the simulation
        string native_file = "Native_in_" + name_arg;

        vessel native_vessel;
        native_vessel.initializeNative(native_file,num_days,step_size);

        if (!vm.count("step")) {
            step_arg = int( num_days / step_size );
        }
        std::cout << "Steps to simulate: " << step_arg << "\n";
        if (vm.count("pressure")){
            native_vessel.P = P_arg;
            std::cout << "Updated pressure: " << P_arg << std::endl;
        }
        //Apparent viscosity vars
        double mu = 0.0;
        if (vm.count("wss")){
            native_vessel.bar_tauw = tauw_arg;
            mu = get_app_visc(&native_vessel, 0);
            native_vessel.Q = native_vessel.bar_tauw / (4 * mu / (3.14159265 * pow(native_vessel.a_h*100, 3)));
            std::cout << "Updated WSS: " << tauw_arg << std::endl;
        }

        //Set mech. induced infl
        if (vm.count("mech_infl_flag")){
            native_vessel.mech_infl_flag = mech_infl_flag;
            std::cout << "Setting mech-induced infl: " << mech_infl_flag << std::endl;
        }

        //Set mech. experiment testing flag
        if (vm.count("mech_exp_flag")){
            native_vessel.mech_exp_flag = mech_exp_flag;
            std::cout << "Setting mech exp flag: " << mech_exp_flag << std::endl;
        }

        //Setup all other output files
        native_vessel.gnr_name = native_vessel.gnr_name + "_" + name_arg;
        native_vessel.equil_gnr_name = native_vessel.equil_gnr_name + "_" + name_arg;
        native_vessel.exp_name = native_vessel.exp_name + "_" + name_arg;
        native_vessel.file_name = native_vessel.file_name + "_" + name_arg;

        //------------------------------------------------------------------------

        //For elastin degradation 
        double Q_p1, s;
        double s_gnd_off = 45, epsilonR_gnd_min = 0.75, k_gnd_h = 0.1;
        //for (int sn = 1; sn < native_vessel.nts; sn++) {
            //Calculate elastin degradation
            //s = sn * native_vessel.dt;
            //Q_p1 = 1.0; // (1 + exp(-k_p1 * zeta_p1)) / (1 + exp(k_p1 * gamma_p_d1 * (s - zeta_p1 * gamma_p_d2)));

            //Q_p1 = (s > s_gnd_off) * ((1 - epsilonR_gnd_min) * exp(-k_gnd_h * (s - s_gnd_off)) + epsilonR_gnd_min)
                //+ (s <= s_gnd_off) * 1;

            //native_vessel.epsilonR_alpha[0 * native_vessel.nts + sn] = Q_p1 * native_vessel.epsilonR_alpha[0 * native_vessel.nts + 0];

            //native_vessel.rhoR_alpha[0 * native_vessel.nts + sn] = native_vessel.epsilonR_alpha[0 * native_vessel.nts + sn] * 
                                                                //native_vessel.rho_hat_alpha_h[0];

        //}

        //PD Tests
        //double P_low = 0.0 * 133.33;
        //double P_high = 10.0 * 133.33;
        //double lambda_z_test = 1.0;
        //int pd_test_check = 0.0;
        //int pd_test_count = 0;
        //vector<int> pd_test_ind = { 1, 11, 91, 151 }; //, 91, 151

        //double gamma_P = 0.5;
        double gamma_Q = 0.0;
        double gamma_lambda_z = 0.0;
        double perturb_time = 30.0;

        if(!restart_arg)
        {
            std::cout << "Initializing new simulation..." << std::endl;
            //Setup file I/O for G&R output
            native_vessel.GnR_out.open(native_vessel.gnr_name);
            native_vessel.Equil_GnR_out.open(native_vessel.equil_gnr_name);
            native_vessel.Exp_out.open(native_vessel.exp_name);

            //Write initial state to file
            int sn = 0;
            native_vessel.printNativeOutputs();
            //Simulate number of timesteps
            //native_vessel.nts = step_arg;

            //Set flag for WSS calculation
            if (!(vm.count("wss"))){
                native_vessel.wss_calc_flag = 1;
            }
            
            //Calculate pressure and flow from pertubations
            if (vm.count("gamma_p")){
                native_vessel.P = (1 + gamma_p) * native_vessel.P_h;
            }
            if (vm.count("gamma_q")){
                native_vessel.Q = (1 + gamma_q) * native_vessel.Q_h;
            }
            if (vm.count("gamma_act")){
                native_vessel.T_act = (1 + gamma_act) * native_vessel.T_act_h;
            }

            //Run the G&R time stepping
            for (int sn = 1; sn < std::min(step_arg,native_vessel.nts); sn++) {
                
                s = native_vessel.dt * sn;

                if(gnr_arg){

                    // if ( s > perturb_time2){
                    //     native_vessel.P = native_vessel.P_h * (1 + gamma_P * (1 - exp(-(s  - perturb_time) / 10))); // + gamma_P2 * (1 - exp(-(s  - perturb_time2) / 10)));
                    //     native_vessel.Q = native_vessel.Q_h * (1 + gamma_Q * (1 - exp(-(s  - perturb_time) / 10))); // + gamma_Q2 * (1 - exp(-(s  - perturb_time2) / 10)));
                    //     native_vessel.lambda_z_curr = native_vessel.lambda_z_h * (1 + gamma_lambda_z * (1 - exp(-(s - perturb_time) / 10))); // + gamma_lambda_z2 * (1 - exp(-(s  - perturb_time2) / 10)));
                    // }

                    native_vessel.s = native_vessel.dt * sn;
                    native_vessel.sn = sn;

                    update_time_step(native_vessel);
                    printf("%s \n", "---------------------------");
                    fflush(stdout);

                    //Store axial stretch history
                    native_vessel.lambda_z_tau[sn] = native_vessel.lambda_z_curr;

                    //Update previous loading state
                    native_vessel.P_prev = native_vessel.P;
                    native_vessel.T_act_prev = native_vessel.T_act;

                }
                // else{
                //     //Print current state
                //     printf("%s %f %s %e %s %e %s %f %s %f %s %f\n", "Time:", native_vessel.s, "a: ", native_vessel.a[native_vessel.sn], "h:", native_vessel.h[native_vessel.sn],
                //    "Cbar_r:", native_vessel.Cbar[0], "Cbar_t:", native_vessel.Cbar[1], "Cbar_z:", native_vessel.Cbar[2]);
                //     printf("%s \n", "xxxxxxxxxxxxxxxxxxxxxxxxxxx");
                //     fflush(stdout); 
                // }

                //To check to run a pressure-diameter test
                //if (pd_test_count < pd_test_ind.size()) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
               native_vessel.printNativeOutputs();
            }

            //Long-term equilibrated solution
            if(gnr_equil_arg){

                sn = native_vessel.nts - 1;
                s = native_vessel.dt * sn;
                
                native_vessel.s = s;
                native_vessel.sn = sn;
            
                native_vessel.P = (1.0 + gamma_p) * native_vessel.P_h;
                native_vessel.Q = (1.0 + gamma_q) * native_vessel.Q_h;
            
                int equil_solve = find_equil_geom(&native_vessel);
                native_vessel.printNativeEquilibratedOutputs();
                //Print equilibrated solution
                printf("%s %f %s %e %s %e %s %f\n", "Time:", s, "a_e: ", native_vessel.a_e, "h_e:", native_vessel.h_e, "mb_equil:", native_vessel.mb_equil_e);
                fflush(stdout);

            }

            //Print vessel to file
            native_vessel.save();

        }
        else
        {
            std::cout << "Continuing simulation from file..." << std::endl;
            //Setup file I/O for G&R output
            native_vessel.GnR_out.open(native_vessel.gnr_name, std::ofstream::out | std::ofstream::app);
            native_vessel.Exp_out.open(native_vessel.exp_name, std::ofstream::out | std::ofstream::app);

            //Read vessel from file
            native_vessel.load();
            //Load pressure and WSS
            if (vm.count("pressure")){
                native_vessel.P = P_arg;
            }
            if (vm.count("wss")){
                native_vessel.bar_tauw = tauw_arg;
                mu = get_app_visc(&native_vessel, native_vessel.sn);
                native_vessel.Q = native_vessel.bar_tauw / (4 * mu / (3.14159265 * pow(native_vessel.a[native_vessel.sn]*100, 3)));
            }

            //Calculate pressure and flow from pertubations
            if (vm.count("gamma_p")){
                native_vessel.P = (1 + gamma_p) * native_vessel.P_h;
            }
            if (vm.count("gamma_q")){
                native_vessel.Q = (1 + gamma_q) * native_vessel.Q_h;
            }
            if (vm.count("gamma_act")){
                native_vessel.T_act = (1 + gamma_act) * native_vessel.T_act_h;
            }

            if (native_vessel.sn == 0){
                //Write initial state to file
                int sn = 0;
                native_vessel.printNativeOutputs();
            }

            if (!(vm.count("wss"))){
                native_vessel.wss_calc_flag = 1;
            }

            int csn = native_vessel.sn+1;
            if (gnr_iter_flag){
                csn = native_vessel.sn;
            }

            //Run the G&R time stepping
            for (int sn = csn; sn < std::min(csn+step_arg,native_vessel.nts); sn++) {
                
                if(gnr_arg){
                    native_vessel.s = native_vessel.dt * sn;
                    native_vessel.sn = sn;
                    update_time_step(native_vessel);
                    printf("%s \n", "---------------------------");
                    fflush(stdout);
                    //Store axial stretch history
                    native_vessel.lambda_z_tau[sn] = native_vessel.lambda_z_curr;
                    native_vessel.P_prev = native_vessel.P;
                    native_vessel.T_act_prev = native_vessel.T_act;
                }
                else{

                }
                
                //if (pd_test_count < pd_test_ind.size()) {
                    //if (sn == pd_test_ind[pd_test_count]) {
                        //pd_test_check = run_pd_test(TEVG, P_low, P_high, lambda_z_test);
                        //pd_test_count++;
                    //}
                //}

                //Write full model outputs
                //if (gnr_iter_flag){
                    //native_vessel.printExpOutputs();
                //}

                if (gnr_out_flag) {
                    native_vessel.printNativeOutputs();
                }
                else{
                    native_vessel.printExpOutputs();
                }

            }

            //Print vessel to file
            native_vessel.save();

            //Long-term equilibrated solution
		    //int sn = native_vessel.nts - 1;
            //double s = native_vessel.dt * sn;
		    //native_vessel.s = s;
		    //native_vessel.sn = sn;
				
            //int equil_solve = find_equil_geom(&native_vessel);
            //native_vessel.printNativeEquilibratedOutputs();
            //Print equilibrated solution
            //printf("%s %f %s %f %s %f %s %f\n", "Time:", s, "a_e: ", native_vessel.a_e, "h_e:", native_vessel.h_e, "mb_equil:", native_vessel.mb_equil_e);
            //fflush(stdout);

        }

        native_vessel.GnR_out.close();
        native_vessel.Equil_GnR_out.close();
        native_vessel.Exp_out.close();

    }
    catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }

    return 0;

}
