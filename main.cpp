/* author: jakob fischer (jakob@automorph.info)
 * date: 22nd March 2013
 * description: 
 * Tool for solving a ODE for a reaction system given as an jrnf-file. The concentration and time
 * are read and written from / to a comma seperated file. Every row / line represents one time step.
 * the first column contains the time of this step, the second the msd per specie number and time
 * step from the previous step...
 * 
 * TODO Implement and comment
 * 
 * TODO Extend so concentration can be taken from external file (allowing for example periodic 
 *      boundary conditions)
 *
 * TODO Increase precission by only calculating effective flow and not forward- and backward-flow 
 *      independently
 * 
 * TODO Use implicit solver from odeint...
 *
 * TODO In connection with function to calculate effective rates, extend file format to contain
 *      concentrations (N values) and effective rates (M values)
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include <algorithm>
#include "tools/cl_para.h"
#include "net_tools/reaction_network.h"
#include "net_tools/reaction_network_fileop.h"
#include "net_tools/network_tools.h"

using namespace std;
using namespace boost::numeric::odeint;


#ifndef GIT_VERSION
#define GIT_VERSION "no version"
#endif



//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
typedef boost::numeric::ublas::matrix< int > matrix_type_int;




class reaction_network_system {
    double beta;
    double initial_t, last_write;
    std::vector<species> sp;
    std::vector<reaction> re;
    vector_type initial_con;
    std::string fn_concentration;

    matrix_type_int N, N_in, N_out;
    vector_type Ea, mu0, e_bs_in, e_bs_out;

public:
    struct stiff_system {
        const reaction_network_system& rns;

        stiff_system(const reaction_network_system& rns_) : rns(rns_)  {}

        /*
         * Prototype for the rhs of ODE
         * TODO: implement + test + speed up
         */

        void operator()( const vector_type &x , vector_type &dxdt , double /* t */ ) {
            vector_type a1(x.size()), a2(x.size()), a3(x.size());
            
            for(size_t k=0; k<x.size(); ++k) {
                a1(k) = exp(-rns.beta*rns.Ea(k));

                a2(k) = 1;
                for(size_t j=0; j<rns.N.size2(); ++j) 
                    if(rns.N_in(j,k) != 0) 
                        a2(k) *= pow(x(j), rns.N_in(j,k));
                    
                a3(k) = 1; 
                for(size_t j=0; j<rns.N.size2(); ++j) 
                    if(rns.N_out(j,k) != 0) 
                        a3(k) *= pow(x(j), rns.N_out(j,k));
            }          

             dxdt = prod(rns.N, element_prod(a1, element_prod(a2, rns.e_bs_in) - element_prod(a3, rns.e_bs_out)));

             for(size_t k=0; k<x.size(); ++k)
                 if(rns.sp[k].is_constant())
                     dxdt(k) = 0;
        }
    };


    struct stiff_system_jacobi {
        const reaction_network_system& rns;

        stiff_system_jacobi(const reaction_network_system& rns_) : rns(rns_)  {}

        void operator()( const vector_type & /* x */ , matrix_type &J , const double & /* t */ , vector_type &dfdt ) {
            J( 0 , 0 ) = -101.0;
            J( 0 , 1 ) = -100.0;
            J( 1 , 0 ) = 1.0;
            J( 1 , 1 ) = 0.0;
            dfdt[0] = 0.0;
            dfdt[1] = 0.0;
        }
    };






    reaction_network_system(const std::string& fn_network, const std::string& fn_concentration_) 
        : fn_concentration(fn_concentration_), beta(1)  {
        
        read_jrnf_reaction_n(fn_network, sp, re);
    
	
	// Remove 1-0 and 0-1 reactions. Such reactions are there to ballance flow through the boundary conditions 
        std::remove_if (re.begin(), re.end(), 
                        [] (reaction& re) -> bool { 
                            return re.get_no_educt() == 1 && re.get_no_product() == 0 || 
                                   re.get_no_educt() == 0 && re.get_no_product() == 1; });
		
	



        //matrix_type_int N, N_in, N_out;
        // TODO: fill these matrices

        


        //vector_type Ea, mu0, e_bs_in, e_bs_out;
        // TODO: fill these vectors



        
        // Read last line of concentration file and write initial concentration to initial_con
        // and initial time to initial_t. After done open the same file for appending...
        initial_con = boost::numeric::ublas::zero_vector<double>(sp.size());

        std::ifstream  data(fn_concentration.c_str());

        if(!data.good()) {
           std::cout << "Could not open concentration file: " << fn_concentration << std::endl;
           return;
        }

        std::string line;
        std::getline(data,line);       // dont want the header
        double last_msd=0;

        // the real 
        while(!std::getline(data,line).eof()) {
            std::stringstream ls(line);
            std::string cell;

            size_t cnt=0;
            while(std::getline(ls,cell,',')) {
                std::stringstream in(cell);       
    
    
                if(cnt == 0) 
                    in >> initial_t ;
                else if(cnt == 1)
                    in >> last_msd;                
                else if(cnt <= sp.size()+1)
                    in >> initial_con(cnt-2);
                else 
                    std::cout << "Error at reading csv / concentration file!" << std::endl;

                ++cnt;
            }
        }

        std::cout << "Simulating file: " << fn_concentration << std::endl;
        std::cout << "Loaded concentration file with starting time " << initial_t << std::endl;
                  
        data.close();        
         
    }



    void run(double Tmax=25000, double deltaT=0.1, double wint=1000, 
             bool write_rates=false, bool solve_implicit=true) {
                
        fstream out(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);
        out.precision(25);

        size_t t0 = time(NULL);

        // Init solver and start
    
        //init_state(x);
        //integrate( next_step , x , initial_t , Tmax , deltaT , write_state );

        vector_type x(initial_con);


        auto write_state = [this]( const vector_type &vec , const double t ) {
            cout << "initial_t=" << initial_t << endl;
            cout << t << "," << vec[0] << "," << vec[1] << endl;
        };


        size_t step_no = 0;

        if(solve_implicit) {
            step_no = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
                                       make_pair( stiff_system(*this) , stiff_system_jacobi(*this) ) ,
                                       x , initial_t , Tmax , deltaT ,
                                       write_state);
        } else {
            step_no = integrate_const( make_dense_output< runge_kutta_dopri5< vector_type > >( 1.0e-6 , 1.0e-6 ) ,
                                       stiff_system(*this) , x , initial_t , Tmax , deltaT ,
                                       write_state);
        }


 
        // 
        write_state(x, Tmax);

        size_t t1 = time(NULL);
        std::cout << "Run took " << t1-t0 << " seconds!" << std::endl;

        out.close();
   }



};













int main(int argc, const char *argv []){
    srand(time(0));
    std::cout << "odeint_rnet version 0x00x04 (commit:" << GIT_VERSION << ")" << std::endl;
    
    cl_para cl(argc, argv);
    
    // Simulates / integrates a reaction network while holding the boundary point 
    // species (constant=true) constant
    // under development...

    if(cl.have_param("simsim_dev")) {
        std::string fn_network, fn_concentration;      

        if(cl.have_param("net")) 
            fn_network=cl.get_param("net");
        else {
            cout << "You have to give the name of reaction network by 'net'!" << endl;  
            return 1;
        }
    
        if(cl.have_param("con")) 
            fn_concentration=cl.get_param("con");
        else {
            cout << "You have to give the name of the concentration file by 'con'!" << endl;   
            return 1;
        }

        double deltaT = cl.have_param("deltaT") ? cl.get_param_d("deltaT") : 0.1;   
        double Tmax = cl.have_param("Tmax") ? cl.get_param_d("Tmax") : 25000;
        double wint = cl.have_param("wint") ? cl.get_param_i("wint") : 1000;
        bool write_rates=cl.have_param("write_rates");
        bool solve_implicit=cl.have_param("solve_implicit");            

            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax << std::endl;

        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration); 
        rns.run(Tmax, deltaT, wint, write_rates, solve_implicit);
    }  
    
    
    if(cl.have_param("help") || cl.have_param("info")) {
        cout << "          odeint_rnet" << endl;  
        cout << "          ===========" << endl;
        cout << " call with parameter 'info' or 'help' for showing this screen" << endl;
        cout << endl;
        cout << "-'simsim_dev' load reaction network 'net' and simulate file 'con'!. Parameters are 'deltaT', and 'Tmax'." << endl;
	cout << endl;
    } 
}
