/* author: jakob fischer (jakob@automorph.info)
 * description: 
 * Tool for solving a ODE for a reaction system given as an jrnf-file. The concentration and time
 * are read and written from / to a comma seperated file. Every row / line represents one time step.
 * the first column contains the time of this step, the second the mean of the quadratic change 
 * (<f(x)> ;dx/dt = f(x)). Then follow the concentrations of all species. If the program is called with
 * the option 'write_rates' this is followed by the effective rates of all reactions.
 *
 * 
 * TODO Extend so boundary concentration can be taken from external file (allowing for example periodic 
 *      boundary conditions)
 *
 * TODO Increase precission by only calculating effective flow and not forward- and backward-flow 
 *      independently...  (Not sure if this is possible)
 *
 *      This version of the program calculates forward and backward reaction rates from formation enthalpies
 *      and activation energies given in the network description. (\beta = 1/(k_b T) = 1)
 * TODO Make this optional + Allow choosing \beta at startup
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

const string odeint_rnet_version_string="0x00x04";



// 
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
typedef boost::numeric::ublas::matrix< int > matrix_type_int;




class reaction_network_system {
    double beta;
    double initial_t, last_write, last_msd;
    std::vector<species> sp;
    std::vector<reaction> re;
    vector_type initial_con, initial_rate, last_flow;
    std::string fn_concentration;
    
    bool write_rate;
    matrix_type_int N, N_in, N_out;
    vector_type Ea, mu0, e_bs_in, e_bs_out, e_m_bEa;

public:
    struct stiff_system {
        reaction_network_system& rns;

        stiff_system(reaction_network_system& rns_) : rns(rns_)  {}

        /*
         * Prototype for the rhs of ODE
         * TODO: test + speed up
         */

        void operator()( const vector_type &x , vector_type &dxdt , double /* t */ ) {
            vector_type a2(rns.re.size()), a3(rns.re.size());
            for(size_t j=0; j<rns.N.size2(); ++j) { 
                a2(j) = 1;
                a3(j) = 1; 

                for(size_t k=0; k<x.size(); ++k) {             
                    if(rns.N_in(k,j) != 0) 
                        a2(j) *= pow(x(k), rns.N_in(k,j));

                    if(rns.N_out(k,j) != 0) 
                        a3(j) *= pow(x(k), rns.N_out(k,j));
                }
            }         



            vector_type k_f=element_prod(rns.e_m_bEa, rns.e_bs_in);
            vector_type k_b=element_prod(rns.e_m_bEa, rns.e_bs_out);

            rns.last_flow = element_prod(rns.e_m_bEa, element_prod(a2, rns.e_bs_in) - element_prod(a3, rns.e_bs_out));
            dxdt = prod(rns.N, rns.last_flow);

            for(size_t k=0; k<x.size(); ++k)
                if(rns.sp[k].is_constant())
                    dxdt(k) = 0;

            // square all elements of dxdt and calculate the mean  (to analyze convergence later)
            rns.last_msd = 0;
            for(size_t i=0; i<dxdt.size(); ++i)
                rns.last_msd += dxdt(i)*dxdt(i);
            rns.last_msd /= dxdt.size();
        }
    };




    struct stiff_system_jacobi {
        const reaction_network_system& rns;

        stiff_system_jacobi(const reaction_network_system& rns_) : rns(rns_)  {}


        /*
         * Prototype for the Jacobi-Matrix of the ODE
         * TODO: test + speed up 
         */

        void operator()( const vector_type & x  , matrix_type &J , const double & /* t */ , vector_type &dfdt ) {
            matrix_type m(x.size(),x.size());
           
            for(size_t i=0; i<x.size(); ++i) {
                for(size_t l=0; l<x.size(); ++l) {
                    double s(0);
                    
                    for(size_t k=0; k<rns.re.size(); ++k) {    
                        double p1(1), p2(1);    // p1 = \prod_{j \neq k} {x_j}^N_{jk}^\mathrm{in}
                                                // p2 = \prod_{j \neq k} {x_j}^N_{jk}^\mathrm{out}

                        for(size_t j=0; j<rns.sp.size(); ++j) {
                            if(j != l && rns.N_in(j,k) != 0)
                                p1 *= x(j);
    
                            if(j != l && rns.N_out(j,k) != 0)
                                p2 *= x(j);
                        } 
                      
                        s += rns.N(i,k)*rns.e_m_bEa(k)*(rns.e_bs_in(k)*rns.N_in(l,k)*pow(x(l), rns.N_in(l,k)-1)*p1 - 
                                                        rns.e_bs_out(k)*rns.N_out(l,k)*pow(x(l),rns.N_out(l,k)-1)*p2);              
                    }                    

                    J(i,l) = s;
                    if(rns.sp[i].is_constant())
                        J(i,l) = 0;
                }
            }                               


            dfdt = boost::numeric::ublas::zero_vector<double>(x.size());
        }
    };






    reaction_network_system(const std::string& fn_network, const std::string& fn_concentration_, bool write_rate_=false) 
        : fn_concentration(fn_concentration_), beta(1), last_msd(0), write_rate(write_rate_)  {
        
        read_jrnf_reaction_n(fn_network, sp, re);
    
	
	// Remove 1-0 and 0-1 reactions. Such reactions are there to ballance flow through the boundary conditions 
        std::remove_if (re.begin(), re.end(), 
                        [] (reaction& re) -> bool { 
                            return re.get_no_educt() == 1 && re.get_no_product() == 0 || 
                                   re.get_no_educt() == 0 && re.get_no_product() == 1; });
		
        // Calculate stoichiometric matrices, activation energies and their logarithms...
        N_in = N_out = boost::numeric::ublas::zero_matrix<double>(sp.size(), re.size());
        Ea =  e_bs_in = e_bs_out = e_m_bEa = boost::numeric::ublas::zero_vector<double>(re.size());
        mu0 = boost::numeric::ublas::zero_vector<double>(sp.size());


        for(size_t i=0; i<sp.size(); ++i)
            mu0(i) = sp[i].get_energy();


        for(size_t i=0; i<re.size(); ++i) {
            double mu0_educts(0), mu0_products(0);

            for(size_t j=0; j<re[i].get_no_educt_s(); ++j) {
                size_t id(re[i].get_educt_id(j)),  mul(re[i].get_educt_mul(j));
                N_in(id,i) += mul;
                mu0_educts += sp[id].get_energy()*mul;
            }

            for(size_t j=0; j<re[i].get_no_product_s(); ++j) {
                size_t id(re[i].get_product_id(j)),  mul(re[i].get_product_mul(j));
                N_out(id,i) += mul;
                mu0_products += sp[id].get_energy()*mul;
            }

            Ea(i) = re[i].get_activation()+max(mu0_educts,mu0_products);
            e_bs_in(i) = exp(beta*mu0_educts);
            e_bs_out(i) = exp(beta*mu0_products);
            e_m_bEa(i) = exp(-beta*Ea(i));
        }
        
        N = N_out-N_in;

        
        // Read last line of concentration file and write initial concentration to initial_con
        // and initial time to initial_t. After done open the same file for appending...
        initial_con = boost::numeric::ublas::zero_vector<double>(sp.size());
        initial_rate = boost::numeric::ublas::zero_vector<double>(re.size());
        last_flow = boost::numeric::ublas::zero_vector<double>(re.size());

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
                else if(cnt <= sp.size()+re.size()+1  && write_rate)
                    in >> initial_rate(cnt-sp.size()-2);
                else
                    std::cout << "Error at reading csv / concentration file!" << std::endl;

                ++cnt;
            }
        }

        std::cout << "Simulating file: " << fn_concentration << std::endl;
        std::cout << "Loaded concentration file with starting time " << initial_t << std::endl;
                  
        data.close();        
         
    }


    void print_rhs(size_t t=0) {
        vector_type x(initial_con);
        vector_type dxdt=boost::numeric::ublas::zero_vector<double>(sp.size());
 
        // calculate rhs
        stiff_system(*this)(x, dxdt, t);

        cout << "righthandside:" << endl;
        for(size_t i=0; i<dxdt.size(); ++i)
            cout << "/  " << dxdt(i) << "  /";
        cout << endl;
    }


    void run(double Tmax=25000, double deltaT=0.1, bool write_rates=false, bool solve_implicit=true) {
                
        fstream out(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);
        out.precision(25);

        size_t t0 = time(NULL);

        // Init solver and start
        vector_type x(initial_con);
        last_flow = initial_rate;
         
        auto write_state = [this, &out, write_rates]( const vector_type &vec , const double t ) {
            out << t << "," << last_msd;

            for(size_t l=0; l<vec.size(); ++l)
                out << "," << vec(l);

            if(write_rates)
                for(size_t l=0; l<last_flow.size(); ++l)
                    out << "," << last_flow(l);
 
            out << std::endl;
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
        std::cout << "Run took " << t1-t0 << " seconds and " << step_no << " steps!" << std::endl;

        out.close();
   }
};



int main(int argc, const char *argv []){
    srand(time(0));
    std::cout << "odeint_rnet version " << odeint_rnet_version_string << " (commit:" 
              << GIT_VERSION << ")" << std::endl;
    
    cl_para cl(argc, argv);
    

    // Prints the right hand side of the ODE (f(x)) for diagnostic purposes.

    if(cl.have_param("print_rhs")) {
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

        // load network and concentration file
        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration); 
        rns.print_rhs(0);  // time is 0 (doesn't matter)
    }  


    // Simulates / integrates a reaction network while holding the boundary point 
    // species (constant=true) constant
    
    if(cl.have_param("simulate")) {
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

        double deltaT = cl.have_param("deltaT") ? cl.get_param_d("deltaT") : 10;   
        double Tmax = cl.have_param("Tmax") ? cl.get_param_d("Tmax") : 25000;
        bool write_rates=cl.have_param("write_rates");
        bool solve_implicit=cl.have_param("solve_implicit");            

            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax << std::endl;
        if(write_rates) 
            cout << "Output contains reactions' effective rates." << endl;

        if(solve_implicit)
            cout << "Stiff solver is used!" << endl;

        // Load reaction network and concentration file
        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration, write_rates); 

        // Simulate ODE (and write results to file)
        rns.run(Tmax, deltaT, write_rates, solve_implicit);
    }  
    
    
    if(cl.have_param("help") || cl.have_param("info")) {
        cout << "          odeint_rnet" << endl;  
        cout << "          ===========" << endl;
        cout << " call with parameter 'info' or 'help' for showing this screen" << endl;
        cout << "current version is " << odeint_rnet_version_string << " (commit:" << GIT_VERSION << ")" << endl;
        cout << endl;
        cout << "-'print_rhs': load reaction network 'net', (last) concentration in 'con'! and" << endl;
        cout << "prints  right hand side of ode for diagnostic purposes." << endl;
	cout << endl;
        cout << "-'simulate': load reaction network 'net' and simulate file 'con'!. Parameters are:" << endl;
        cout << "   x'deltaT': Interval in which concentrations are saved to file 'con'" << endl;
        cout << "   x'Tmax': Time up to which the system is simulated" << endl;
        cout << "   x'write_rates': 'con' file contains effective reaction rates" << endl;
        cout << "   x'solve_implicit': use stiff solver" << endl;
	cout << endl;
    } 
}
