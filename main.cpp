/* author: jakob fischer (jakob@automorph.info)
 * description: 
 * Tool for solving a ODE for a reaction system given as an jrnf-file. The concentration and time
 * are read and written from / to a comma seperated file. Every row / line represents one time step.
 * the first column contains the time of this step, the second the mean of the quadratic change 
 * (<f(x)> ;dx/dt = f(x)). Then follow the concentrations of all species. If the program is called with
 * the option 'write_rates' this is followed by the effective rates of all reactions.
 *
 * One version of the program (option "simulate" calculates effective reaction rates directly
 * from formation enthalpies and activation energies given in the network description. 
 * (\beta = 1/(k_b T) = 1). This allows higher precision and usage of an explicit solver which
 * is over all faster for small networks (< 50 species). 
 *
 * TODO There seems to be a problem with exceptionally short steps (distances in the
 *      <times> vector given to integrate_times) as well as with very wide steps. This
 *      occurs when the "write_log" option lead to output step sizes below 1e-7 or 
 *      above 1e7. The solution for now is to first use linear stepping, then 
 *      "write_log_abs" and at larger timescale linear stepping again. 
 *      Giving a minimal and maximal step size as parameter might be a good 
 *      improvement for future versions! 
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

// version string increment if mayor changes happen
const string odeint_rnet_version_string="0x00x04";

// various ublas types needed for implicit solver
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
typedef boost::numeric::ublas::matrix< int > matrix_type_int;


/*
 *
 *
 */

class reaction_network_system {
    double beta;                // beta = 1/(k_b T)
    std::vector<species> sp;    // The reaction
    std::vector<reaction> re;   // network
    // initial time, concentration and rates loaded with concentration file
    double initial_t;    
    vector_type initial_con, initial_rate;
    std::string fn_concentration;      // filename
    
    bool write_rate;                   // do write rate?
    matrix_type_int N, N_in, N_out;    // stoichiometric matrices
    // List that contains for each column in N_in / N_out (= for each reaction)
    // a vector of those rows (species) that are nonzero.
    std::vector< std::vector<size_t> > N_in_list, N_out_list;
    // 
    vector_type Ea, mu0, e_bs_in, e_bs_out, e_m_bEa;

public:
    struct fast_system {
        reaction_network_system& rns;

        fast_system(reaction_network_system& rns_) : rns(rns_)  {}

        /*
         * Prototype for the rhs of ODE
         * 
         */

        void operator()( const vector_type &x , vector_type &dxdt , double /* t */ ) {

            for(size_t i=0; i<rns.sp.size(); ++i) 
                dxdt[i]=0;
    
            for(size_t i=0; i<rns.re.size(); ++i) {
                double rate_f=rns.re[i].get_k();
	        double rate_b=rns.re[i].get_k_b();
	 
                for(size_t j=0; j<rns.re[i].get_no_educt_s(); ++j)
                    rate_f *= x[rns.re[i].get_educt_id(j)];

                for(size_t j=0; j<rns.re[i].get_no_product_s(); ++j)
                    rate_b *= x[rns.re[i].get_product_id(j)];
	   
                for(size_t j=0; j<rns.re[i].get_no_educt_s(); ++j)
                    dxdt[rns.re[i].get_educt_id(j)] -= (rate_f - rate_b);

                for(size_t j=0; j<rns.re[i].get_no_product_s(); ++j)
                    dxdt[rns.re[i].get_product_id(j)] += (rate_f - rate_b);
             }

            for(size_t i=0; i<rns.sp.size(); ++i)
                if(rns.sp[i].is_constant())
	            dxdt[i]=0;
         }
    };

    struct stiff_system {
        reaction_network_system& rns;

        stiff_system(reaction_network_system& rns_) : rns(rns_)  {}


        /*
         * 
         */

        void calculate_rates(const vector_type &x, vector_type &rates) {
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

            rates = element_prod(rns.e_m_bEa, element_prod(a2, rns.e_bs_in) - element_prod(a3, rns.e_bs_out));
        }

        /*
         * Prototype for the rhs of ODE
         * 
         */

        void operator()( const vector_type &x , vector_type &dxdt , double /* t */ ) {
            vector_type a2(rns.re.size()), a3(rns.re.size());
            for(size_t j=0; j<rns.N.size2(); ++j) { 
                a2(j) = 1;
                a3(j) = 1; 

                for(size_t l=0; l<rns.N_in_list[j].size(); ++l) 
                        a2(j) *= pow(x(rns.N_in_list[j][l]), rns.N_in(rns.N_in_list[j][l],j));

                for(size_t l=0; l<rns.N_out_list[j].size(); ++l) 
                        a3(j) *= pow(x(rns.N_out_list[j][l]), rns.N_out(rns.N_out_list[j][l],j));

            }         

            dxdt = prod(rns.N, element_prod(rns.e_m_bEa, element_prod(a2, rns.e_bs_in) - element_prod(a3, rns.e_bs_out)));

            for(size_t k=0; k<x.size(); ++k)
                if(rns.sp[k].is_constant())
                    dxdt(k) = 0;
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
            //matrix_type m(x.size(),x.size());
           
            for(size_t i=0; i<x.size(); ++i) {
                for(size_t l=0; l<x.size(); ++l) {
                    double s(0);
                    
                    for(size_t k=0; k<rns.re.size(); ++k) {    
                        double p1(1), p2(1);    // p1 = \prod_{j \neq k} {x_j}^N_{jk}^\mathrm{in}
                                                // p2 = \prod_{j \neq k} {x_j}^N_{jk}^\mathrm{out}

                        for(size_t m=0; m<rns.N_in_list[k].size(); ++m) {
                            size_t j=rns.N_in_list[k][m];
                            
                            if(j != l)
                                p1 *= x(j);
                        }
    
                        for(size_t m=0; m<rns.N_out_list[k].size(); ++m) {
                            size_t j=rns.N_out_list[k][m];

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




    /*
     *
     */

    reaction_network_system(const std::string& fn_network, const std::string& fn_concentration_, bool write_rate_=false) 
        : fn_concentration(fn_concentration_), beta(1), write_rate(write_rate_)  {
        
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

        /* Build lists (for speedup)
         */
        
        for(size_t i=0; i<re.size(); ++i) {
            N_in_list.push_back(std::vector<size_t>());
            N_out_list.push_back(std::vector<size_t>());

            for(size_t j=0; j<sp.size(); ++j) {
                if(N_in(j,i) != 0)
                    N_in_list[i].push_back(j);

                if(N_out(j,i) != 0)
                    N_out_list[i].push_back(j);
            }
        }
        
        // Read last line of concentration file and write initial concentration to initial_con
        // and initial time to initial_t. When done, open the same file for appending...
        initial_con = boost::numeric::ublas::zero_vector<double>(sp.size());
        initial_rate = boost::numeric::ublas::zero_vector<double>(re.size());

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


    /*
     * To check some differential equations or the algorithm it might be usefull
     * to calculate the right hand side of the differential equation for a given
     * concentration vector.
     */

    void print_rhs() {
        vector_type x(initial_con);
        vector_type dxdt=boost::numeric::ublas::zero_vector<double>(sp.size());
 
        // calculate rhs - operator saves it do dxdt
        stiff_system(*this)(x, dxdt, 0);

        cout << "righthandside:" << endl;
        for(size_t i=0; i<dxdt.size(); ++i)
            cout << "/  " << dxdt(i) << "  /";
        cout << endl;
    }


    /*
     *
     *
     */

    void run(double Tmax=25000, double deltaT=0.1, bool write_rates=false, 
             bool solve_implicit=true, size_t wint=500, size_t write_log=0) {
                
        fstream out(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);
        out.precision(25);

        size_t t0 = time(NULL);

        // Init solver and start
        vector_type x(initial_con), last_con(initial_con);
        double last_write(initial_t);
         
        auto write_state = [this, &out, t0, wint, deltaT, Tmax, &last_con, &last_write]( const vector_type &vec , const double t ) {
            double last_msd=0;
            for(size_t i=0; i<vec.size(); ++i)
                 last_msd += (vec[i]-last_con[i])*(vec[i]-last_con[i]);

            if(t != last_write)
                last_msd /= (vec.size()*(t-last_write));
            last_write=t;
            last_con=vec;

            out << t << "," << last_msd;

            for(size_t l=0; l<vec.size(); ++l)
                out << "," << vec(l);

            out << std::endl;
        };


        size_t step_no = 0;  // Number of steps done by integrator (for diagnostics)
        std::vector<double> times( wint );
        for( size_t i=0 ; i<wint ; ++i ) 
            if(write_log == 0) 
                times[i] = initial_t + double(i+1)/(wint)*(Tmax-initial_t);
            else if(write_log == 1 || initial_t == 0) 
                times[i] = initial_t + (exp(double(i+1))-1)/(exp(double(wint))-1)*(Tmax-initial_t);
            else 
                times[i] = initial_t*pow(exp(1/double(wint)*log(Tmax/initial_t)),double(i+1));

        step_no = solve_implicit ?
                      integrate_times( make_controlled< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
                                       make_pair( stiff_system(*this) , stiff_system_jacobi(*this) ) ,
                                       x , times, deltaT, write_state) :
                      integrate_times( make_controlled< runge_kutta_dopri5< vector_type > >( 1.0e-6 , 1.0e-6 ) ,
                                       stiff_system(*this) , x, times, deltaT, write_state);


        size_t t1 = time(NULL);
        std::cout << "Run took " << t1-t0 << " seconds and " << step_no << " steps!" << std::endl;

        out.close();
   }


   /* 
    * The fastint integration (calculating forward and backward rates directly 
    * through reaction constants) needs the networks to be exanded what this 
    * method ensures (multiple occurance of same educt product in same reaction
    * has to be explicit - "A + A" instead of "2 A").
    */

    void initialize_fastint() {
        for(size_t i=0; i<re.size(); ++i) {
            for(size_t j=0; j<re[i].get_no_educt_s(); ++j) { 
                if(re[i].get_educt_mul(j) > 1.1) {
                    re[i].add_educt_s(re[i].get_educt_id(j), re[i].get_educt_mul(j)-1);
                    re[i].set_educt_mul(j,1.0);
                }
            }

            for(size_t j=0; j<re[i].get_no_product_s(); ++j) { 
                if(re[i].get_product_mul(j) > 1.1) {
                    re[i].add_product_s(re[i].get_product_id(j), re[i].get_product_mul(j)-1);
                    re[i].set_product_mul(j,1.0);
                }
            }
        }
    }

    void run_fastint(double Tmax=25000, double deltaT=0.1, size_t wint=500, size_t write_log=0) {
                
        fstream out(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);
        out.precision(25);

        size_t t0 = time(NULL);

        vector_type x(initial_con), last_con(initial_con);
        double last_write(initial_t);
         
        auto write_state = [this, &out, t0, wint, deltaT, Tmax, &last_con, &last_write]( const vector_type &vec , const double t ) {
            // calculate meas square distance from last write
            double last_msd=0;
            for(size_t i=0; i<vec.size(); ++i)
                 last_msd += (vec[i]-last_con[i])*(vec[i]-last_con[i]);
            if(t != last_write)
                last_msd /= (vec.size()*(t-last_write));
            last_write=t;
            last_con=vec;

            // write time + msd...
            out << t << "," << last_msd;
            // ...and all concentrations
            for(size_t l=0; l<vec.size(); ++l)
                out << "," << vec(l);

            // If 
            if(last_msd > 1e-44 && last_msd < 1e-20) {
                std::cout << "Reached msd < 1e-20 condition - exiting early!" << std::endl;
                size_t t1 = time(NULL);
                std::cout << "Run took " << t1-t0 << " seconds!" << std::endl;
                exit(0);
            }

            out << std::endl;
        };

        // calculate vector of time points at which "write_state" will be called
        std::vector<double> times( wint );
        for( size_t i=0 ; i<wint ; ++i ) 
            if(write_log == 0) 
                times[i] = initial_t + double(i+1)/(wint)*(Tmax-initial_t);
            else if(write_log == 1 || initial_t == 0) 
                times[i] = initial_t + (exp(double(i+1))-1)/(exp(double(wint))-1)*(Tmax-initial_t);
            else 
                times[i] = initial_t*pow(exp(1/double(wint)*log(Tmax/initial_t)),double(i+1));

        // Call integrator - returns number of steps
        size_t step_no = integrate_times( make_dense_output< runge_kutta_dopri5< vector_type > >( 1.0e-6 , 1.0e-6 ) ,
                         fast_system(*this) , x, times, deltaT, write_state);

        // print time + steps needed and close file
        size_t t1 = time(NULL);
        std::cout << "Run took " << t1-t0 << " seconds and " << step_no << " steps!" << std::endl;
        out.close();
   }
};


/*
 * Main-Function. Manages commandline parameter and prints help screen.
 * All the real work is done through using the reaction_network_system 
 * class (above).
 */

int main(int argc, const char *argv []){
    srand(time(0));
    std::cout << "odeint_rnet version " << odeint_rnet_version_string << " (commit:" 
              << GIT_VERSION << ")" << std::endl;
    
    cl_para cl(argc, argv);  // create class to query command line parameters
    

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

        // load network and concentration file and call function to print right hand side
        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration); 
        rns.print_rhs();  
    }  


    // Simulates / integrates a reaction network while holding the boundary point 
    // species (constant=true) constant. Effective reaction rates are calculated 
    // directly from chemical energies. Stiff solver is available with option
    // "solve_implicit". Very slow for larger networks (>50 species).
    
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

        // Initial time-step for integrator
        double deltaT = cl.have_param("deltaT") ? cl.get_param_d("deltaT") : 0.1;   
        // ODE is integrated up to <Tmax>
        double Tmax = cl.have_param("Tmax") ? cl.get_param_d("Tmax") : 25000;  
        // Number of times the output is written between initial time and <Tmax> 
        double wint = cl.have_param("wint") ? cl.get_param_d("wint") : 500;
        // Also rates are given as output (to the file fn_concentration)
        bool write_rates=cl.have_param("write_rates");
        // Use implicit solver?
        bool solve_implicit=cl.have_param("solve_implicit");    
        // Is output given logarithmically spaced? (Or linearly?)        
        size_t write_log=0;
        if(cl.have_param("write_log"))
            write_log=1;
        if(cl.have_param("write_log_abs"))
            write_log=2;
            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax
                  << "   wint=" << wint << std::endl;
        if(write_rates) 
            cout << "Output contains reactions' effective rates." << endl;

        if(write_log == 1) 
            cout << "Period of output will be equidistant on logscale (from t=t_0)." << endl;

        if(write_log == 2) 
            cout << "Period of output will be equidistant on logscale (from t=0)." << endl;

        if(solve_implicit)
            cout << "Stiff solver is used!" << endl;

        // Load reaction network and concentration file
        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration, write_rates); 

        // Simulate ODE (and write results to file)
        // The method gives diagnostic feedback to the user through console output
        rns.run(Tmax, deltaT, write_rates, solve_implicit, wint, write_log);
    }  


    // Simulates / integrates a reaction network while holding the boundary point 
    // species (constant=true) constant. Official parameter to call this is 
    // "fastint", the option to call using "simsim" is maintained for backward
    // compatibility. This integrator is faster as the one above, especially for
    // larger (>50 species) networks. It uses the reaction constants present in 
    // the network and no stiff solver. 

    if(cl.have_param("fastint") || cl.have_param("simsim")) {
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

        // Initial time-step for integrator
        double deltaT = cl.have_param("deltaT") ? cl.get_param_d("deltaT") : 0.1;   
        // ODE is integrated up to <Tmax>
        double Tmax = cl.have_param("Tmax") ? cl.get_param_d("Tmax") : 25000;  
        // Number of times the output is written between initial time and <Tmax> 
        double wint = cl.have_param("wint") ? cl.get_param_d("wint") : 500;
        // Is output given logarithmically spaced? (Or linearly?)        
        size_t write_log=0;
        if(cl.have_param("write_log"))
            write_log=1;
        if(cl.have_param("write_log_abs"))
            write_log=2;
            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax
                  << "   wint=" << wint << std::endl;

        if(write_log == 1) 
            cout << "Period of output will be equidistant on logscale (from t=t_0)." << endl;

        if(write_log == 2) 
            cout << "Period of output will be equidistant on logscale (from t=0)." << endl;


        // Load reaction network and concentration file
        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration, false); 

        // Simulate ODE (and write results to file)
        // The method gives diagnostic feedback to the user through console output
        rns.initialize_fastint();
        rns.run_fastint(Tmax, deltaT, wint, write_log);
    }  

    
    // User interface / help dialogue
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
        cout << "   x'deltaT': Integration interval (relevant for integrator, not for output!)" << endl;
        cout << "   x'Tmax': Time up to which the system is simulated" << endl;
        cout << "   x'write_rates': 'con' file contains effective reaction rates" << endl;
        cout << "   x'wint': number of output times between Tstart and Tmax" << endl;
        cout << "   x'write_log': write output logarithmically spaced (from t=t_0)" << endl;
        cout << "   x'write_log_abs': write output logarithmically spaced from (t=0)" << endl;
        cout << "   x'solve_implicit': use stiff solver (might be slower for big nets!)" << endl;
        cout << endl;
        cout << "-'fastint': load reaction network 'net' and simulate file 'con'!. Other than above" << endl;
        cout << "   this integrator does not offer a stiff solver. It also does not use energies to" << endl;
        cout << "   calculate effective rates directly but uses reaction constants. Parameters are:" << endl;
        cout << "   x'deltaT': Integration interval (relevant for integrator, not for output!)" << endl;
        cout << "   x'Tmax': Time up to which the system is simulated" << endl;
        cout << "   x'wint': number of output times between Tstart and Tmax" << endl;    
        cout << "   x'write_log': write output logarithmically spaced (from t=t_0)" << endl;
        cout << "   x'write_log_abs': write output logarithmically spaced from (t=0)" << endl;
	cout << endl;
    } 
}
