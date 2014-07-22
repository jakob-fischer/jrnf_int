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



typedef std::vector< double > state_type;

std::vector<species> sp;
std::vector<reaction> re;
std::vector<double> initial_con, last_con;
double initial_t(0.0), last_write(0.0), deltaT(0.1), Tmax(25000);
size_t wint=10000;
size_t t0=0;


std::vector<bool> bi_reaction;
std::vector<bool> const_vec;
ofstream out;



void init_state(std::vector<double>& vec) {
    for(size_t i=0; i<sp.size(); ++i) {
        vec.push_back(initial_con[i]);
        last_con.push_back(initial_con[i]);
	const_vec.push_back(sp[i].is_constant());
    }
}


//
//
//

void next_step( const state_type &vec , state_type &dxdt , double t ) {
    for(size_t i=0; i<sp.size(); ++i) 
        dxdt[i]=0;
    
    for(size_t i=0; i<re.size(); ++i) {
        double rate_f=re[i].get_k()*vec[re[i].get_educt_id(0)];
	double rate_b=re[i].get_k_b()*vec[re[i].get_product_id(0)];
	
	if(bi_reaction[i]) {
	    rate_f*=vec[re[i].get_educt_id(1)];
	    rate_b*=vec[re[i].get_product_id(1)];  
	}
	   
	// Educt 0
	dxdt[re[i].get_educt_id(0)] -= (rate_f - rate_b);
	
	// Product 1
	dxdt[re[i].get_product_id(0)] += (rate_f - rate_b);
	
	if(bi_reaction[i]) {
	    // Educt 1
	    dxdt[re[i].get_educt_id(1)] -= (rate_f - rate_b); 
	  
	    // Product 1
	    dxdt[re[i].get_product_id(1)] += (rate_f - rate_b);
	}
    }

    for(size_t i=0; i<const_vec.size(); ++i)
        if(const_vec[i])
	    dxdt[i]=0;
}



void next_step_ud_c( const state_type &vec , state_type &dxdt , double t ) {
    for(size_t i=0; i<sp.size(); ++i) 
        dxdt[i]=0;

    //cout << "re.size()=" << re.size() << endl;
    
    for(size_t i=0; i<re.size(); ++i) {
        double rate_f=re[i].get_k()*pow(vec[re[i].get_educt_id(0)], re[i].get_educt_mul(0));
	
	for(size_t j=1; j<re[i].get_no_educt_s(); ++j)
	    rate_f*=pow(vec[re[i].get_educt_id(j)], re[i].get_educt_mul(j));
	   
	for(size_t k=0; k<re[i].get_no_educt_s(); ++k) 
	    dxdt[re[i].get_educt_id(k)] -= (rate_f*re[i].get_educt_mul(k));   
	
	for(size_t k=0; k<re[i].get_no_product_s(); ++k) 
	    dxdt[re[i].get_product_id(k)] += (rate_f*re[i].get_product_mul(k));
    }

    for(size_t i=0; i<const_vec.size(); ++i)
        if(const_vec[i])
	    dxdt[i]=0;
}




void do_write_state( const state_type &vec , const double t, size_t count=0 ) {
    // calculate qdiff
    double qdiff=0;
    for(size_t i=0; i<vec.size(); ++i)
        qdiff += (last_con[i] - vec[i])*(last_con[i] - vec[i]);

    if(t - last_write != 0)
        qdiff /= ((t-last_write)*vec.size());

    // update last_con
    for(size_t i=0; i<vec.size(); ++i) 
        last_con[i] = vec[i];


    bool abort_early=qdiff < 1e-20 && qdiff > 1e-44;


    // output if sufficient time has passed
    if((t-last_write) > 2.0*(last_write-initial_t) || t > Tmax-deltaT ||
       abort_early) {
        if(count == 0)
            out << std::endl;

        out << t << "," << qdiff;
        
        for(size_t i=0; i<sp.size(); ++i) 
            out << "," << vec[i];  

        last_write = t;

        out << std::endl;
    }

    if(qdiff < 1e-20 && qdiff > 1e-44) {
        size_t t1 = time(NULL);
        std::cout << "Run took " << t1-t0 << " seconds (finishing with <e-20 cond)!" << std::endl;

        exit(0);
    }
}


void write_state( const state_type &vec , const double t ) {
    static size_t count=0;
    
    if(count % wint == 0) 
        do_write_state(vec, t, count);        

    ++count;
}



//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


bool is_10or01_reaction(reaction& re) { 
    return re.get_no_educt() == 1 && re.get_no_product() == 0 || re.get_no_educt() == 0 && re.get_no_product() == 1;
}




class reaction_network_system {
    double initial_t, last_write;
    std::vector<species> sp;
    std::vector<reaction> re;
    //std::vector<double> initial_con, last_con;
    std::string fn_concentration;

public:
    struct stiff_system {
        const reaction_network_system& rns;

        stiff_system(const reaction_network_system& rns_) : rns(rns_)  {}

        void operator()( const vector_type &x , vector_type &dxdt , double /* t */ ) {
            dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ];
            dxdt[ 1 ] = x[ 0 ];
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
        : fn_concentration(fn_concentration_)  {
        
        read_jrnf_reaction_n(fn_network, sp, re);
    
	
	// Remove 1-0 and 0-1 reactions. Such reactions are there to ballance flow through the boundary conditions 
        std::remove_if (re.begin(), re.end(), 
                        [] (reaction& re) -> bool { 
                            return re.get_no_educt() == 1 && re.get_no_product() == 0 || 
                                   re.get_no_educt() == 0 && re.get_no_product() == 1; });
		
	
        // Checking that there are only 2-2 and 1-1 reactions. Also bring all 2-2 reactions
	// into a normal form.

/*

        for(size_t i=0; i<re.size(); ++i) {
            if(re[i].get_no_educt() == 2 && re[i].get_no_product() == 2) {
	            bi_reaction.push_back(true);  
	        
		        if(re[i].get_no_educt_s() == 1) {
	                re[i].add_educt_s(re[i].get_educt_id(0));
		            re[i].set_educt_mul(0,1.0);
		        }
	    
	            if(re[i].get_no_product_s() == 1) {
	                re[i].add_product_s(re[i].get_product_id(0));	   
		            re[i].set_product_mul(0,1.0);
	            }
	        } else if(re[i].get_no_educt() == 1 && re[i].get_no_product() == 1) {
	            bi_reaction.push_back(false);
	        } else {
	            cout << "Reaction " << i << " : invalid reaction. Only 1-1 and 2-2 r. allowed!" << endl;
		    cout << re[i].get_string() << endl;
	        }
            }


*/       




        
        // Read last line of concentration file and write initial concentration to initial_con
        // and initial time to initial_t. After done open the same file for appending...
        for(size_t i=0; i<sp.size(); ++i) 
            initial_con.push_back(0.0);

        std::ifstream  data(fn_concentration.c_str());

        if(!data.good()) {
           std::cout << "Could not open concentration file: " << fn_concentration << std::endl;
        }

        std::string line;
        std::getline(data,line);       // dont want the header
        double last_msd=0;

        // the real 
        while(!std::getline(data,line).eof()) {
            std::stringstream ls(line);
            std::string cell;
                 
            //cout << "line: " << line << std::endl;           
   
            size_t cnt=0;
            while(std::getline(ls,cell,',')) {
                std::stringstream in(cell);       
    
    
                if(cnt == 0) 
                    in >> initial_t ;
                else if(cnt == 1)
                    in >> last_msd;                
                else if(cnt <= sp.size()+1)
                    in >> initial_con[cnt-2];
                else 
                    std::cout << "Error at reading csv / concentration file!" << std::endl;

                ++cnt;
            }
        }

        std::cout << "Simulating file: " << fn_concentration << std::endl;
        std::cout << "Loaded concentration file with starting time " << initial_t << std::endl;
        
        if(last_msd > 1e-44 & last_msd < 1e-20) {
            std::cout << "1e-20-condition fulfilled: system already relaxed sufficiently!" << std::endl;
        }

          
        data.close();        
         
    }



    void run(double Tmax_=25000, double deltaT_=0.1, size_t wint_=1000) {
       out.open(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);
        out.precision(25);

        t0 = time(NULL);

        // Init solver and start
        state_type x;
        init_state(x);
        integrate( next_step , x , initial_t , Tmax , deltaT , write_state );
        //runge_kutta_dopri5< state_type > rk;

        //integrate_adaptive(rk, next_step , x , initial_t , Tmax , deltaT , write_state );
        //runge_kutta4

         do_write_state(x, Tmax);

        size_t t1 = time(NULL);
        std::cout << "Run took " << t1-t0 << " seconds!" << std::endl;


   }



};








//void write_state( const vector_type &vec , const double t ) {
//    cout << t << "," << vec[0] << "," << vec[1] << endl;
//}





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

        deltaT = cl.have_param("deltaT") ? cl.get_param_d("deltaT") : 0.1;   
        Tmax = cl.have_param("Tmax") ? cl.get_param_d("Tmax") : 25000;
        wint = cl.have_param("wint") ? cl.get_param_i("wint") : 1000;
        bool write_rates=cl.have_param("write_rates");
        bool solve_implicit=cl.have_param("solve_implicit");            

            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax << std::endl;

        reaction_network_system rns = reaction_network_system(fn_network, fn_concentration); 
        rns.run(Tmax, deltaT, wint);
    }  
    
    
    if(cl.have_param("help") || cl.have_param("info")) {
        cout << "          odeint_rnet" << endl;  
        cout << "          ===========" << endl;
        cout << " call with parameter 'info' or 'help' for showing this screen" << endl;
        cout << endl;
        cout << "-'simsim' load reaction network 'net' and simulate file 'con'!. Parameters are 'deltaT', and 'Tmax'." << endl;
	cout << endl;
    } 
}
