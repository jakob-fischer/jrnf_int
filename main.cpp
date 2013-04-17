/* author: jakob fischer (jakob@automorph.info)
 * date: 22nd March 2013
 * description: 
 * Tool for solving a ODE for a reaction system given as an jrnf-file. The concentration and time
 * are read and written from / to a comma seperated file. Every row / line represents one time step.
 * the first column contains the time of this step.
 * 
 * TODO Implement and comment
 * 
 * TODO Extend so concentration can be taken from external file (allowing for example periodic 
 *      boundary conditions)
 */

#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <algorithm>
#include "tools/cl_para.h"
#include "net_tools/reaction_network.h"
#include "net_tools/reaction_network_fileop.h"
#include "net_tools/network_tools.h"

using namespace std;
using namespace boost::numeric::odeint;



typedef std::vector< double > state_type;

std::vector<species> sp;
std::vector<reaction> re;
std::vector<double> initial_con, last_con;
double initial_t(0.0), last_write(0.0), deltaT(0.1), Tmax(25000);
size_t wint=10000;

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


void write_state( const state_type &vec , const double t ) {
    static size_t count=0;
    
    if(count % wint == 0) {
        // calculate qdiff
        double qdiff=0;
        for(size_t i=0; i<vec.size(); ++i)
            qdiff += (last_con[i] - vec[i])*(last_con[i] - vec[i]);

        if(t - last_write != 0)
            qdiff /= ((t-last_write)*vec.size());

        // update last_con
        for(size_t i=0; i<vec.size(); ++i) 
            last_con[i] = vec[i];


        // output if sufficient time has passed
        if((t-last_write) > 2.0*(last_write-initial_t)) {
            if(count == 0)
                out << std::endl;

            out << t << "," << qdiff;
        
            for(size_t i=0; i<sp.size(); ++i) 
	        out << "," << vec[i];  

            last_write = t;

            out << std::endl;
        }
        
    }

    ++count;
}


bool is_10or01_reaction(reaction& re) { 
    return re.get_no_educt() == 1 && re.get_no_product() == 0 || re.get_no_educt() == 0 && re.get_no_product() == 1;
}


int main(int argc, const char *argv []){
    srand(time(0));
    
    cl_para cl(argc, argv);
    
    if(cl.have_param("simsim")) {
        // Simulates / integrates a reaction network while holding the boundary point 
        // species (constant=true) constant
      
        if(!cl.have_param("net")) {
            cout << "You have to give the name of reaction network by 'net'!" << endl;  
            return 1;
        }
    
        if(!cl.have_param("con")) {
            cout << "You have to give the name of the concentration file by 'con'!" << endl;   
            return 1;
        }

        if(cl.have_param("deltaT")) 
            deltaT = cl.get_param_d("deltaT");
    
        if(cl.have_param("Tmax")) 
            Tmax = cl.get_param_d("Tmax");

        if(cl.have_param("wint"))
            wint = cl.get_param_i("wint");

        std::string fn_network=cl.get_param("net");
        std::string fn_concentration=cl.get_param("con");
            
        std::cout << "Parameters are deltaT=" << deltaT << "  and Tmax=" << Tmax << std::endl;


        read_jrnf_reaction_n(fn_network, sp, re);
    
	
	// Remove 1-0 and 1-0 reactions which are with constant species. Such reactions are
	// there to ballance flow over the boundary conditions
        std::remove_if (re.begin(), re.end(), is_10or01_reaction);
		
	
        // Checking that there are only 2-2 and 1-1 reactions. Also bring all 2-2 reactions
	// into a normal form.



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
	            return 1;
	        }
        }

        
        // Read last line of concentration file and write initial concentration to initial_con
        // and initial time to initial_t. After done open the same file for appending...
        for(size_t i=0; i<sp.size(); ++i) 
            initial_con.push_back(0.0);

        std::ifstream  data(fn_concentration.c_str());

        if(!data.good()) {
           std::cout << "Could not open concentration file: " << fn_concentration << std::endl;
           return 0;
        }

        std::string line;
        std::getline(data,line);       // dont want the header

        // the real 
        while(!std::getline(data,line).eof()) {
            std::stringstream ls(line);
            std::string cell;
                 
            cout << "line: " << line << std::endl;           
   
            size_t cnt=0;
            while(std::getline(ls,cell,',')) {
                std::stringstream in(cell);       
                double last_msd=0;
     
    
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

          
        data.close();
        out.open(fn_concentration.c_str(), std::ios_base::out | std::ios_base::app);


        // Init solver and start
        state_type x;
        init_state(x);
        integrate( next_step , x , initial_t , Tmax , deltaT , write_state );
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
