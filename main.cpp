#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "cl_para.h"
#include "reaction_network.h"
#include "reaction_network_fileop.h"
#include "network_tools.h"


using namespace std;
using namespace boost::numeric::odeint;



typedef std::vector< double > state_type;

std::vector<species> sp;
std::vector<reaction> re;
size_t fixed_1=0; 
double value_1=2.0;
double value_2=1.0;
size_t fixed_2=0;
std::vector<bool> bi_reaction;
ofstream out;

std::vector<bool> const_vec;


void init_state(std::vector<double>& vec) {
    for(size_t i=0; i<sp.size(); ++i) {
	    vec.push_back(1.0);
    
	    if(fixed_1 == i || fixed_2 == i)
	        const_vec.push_back(true);
	    else 
	        const_vec.push_back(false);
    }
    
    vec[fixed_1]=value_1;
    vec[fixed_2]=value_2;
}



void init_state(std::vector<double>& vec, const std::string& fixed_ic) {
    ifstream fx_ic(fixed_ic.c_str());
  
    for(size_t i=0; i<sp.size(); ++i) {
        double ic=-1;
	bool fx=false;
	fx_ic >> ic >> fx;
	    
	if(ic < 0) {
	    cout << "ERROR: with reading 'fixed'-file!" << endl;
	    return;  
	}
	    
	vec.push_back(ic);
	const_vec.push_back(fx);
    }
}




void next_step( const state_type &vec , state_type &dxdt , double t ) {
    for(size_t i=0; i<sp.size(); ++i) 
        dxdt[i]=0;

    //cout << "re.size()=" << re.size() << endl;
    
    for(size_t i=0; i<re.size(); ++i) {
        double rate_f=re[i].get_k()*vec[re[i].get_educt_id(0)];
	double rate_b=re[i].get_k_b()*vec[re[i].get_product_id(0)];
	
	if(bi_reaction[i]) {
	    rate_f*=vec[re[i].get_educt_id(1)];
	    rate_b*=vec[re[i].get_product_id(1)];  
	}
	   
        //cout << "rate_f=" << rate_f << "   rate_b=" << rate_b << endl;
	   
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
    if(count % 100 == 0) {
        out << t;
        for(size_t i=0; i<sp.size(); ++i) {
	    out << "   " << vec[i];  
	}
	out << endl;
    }
}

int main(int argc, const char *argv []){
    srand(time(0));
    
    cl_para cl(argc, argv);
    
    if(cl.have_param("sim_sim_1")) {
        // Simulates / integrates a reaction network while houlding two boundary points constant
        // 
      
        if(!cl.have_param("in")) {
            cout << "You have to give the name of input by 'in'!" << endl;  
            return 1;
        }
    
        if(!cl.have_param("out")) {
            cout << "You have to give the name of the output file!" << endl;   
            return 1;
        }
    
        std::string fn_network=cl.get_param("in");
        std::string fn_output=cl.get_param("out");
        out.open(fn_output.c_str());
    
        if(cl.have_param("f1"))
          fixed_1=cl.get_param_i("f1");
    
        if(cl.have_param("v1"))
          value_1=cl.get_param_d("v1");
    
        if(cl.have_param("f2"))
          fixed_2=cl.get_param_i("f2");
    
        if(cl.have_param("v2"))
          value_2=cl.get_param_d("v2");
    
        read_jrnf_reaction_n(fn_network, sp, re);
    
        // Checking that there are only 2-2 and 1-1 reactions 
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
	        bi_reaction.push_back(false);        // Checking that there are only 2-2 and 1-1 reactions 
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
    
        // Calculate rate coefficients for forward and backward reaction
        for(size_t i=0; i<re.size(); ++i) {
            double e_E = sp[re[i].get_educt_id(0)].get_energy();
            double e_P = sp[re[i].get_product_id(0)].get_energy();
	    double e_A = re[i].get_activation();
      
            if(bi_reaction[i]) {
	        e_E += sp[re[i].get_educt_id(1)].get_energy();
                e_P += sp[re[i].get_product_id(1)].get_energy();	  
	    }
	
	    re[i].set_k(exp(-(e_A-e_E)));
	    re[i].set_k_b(exp(-(e_A-e_P)));
        }
	    } else {
	        cout << "Reaction " << i << " : invalid reaction. Only 1-1 and 2-2 r. allowed!" << endl;
		cout << re[i].get_string() << endl;
	        return 1;
	    }
        }
    
        // Calculate rate coefficients for forward and backward reaction
        for(size_t i=0; i<re.size(); ++i) {
            double e_E = sp[re[i].get_educt_id(0)].get_energy();
            double e_P = sp[re[i].get_product_id(0)].get_energy();
	    double e_A = re[i].get_activation();
      
            if(bi_reaction[i]) {
	        e_E += sp[re[i].get_educt_id(1)].get_energy();
                e_P += sp[re[i].get_product_id(1)].get_energy();	  
	    }
	
	    re[i].set_k(exp(-(e_A-e_E)));
	    re[i].set_k_b(exp(-(e_A-e_P)));
        }
    
        state_type x;
        init_state(x);
        integrate( next_step , x , 0.0 , 2500.0 , 0.01 , write_state );
    }

    
    
    if(cl.have_param("sim_sim_2")) {
        // Simulates / integrates a reaction network while houlding two boundary points constant
        // 
      
        if(!cl.have_param("rnet")) {
            cout << "You have to give the name of input by 'rnet'!" << endl;  
            return 1;
        }
        
        if(!cl.have_param("initial")) {
            cout << "You have to give the name of input by 'initial'!" << endl;  
            return 1;
        }
    
        if(!cl.have_param("out")) {
            cout << "You have to give the name of the output file!" << endl;   
            return 1;
        }
    
        if(cl.have_param("do_pf")) {
	    if(!cl.have_param("pf_var") || 
	       !cl.have_param("pf_per")) {
                 cout << "Periodic forcing enabled. You have to give 'pf_var' and 'pf_per' parameter!" << endl;   
                 return 1;	      
	    } 
	}
    
        std::string fn_network=cl.get_param("rnet");
        std::string fn_output=cl.get_param("out");
	std::string fn_initial=cl.get_param("initial");
	
        out.open(fn_output.c_str());
    
	bool do_pf=cl.have_param("do_pf");
	size_t pf_var=cl.get_param_i("pf_var");
	double pf_per=cl.get_param_d("pf_per");
    
        read_jrnf_reaction_n(fn_network, sp, re);
    
    
        state_type x;
	cout << "building initial state." << endl;
        init_state(x, fn_initial);
	cout << "runing simulation." << endl;
        integrate( next_step_ud_c , x , 0.0 , 2500.0 , 0.01 , write_state );
    }

    
    
    
    if(cl.have_param("help") || cl.have_param("info")) {
        cout << "          odeint_rnet" << endl;  
        cout << "          ===========" << endl;
        cout << " call with parameter 'info' or 'help' for showing this screen" << endl;
        cout << endl;
        cout << "not implemented yet" << endl;
	cout << endl;
    } 
}