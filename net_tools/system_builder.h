#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include <string>
#include <vector>


/*
 * Builds a reaction network from a simple reaction and a connective structure
 * of compartments...
 */

class system_builder {
    std::vector<species> &sp_ref;
    std::vector<reaction> &re_ref;
    std::vector<compartment> &cp_ref;

    std::vector<species> sp;
    std::vector<reaction> re;
    
  public:
    std::vector<species>& get_species_v() {
        return sp;  
    }
    
    
    std::vector<reaction>& get_reactions_v() {
        return re;  
    }
    
    
    size_t translate_sp(size_t sp, size_t comp) {
        return comp*sp_ref.size() + sp;  
    }
    
    

    system_builder(std::vector<species> &s, std::vector<reaction> &r, std::vector<compartment> &c) 
	  : sp_ref(s), re_ref(r), cp_ref(c) {
        // Create a new species in sp for every combination of elements from 
	// sp_ref and cp_ref. Even if they aren't all needed
	for(size_t i=0; i<cp_ref.size(); ++i) {
	    for(size_t j=0; j<sp_ref.size(); ++j) {
	        sp.push_back(sp_ref[j]);
		sp.back().set_id(sp.size()-1);
		sp.back().set_name(sp_ref[j].get_name() + "_" + cp_ref[i].get_name());
	    }
	}
    }
    
    
    system_builder& add_reaction(size_t r_id, size_t cp_id) {
        // first copy reaction
        reaction rn(re_ref[r_id]); 
	
	// and then translate all educt and product ids
        for(size_t i=0; i<rn.get_no_educt_s(); ++i)
	    rn.set_educt_id(i, translate_sp(rn.get_educt_id(i), cp_id));
	
        for(size_t i=0; i<rn.get_no_product_s(); ++i)
	    rn.set_product_id(i, translate_sp(rn.get_product_id(i), cp_id));
	
	re.push_back(rn);
	  
        return *this;  
    }
    
        
    system_builder& add_reaction_c_by_T(size_t r_id, size_t cp_id) {
        // first copy reaction
        reaction rn(re_ref[r_id]); 
	double E_educt=0;
	
	// and then translate all educt and product ids
        for(size_t i=0; i<rn.get_no_educt_s(); ++i) 
	    rn.set_educt_id(i, translate_sp(rn.get_educt_id(i), cp_id));
	    
        for(size_t i=0; i<rn.get_no_product_s(); ++i)
	    rn.set_product_id(i, translate_sp(rn.get_product_id(i), cp_id));
	
	rn.set_c_by_T(sp, cp_ref[cp_id].get_temperature());
	
	
	re.push_back(rn);
	  
        return *this;  
    }
    
    
    

    system_builder& add_flow(size_t s_id, size_t cp1_id, size_t cp2_id, double rate) {
        reaction rn;
	rn.set_activation(0.0);
	rn.set_c(rate).set_k(rate);
	rn.add_educt(translate_sp(s_id, cp1_id), 1);
	rn.add_product(translate_sp(s_id, cp2_id), 1);
	
	re.push_back(rn);
      
        return *this;  
    }


};


#endif // SYSTEM_BUILDER_H