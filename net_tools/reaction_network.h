#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <cmath>
#include <map>



/*
 * Class 'species' represents one species in a chemical reaction or in a  
 * reaction-network. Species are identified by a name and an integer id
 * which corresponds to the indice in std::vector<species> ...
 *
 * Every species has a type, Normal ones is set to '0'. Other values imply
 * that the species has a none standard dynamics or is only used for monitoring
 * flux.
 */

class species {
    // Integer id
    size_t id;
    
    // 0 == standard   //  1 == monitoring heat species
    size_t type;    
    
    // Is this species going to stay constant 
    bool constant;
    
    // Bind-Enthalpy of this species
    double energy;
    
    // Name / Identifier of the species
    std::string name;
  
  public:
    
   
    
    /*
     * Constructor. Set the name of species as you want
     */
    
    species(size_t i, const std::string& n, bool c=false, size_t t=0) {
	id=i;  name=n;  type=t;    constant=c;    
    }

    
    /*
     * Set name of species
     */
    
    species& set_name(const std::string& n) {
	name=n;
	return *this;
    }
    
    
    /* 
     * Get name of species
     */
    
    const std::string& get_name() const {
	return name;
    }
    
    /*
     * Get the id
     */
    
    size_t get_id() const {
        return id;
    }
    
    /*
     * Change the id
     */
    
    species& set_id(size_t i) {
        id=i;
	return *this;
    } 
    
    
    /*
     * Get species type
     */
    
    const size_t get_type() const {
        return type;  
    }
    
    
    /*
     * Change species type
     */
    
    species& set_type(size_t t) {
        type=t;
	return *this;
    }
    
    /* 
     * Get species energy
     */
    
    double get_energy() const {
        return energy;  
    }
    
    /*
     * Change species energy
     */
    
    species& set_energy(double e) {
        energy=e;  
	return *this;
    }
    
    /*
     * Is this species kept constant in simulation
     */
    
    bool is_constant() const {
        return constant;  
    }
    
    /*
     * Set whether this species should be kept constant
     */
    
    void set_constant(bool c) {
        constant = c;  
    }
    
    // Friend declaration of output function / operator
    friend std::ostream& operator<<(std::ostream& os, const species& sp);
    
    // Friend declaration of input function / operator
    friend std::istream& operator>>(std::istream& is, species& sp);
};


/*
 * Operator for output of a specis object. Output can be understood by humans
 * but also be done to a file and reread into a species object later...
 */

std::ostream& operator<<(std::ostream& os, const species& sp) {
    os << sp.type << " " << sp.name << " " << sp.energy << " " << sp.constant;
        
    return os;
} 


/*
 * Operator for input of a species object
 */

std::istream& operator>>(std::istream& is, species& sp) {
    is >> sp.type >> sp.name >> sp.energy >> sp.constant;
    
    return is;
} 


/*
 * Class 'compartment'
 */

class compartment {
    // Integer id
    size_t id;
    
    // Volume
    double volume;
    
    // Temperature
    double temperature;
    
    // Name / Identifier of the compartment
    std::string name;
  
  public:
    
    /*
     * Constructor. Name is derived from id... 
     */
    
    compartment(size_t i) : temperature(1.0) {
        id=i;
	
	std::stringstream ss;
	ss << i;
	name = ss.str();
    }
    
    
    /*
     * Constructor. Set the name of compartment as you want
     */
    
    compartment(size_t i, const std::string& n) : temperature(1.0) {
	id=i;  name=n;      
    }

    
    /*
     * Set name of compartment
     */
    
    void set_name(const std::string& n) {
	name=n;
    }
    
    
    /* 
     * Get name of compartment
     */
    
    const std::string& get_name() const {
	return name;
    }
    
    /*
     * Get the id
     */
    
    size_t get_id() const {
        return id;
    }    
    
    /*
     * Change the id
     */
    
    compartment& set_id(size_t i) {
        id=i;
	return *this;
    } 
    
    /*
     * Get volume
     */
    
    double get_volume() const {
        return volume;
    }    
    
    /*
     * Change volume
     */
    
    compartment& set_volume(double v) {
        volume=v;
	return *this;
    } 
    
    /*
     *
     *
     */
    
    double get_temperature() const {
        return temperature;  
    }
    
    compartment& set_temperature(double t) {
        temperature=t;
	
	return *this;
    }
    
    

    friend std::ostream& operator<<(std::ostream& os, const compartment& cp);
    
    friend std::istream& operator>>(std::istream& is, compartment& cp);
};


std::ostream& operator<<(std::ostream& os, const compartment& cp) {
    os << cp.name << " " << cp.volume;
    
    return os;
} 


std::istream& operator>>(std::istream& is, compartment& cp) {
    is >> cp.name >> cp.volume;
    
    return is;
} 




/*
 * The 'reaction' class represents one single reaction in the reaction network
 * The reaction is initialiced with reaction propensity 1 and no species as 
 * educts and products.
 */

class reaction {
    double c, k, k_b, activation;
    bool reversible;
    std::vector< std::pair<size_t,double> > educts, products;
  
  
  public:
    reaction() {
       c=1.0;
       k=1.0;
       k_b=1.0;
       activation=0.0;
       reversible=false;
    }
   
    bool is_reversible() const {
        return reversible;  
    }
   
    reaction& set_reversible(bool rev) {
        reversible = rev;
	return *this;      
    }
   
    double get_c() const {
        return c;
    }
   
    reaction& set_c(double c_) {
	c=c_;
	return *this;
    }
    
    reaction& set_c_by_T(const std::vector<species>& sp, double T) {
      	double E_educt=0;
	
	// and then translate all educt and product ids
        for(size_t i=0; i<educts.size(); ++i) 
	    E_educt += sp[educts[i].first].get_energy()*educts[i].second;
	
	c = exp(-(activation-E_educt)/T);
      
      
        return *this;  
    }
    
    double get_k() const {
        return k;
    }
   
    reaction& set_k(double k_) {
	k=k_;
	return *this;
    }
    
    double get_k_b() const {
        return k_b;  
    }
    
    reaction& set_k_b(double k_b_) {
        k_b=k_b_;  
        return *this;
    }
    
    double get_activation() const {
        return activation;
    }
   
    reaction& set_activation(double activation_) {
	activation=activation_;
	return *this;
    }
    
    
   
  
    size_t get_no_educt_s() const {
	return educts.size();
    }
    
    size_t get_no_educt() const {
        size_t a=0;
	
	for(size_t i=0; i<educts.size(); ++i)
	    a += educts[i].second;
	
	return a;
    }
    
    size_t get_educt_id(size_t i) const {  return educts[i].first;  }
    
    size_t set_educt_id(size_t i, size_t n) {  educts[i].first = n;  }
    
    double get_educt_mul(size_t i) const {  return educts[i].second;  }
    
    reaction& set_educt_mul(size_t i, double m) {  educts[i].second=m;  return *this;  }

    reaction& set_product_mul(size_t i, double m) {  products[i].second=m;  return *this;  }
    
    reaction& add_educt_s(size_t id, double m=1.0) {
       educts.push_back( std::pair<size_t,size_t>(id, m) );
       
       return *this;
    }
    
    reaction& add_educt(size_t id, double m=1.0) {
       for(size_t i=0; i<educts.size(); ++i) 
	   if(educts[i].first == id) {
	      educts[i].second += m;
	      return *this;
	   }
	   
        educts.push_back( std::pair<size_t,size_t>(id, m) );
       
       return *this;
    }


    size_t get_no_product_s() const {
	return products.size();
    }
    
    size_t get_no_product() const {
        size_t a=0;
	
	for(size_t i=0; i<products.size(); ++i)
	    a += products[i].second;
	
	return a;
    }
    
    size_t get_product_id(size_t i) const {  return products[i].first;  }
    
    size_t set_product_id(size_t i, size_t n) {  products[i].first = n;  }
    
    
    double get_product_mul(size_t i) const {  return products[i].second;  }
    
    reaction& add_product_s(size_t id, double m=1.0) {
	products.push_back( std::pair<size_t,double>(id, m) );
       
       return *this;
    }
    
    reaction& add_product(size_t id, double m=1.0) {
       for(size_t i=0; i<products.size(); ++i) 
	   if(products[i].first == id) {
	      products[i].second += m;
	      return *this;
	   }        
	   
	products.push_back( std::pair<size_t,double>(id, m) );
       
       return *this;
    }
    
    std::string get_string() const {
        std::stringstream ss;
	
	for(size_t i=0; i<educts.size(); ++i) {
	    ss << educts[i].second << " ";
	    ss << "<" << educts[i].first << ">";
	    
	    if(i+1!=educts.size())
	      ss << " + ";
	}
	
	if(reversible)
	    ss << " <--" << c << "-> ";  
        else
	    ss << " --" << c << "-> ";
	
	for(size_t i=0; i<products.size(); ++i) {
	    ss << products[i].second << " ";
	    ss << "<" << products[i].first << "> ";
	    
	    if(i+1!=products.size())
	      ss << " + ";
	}	
      
        return ss.str();
   }
   
   
    std::string get_string_b(const std::vector<species> &sp) const {
        std::stringstream ss;
	
	if(!reversible)
	    return std::string("");
	
	double e_e(0.0), e_p(0.0);
	
	for(size_t i=0; i<educts.size(); ++i)
	    e_e += sp[educts[i].first].get_energy()*educts[i].second;
	
	for(size_t i=0; i<products.size(); ++i)
	    e_p += sp[products[i].first].get_energy()*products[i].second;    
	
	
	for(size_t i=0; i<educts.size(); ++i) {
	    ss << educts[i].second << " ";
	    ss << sp[educts[i].first].get_name();
	    
	    if(i+1!=educts.size())
	      ss << " + ";
	}
	
	
	ss << " <-- " << e_e << " --- " << activation << " --- " << e_p << " --> ";  
        
	for(size_t i=0; i<products.size(); ++i) {
	    ss << products[i].second << " ";
	    ss << sp[products[i].first].get_name();
	    
	    if(i+1!=products.size())
	      ss << " + ";
	}	
      
        return ss.str();
    }
    
    
        std::string get_string(const std::vector<species> &sp) const {
        std::stringstream ss;
	
	for(size_t i=0; i<educts.size(); ++i) {
	    ss << educts[i].second << " ";
	    ss << sp[educts[i].first].get_name();
	    
	    if(i+1!=educts.size())
	      ss << " + ";
	}
	
	if(reversible)
	    ss << " <--" << c << "-> ";  
        else
	    ss << " --" << c << "-> ";
	
	for(size_t i=0; i<products.size(); ++i) {
	    ss << products[i].second << " ";
	    ss << sp[products[i].first].get_name();
	    
	    if(i+1!=products.size())
	      ss << " + ";
	}	
      
        return ss.str();
    }
    
    bool has_species(size_t sp) {
        for(size_t i=0; i<educts.size(); ++i)
	    if(educts[i].first == sp)
	        return true;

	for(size_t i=0; i<products.size(); ++i)
	    if(products[i].first == sp)
	        return true;
      
        return false;  
    }
    
    void dump_species(std::vector<bool> &dump) {
        for(int i=educts.size()-1; i>=0; --i) { 
	    if(dump[educts[i].first])
	        educts.erase(educts.begin()+i);
	}
	
        for(int i=products.size()-1; i>=0; --i) 
	    if(dump[products[i].first])
	        products.erase(products.begin()+i);
    }
   
   
    void transform_species(std::map<size_t, size_t> &old_new) {
        for(size_t i=0; i<educts.size(); ++i)
	    educts[i].first = old_new[educts[i].first];
	
        for(size_t i=0; i<products.size(); ++i)
	    products[i].first = old_new[products[i].first];
    }
    
    bool operator== (reaction& rhs) {
        std::map<size_t, double> m1, m2;
	
	if(get_activation() != rhs.get_activation() ||
	   get_c() != rhs.get_c() ||
	   get_k() != rhs.get_k())
	    return false;
	    
	for(size_t i=0; i<educts.size(); ++i)
	    m1[educts[i].first] += educts[i].second;
	
	for(size_t i=0; i<rhs.educts.size(); ++i)
	    m2[rhs.educts[i].first] += rhs.educts[i].second;
	
	if(m1.size() != m2.size())
	    return false;
	
	for(std::map<size_t, double>::iterator it=m1.begin(); it != m1.end(); ++it)
	    if(m2[it->first] != it->second)
	        return false;
	
	
	m1.clear();
	m2.clear();
	
	for(size_t i=0; i<products.size(); ++i)
	    m1[products[i].first] += products[i].second;
	
	for(size_t i=0; i<rhs.educts.size(); ++i)
	    m2[rhs.products[i].first] += rhs.products[i].second;
	
	
	if(m1.size() != m2.size())
	    return false;
	
        for(std::map<size_t, double>::iterator it=m1.begin(); it != m1.end(); ++it)
	    if(m2[it->first] != it->second)
	        return false;
	
	
	
        return true;  
    }


    friend std::ostream& operator<<(std::ostream& os, const reaction& re);
    
    friend std::istream& operator>>(std::istream& is, reaction& re);
    
    void photodiss(double a, double b, double c) {}
    
    void sorder(double a, double b, double c, double d, double e)  {}
    
    void assrea(double a, double b, double c, double d) {}
    
    void eqrea(double a, double b, double c, double d) {}
};


std::ostream& operator<<(std::ostream& os, const reaction& re) {
    os << re.reversible << " " << re.c << " " << re.k << " "<< re.k_b << " " << re.activation << " " << re.educts.size() << " " << re.products.size();
    
    for(size_t i=0; i<re.educts.size(); ++i)
        os << " " << re.educts[i].first << " " << re.educts[i].second;
    
    
    for(size_t i=0; i<re.products.size(); ++i)
        os << " " << re.products[i].first << " " << re.products[i].second;
    
    return os;
} 


std::istream& operator>>(std::istream& is, reaction& re) {
  size_t educt_no, product_no;
  is >> re.reversible >> re.c >> re.k >> re.k_b >> re.activation >> educt_no >> product_no;
  
  for(size_t i=0; i<educt_no; ++i) {
      size_t s; double m;
      is >> s >> m;
      re.add_educt(s, m);    
  }
  
  
  for(size_t i=0; i<product_no; ++i) {
      size_t s; double m;
      is >> s >> m;
      re.add_product(s, m);    
  }
  
    //is >> cp.name >> cp.volume;
    
    return is;
} 



/*
 * Function for merging two reaction networks
 *
 * TODO add function for removing reaction duplicates
 *         -> done, but has to be tested
 */

void combine_r_networks(std::vector<species> &sp_1, std::vector<reaction> &re_1, 
			std::vector<species> &sp_2, std::vector<reaction> &re_2, 
			std::vector<species> &sp, std::vector<reaction> &re) {
    std::map<size_t, size_t> old_new;
  
    // Put species and reaction of first network into output network
    for(size_t i=0; i<sp_1.size(); ++i) 
        sp.push_back(sp_1[i]);
    
    for(size_t i=0; i<re_1.size(); ++i) 
        re.push_back(re_1[i]);
  
    // Add all species of the second reaction network. They are not added if there
    // is already one species of the first network with the same name...   
    for(size_t i=0; i<sp_2.size(); ++i) {
        bool found_dup=false;
        for(size_t j=0; j<sp.size() && !found_dup; ++j) {
	    if(sp[j].get_name() == sp_2[i].get_name()) {
	        std::cout << "found duplicat: " << sp[j].get_name() << std::endl;
		old_new[i] = j;
		found_dup=true;
	    }
	}
      
        // Not found a duplicate so we have to add
        if(!found_dup) {
	    old_new[i] = sp.size();
	    sp.push_back(sp_2[i]);
	}
    }
    
    // Now add reactions of second network
    for(size_t i=0; i<re_2.size(); ++i) {
        bool found_dup=false;
	
	for(size_t j=0; j<re.size(); ++j)
	    if(re[j] == re_2[i])
	        found_dup=true;
	
	if(!found_dup) {
            re.push_back(re_2[i]);
            re.back().transform_species(old_new);
	}
    }
}

// TODO test and comment
//

void filter_r_network(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		      std::vector<species> &sp, std::vector<reaction> &re,
		      std::vector<bool> &dump) {
    std::map<size_t, size_t> old_new;
    
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(!dump[i]) {
	    old_new[i] = sp.size();
	    sp.push_back(sp_in[i]);
	}
    }
    
    for(size_t i=0; i<re_in.size(); ++i) {	
        re.push_back(re_in[i]);
        
        re.back().dump_species(dump);
	re.back().transform_species(old_new);
    }
    
    for(int i=re.size()-1; i>=0; --i) {
        if(re[i].get_no_educt() == 0 &&
    	   re[i].get_no_product() == 0)
    	    re.erase(re.begin()+i);
    }
}

// TODO test + comment
//
void filter_r_network_dr(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		         std::vector<species> &sp_out, std::vector<reaction> &re_out,
		         std::vector<bool> &dump) {
    std::map<size_t, size_t> old_new;
    std::vector<species> sp;
    std::vector<reaction> re;
    
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(!dump[i]) {
	    old_new[i] = sp.size();
	    sp.push_back(sp_in[i]);
	}
    }
    
    for(size_t i=0; i<re_in.size(); ++i) {
        //
        bool all_sp_av=true;
        
        for(size_t j=0; j<dump.size(); ++j) {
	    if(dump[j]) {      
	        for(size_t k=0; k<re_in[i].get_no_educt_s(); ++k)
		    if(re_in[i].get_educt_id(k) == j)
		        all_sp_av = false;

	        for(size_t k=0; k<re_in[i].get_no_product_s(); ++k)
		    if(re_in[i].get_product_id(k) == j)
		        all_sp_av = false;
	    }  
	}
      
        if(all_sp_av) {
            re.push_back(re_in[i]);
            re.back().transform_species(old_new);
	}
    }
    
    for(int i=re.size(); i>=0; --i) {
        if(re[i].get_no_educt() == 0 &&
	   re[i].get_no_product() == 0)
	    re.erase(re.begin()+i);
    }
    
    
    std::vector<bool> dump_2;
    for(size_t i=0; i<sp.size(); ++i)
        dump_2.push_back(true);
    
    for(size_t i=0; i<re.size(); ++i) {
        for(size_t j=0; j<re[i].get_no_educt_s(); ++j)
	    dump_2[re[i].get_educt_id(j)] = false;
      
        for(size_t j=0; j<re[i].get_no_product_s(); ++j)
	    dump_2[re[i].get_product_id(j)] = false;
    }
    
    filter_r_network(sp, re, sp_out, re_out, dump_2);
}
      

      
/*
 * This macro function filters a network removing the species with name 'dump' 
 */
      
void filter_r_network_s(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		        std::vector<species> &sp, std::vector<reaction> &re,
		        std::string &dump) {
    // Which species to dump
    std::vector<bool> dump_f;
  
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(sp_in[i].get_name() == dump) {
	    dump_f.push_back(true);
	    std::cout << "Dumping species with id " << i << "!" << std::endl;    
	}
	else 
	    dump_f.push_back(false);
    }	
	
    // Filter out the species
    filter_r_network(sp_in, re_in, sp, re, dump_f);
}


/*
 * This macro function filters a network removing the species with name 'dump' 
 * TODO This function should dump all reactions containing species dump
 */
      
void filter_r_network_r(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		        std::vector<species> &sp, std::vector<reaction> &re,
		        std::string &dump) {
    // Which species to dump
    std::vector<bool> dump_f;
    size_t dump_s=0;
    
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(sp_in[i].get_name() == dump) {
	    dump_f.push_back(true);
	    std::cout << "Dumping species with id " << i << "!" << std::endl;
	    dump_s=i;
	}
	else 
	    dump_f.push_back(false);
    }	
    
    for(int i=re_in.size()-1; i>=0; --i) 
        if(re_in[i].has_species(dump_s)) 
	    re_in.erase(re_in.begin()+i);
	
    // Filter out the species
    filter_r_network(sp_in, re_in, sp, re, dump_f);
}




#endif // REACTION_NETWORK_H