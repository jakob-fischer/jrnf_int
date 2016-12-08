/* author: jakob fischer (jakob@automorph.info)
 * date: 9th February 2013
 * description: 
 * Simple class representing species and reactions of reaction networks.
 */

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
     * Constructor. Parameters:
     * i - id of species / position in array / in network
     * n - name string
     * c - kept constant in network / model
     * t - type, allways zero for now
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
     * Return temperature of compartement
     */
    
    double get_temperature() const {
        return temperature;  
    }
    
    
    /*
     * Set temperature of compartement
     */
    
    compartment& set_temperature(double t) {
        temperature=t;
	    return *this;
    }
    
    // Friend declaration of output function / operator
    friend std::ostream& operator<<(std::ostream& os, const compartment& cp);
    
    // Friend declaration of input function / operator
    friend std::istream& operator>>(std::istream& is, compartment& cp);
};


/*
 * Operator for writing compartement class to ostream
 */

std::ostream& operator<<(std::ostream& os, const compartment& cp) {
    os << cp.name << " " << cp.volume;
    
    return os;
} 


/*
 * Operator for reading compartement class from istream 
 */
 
std::istream& operator>>(std::istream& is, compartment& cp) {
    is >> cp.name >> cp.volume;
    
    return is;
} 



/*
 * The 'reaction' class represents one single reaction in the reaction network
 * The reaction is initialiced with reaction propensity 1 and no species as 
 * educts and products.
 * 
 * TODO This class does encapsulate bad because the ordering of reactants
 * inside in the class is visible to the outside. The methods that are
 * affected by this are marked as DEPRECATED.
 */

class reaction {
	// Reaction konstants, and activation energy
	// activation is defined as the minimum a. e. of forward and backward r.
    double c, k, k_b, activation;
    
    // True for reversible reactions
    bool reversible;
    
    // Vector containing educts and products of the network
    std::vector< std::pair<size_t,double> > educts, products;
  
  public:
  
    /*
     * Standard c. Constants are set 1.0, a. e. 0, and reaction irrevers.
     */
      
    reaction() {
       c=1.0;
       k=1.0;
       k_b=1.0;
       activation=0.0;
       reversible=false;
    }
  

    /*
     * Is reversible?
     */  
   
    bool is_reversible() const {
        return reversible;  
    }
   
   
    /*
     * Set reversibility of reaction
     */
   
    reaction& set_reversible(bool rev) {
        reversible = rev;
	return *this;      
    }
   
    
    /*
     * Get constant `c` (was used in stoch. simulation)
     */
   
    double get_c() const {
        return c;
    }
   
    
    /*
     * Set constant `c` (was used in stoch. simulation)
     */
   
    reaction& set_c(double c_) {
	    c=c_;
	    return *this;
    }
    
    /*
        Outdated - TODO replace through similar set_k_by_T?
    reaction& set_c_by_T(const std::vector<species>& sp, double T) {
      	double E_educt=0;
	
	// and then translate all educt and product ids
        for(size_t i=0; i<educts.size(); ++i) 
	    E_educt += sp[educts[i].first].get_energy()*educts[i].second;
	
	c = exp(-(activation-E_educt)/T);
      
      
        return *this;  
    }*/
    
    
    /*
     * Get forward reaction constant 
     */
     
    double get_k() const {
        return k;
    }
   
   
    /*
     * Set forward reaction constant
     */
   
    reaction& set_k(double k_) {
	    k=k_;
	    return *this;
    }
    
    
    /*
     * Get backward reaction constant
     */
    
    double get_k_b() const {
        return k_b;  
    }
    
    
    /*
     * Set backward reaction constant
     */
    
    reaction& set_k_b(double k_b_) {
        k_b=k_b_;  
        return *this;
    }
    
    
    /*
     * Get activation energy
     */    
    
    double get_activation() const {
        return activation;
    }
    
    
    /*
     * Set activation energy
     */
   
    reaction& set_activation(double activation_) {
	    activation=activation_;
	    return *this;
    }
    
    
    /*
     * Returns the size of the educts vector in this reaction 
     * DEPRECATED
     */   
  
    size_t get_no_educt_s() const {
	    return educts.size();
    }
    
    
    /*
     * Returns the number of different educts in the reaction 
     */
    
    size_t get_no_educt() const {
        size_t a=0;
	
	    for(size_t i=0; i<educts.size(); ++i)
	        a += educts[i].second;
	
	    return a;
    }
    
    
    /*
     * Get id of ith educt in educts vector in this reaction 
     * DEPRECATED
     */
    
    size_t get_educt_id(size_t i) const {  return educts[i].first;  }
    
    
    /*
     * Set id of ith educt in educts vector in this reaction 
     * DEPRECATED
     */
    
    size_t set_educt_id(size_t i, size_t n) {  educts[i].first = n;  }


    /*
     * Get stoichiometric coefficient of ith educt in educts vector  
     * DEPRECATED
     */
    
    double get_educt_mul(size_t i) const {  return educts[i].second;  }


    /*
     * Set stoichiometric coefficient of ith educt in educts vector  
     * DEPRECATED
     */
        
    reaction& set_educt_mul(size_t i, double m) {  educts[i].second=m;  return *this;  }


    /*
     * Set stoichiometric coefficient of ith product in products vector  
     * DEPRECATED
     */

    reaction& set_product_mul(size_t i, double m) {  products[i].second=m;  return *this;  }
    
    
    /*
     * Add the species `id` as educt with the stoichiometric
     * coficcent `m` to the reaction. 
     */    
     
    reaction& add_educt_s(size_t id, double m=1.0) {
       educts.push_back( std::pair<size_t,size_t>(id, m) );
       
       return *this;
    }


    /*
     * Add the species `id` as educt with the stoichiometric
     * coficcent `m` to the reaction. If `id` is already a educt its
     * stoichiometric coefficient is increased by `m`
     */
    
    reaction& add_educt(size_t id, double m=1.0) {
       for(size_t i=0; i<educts.size(); ++i) 
	   if(educts[i].first == id) {
	      educts[i].second += m;
	      return *this;
	   }
	   
        educts.push_back( std::pair<size_t,size_t>(id, m) );
       
       return *this;
    }


    /*
     * Returns the size of the products vector in this reaction 
     * DEPRECATED
     */
     
    size_t get_no_product_s() const {
	    return products.size();
    }
    
    /*
     * Returns the number of different products in the reaction 
     */
     
    size_t get_no_product() const {
        size_t a=0;
	
	    for(size_t i=0; i<products.size(); ++i)
	        a += products[i].second;
	
	    return a;
    }
    
    
    /*
     * Get id of ith product in products vector in this reaction 
     * DEPRECATED
     */
    
    size_t get_product_id(size_t i) const {  return products[i].first;  }
    
    /*
     * Set id of ith product in products vector in this reaction 
     * DEPRECATED
     */
    
    size_t set_product_id(size_t i, size_t n) {  products[i].first = n;  }
        
        
    /*
     * Get stoichiometric coefficient of ith product in products vector 
     * in this reaction 
     * DEPRECATED
     */

    double get_product_mul(size_t i) const {  return products[i].second;  }


    /*
     * Add the species `id` as product with the stoichiometric
     * coficcent `m` to the reaction. 
     */    
        
    reaction& add_product_s(size_t id, double m=1.0) {
	    products.push_back( std::pair<size_t,double>(id, m) );
       
        return *this;
    }
    
    
    /*
     * Add the species `id` as product with the stoichiometric
     * coficcent `m` to the reaction. If `id` is already a product its
     * stoichiometric coefficient is increased by `m`
     */
    
    reaction& add_product(size_t id, double m=1.0) {
        for(size_t i=0; i<products.size(); ++i) 
	        if(products[i].first == id) {
	            products[i].second += m;
	            return *this;
	        }        
	   
	    products.push_back( std::pair<size_t,double>(id, m) );
       
       return *this;
    }


    /*
     * Calculate forward reaction rate for this reaction
     */
  
    double calculate_rate_f(double Ea, const std::vector<double>& species_E, bool logscale) {
        if(logscale) {
            double Ee=1;

            for(size_t i=0; i<educts.size(); ++i)
                Ee *= pow(species_E[educts[i].first], educts[i].second);           

            return Ee/Ea;
        } else { 
            double Ee=0;

            for(size_t i=0; i<educts.size(); ++i)
                Ee += species_E[educts[i].first]*educts[i].second;

            return exp(-(Ea-Ee));
        }
    }


    /*
     * Calculate backward reaction rate for this reaction
     */

    double calculate_rate_b(double Ea, const std::vector<double>& species_E, bool logscale) {
        if(logscale) {
            double Ep=1;

            for(size_t i=0; i<products.size(); ++i)
                Ep *= pow(species_E[products[i].first], products[i].second);           

            return Ep/Ea;
        } else { 
            double Ep=0;

            for(size_t i=0; i<products.size(); ++i)
                Ep += species_E[products[i].first]*products[i].second;

            return exp(-(Ea-Ep));
        }
    }


    /*
     * Calculate the effective rate of the reaction from forward and backward reaction rate
     */

     double calculate_rate(double Ea, const std::vector<double>& species_E, bool logscale) {
         return (calculate_rate_f(Ea, species_E, logscale) - calculate_rate_b(Ea, species_E, logscale));
     }
    
    
    /*
     * Converts the reaction into a C++ string using the ids of species
     * as their names
     */
    
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


    /*
     * Returns the reaction as string representation using the name of 
     * the species in `sp`. Works only for reversible reactions for which
     * it also print the educts, products and activation-energy.
     */   
   
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
    
    
    /*
     * Returns the reaction as string representation using the name of 
     * the species in `sp`
     */
    
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
    
    
    /*
     * Checks if reaction has species with id `sp`
     */
    
    bool has_species(size_t sp) {
        for(size_t i=0; i<educts.size(); ++i)
	        if(educts[i].first == sp)
	            return true;

	    for(size_t i=0; i<products.size(); ++i)
	        if(products[i].first == sp)
	            return true;
      
        return false;  
    }
    
    
    /*
     * Dumps the species hilighted by `dump` in this reaction 
     */
     
    void dump_species(std::vector<bool> &dump) {
        for(int i=educts.size()-1; i>=0; --i) { 
	        if(dump[educts[i].first])
	            educts.erase(educts.begin()+i);
	    }
	
        for(int i=products.size()-1; i>=0; --i) 
	        if(dump[products[i].first])
	            products.erase(products.begin()+i);
    }
   
   
    /*
     * Transforms all species of this reaction my the map given
     */
   
    void transform_species(std::map<size_t, size_t> &old_new) {
        for(size_t i=0; i<educts.size(); ++i)
	        educts[i].first = old_new[educts[i].first];
	
        for(size_t i=0; i<products.size(); ++i)
	        products[i].first = old_new[products[i].first];
    }
    
    
    /*
     * Operator for comparing if two reactions are equal
     */
    
    bool operator== (reaction& rhs) {
        std::map<size_t, double> m1, m2;
	
	    // First compare reaction konstants
	    if(get_activation() != rhs.get_activation() ||
	       get_c() != rhs.get_c() ||
	       get_k() != rhs.get_k() ||
	       get_k_b() != rhs.get_k_b())
	        return false;
	    
	    // Put educts into map
	    for(size_t i=0; i<educts.size(); ++i)
	        m1[educts[i].first] += educts[i].second;
	
	    for(size_t i=0; i<rhs.educts.size(); ++i)
	        m2[rhs.educts[i].first] += rhs.educts[i].second;
	
	    // Compare maps size
	    if(m1.size() != m2.size())
	        return false;
	
	    // Compare maps
	    for(std::map<size_t, double>::iterator it=m1.begin(); it != m1.end(); ++it)
	        if(m2[it->first] != it->second)
	            return false;
	
	    // Clear maps and redo for educts
	    m1.clear();
	    m2.clear();
	
	    for(size_t i=0; i<products.size(); ++i)
	        m1[products[i].first] += products[i].second;
	
	    for(size_t i=0; i<rhs.educts.size(); ++i)
	        m2[rhs.products[i].first] += rhs.products[i].second;
		
		// Size
	    if(m1.size() != m2.size())
	        return false;
	
	    // Contents
        for(std::map<size_t, double>::iterator it=m1.begin(); it != m1.end(); ++it)
	        if(m2[it->first] != it->second)
	            return false;
	
	    // If no mismatch was found reactions are equal
        return true;  
    }


    // Friend declaration of output function / operator
    friend std::ostream& operator<<(std::ostream& os, const reaction& re);
    
    // Friend declaration of input function / operator
    friend std::istream& operator>>(std::istream& is, reaction& re);
};


/*
 * Operator for writing to ostream
 */

std::ostream& operator<<(std::ostream& os, const reaction& re) {
    os << re.reversible << " " << re.c << " " << re.k << " "<< re.k_b << " " << re.activation << " " << re.educts.size() << " " << re.products.size();
    
    for(size_t i=0; i<re.educts.size(); ++i)
        os << " " << re.educts[i].first << " " << re.educts[i].second;
    
    
    for(size_t i=0; i<re.products.size(); ++i)
        os << " " << re.products[i].first << " " << re.products[i].second;
    
    return os;
} 


/*
 * Operator for reading a reaction from istream
 */

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
  
    return is;
} 



/*
 * Function for merging two reaction networks
 *
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



/*
 * Copy the reaction network while removing all species that are hilighted
 * in the `dump` vector. 
 */
 
void filter_r_network(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		      std::vector<species> &sp, std::vector<reaction> &re,
		      std::vector<bool> &dump) {
    std::map<size_t, size_t> old_new;
    
    // Copy the species that are kept
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(!dump[i]) {
	        old_new[i] = sp.size();
	        sp.push_back(sp_in[i]);
	    }
    }
    
    // Copy reactions while removing the species to dump
    for(size_t i=0; i<re_in.size(); ++i) {	
        re.push_back(re_in[i]);
        
        re.back().dump_species(dump);
	    re.back().transform_species(old_new);
    }
    
    // Erase reactions that are empty
    for(int i=re.size()-1; i>=0; --i) {
        if(re[i].get_no_educt() == 0 &&
    	   re[i].get_no_product() == 0)
    	    re.erase(re.begin()+i);
    }
}


/*
 * Function copys a reaction network while removing all reactions that
 * contain species hilighted in dump.
 */
 
void filter_r_network_dr(std::vector<species> &sp_in, std::vector<reaction> &re_in,  
		         std::vector<species> &sp_out, std::vector<reaction> &re_out,
		         std::vector<bool> &dump) {
    std::map<size_t, size_t> old_new;
    
    // Copy those species that are not dumped
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(!dump[i]) {
	        old_new[i] = sp_out.size();
	        sp_out.push_back(sp_in[i]);
	    }
    }
    
    // Tests for every reaction if all species are still available after
    // transformation / not dumped
    for(size_t i=0; i<re_in.size(); ++i) {
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
            re_out.push_back(re_in[i]);
            re_out.back().transform_species(old_new);
	    }
    }
    
    // Erase those reactions that are empty 
    for(int i=re_out.size(); i>=0; --i) {
        if(re_out[i].get_no_educt() == 0 &&
	       re_out[i].get_no_product() == 0)
	        re_out.erase(re_out.begin()+i);
    }
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
	    } else 
	        dump_f.push_back(false);
    }	
	
    // Filter out the species
    filter_r_network(sp_in, re_in, sp, re, dump_f);
}


/*
 * This macro function filters a network removing the species with name 'dump'.
 * This function does not delete the species from the reactions containing
 * it but deletes all those reactions all together.
 */
      
void filter_r_network_r(std::vector<species> &sp_in, std::vector<reaction> &re_in_,  
		        std::vector<species> &sp, std::vector<reaction> &re,
		        std::string &dump) {
	std::vector<reaction> re_in(re_in_);
					
    // Which species to dump
    std::vector<bool> dump_f;
    size_t dump_s=0;
    
    for(size_t i=0; i<sp_in.size(); ++i) {
        if(sp_in[i].get_name() == dump) {
	        dump_f.push_back(true);
	        std::cout << "Dumping species with id " << i << "!" << std::endl;
	        dump_s=i;
	    } else
	        dump_f.push_back(false);
    }	
    
    // Delete all reactions containing the species
    for(int i=re_in.size()-1; i>=0; --i) 
        if(re_in[i].has_species(dump_s)) 
	        re_in.erase(re_in.begin()+i);
	
    // Filter out the species
    filter_r_network(sp_in, re_in, sp, re, dump_f);
}




#endif // REACTION_NETWORK_H
