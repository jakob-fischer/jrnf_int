/* author: jakob fischer (jakob@automorph.info)
 * date: 9th February 2013
 * description: 
 * Contains methods for writing reaction networks to files and to the 
 * shell.
 */

#ifndef REACTION_NETWORK_FILEOP_H
#define REACTION_NETWORK_FILEOP_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


/*
 * Writes the reaction network given by the species 'sp' and the reaction
 * 're' in an textual representation of reactions to the file named
 * 'filename'
 */

void write_text_reaction_n_b(const std::string& filename, 
                            const std::vector<species>& sp,
			                const std::vector<reaction>& re) {
    std::ofstream out(filename.c_str());
    
    for(size_t i=0; i<re.size(); ++i)
        out << re[i].get_string_b(sp) << std::endl;
}


/*
 * Writes a reaction network given by species-vector 'sp' and reaction-
 * vector 're' to the file with filename 'filename' in jrnf format.
 */

int write_jrnf_reaction_n(const std::string& filename, const std::vector<species>& sp,
			  const std::vector<reaction>& re) {
    std::ofstream out(filename.c_str());
    
    // Writing header
    out << "jrnf0003" << std::endl;
    out << sp.size() << " " << re.size() << std::endl;
    
    // Writing the species
    for(size_t i=0; i<sp.size(); ++i)
        out << sp[i] << std::endl;
    
    // and the reactions
    for(size_t j=0; j<re.size(); ++j)
        out << re[j] << std::endl;

    return 0;
}


/*
 * Reads the reaction network in 'filename' to 'sp' and 're'.
 */

int read_jrnf_reaction_n(const std::string& filename, std::vector<species>& sp,
			  std::vector<reaction>& re) {
    std::ifstream in(filename.c_str());
    
    // reading and validating header
    std::string header;
    in >> header;
    
    if(header != "jrnf0003") {
        std::cout << "File is not in jrnf-format!" << std::endl;
        return 1;  
    }
    
    // Species number and reaction number are in second row
    size_t sp_no, re_no;
    in >> sp_no >> re_no;
    
    // Reading species
    for(size_t i=0; i<sp_no; ++i) {
        species s(0,"");
	    in >> s;
	    sp.push_back(s);      
    }
    
    // Reading reactions
    for(size_t j=0; j<re_no; ++j) {
        reaction r;
	    in >> r;
	    re.push_back(r);
    }
    
    return 0;
}


/*
 * Function for writing the reaction network 'sp', 're' to the file
 * 'filename' in sbml-format. Not all (thermo)dynamic information is
 * maintained. Thus the file is only usable for topologic analysis.
 */

int write_sbml_reaction_n(const std::string& filename, 
                          const std::vector<species>& sp,
			              const std::vector<reaction>& re) {
    std::ofstream out(filename.c_str());
   
    // Writing header     
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    out << "<sbml xmlns=\"http://www.sbml.org/sbml/level2\" xmlns:html=\"http://www.w3.org/1999/xhtml\" ";
    out << "xmlns:jigcell=\"http://www.sbml.org/2001/ns/jigcell\" xmlns:math=";
    out << "\"http://www.w3.org/1998/Math/MathML\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#";
    out << "\" xmlns:sbml=\"http://www.sbml.org/sbml/level2\" xmlns:xlink=\"http://www.w3.org/1999/xlink";
    out << "\" level=\"2\" version=\"1\">" << std::endl << std::endl;
    
    // 'Standard' information on model
    out << "<model>" << std::endl;
    out << "  <listOfCompartments>" << std::endl;
    out << "    <compartment id=\"cell_a\" name=\"cell_a\" />" << std::endl;
    out << "  </listOfCompartments>" << std::endl;
    
    // Writing list of species
    out << "  <listOfSpecies>" << std::endl;     
    for(size_t i=0; i<sp.size(); ++i) {
	    if(sp[i].get_type() == 0) {
	        out << "    <species id=\"" << sp[i].get_name();
	        out << "\" name=\"" << sp[i].get_name();
	        out << "\" compartment=\"cell_a\" />" << std::endl;
	    }      
    }
     
    out << "  </listOfSpecies>" << std::endl;  
    out << std::endl;
    
    // Writing list of reactions
    out << "  <listOfReactions>" << std::endl;

    for(size_t i=0; i<re.size(); ++i) {
        out << "      <reaction id=\"R" << i+1;
	    if(re[i].is_reversible())
	        out << "\" reversible=\"true\">" << std::endl;  
	    else
	        out << "\" reversible=\"false\">" << std::endl;
	    out << "        <listOfReactants>" << std::endl;
	
	    for(size_t j=0; j<re[i].get_no_educt_s(); ++j) {
	        if(sp[re[i].get_educt_id(j)].get_type() == 0) {
	            out << "          <speciesReference species=\"" << sp[re[i].get_educt_id(j)].get_name();
		       out << "\" stoichiometry=\"" << re[i].get_educt_mul(j) << "\" />" << std::endl;  
	        }
	    }
	
	    out << "        </listOfReactants>" << std::endl;
        out << "        <listOfProducts>" << std::endl;
	
        for(size_t j=0; j<re[i].get_no_product_s(); ++j) {
	        if(sp[re[i].get_product_id(j)].get_type() == 0) {
	            out << "          <speciesReference species=\"" << sp[re[i].get_product_id(j)].get_name();
		        out << "\" stoichiometry=\"" << re[i].get_product_mul(j) << "\" />" << std::endl;  
	        }
        }
	
	    out << "        </listOfProducts>" << std::endl;
	    out << std::endl;
        out << "        <kineticLaw>" << std::endl;
        out << "          <listOfParameters>" << std::endl;
        out << "            <parameter id=\"c\"  value=\"" << re[i].get_c() << "\" />" << std::endl;
	    out << "          </listOfParameters>" << std::endl;
        out << "	</kineticLaw>" << std::endl;
        out << "      </reaction>" << std::endl;
      
      
    }
    
    out << "  </listOfReactions>" << std::endl;
    
    // Closing tags
    out << "</model>" << std::endl;
    out << "</sbml>" << std::endl;
    
    return 0;
}



#endif // REACTION_NETWORK_FILEOP_H
