/* author: jakob fischer (jakob@automorph.info)
 * date: 9th February 2013
 * description: 
 * Includes class for handling C++ console parameters.
 */

#ifndef CL_PARA_H
#define CL_PARA_H

#include <string>
#include <sstream>
#include <map>
#include <iostream>


/* 
 * Class for handling C++ console parameters.
 */

class cl_para {
	// map that asocciates console parameters with its values
	// parameters given without value get the value "empty"
    std::map<std::string, std::string> par;
  
public:

    /* Constructor gets the parameters exactly as given to the C(++)
     * main function.
     */

    cl_para(int argc, const char* argv[]) {
        for(size_t i=1; i<argc; ++i) { 
	        std::string str(argv[i]);
	        size_t fnd=str.find('=');
	    
	        if(fnd == str.npos) 
	            par[str] = std::string("empty");  
	        else 
	            par[str.substr(0, fnd)] = str.substr(fnd+1); 
	    }
    }
    
    
    /*
     * Is parameter str in the parameter list?
     */
    
    bool have_param(const std::string& str) {
        if(par[str].length() > 0) 
	        return true;
	    else
	        return false;
    }
    
    
    /*
     * Tries to convert the value of parameter `str` to double and
     * returns the result.
     */
    
    double get_param_d(const std::string& str) {
        double d;
      
	    std::istringstream ss(par[str]);
	    ss >> d;
      
        return d;  
    }
    
    
    /*
     * Returns the value of parameter `str` as string.
     */
    
    std::string get_param(const std::string& str) {
        return par[str];  
    }
    

    /*
     * Tries to convert the value of parameter `str` as integer and
     * returns the result.
     */
    
    int get_param_i(const std::string& str) {
        int i;
      
	    std::istringstream ss(par[str]);
	    ss >> i;
      
        return i;  
    }

};





#endif // CL_PARA_H
