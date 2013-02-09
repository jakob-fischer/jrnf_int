#ifndef CL_PARA_H
#define CL_PARA_H

#include <string>
#include <sstream>
#include <map>
#include <iostream>


class cl_para {
    std::map<std::string, std::string> par;
  
public:
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
    
    
    bool have_param(const std::string& str) {
        if(par[str].length() > 0) 
	    return true;
	else
	    return false;
    }
    
    
    double get_param_d(const std::string& str) {
        double d;
      
	std::istringstream ss(par[str]);
	ss >> d;
      
        return d;  
    }
    
    
    std::string get_param(const std::string& str) {
        return par[str];  
    }
    
    
    int get_param_i(const std::string& str) {
        int i;
      
	std::istringstream ss(par[str]);
	ss >> i;
      
        return i;  
    }
  
  
};





#endif // CL_PARA_H