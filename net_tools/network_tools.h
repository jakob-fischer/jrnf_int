#ifndef NETWORK_TOOLS_H
#define NETWORK_TOOLS_H

#include <vector>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <iml.h>
#include <iostream>
#include <numeric>
#include <algorithm>

/*
 * Creates a modular hierarchic network following the model of Pan and Sinha
 * the edges are inserted with random / stochastic probability
 * 
 * TODO To be implemented
 */

void build_pan_sinha_stoch(std::vector< std::pair<size_t, size_t> >& edges,
		     size_t h, size_t n, size_t m, double r, double rho_0, bool self_loop=false, bool directed=false) { 
}


/*
 * Helper for "build_pan_sinha".
 * Calculates the size (number of nodes) of the module on level 'level' with the
 * index 'element. It needs a vector 'elem_c' with sizes of elementary modules (level 0)
 * and the parameter 'm' 
 */

size_t bps_module_size(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    if(level == 0) 
        return elem_c[element];
    
    if(level == 1) {
        size_t sum=0;
	
	for(size_t i=0; i<m; ++i) 
	    sum += elem_c[element*m+i];
	
	return sum;
    }
    
    return bps_module_size(elem_c, m, level-1, 2*element) +
           bps_module_size(elem_c, m, level-1, 2*element+1);
  
}


/*
 * Similar to previous function, but calculates the offset. 
 * (between first element of the network and first one in the module)
 */

size_t bps_module_offset(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    if(level == 0) {
        size_t sum=0;
	
        for(size_t i=0; i<element; ++i)
            sum += elem_c[i];
	  
        return sum;
    }
	
    if(level == 1)
	return bps_module_offset(elem_c, m, 0, m*element);
    
    return bps_module_offset(elem_c, m, level-1, 2*element);
}


/*
 * Helper for 'create_pan_sinha' that checks if two nodes 'n1' and 'n2' are in 
 * the same level-0 module.
 */

bool bps_same_emodule(const std::vector<size_t> & elem_c, size_t n1, size_t n2) {
    // ensure n1 <= n2
    
    size_t first=0;
    size_t sum=0;
    while(sum + elem_c[first] < n1) 
        sum += elem_c[first++];  

    size_t second=0;
    sum=0;
    while(sum + elem_c[second] < n2) 
        sum += elem_c[second++];  
 
    return first==second;  
}

/*
 * Calculates the maximal number of edges for networks with 'nodes' nodes.
 * Depending on if 'self_loop's are allowed and if the network is 'directed'.
 * Multiple edges are not taken into account.
 */

size_t bps_max_edges(size_t nodes, bool self_loop, bool directed) {
    if(self_loop && directed)
        return (nodes*nodes);
    else if(self_loop && !directed)
        return ((nodes*(nodes-1))/2+nodes);
    else if(!self_loop && directed)
        return nodes*nodes-nodes;
    else if(!self_loop && !directed)
        return (nodes*(nodes-1))/2;
}



/*
 * The function creates a hierarchical-modular network using the algorithm of pan-sinha
 * 
 * 
 */

void create_pan_sinha(std::vector< std::pair<size_t, size_t> >& edges, size_t N, size_t M,
		     size_t h, size_t m, double r, bool multiple =false, bool self_loop=false, bool directed=false) {
    size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor
    {                                // and adding 1 until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	    N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }
    

    // Calculate maximum number of edges for each level
    std::vector<double> max_edges_level, max_edges_level_r;
    std::vector<size_t> possible_base;    
    
    {   // lowest level
        max_edges_level.push_back(0);
        for(size_t k=0; k<n_el_mod; ++k) {
            size_t nodes=bps_module_size(elem_c, m, 0, k);
            max_edges_level[0] += bps_max_edges(nodes, self_loop, directed);
        }
    }     

    {   // medium level
        max_edges_level.push_back(0);
        for(size_t k=0; k<n_el_mod/m; ++k) {
            size_t nodes=bps_module_size(elem_c, m, 1, k);
            size_t d=bps_max_edges(nodes, self_loop, directed);

            for(size_t j=0; j<m; ++j)
                d -= bps_max_edges(elem_c[m*k+j],self_loop, directed);

            max_edges_level[1] += d;
        }
    }

    // the upper levels
    for(size_t i=1; i<=h; ++i) {
        max_edges_level.push_back(0);
        for(size_t k=0; k<pow(2,h-i); ++k) {
            size_t nodes1=bps_module_size(elem_c, m, i, k*2);
            size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);

            if(directed)
                max_edges_level[i+1] += 2*nodes1*nodes2;
            else 
                max_edges_level[i+1] += nodes1*nodes2;
        }
    }

    // rescale in terms of r
    for(size_t i=0; i<max_edges_level.size(); ++i)
        max_edges_level_r.push_back(max_edges_level[i]*pow(r, i));
    

    // Distribute M edges on all levels according to max_edges_level-distribution
    double sum_mel=std::accumulate(max_edges_level_r.begin(),max_edges_level_r.end(),0);
    size_t M_tmp=M;
    std::vector<size_t> edges_level_partition;
    std::vector<double> edges_level_partition_dif;

    for(size_t i=0; i<max_edges_level.size(); ++i) {
        size_t tmp=floor(M*max_edges_level_r[i]/sum_mel);    
        M_tmp -= tmp;
        edges_level_partition.push_back(tmp);
        edges_level_partition_dif.push_back(M*max_edges_level_r[i]/sum_mel-tmp);
    }

    while(M_tmp > 0) {
        std::vector<double>::iterator it= 
              std::max_element(edges_level_partition_dif.begin(), 
                               edges_level_partition_dif.end());
        size_t next=it-edges_level_partition_dif.begin();
        ++edges_level_partition[next];
        edges_level_partition_dif[next]=0;
        --M_tmp;
    }
       
    std::cout << "Show edge -> level assignment:" << std::endl;
    for(size_t i=0; i<edges_level_partition.size(); ++i)
        std::cout << i << "  :  " << edges_level_partition[i] << "   (" << max_edges_level[i]/pow(r,i) << ")" << std::endl;

    
    // Now generate the edges
    
    {   // lowest level - distribute 
        std::vector<size_t> dist_low_level_edges;
	std::vector<double> dist_low_level_edges_dif;
	size_t elp_tmp=edges_level_partition[0];
	
	for(size_t i=0; i<n_el_mod; ++i) {
	    //std::cout << "elem_c[" << i << "]=" << elem_c[i] << "     bpsm=" << bps_max_edges(elem_c[i], multiple, self_loop, directed) <<  std::endl;
	  
	  
	    double rel=double(bps_max_edges(elem_c[i], self_loop, directed)/max_edges_level[0]*edges_level_partition[0]);
	    dist_low_level_edges.push_back(floor(rel));
	    dist_low_level_edges_dif.push_back(rel-floor(rel));
	    elp_tmp -= floor(rel);
	}
	
	while(elp_tmp > 0) {
            std::vector<double>::iterator it= 
                  std::max_element(dist_low_level_edges_dif.begin(), 
                                   dist_low_level_edges_dif.end());
            size_t next=it-dist_low_level_edges_dif.begin();
            ++dist_low_level_edges[next];
            dist_low_level_edges_dif[next]=0;
            --elp_tmp;
	}
	
	//std::cout << "elp_tmp=" << elp_tmp << std::endl;
	
	//std::cout << "LEVEL ZERO: ";
	//for(size_t i=0; i<dist_low_level_edges.size(); ++i) 
	//    std::cout << " " << dist_low_level_edges[i] << " ";
	//std::cout << std::endl;
	 
	//  DO IT
	for(size_t i=0; i<n_el_mod; ++i) {
	    size_t offset = bps_module_offset(elem_c, m, 0, i);
	    
	    //std::cout << "i=" << i << std::endl;
	    //    std::cout << "dist_low_level_edges[" << i << "] = " << dist_low_level_edges[i] << std::endl;
	    
	    while(dist_low_level_edges[i] > 0) {
	        size_t n1=offset+rand()%elem_c[i];
	        size_t n2=offset+rand()%elem_c[i];
	        
		if(n1 != n2 || self_loop) {
		    if(directed) {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_low_level_edges[i];
		        }	 
		    } else {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() &&
			   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_low_level_edges[i];
		        }
		    } 
		}
	    }
	}
    }      
    
    //std::cout << "DONE" << std::endl;
    
    {   // medium level - the most complicated ;)
        std::vector<size_t> dist_med_level_edges;
	std::vector<double> dist_med_level_edges_dif;
	size_t elp_tmp=edges_level_partition[1];
	
	for(size_t k=0; k<n_el_mod/m; ++k) {
	    size_t nodes=bps_module_size(elem_c, m, 1, k);
            size_t d=bps_max_edges(nodes, self_loop, directed);

            for(size_t j=0; j<m; ++j)
                d -= bps_max_edges(elem_c[m*k+j] ,self_loop, directed);

	    double rel=double(d/max_edges_level[1]*edges_level_partition[1]);
	    
	    dist_med_level_edges.push_back(floor(rel));
 	    dist_med_level_edges_dif.push_back(rel-floor(rel));
	    elp_tmp -= floor(rel);
	}
	
	while(elp_tmp > 0) {
            std::vector<double>::iterator it= 
                  std::max_element(dist_med_level_edges_dif.begin(), 
                                   dist_med_level_edges_dif.end());
            size_t next=it-dist_med_level_edges_dif.begin();
            ++dist_med_level_edges[next];
            dist_med_level_edges_dif[next]=0;
            --elp_tmp;
	}
	
	//std::cout << "elp_tmp=" << elp_tmp << std::endl;
	
	//std::cout << "LEVEL ONE: ";
	//for(size_t i=0; i<dist_med_level_edges.size(); ++i) 
	//    std::cout << " " << dist_med_level_edges[i] << " ";
	//std::cout << std::endl;
	 
	// DO IT
	 
        for(size_t k=0; k<n_el_mod/m; ++k) {
	    size_t offset = bps_module_offset(elem_c, m, 1, k);	
	    size_t nodes=bps_module_size(elem_c, m, 1, k);
	    
	    while(dist_med_level_edges[k] > 0) {
	        size_t n1=offset+rand()%nodes;  
	        size_t n2=offset+rand()%nodes; 
	      
		if(!bps_same_emodule(elem_c, n1, n2)) {
                    if(directed) {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_med_level_edges[k];
		        }	 
		    } else {
	                if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() ||
			   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_med_level_edges[k];
		        }
		    } 		  
		  
		  
		  
		}
	    }
	}
    }
    
    {    // connect the higher levels 
        for(size_t i=1; i<=h; ++i) {
            std::vector<size_t> dist_x_level_edges;
	    std::vector<double> dist_x_level_edges_dif;
	    size_t elp_tmp = edges_level_partition[i+1];
	  
            for(size_t k=0; k<pow(2,h-i); ++k) {
                size_t nodes1=bps_module_size(elem_c, m, i, k*2);
                size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);
                size_t pot_edges=nodes1*nodes2;
                if(directed)
                    pot_edges *= 2;
	
		double rel=double(pot_edges/max_edges_level[i+1]*edges_level_partition[i+1]);
	        dist_x_level_edges.push_back(floor(rel));
	        dist_x_level_edges_dif.push_back(rel-floor(rel));
	        elp_tmp -= floor(rel);
            }
             
            while(elp_tmp > 0) {
                std::vector<double>::iterator it= 
                      std::max_element(dist_x_level_edges_dif.begin(), 
                                       dist_x_level_edges_dif.end());
                size_t next=it-dist_x_level_edges_dif.begin();
                ++dist_x_level_edges[next];
                dist_x_level_edges_dif[next]=0;
                --elp_tmp;
	    }
	    
	   // std::cout << "LEVEL " << i+1 << " : ";
	   // for(size_t k=0; k<dist_x_level_edges.size(); ++k) 
	   //     std::cout << " " << dist_x_level_edges[k] << " ";
	   // std::cout << std::endl;
	    
	    
	    // NOW we are doing the connecting
            for(size_t k=0; k<pow(2,h-i); ++k) {
                size_t nodes1=bps_module_size(elem_c, m, i, k*2);
                size_t nodes2=bps_module_size(elem_c, m, i, k*2+1);
		size_t offset=bps_module_offset(elem_c, m, i+1, k);
		
		while(dist_x_level_edges[k] > 0) {
		    size_t n1=offset+rand()%nodes1;
		    size_t n2=offset+nodes1+rand()%nodes2;
		    
		    if(directed) {
		        if(rand()%2) {
			    if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		                --dist_x_level_edges[k];
		            }	  
			} else {
			    if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n2, n1)) == edges.end() || multiple) {
		                edges.push_back(std::pair<size_t, size_t> (n2, n1) );
		                --dist_x_level_edges[k];
		            }	
			}
		    }else {
		        if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(n1, n2)) == edges.end() ||
			   multiple) {
		            edges.push_back(std::pair<size_t, size_t> (n1, n2) );
		            --dist_x_level_edges[k];
		        }	  
		    }
		  
		}
            }
        }    
    }
}




void create_barabasi_albert_old(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t m, size_t m0, 
			    bool allow_multiple=false) {
    size_t no_nodes=0;	    
    std::vector<size_t> degree;
    std::vector<bool> available;
    for(size_t i=0; i<N; ++i) {
        degree.push_back(0);
        available.push_back(true);
    }
	    
    for(; no_nodes<m0; ++no_nodes) {
	edges.push_back(std::pair<size_t, size_t>(no_nodes, no_nodes+1));
	++degree[no_nodes];
	++degree[no_nodes+1];
    }
    ++no_nodes;
	        
    for(; no_nodes<N; ++no_nodes) {
        size_t sum=0;
	for(size_t i=0; i<no_nodes; ++i) {
	    sum += degree[i];
	    available[i] = true;    
	}
	
	for(size_t j=0; j<m; ++j) {
	     double rnd=double(rand())/RAND_MAX;
	     size_t next=0;
	     
	     while(rnd > double(degree[next])/ sum || !available[next]) {		
	         if(available[next])
		     rnd -= double(degree[next])/sum;
		 
	         ++next;
		 
		 if(next > N) exit(0);
	     }
	     
	     available[next]=false;
	     sum -= degree[next];
	     edges.push_back(std::pair<size_t, size_t>(next, no_nodes));
             ++degree[next];
             ++degree[no_nodes];
	}	      
    } 
	    
    if(no_nodes != N) {
        std::cout << "ERROR: Something went wrong!" << std::endl;  
    }
}




void create_barabasi_albert(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t M, 
			    bool multiple=false, bool self_loop=false, bool directed=true) {
    std::vector<size_t> fun;
    size_t fun_sum=0;
    for(size_t i=0; i<N; ++i)
        fun.push_back(0);
  
    // Implement first edge between node 0 and 1
    edges.push_back( std::pair<size_t, size_t>(0,1) );
    ++fun[0];
    ++fun[1];
    fun_sum += 2;
    M -= 1;
    
    // Proceeding with node 2 and higher
    for(size_t i=2; i<N; ++i) {
        size_t no_n_edges = M/(N-i) > i ? i : M/(N-i);
	
	//std::cout << "i=" << i << "  no_n_edges=" << no_n_edges << std::endl;
	
	for(size_t k=0; k<no_n_edges; ++k) {
	    //std::cout << "K=" << k << "   fun_sum=" << fun_sum << std::endl;
	    size_t l=0;
	    bool done=false;
	    
	    do {
	        l=0;
	        long r=rand()%fun_sum;
		while(r >= fun[l]) {
		    r -= fun[l];
		    ++l;
		}
	      
	        if(l != i || self_loop) {
	            if(directed) {
	                if(find(edges.begin(), edges.end(), 
			        std::pair<size_t, size_t>(i, l)) == edges.end() ||
		           multiple) {
			    if(rand()%2)
		                edges.push_back(std::pair<size_t, size_t> (i, l) );
			    else
			        edges.push_back(std::pair<size_t, size_t> (l, i) );
			    
			    ++fun[i];
			    ++fun[l];
			    fun_sum += 2;
			    
		            done = true;
		        }		  
		    } else {
	                if(find(edges.begin(), edges.end(), 
			        std::pair<size_t, size_t>(i, l)) == edges.end() ||
		           multiple) {
		            edges.push_back(std::pair<size_t, size_t> (i, l) );
			    ++fun[i];
			    ++fun[l];
			    fun_sum += 2;
		            done = true;
		        }
		    }
                }
	    } while(!done);
	}
	
	M -= no_n_edges;
    }
}







/*
 * Creates a simple erdos-renyi-network with n nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero. If 
 * allow_multiple is set true self loops and the multiple occurance of the same 
 * link is allowed.
 */

void create_erdos_renyi(std::vector< std::pair<size_t, size_t> > &edges, size_t n, size_t M, 
			bool allow_multiple=false, bool self_loop=false, bool directed = false) {
    for(size_t j=edges.size(); j<M; ++j) {
        bool found=false;
	
	do {  
	    size_t first=rand()%n;  // Create node indices for link j randomly
            size_t second=rand()%n;
		
	    // If applicable check for self loops and multiple occorance
	    if(first != second || self_loop) {
	        if(allow_multiple) {
	            found = true;	 
		    edges.push_back(std::pair<size_t, size_t> (first, second) );
	        } else if(directed) {

		} else {
		    if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() &&
		       find(edges.begin(), edges.end(), std::pair<size_t, size_t>(second, first)) == edges.end()) {
		        edges.push_back(std::pair<size_t, size_t> (first, second) );
		        found = true;
		    }		  
		}  
	    }
	} while (!found);
    }
}



/*
 * Creates an Erdos Renyi network with an stochastical algorithm with N nodes
 * and on average M edges. From M and the maximal number of edges a global 
 * probability for creating edges is calculated.
 */

void create_erdos_renyi_stoch(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t M, 
			bool self_loop=false, bool directed = false) {
    if(!self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5);  
      
        for(size_t i=0; i<N; ++i) 
	    for(size_t j=i+1; j<N; ++j)
	        if(double(rand())/RAND_MAX < d)
		    edges.push_back(std::pair<size_t, size_t> (i, j) );        
    } else if (self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5+N);  
      
        for(size_t i=0; i<N; ++i) 
	    for(size_t j=i+1; j<N; ++j)
	        if(double(rand())/RAND_MAX < d)
		    edges.push_back(std::pair<size_t, size_t> (i, j) );        
    } else if (!self_loop && directed) {
        double d=double(M)/(N*N-N); 
      
        for(size_t i=0; i<N; ++i) 
	    for(size_t j=0; j<N; ++j)
	        if(double(rand())/RAND_MAX < d && i!=j)
		    edges.push_back(std::pair<size_t, size_t> (i, j) );      
      
    } else if (self_loop && directed) {
        double d=double(M)/(N*N);  
      
        for(size_t i=0; i<N; ++i) 
	    for(size_t j=0; j<N; ++j)
	        if(double(rand())/RAND_MAX < d)
		    edges.push_back(std::pair<size_t, size_t> (i, j) );
    }
}


/*
 * An older version - implementation of the Watts Strogatz algorithm.
 */

void create_watts_strogatz_old(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t K, double beta, bool allow_multiple=false) {
    if(allow_multiple)
        std::cerr << "Warning: create_watts_strogatz called with allow_multiple. Not implemented yet!" << std::endl;
  
    // create regular lattice
    for(size_t i=1; i<=K/2; ++i) 
        for(size_t j=0; j<N; ++j) 
	    edges.push_back( std::pair<size_t, size_t>(j, (j+i)%N) );    
	
 
    for(size_t i=0; i<edges.size(); ++i) {
        if(double(rand())/RAND_MAX < beta) {
	    size_t first=edges[i].first;
	    bool found=false;
	    
	    do {
	        size_t second=rand()%N;
		
		if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() &&
		   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(second, first)) == edges.end() &&
		   second != first) {
		    edges[i].second = second;
		    found = true;
	        }
	    } while (!found);
	}
    }
}






void create_watts_strogatz(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t M, double beta, 
			   bool multiple=false, bool self_loop=false, bool directed=true) {
    // PARAMETERS ARE ONLY CHECKED FOR REWIRING
  
    size_t M_tmp=M;
  
    for(size_t dif=1; M_tmp > 0; ++dif) {
        for(size_t base=0; base<N && M_tmp > 0; ++base) {
	    size_t n1=base;
	    size_t n2=(base+dif)%N;
	  
	    if(rand()%2) {
	        edges.push_back( std::pair<size_t, size_t> (n1, n2) );  
	    } else {
	        edges.push_back( std::pair<size_t, size_t> (n2, n1) );
	    }
	    
	    --M_tmp;
	}   
    }
 
    // Shuffle
    for(size_t i=0; i<edges.size(); ++i) {
        if(double(rand())/RAND_MAX < beta) {
	    size_t first=edges[i].first;
	    bool found=false;
	    
	    do {
	        size_t second=rand()%N;
		
		if(second != first || self_loop) {
		    if(directed) {
		        if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() ||
		           multiple) {
		            edges[i].second = second;
		            found = true;
			}		      
		    } else {
		        if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() &&
		           find(edges.begin(), edges.end(), std::pair<size_t, size_t>(second, first)) == edges.end() ||
		           multiple) {
		            edges[i].second = second;
		            found = true;
			}
	            }
		}
	    } while (!found);
	}
    }
}





/* bad written class for simple bidirectional network analysis. not developed at the 
 * moment because R and igraph can do better.
 */


class bd_network {
    size_t nodes, edges;
    bool* ad_m;
    size_t* sh_path;
    
    std::vector<size_t> degrees, deg_hist, sh_path_hist, cluster_id, cluster_size;
    std::vector<double> deg_dist, sh_path_dist, clustering_coeff;
    size_t min_deg, max_deg, max_sh_path, smallest_cluster, biggest_cluster;
    double avg_deg, avg_sh_path, biggest_cluster_ratio, avg_clustering_coeff;
    
  
public:
    /*
     * INPUT OUTPUT CODE
     *
     */
    void add_edge(size_t a, size_t b) {
        if(ad_m[a*nodes+b])
            return;
	  
	++edges;
        ad_m[a*nodes+b]=true;
	ad_m[b*nodes+a]=true;
    }
    
  
    bd_network(size_t nd) : nodes(nd), edges(0), avg_deg(0), avg_sh_path(0), 
	    biggest_cluster_ratio(0), min_deg(0), max_deg(0), max_sh_path(0), 
	    smallest_cluster(0), biggest_cluster(0), avg_clustering_coeff(0) {
        ad_m = new bool[nodes*nodes];
	sh_path = new size_t[nodes*nodes];
	
	for(size_t i=0; i<nodes*nodes; ++i) {
	    ad_m[i]=false;
	    sh_path[i]=0;    
	}
    }
    
    
    bd_network(size_t N, std::vector< std::pair<size_t, size_t> > &edg) : nodes(N), 
            edges(0), avg_deg(0), avg_sh_path(0), biggest_cluster_ratio(0), min_deg(0), 
            max_deg(0), max_sh_path(0), smallest_cluster(0), biggest_cluster(0), avg_clustering_coeff(0)  {
        ad_m = new bool[nodes*nodes];
	sh_path = new size_t[nodes*nodes];
	
	for(size_t i=0; i<nodes*nodes; ++i) {
	    ad_m[i]=false;
	    sh_path[i]=0;    
	}           
      
        for(size_t i=0; i<edg.size(); ++i) 
	    add_edge(edg[i].first, edg[i].second);
    }
    
    
    ~bd_network() {
        delete[] ad_m; 
	delete[] sh_path;
    } 
    
    
    bool is_connected() {
        for(size_t i=0; i<nodes*nodes; ++i)
	    if(sh_path[i] > edges)
	        return false;
	    
	return true;
    }
  
    std::vector<size_t>& get_degrees() {  return degrees;  }
    
    std::vector<size_t>& get_degree_hist() {  return deg_hist;  }
    
    std::vector<double>& get_degree_dist() {  return deg_dist;  }
    
    size_t get_shortest_path(size_t a, size_t b) {  return sh_path[b*nodes+a];  }
    
    std::vector<size_t>& get_path_hist() {   return sh_path_hist;  }
    
    std::vector<double>& get_path_dist()  {  return sh_path_dist;  }
    
    size_t get_min_degree() {  return min_deg;  }
    
    size_t get_max_degree() {  return max_deg;  }
    
    double get_avg_degree() {  return avg_deg;  }
    
    size_t get_max_path() {  return max_sh_path;  }
    
    size_t get_diameter() {  return max_sh_path;  }
    
    double get_avg_path() {  return avg_sh_path;  }
    
    double get_biggest_cluster_ratio() {  return biggest_cluster_ratio;  }
    
    size_t get_smallest_cluster() {  return smallest_cluster;  }
    
    size_t get_biggest_cluster() {  return biggest_cluster;  }
    
    size_t get_cluster_number() {  return cluster_id.size(); }
        
    std::vector<double>& get_clustering_coeff() {  return clustering_coeff;  } 
    
    double get_avg_clustering_coeff() {  return avg_clustering_coeff;  }
    
  
    /*
     * ANALYSIS CODE
     *
     */
    
    
    double get_steady_state_A(size_t a, size_t b, std::vector<double>& concentration) {
        long i, j, bd, s;
        mpz_t *mp_B, *mp_N;

        /* generate a n x m random left hand side matrix A */
        size_t n = nodes;
        size_t m = nodes+1;
        long* A = new long[n*m];    

	// initializing
	for(size_t i=0; i<n*m; ++i)
	    A[i] = 0;
	
	// Creating stoichiometry matrix for (inner) network interaction
	for(size_t i=0; i<nodes; ++i) {
	    for(size_t j=i+1; j<nodes; ++j) {
	        if(ad_m[i*nodes+j]) {	  
	            A[i*m+j] += 1;
	            A[i*m+i] -= 1;
	            A[j*m+i] += 1;
	            A[j*m+j] -= 1;
		}
	    }
	}
	
	// Adding boundary condition / flow
	A[a*m+m-1]=1;
        A[b*m+m-1]=-1;
	
	

        s = nullspaceLong (n, m, A, &mp_N);
        

    
	mpz_t C_11, C_12, C_21, C_22, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, det;
	mpz_init(C_11);   mpz_init(C_12);   mpz_init(C_21);   mpz_init(C_22);   mpz_init(tmp1); 
	mpz_init(tmp2);   mpz_init(tmp3);   mpz_init(tmp4);   mpz_init(tmp5);   mpz_init(tmp6);   mpz_init(det);
	
	mpz_set(C_11, mp_N[a * s + 0]);
	mpz_set(C_12, mp_N[a * s + 1]);
	mpz_set(C_21, mp_N[b * s + 0]);
	mpz_set(C_22, mp_N[b * s + 1]);
	
        delete [] A;
  

	mpz_mul(tmp1, C_11, C_22);
	mpz_mul(tmp2, C_21, C_12);
	mpz_sub(det, tmp1, tmp2);
	
	
	mpz_set(tmp3, mp_N[nodes * s + 0]);
	mpz_set(tmp4, mp_N[nodes * s + 1]);
	mpz_mul(tmp1, C_22, tmp3);
	mpz_mul(tmp2, C_21, tmp4);
	mpz_sub(tmp3, tmp2, tmp1);
	
	mpq_t v;
	mpq_init(v);
        mpq_set_num(v, tmp3);
	mpq_set_den(v, det);
	double rs=mpq_get_d(v);
	
	{   // Calculate the solution vector
	    concentration.clear();
	    
	    for(size_t i=0; i<=nodes; ++i) {
	        double d=0.0;
	      	
		mpz_set(tmp3, mp_N[i * s + 0]);
	        mpz_set(tmp4, mp_N[i * s + 1]);
	        mpz_mul(tmp1, C_22, tmp3);
	        mpz_mul(tmp2, C_21, tmp4);
	        mpz_sub(tmp3, tmp2, tmp1);
	
                mpq_set_num(v, tmp3);
	        mpq_set_den(v, det);
	        d=mpq_get_d(v);
		
	        concentration.push_back(d);
	    }
	}
	
        for (i = 0; i < m * s; i++)
            mpz_clear (mp_N[i]);
        free (mp_N);
	
	// mpz cleanup
	mpz_clear(C_11);   mpz_clear(C_12);   mpz_clear(C_21);   mpz_clear(C_22);
	mpz_clear(tmp1);   mpz_clear(tmp2);   mpz_clear(tmp3);   mpz_clear(tmp4);	
	mpz_clear(tmp5);   mpz_clear(tmp6);   mpz_clear(det);    mpq_clear(v);
	
	return rs;
    }
    
    
   
    
    double get_steady_state_B(size_t a, size_t b, std::vector<double>& concentration) {
        double alpha=get_steady_state_A(a, b, concentration);
        
	for(size_t i=0; i<concentration.size(); ++i)
	    concentration[i] *= (alpha)/(1+2*alpha);
	  
      
	return (alpha)/(1+2*alpha);
    }
    
  
    void calculate_degrees() {  
        degrees.clear();
	deg_hist.clear();
	deg_dist.clear();
      
        // Calculate degree of every node
        for(size_t i=0; i<nodes; ++i) {
	    degrees.push_back(0);  
	    
	    for(size_t j=0; j<nodes; ++j)
	        if(ad_m[i*nodes+j])
		    ++degrees[i];
	}      
	
	// Calculate minimal, maximal and average degree
	avg_deg=0;
	min_deg=max_deg=degrees[0];
	
	for(size_t i=0; i<nodes; ++i) {
	    avg_deg += degrees[i];  
	    
	    if(degrees[i] < min_deg) min_deg=degrees[i];
	    if(degrees[i] > max_deg) max_deg=degrees[i];
	}
	
	avg_deg /= nodes;
	
	// Calculate degree distribution (first histogramm than norm. dist.)	
	for(size_t i=0; i<=max_deg; ++i) {
	    deg_hist.push_back(0);
	    
	    for(size_t j=0; j<nodes; ++j) 
	        if(degrees[j] == i)
		    ++deg_hist[i];
		
	    deg_dist.push_back(double(deg_hist[i])/nodes);
	}
      
    }
  
 
    /*
     * Calculating the shortest paths between all pairs of network nodes...
     */

    void calculate_sh_paths() {
        sh_path_dist.clear();
	sh_path_hist.clear();
      
        // Calculate length of paths	
        for(size_t i=0; i<nodes; ++i) {
	    // The row is inicialiced
            for(size_t j=0; j<nodes; ++j) 
	        sh_path[i*nodes+j]=nodes+2;  // Marks a not yet known distance
	    
	    sh_path[i*nodes+i]=0;  // initial condition
	    
	    // now the distance is calculated recursively
	    //calculate_dmap(i, i);
	    
	    size_t current=0;
	    size_t highest=nodes+2;
	    while(highest==nodes+2) {
	        highest = 0;
		
		//std::cout << "current is " << current << std::endl;
		
		for(size_t j=0; j<nodes; ++j) {
		    if(sh_path[i*nodes+j] == current) {		      
		        for(size_t k=0; k<nodes; ++k) 
			    if(ad_m[j*nodes+k] && sh_path[i*nodes+k] > current)
			        sh_path[i*nodes+k] = current+1;
		    } else if(sh_path[i*nodes+j] > highest) {
		        highest = sh_path[i*nodes+j];
		    }
		}
	        
	        ++current;
	    }
	}
		
	
	// Calculate shortest path length distribution / average / maximum 
	max_sh_path=0;
	size_t count=0;
	
	// Iterating all pairs of nodes
        for(size_t i=0; i<nodes; ++i) {
            for(size_t j=0; j<nodes; ++j) {
	        if(sh_path[i*nodes+j] != 0 &&        // Path length has to be finite
		   sh_path[i*nodes+j] <= edges) {
		   
		   avg_sh_path += sh_path[i*nodes+j];
		   ++count;
		
		   // updating maximum shortest path
		   if(max_sh_path < sh_path[i*nodes+j])
		       max_sh_path = sh_path[i*nodes+j];
		}
	    }
	} 
	
	// Initializing path histogramm vector
	for(size_t i=0; i<=max_sh_path; ++i)
	    sh_path_hist.push_back(0);
	
	// Calculating histogramm
	for(size_t i=0; i<nodes; ++i) 
            for(size_t j=0; j<nodes; ++j) 
	        if(sh_path[i*nodes+j] != 0 &&
		   sh_path[i*nodes+j] <= edges) 
		    ++sh_path_hist[sh_path[i*nodes+j]];  
	
	// .. average path length..
	avg_sh_path /= count;
	
	// ..and shortest path length distribution.
	for(size_t i=0; i<=max_sh_path; ++i)
	    sh_path_dist.push_back(double(sh_path_hist[i])/count);
    }
  
  
    /*
     * Method calculates the clustering coefficient for each node and the average
     * clustering coefficient.
     *
     * Degrees have to be calculated first!
     */
     
    void calculate_clustering_coeff() {
        clustering_coeff.clear();
        avg_clustering_coeff=0;
	
	
	for(size_t i=0; i<nodes; ++i)
	    clustering_coeff.push_back(0);
	    
	for(size_t i=0; i<nodes; ++i) 
	    for(size_t j=0; j<nodes; ++j) 
	        if(ad_m[i*nodes+j]) 
		    // Found pair of connected nodes. Now look for nodes that are 
		    // connected to both....  
		    for(size_t k=0; k<nodes; ++k)
		        if(ad_m[k*nodes+i] && ad_m[k*nodes+j]) 
			    clustering_coeff[k] += 1.0;
	
	// Normalizing the counted size...
	for(size_t i=0; i<nodes; ++i) 
	    if(degrees[i] > 1)
	        clustering_coeff[i] /= (degrees[i]*(degrees[i]-1));
	    else 
	        clustering_coeff[i] = -1;
	
	
	// calculate average clustering coefficient
	avg_clustering_coeff=0;
	size_t cc_count=0;
	for(size_t i=0; i<clustering_coeff.size(); ++i)
	    if(clustering_coeff[i] >= 0) { 
	        avg_clustering_coeff += clustering_coeff[i];
	        ++cc_count;
	    }
	
	avg_clustering_coeff /= cc_count;
    }
  
  
  
    void m_draw_clusters_rec(size_t node, size_t next_color=0) {
        // next_color==0 bedeutet es sollen alle cluster bestimmt werden...
        if(next_color==0) {
            for(size_t i=0; i<nodes; ++i)
	        if(cluster_id[i] == 0)
	            m_draw_clusters_rec(i, ++next_color);
      
	     return;
        }
  
        // Ist der aktuelle Node schon gefärbt wird die Rekursion abgebrochen
        if(cluster_id[node] != 0)
            return;
  
        // Setze die Farbe des aktuellen Nodes (node) auf next_color
        cluster_id[node]=next_color;
  
        // färbe rekursiv alle mit node verbundenen nodes
        for(size_t i=0; i<nodes; ++i) 
	  if(ad_m[node*nodes+i])
	        m_draw_clusters_rec(i, next_color);
      
    }
  
    
    void calculate_con_clusters() {
        cluster_id.clear();
	cluster_size.clear();
	
	for(size_t i=0; i<nodes; ++i)
	    cluster_id.push_back(0);
      
	m_draw_clusters_rec(0, 0);
	
	// calculate number of clusters
	size_t no_clusters=0;
	for(size_t i=0; i<nodes; ++i)
	    if(cluster_id[i] > no_clusters)
	        no_clusters = cluster_id[i];
	
	// calculating cluster sizes
	for(size_t i=0; i<no_clusters; ++i)
	    cluster_size.push_back(0);
	
	for(size_t i=0; i<nodes; ++i) 
	    ++cluster_size[cluster_id[i]-1];
	
	// calculating maximum and minimum cluster size
        smallest_cluster=edges;
	biggest_cluster=0;
	
	for(size_t i=0; i<cluster_size.size(); ++i) {
	    if(cluster_size[i]<smallest_cluster)
	        smallest_cluster=cluster_size[i];
	    
	    if(cluster_size[i]>biggest_cluster)
	        biggest_cluster=cluster_size[i]; 
	}
	
	biggest_cluster_ratio = double(biggest_cluster)/nodes; 
    }
  
  
    void do_analyis() {
        std::cout << "deg..";  	std::cout.flush();
	calculate_degrees();
        std::cout << "sp..";	std::cout.flush();
	calculate_sh_paths();
        std::cout << "cc..";	std::cout.flush();
	calculate_clustering_coeff();
        std::cout << "con..";	std::cout.flush();
	calculate_con_clusters();
	std::cout << std::endl;
     } 
};




















#endif // NETWORK_TOOLS_H