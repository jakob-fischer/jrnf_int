/* author: jakob fischer (jakob@automorph.info)
 * date: 9th February 2013
 * description: 
 * Tools for creating different types of complex networks and for
 * analyzing them.
 *
 * TODO: In couple-BA  how should probability of couples be calculated?
 */

#ifndef NETWORK_TOOLS_H
#define NETWORK_TOOLS_H

#include <vector>
#include <set>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <iostream>
#include <numeric>
#include <algorithm>



/*
 * Class that allows it to sample from a huge distribution given 
 * by (unnormalized) integer values. Using a binary tree makes this
 * class very fast and also allows it to change parts of the 
 * distribution on log(N) time.
 */

class bt_draw {
    // Maximum number of generations of the binary tree
    static const size_t e_max=40; 
	
    /*
     * Contains the probabilities of the binary tree. probs[0] is the
     * vector with elementary probabilities, prob[1], prob[2], ... are
     * the higher / earlier generations. 'last' contains the number of
     * generations.
     */
    size_t last;
    std::vector<size_t> probs[e_max];
    
public:	


    /*    
     * Standard constructor
     */

    bt_draw() : last(0) {}
	
	
    /*
     * Function for changing the probability of the 'value' 'id' to 'v'.
     * The method is changing all summed probabilities higher in the 
     * tree.
     */	
	
    bt_draw& set(size_t id, size_t v) {
        size_t old=probs[0][id];

	for(size_t i=0; i<e_max && probs[i].size() != 0; ++i) {
	    probs[i][id/pow(2,i)]  -= old;
	    probs[i][id/pow(2,i)]  += v;	    	    
	}
	    
	return (*this);
    }
	

    /*
     * Extends the range of the random value by one. This outcome is 
     * occuring with (relative) probability t 
     */
	
    bt_draw& add(size_t t) {
        // create new entry on lower level
        probs[0].push_back(0);
	
        // Change tree accordingly - As long as one level higher there
        // is half the number of nodes we have to add one to connect to    
	size_t l=0;
	while(probs[l].size() != 2*probs[l+1].size() &&
	      probs[l].size() != 1) {
	    probs[l+1].push_back(0);
	    	    
	    if(probs[l+1].size() == 1) {
	        ++last;
	        probs[l+1][0] = probs[l][0];
	    }
	        
	    ++l;
	}
	    
        // use set-function to change entry+summation in tree
	return set(probs[0].size()-1, t);	
    }
	
	
    /*
     * Returns the probability of the value with id 'i'
     */

    size_t get(size_t i) {
	return probs[0][i];
    }
	

    /* 
     * Draws one number from the distribution. Ideally 
     * srand should have been called somewhere before.
     */
	
    size_t draw() const {
        // Object is empty. Report error
        if(last == 0) {
	    if(probs[0].size() == 0)
                std::cout << "bt_draw::draw() - 'last' is zero!" << std::endl;

	    return 0;
	}
			
        // Draw random number modulo sum of all probabilities
	size_t rnd=rand()%probs[last][0];
	    
        size_t z=0; // Selected element on the level 'i'

        // Going down the tree from top to bottom
	for(size_t i=last; i>0; --i) {
	    // Either not at the right flank or there, but having two childs
	    if(z != probs[i].size()-1 || 2*probs[i].size() == probs[i-1].size()) {
	        if(probs[i-1][z*2] > rnd) { // going down on the left
	            z=z*2;
	        } else {                    // going down on the right
                    rnd = rnd - probs[i-1][z*2];
	            z=z*2+1;
	        }
	    } else { // Just having one child
	        z=z*2;
	    }    	    
	}

	return z;
    }

      
    /*
     * Returns size of the distribution
     */
	
    size_t size() const {
	return probs[0].size();	
    }


    /*
     * Copy constructor generating bt_draw from vector
     */
  
    bt_draw(const std::vector<size_t>& p) : last(0)  {        	
	for(size_t i=0; i<p.size(); ++i)
	    add(p[i]);
    }

    /*
     * Copy constructor, constructing bt_draw from bt_draw
     */

    bt_draw(const bt_draw& ref) : last(0) {
	for(size_t i=0; i<ref.probs[0].size(); ++i)
	    add(ref.probs[0][i]);
    }
};






/*
 * Helper for "create_pan_sinha".
 * Calculates the size (number of nodes) of the module on level 'level' with the
 * index 'element. It needs a vector 'elem_c' with sizes of elementary modules (level 0)
 * and the parameter 'm' 
 */

size_t bps_module_size(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    // Level zero is easy
    if(level == 0) 
        return elem_c[element];
    
    // Sum over m elementary modules for level = 1
    if(level == 1) {
        size_t sum=0;
	
	for(size_t i=0; i<m; ++i) 
	    sum += elem_c[element*m+i];
	
	return sum;
    }
    
    // Recursive addition on higher level
    return bps_module_size(elem_c, m, level-1, 2*element) +
           bps_module_size(elem_c, m, level-1, 2*element+1);  
}


/*
 * Similar to previous function, but calculates the offset. 
 * (between first element of the network and first one in the module)
 */

size_t bps_module_offset(const std::vector<size_t> &elem_c, size_t m,
			      size_t level, size_t element) {
    // On level 0 sum over elementary modules before 'element'
    if(level == 0) {
        size_t sum=0;
	
        for(size_t i=0; i<element; ++i)
            sum += elem_c[i];
	  
        return sum;
    }
	
    // On level 1, take the offset of the first module one level lower
    if(level == 1)
	    return bps_module_offset(elem_c, m, 0, m*element);
    
    // For all other do the same, but here replace 'm' with 2
    return bps_module_offset(elem_c, m, level-1, 2*element);
}


/*
 * Find the number (0-indexed) of the module on level "level" in which
 * the node "node" is located
 */

size_t bps_module_no(const std::vector<size_t> &elem_c, size_t m,
	             size_t level, size_t node) {
    size_t max_modules=0;

    // Get number of modules
    if(level==0) 
	max_modules=elem_c.size();
    else 
        max_modules=elem_c.size()/(m*pow(2,level-1)); 
		
    // Iterate through all modules		  
    for(size_t i=0; i<max_modules; ++i) {
	size_t offset=bps_module_offset(elem_c, m, level, i);
        size_t size=bps_module_size(elem_c, m, level, i);
        
	if(node >= offset  &&
           node < offset+size)
            return i; 
    }			  
	
    std::cout << "Error: bps_module_no - did not find module!" << std::endl;
    return 0;
}


/*
 * Calculates the maximal number of edges for networks with 'nodes' nodes.
 * Depending on if 'self_loop's are allowed and if the network is 'directed'.
 * Multiple edges are not taken into account, because if they are allowed 
 * there is obviously no maximum
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
 * This function returns the (lowest) level on which the two nodes 'a' and 'b'
 * are in the same module. This allows to calculate the probability the 
 * edge betwen them exists.
 */
 
 size_t bps_nodes_level_comp(size_t a, size_t b, 
                             std::vector<size_t>& elem_c, size_t m) {
    size_t lvl=0;
    while(bps_module_no(elem_c, m, lvl, a) !=	
          bps_module_no(elem_c, m, lvl, b))
        ++lvl;							 
					
    return lvl;							 
}


/*
 * Macro function that does the same as bps_nodes_level_comp but uses a lookup
 * array in which for each node and level the id / number of the module is stored.
 */

size_t bps_nodes_level_comp_lu(const std::vector<size_t> &modules_no, size_t N,
                               size_t a, size_t b) {
    size_t lvl=0;
    while(modules_no[N*lvl+a] != modules_no[N*lvl+b])
        ++lvl;							 
					
    return lvl;										   
}


/*
 * Function which is given two nodes <a> and <b>  
 */

bool bps_same_module(const std::vector<size_t> modules_size, size_t a, size_t b) {
    bool same_module=true;

    size_t border=0;
    for(size_t k=0; k<modules_size.size(); ++k) {
        border += modules_size[k];
        if(a < border && b >= border || a >= border && b < border)
            same_module=false; 
    }

    return same_module;
}



/*
 * Creates complex network following the model of Pan & Sinha. In this model a
 * network with hierarchical modular structure is created by drawing links from
 * the set of possible links (NxN) and choosing the linking probably accordingly.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * h               - number of binary levels
 * m               - Number of modules on level 0 that are put in one module on
 *                   level 1 (on higher levels it is allways to / binary tree)
 * r               - Change linking probability by going one level up (linking
 *                   probability between two nodes is choosen from the lowest 
 *                   level on which they are in the same module)
 */ 

void create_pan_sinha(std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t N, size_t M, size_t h, size_t m, double r, 
                      bool multiple =false, bool self_loop=false, 
                      bool directed=false) {
    size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 to each until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) {  
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }

    std::cout << "Elementary size is " << elem_c[0] << std::endl;
    
    // Creating lookup table for asociating every node on every level with one module
    // starting with the N entries for level 0, than N entries for level 1, then 2...
    std::vector<size_t> module_no;
    size_t module_no_lvl=h+2;
    for(size_t i=0; i<module_no_lvl; ++i) 
        for(size_t j=0; j<N; ++j) 
	    module_no.push_back(bps_module_no(elem_c, m, i, j));
    
    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<N*N; ++count) {    
        // translate to nodes
	size_t i=count/N;
	size_t j=count%N;
		
        // If non-directed than we require i <= j (more efficient and leads to i==j 
        // having same probability than for example i-2 == j
	if((directed || i <= j) && (self_loop  || i != j))  
	    nw_prob.add( pow(r, bps_nodes_level_comp_lu(module_no, N, i, j))*100000 );
	else 
	    nw_prob.add(0);
    }
    
    // Now create network
    for(size_t a=0; a<M; ++a) {  	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
	// translate to edge
	size_t i=count/N;
	size_t j=count%N;
			
        edges.push_back(std::pair<size_t, size_t> (i, j) );
            
        if(!multiple)  // if no multiples allowed set prob. for next run to 0
            nw_prob.set(count, 0);	
    }
}



/*
 * Creates a simple barabasi-albert-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */ 
 
void create_barabasi_albert(std::vector< std::pair<size_t, size_t> > &edges, 
                            size_t N, size_t M, bool multiple=false, 
                            bool self_loop=false, bool directed=true) {

    bt_draw fun;
 
    // Implement first edge between node 0 and 1
    edges.push_back( std::pair<size_t, size_t>(0,1) );
    fun.add(1);
    fun.add(1);
    M -= 1;
    
    // Proceeding with node 2 and higher
    for(size_t i=2; i<N; ++i) {
        // - distribute number of edges that remain to be placed on number
	// of nodes that have to be "created"
	// - maximum number of edges that are going to be connected to the 
        // current node is the number of already existing nodes
        size_t no_n_edges = M/(N-i) > i ? i : M/(N-i);
        
        fun.add(0);
	
        for(size_t k=0; k<no_n_edges; ++k) {
	    // try creating edge number k (until done is true)
	    bool done=false;
	            
	    do {
		// select node randomly with probability of their functionality
	        size_t l=fun.draw();
                 
                //std::cout << "   drawn  l=" << l << std::endl;

	        // proceed if no self loop or self loop allowed
	        if(l != i || self_loop) {
	            if(directed) {
	                if(find(edges.begin(), edges.end(), 
		                std::pair<size_t, size_t>(i, l)) == edges.end() ||
		           multiple) {
							   
			    // directed edges are directed randomly   
			    if(rand()%2)
		                edges.push_back(std::pair<size_t, size_t> (i, l) );
			    else
			        edges.push_back(std::pair<size_t, size_t> (l, i) );
			    
                            fun.set(i, fun.get(i)+1);
                            fun.set(l, fun.get(l)+1);
			    
		            done = true;
		        }		  
		    } else {
	                if(find(edges.begin(), edges.end(), 
			        std::pair<size_t, size_t>(i, l)) == edges.end() ||
		           multiple) {
	                    // for non directed the higher numbered node is always
	                    // put in first
		            edges.push_back(std::pair<size_t, size_t> (i, l) );

                            fun.set(i, fun.get(i)+1);
                            fun.set(l, fun.get(l)+1);

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
 * Creates a simple watts-strogatz-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * beta            - parameter governing what fraction of nodes is shuffled 
 * multiple        - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void create_watts_strogatz(std::vector< std::pair<size_t, size_t> > &edges, 
                           size_t N, size_t M, double beta, bool multiple=false, 
                           bool self_loop=false, bool directed=true) {
    // the parameters self_loop, directed and multiple are only checked
    // for the shuffling
    size_t M_tmp=M;     // counts the remaining number of links to be created
    
    for(size_t dif=1; M_tmp > 0; ++dif) {
        for(size_t base=0; base<N && M_tmp > 0; ++base) {
            size_t n1=base;
	    size_t n2=(base+dif)%N;
	  
	    // Random direction if directed, if not there is no harm done
	    if(rand()%2) {
	         edges.push_back( std::pair<size_t, size_t> (n1, n2) );  
	    } else {
	        edges.push_back( std::pair<size_t, size_t> (n2, n1) );
	    }
            --M_tmp;
	}   
    }
 
    // Shuffle second node of edge 'i' with probability beta
    for(size_t i=0; i<edges.size(); ++i) {
        if(double(rand())/RAND_MAX < beta) {
	    size_t first=edges[i].first;
	    bool found=false;
	    
	    // draw random numbers for second node until you find one 
	    // that matches the constraints
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



/*
 * Creates a simple erdos-renyi-network with N nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void create_erdos_renyi(std::vector< std::pair<size_t, size_t> > &edges, 
                        size_t N, size_t M, bool allow_multiple=false, 
                        bool self_loop=false, bool directed = false) {
    for(size_t j=edges.size(); j<M; ++j) {
        // Found a link that is allowed given the existing edges?
        bool found=false;
	
	do {  
	    size_t first=rand()%N;  
            size_t second=rand()%N;
		
	    // If applicable check for self loops and multiple occorance
	    if(first != second || self_loop) {
	        if(allow_multiple) { // If multiple links is allowed no reason to check for direction
	            found = true;	 
	            edges.push_back(std::pair<size_t, size_t> (first, second) );
	        } else if(directed) { // If directed and no multiples allowed check once
                    if(find(edges.begin(), edges.end(), 
		            std::pair<size_t, size_t>(first, second)) == edges.end()) {
		        edges.push_back(std::pair<size_t, size_t> (first, second) );
		        found = true;
		    }	
                } else {              // Else link can have two forms
                    if(find(edges.begin(), edges.end(), 
		            std::pair<size_t, size_t>(first, second)) == edges.end() &&
		       find(edges.begin(), edges.end(), 
		            std::pair<size_t, size_t>(second, first)) == edges.end()) {
		        edges.push_back(std::pair<size_t, size_t> (first, second) );
		        found = true;
		    }		  
		}  
	    }
	} while (!found);
    }
}



/*
 * Creates a simple modular network. In this model a network with hierarchical modular structure is created by drawing links from
 * the set of possible links (NxN) and choosing the linking probably accordingly.
 * 
 * edges           - reference to edgeliste where links are written to (0-based)
 * N               - number of nodes in network
 * M               - number of edges / links in network
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * m                - Number of modules 
 * r               - Change linking probability by going one level up (linking
 *                   probability between two nodes is choosen from the lowest 
 *                   level on which they are in the same module)
 */ 

void create_simple_modular(std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t N, size_t M, size_t m, double r, 
                      bool multiple =false, bool self_loop=false, 
                      bool directed=false) {
    size_t n_el_mod=m;                             // number of modules 
    size_t n_floor = floor(double(N)/n_el_mod);    // min no nodes in modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 to each until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) {  
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }

    std::cout << "Elementary size is " << elem_c[0] << std::endl;
    
    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<N*N; ++count) {    
        // translate to nodes
	size_t i=count/N;
	size_t j=count%N;
		
        // If non-directed than we require i <= j (more efficient and leads to i==j 
        // having same probability than for example i-2 == j
	if((directed || i <= j) && (self_loop  || i != j)) { 
            if(bps_same_module(elem_c,i,j))
	        nw_prob.add( 1*100000 );
            else
                nw_prob.add( r*100000 );
	} else { 
	    nw_prob.add(0);
        }
    }
    
    // Now create network
    for(size_t a=0; a<M; ++a) {  	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
	// translate to edge
	size_t i=count/N;
	size_t j=count%N;
			
        edges.push_back(std::pair<size_t, size_t> (i, j) );
            
        if(!multiple)  // if no multiples allowed set prob. for next run to 0
            nw_prob.set(count, 0);	
    }
}



/*
 * For an Erdos-Renyi-Network this function selects edges to couple.
 * This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Erdos-Renyi statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 * 
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void couple_erdos_renyi(std::vector< std::pair<size_t, size_t> > &couples,
                        size_t C,
                        std::vector< std::pair<size_t, size_t> > &edges, 
                        bool maintain_p=false, bool allow_multiple=false, 
                        bool self_loop=false, bool directed = false) {
    // Calculate N and M 
    size_t M=edges.size();
    size_t N=0;
    for(size_t i=0; i<M; ++i) {
        if(edges[i].first > N)
            N = edges[i].first+1;

        if(edges[i].second > N)
            N = edges[i].second+1;
    }

    
    std::vector<bool> taken_f;  // which links are already part of a link pair
    for(size_t i=0; i<edges.size(); ++i)
        taken_f.push_back(false);

    std::vector<bool> edges_m; // edges matrix (necessary if you want to maintain (not) allow multiple) 
    if(maintain_p && !allow_multiple) {
        for(size_t i=0; i<N*N; ++i) 
            edges_m.push_back(false);

        for(size_t i=0; i<edges.size(); ++i)
            edges_m[edges[i].first*N+edges[i].second]=true;
    }


    while(C > 0) {  // C pairs have to be found
        size_t first=rand()%edges.size();	
	size_t second=rand()%edges.size();	
        
        bool legal=true;  // is this pair legal with all boundary conditions	



        // both partners have to be different and free...    
	if(first == second || taken_f[first] || taken_f[second])
            legal=false;

        // If properties from network creation are to be maintained special
        // checks are necessary
        if(maintain_p) {
            size_t ff(edges[first].first), fs(edges[first].second);
            size_t sf(edges[second].first), ss(edges[second].second);
                
            if(!self_loop && (ff == ss || fs == sf))
                legal=false;

            if(!allow_multiple) {
                if(edges_m[ff*N+ss] || edges_m[sf*N+fs]) 
                    legal=false;

                if(!directed)
                    if(edges_m[ss*N+ff] || edges_m[fs*N+sf]) 
                        legal=false;
            }
        } 


        // if still legal: do it!
        if(legal) {
	    couples.push_back(std::pair<size_t, size_t> (first, second));
	    taken_f[first]=true;
	    taken_f[second]=true;
	    --C;	
            
            size_t ff(edges[first].first), fs(edges[first].second);
            size_t sf(edges[second].first), ss(edges[second].second);

            if(maintain_p && !allow_multiple) 
                edges_m[ff*N+ss]=edges_m[sf*N+fs]=true;
	}
    }
}


// Macro helper to determine "distance" 
size_t cws_links_cy_dist(int a, int b, size_t N) {
    size_t a_amb = std::abs(a-b);
    size_t a_amb2 = std::abs(N-a_amb);

    return std::min(a_amb, a_amb2);
}


/*
 * For a Watts-Strogatz-Network this function this function selects edges to 
 * couple edges. This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Watts-Strogatz statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 *
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */
 
void couple_watts_strogatz(std::vector< std::pair<size_t, size_t> > &couples,
                           size_t C, std::vector< std::pair<size_t, size_t> > &edges,
                           double beta, bool maintain_p=false, bool allow_multiple=false, 
                           bool self_loop=false, bool directed = false) {
    // Calculate N and M 
    size_t M=edges.size();
    size_t N=0;
    for(size_t i=0; i<M; ++i) {
        if(edges[i].first+1 > N)
            N = edges[i].first+1;

        if(edges[i].second+1 > N)
            N = edges[i].second+1;
    }


    // Which edges are still available vor coupling?
    std::set<size_t> edge_available;
    
    for(size_t j=0; j<M; ++j)
	edge_available.insert(j);


    size_t close_v = ceil(double(M)/N); // needed to count which additional "connections" are far
    std::cout << "close_v=" << close_v << std::endl;

    std::vector<bool> edges_m; // edges matrix (necessary if you want to maintain (not) allow multiple) 
    std::vector<std::set<size_t> > edges_on_node;    // list for every node (to which edges is it connected?)
    for(size_t i=0; i<N; ++i) 
        edges_on_node.push_back(std::set<size_t>());

    if(maintain_p && !allow_multiple) {
        for(size_t i=0; i<N*N; ++i) 
            edges_m.push_back(false);

        for(size_t i=0; i<edges.size(); ++i) {
            edges_m[edges[i].first*N+edges[i].second]=true;

            edges_on_node[edges[i].first].insert(i);
            edges_on_node[edges[i].second].insert(i);
        }
    }

    // count number of far links (regarding cyclic structure further appart edges than close_v)
    //std::cout << "N=" << N << std::endl;
    size_t far_links = 0;
    for(size_t i=0; i<edges.size(); ++i) {
        if(cws_links_cy_dist(edges[i].first, edges[i].second, N) > close_v)
            ++far_links;

    }

    //double p=double(far_links)/M;
    //std::cout << "p=" << p << "   beta=" << beta << std::endl;
    
    size_t prob_far_far=1; //1000000.0*beta*beta;
    size_t prob_clo_clo=1000000.0;
    size_t prob_far_clo=10; //1000000.0*beta*(1-beta);

    std::cout << "far_links=" << far_links << "  N=" << N << "   M=" << M << std::endl;
    std::cout << "prob_ff=" << prob_far_far << "   prob_cc=" << prob_clo_clo << "   prob_fc=" << prob_far_clo << std::endl;

    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<M*(M-1)/2; ++count) {    
        // translate to edge
	size_t i = (size_t) (sqrt(2.0*count+0.5)+0.5);
	size_t j = count - i*(i-1)/2;
        
        if(j == 0 && i%1000 == 0) 
            std::cout << i << " / " << M << std::endl;

        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);  
        size_t closenes=size_t(cws_links_cy_dist(ff, ss, N) > close_v) + 
                        size_t(cws_links_cy_dist(sf, fs, N) > close_v);

        bool legal=true;  // is this pair legal with all boundary conditions	



        // both partners have to be different    
	//if(i == j)
        //    legal=false;

        // If properties from network creation are to be maintained special
        // checks are necessary
        if(maintain_p) {                
            if(!self_loop && (ff == ss || fs == sf))
                legal=false;

            if(!allow_multiple) {
                if(edges_m[ff*N+ss] || edges_m[sf*N+fs]) 
                    legal=false;

                if(!directed)
                    if(edges_m[ss*N+ff] || edges_m[fs*N+sf]) 
                        legal=false;
            }
        }

        if(!legal)
            nw_prob.add(0);
        else if(closenes == 0) 
            nw_prob.add(prob_clo_clo);
        else if(closenes == 1) 
            nw_prob.add(prob_far_clo);
        else
            nw_prob.add(prob_far_far);
    }


    for(size_t a=0; a<C; ++a) {  
        if(a%100 == 0) 
            std::cout << a << " / " << C << std::endl;
	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();


	// translate to edges
	size_t i=(size_t) (sqrt(2.0*count+0.5)+0.5);
	size_t j=count - i*(i-1)/2;
			
        couples.push_back(std::pair<size_t, size_t> (i, j) );
            
        // UPDATE
        nw_prob.set(count, 0); // Don't want to see this again
        
        // If we don't allow multiple we have to set all probabilities of couples
        // coinciding with the newly created couple to zero (may take a while)
        if(maintain_p && !allow_multiple) {
            size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);            
            
            // TODO comment
            std::set<size_t> sv;   // potential edges...
            sv.insert(edges_on_node[ff].begin(), edges_on_node[ff].end());
            sv.insert(edges_on_node[fs].begin(), edges_on_node[fs].end());
            sv.insert(edges_on_node[sf].begin(), edges_on_node[sf].end());
            sv.insert(edges_on_node[ss].begin(), edges_on_node[ss].end());


            std::set<size_t> su;   // potential couples  
            for(std::set<size_t>::iterator it1=sv.begin(); it1!=sv.end(); ++it1) {
                for(std::set<size_t>::iterator it2=edge_available.begin(); it2!=edge_available.end(); ++it2) {
                    size_t i1=(*it1);
                    size_t i2=(*it2);

                    if(i1 < i2)
                        su.insert(i1+i2*(i2-1)/2);
                    
                    if(i2 > i1)
                        su.insert(i2+i1*(i1-1)/2);
                }
            }

            

            // Now iterate through all couples which probability (can have been) changed
            for(std::set<size_t>::iterator it=su.begin(); it!=su.end(); ++it) {

	        size_t e1 = (size_t) (sqrt(2.0*(*it)+0.5)+0.5);
                size_t e2 = (*it) - e1*(e1-1)/2;
                size_t ff_(edges[e1].first), fs_(edges[e1].second), sf_(edges[e2].first), ss_(edges[e2].second);  
                     
                if(ff_ == ff && ss_ == ss || ff_ == sf && ss_ == fs ||
                   sf_ == sf && fs_ == fs || sf_ == ff && fs_ == ss)
                    nw_prob.set((*it), 0);

                if(!directed)   // 4 more cases
                    if(ff_ == ss && ss_ == ff || ff_ == fs && ss_ == sf ||
                       sf_ == fs && fs_ == sf || sf_ == ss && fs_ == ff)
                        nw_prob.set((*it), 0);
                
            }
        }

        for(size_t k=0; k<N; ++k) {
            edges_on_node[k].erase(i);
            edges_on_node[k].erase(j);
        }

    }
}


/*
 * For a Barabasi-Albert-Network this function this function selects edges to 
 * couple edges. This can be used when a mixed 1 -> 1 / 2->2 - reaction
 * network should be created with an Barabasi Albert statistics of the reactants 
 * graph. There are just 'C' pairs of links selected (without replacement).
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 *        
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * 
 * TODO: How should the probability of adding a certain couple be estimated.
 *       maximum product of "new edges" functionality or sum of products? 
 */
 
void couple_barabasi_albert(std::vector< std::pair<size_t, size_t> > &couples,
                            size_t C, std::vector< std::pair<size_t, size_t> > &edges,
                            bool maintain_p=false, bool allow_multiple=false, 
                           bool self_loop=false, bool directed = false) {
    // Calculate N and M 
    size_t M=edges.size();
    size_t N=0;

    for(size_t i=0; i<M; ++i) {
        if(edges[i].first+1 > N)
            N = edges[i].first+1;

        if(edges[i].second+1 > N)
            N = edges[i].second+1;
    }

    std::vector<bool> edges_m; // edges matrix (necessary if you want to maintain (not) allow multiple) 
    std::vector<std::set<size_t> > edges_on_node;    // list for every node (to which edges is it connected?)
    for(size_t i=0; i<N; ++i) 
        edges_on_node.push_back(std::set<size_t>());

    if(maintain_p) {
        for(size_t i=0; i<N*N; ++i) 
            edges_m.push_back(false);

        for(size_t i=0; i<edges.size(); ++i) {
            edges_m[edges[i].first*N+edges[i].second]=true;

            edges_on_node[edges[i].first].insert(i);
            edges_on_node[edges[i].second].insert(i);
        }
    }

    // Vector for functionality of nodes
    std::vector<size_t> fun_n;
    size_t fun_sum=0;

    for(size_t i=0; i<N; ++i)
        fun_n.push_back(0);

    // Which edges are still available vor coupling?
    std::set<size_t> edge_available;
    
    for(size_t j=0; j<M; ++j) {
        ++fun_n[edges[j].first];	
	++fun_n[edges[j].second];
	fun_sum += 2;
	edge_available.insert(j);
    }


    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<M*(M-1)/2; ++count) {    
        // translate to edge
	size_t i = (size_t) (sqrt(2.0*count+0.5)+0.5);
	size_t j = count - i*(i-1)/2;

        if(j == 0 && i%1000 == 0) 
            std::cout << i << " / " << M << std::endl;
        
        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);  

        bool legal=true;  // is this pair legal with all boundary conditions	



        // both partners have to be different     
	//if(i >= j)
        //    legal=false;

        // If properties from network creation are to be maintained special
        // checks are necessary
        if(maintain_p) {                
            if(!self_loop && (ff == ss || fs == sf))
                legal=false;

            if(!allow_multiple) {
                if(edges_m[ff*N+ss] || edges_m[sf*N+fs]) 
                    legal=false;

                if(!directed)
                    if(edges_m[ss*N+ff] || edges_m[fs*N+sf]) 
                        legal=false;
            }
        }

        if(!legal)
            nw_prob.add(0);
        else 
            nw_prob.add(fun_n[fs]*fun_n[ss]);
        
    }


    for(size_t a=0; a<C; ++a) {  	
        if(a%100 == 0) 
            std::cout << a << " / " << C << std::endl;

	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
	// translate to edges and nodes
	size_t i = (size_t) (sqrt(2.0*count+0.5)+0.5);
	size_t j = count - i*(i-1)/2;
        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);   
			
        couples.push_back(std::pair<size_t, size_t> (i, j) );

        // UPDATE
        nw_prob.set(count, 0); // Don't want to see this again


        // If we want to maintain probabilities for substrate graphs we have to update
        // "virtual" functionalities und change some couples probabilities...     
        if(maintain_p) {
            ++fun_n[ff];	
	    ++fun_n[fs];	
	    ++fun_n[sf];	
	    ++fun_n[ss];
	    fun_sum += 4;

            // TODO comment
            std::set<size_t> sv;   // potential edges...
            sv.insert(edges_on_node[ff].begin(), edges_on_node[ff].end());
            sv.insert(edges_on_node[fs].begin(), edges_on_node[fs].end());
            sv.insert(edges_on_node[sf].begin(), edges_on_node[sf].end());
            sv.insert(edges_on_node[ss].begin(), edges_on_node[ss].end());


            std::set<size_t> su;   // potential couples  
            for(std::set<size_t>::iterator it1=sv.begin(); it1!=sv.end(); ++it1) {
                for(std::set<size_t>::iterator it2=edge_available.begin(); it2!=edge_available.end(); ++it2) {
                    size_t i1=(*it1);
                    size_t i2=(*it2);

                    if(i1 < i2)
                        su.insert(i1+i2*(i2-1)/2);
                    
                    if(i2 > i1)
                        su.insert(i2+i1*(i1-1)/2);
                }
            }

            

            // Now iterate through all couples which probability (can have been) changed
            for(std::set<size_t>::iterator it=su.begin(); it!=su.end(); ++it) {

	        size_t k = (size_t) (sqrt(2.0*(*it)+0.5)+0.5);
                size_t l = (*it) - k*(k-1)/2;

                size_t ff_(edges[k].first), fs_(edges[k].second);
                size_t sf_(edges[l].first), ss_(edges[l].second);
                
                nw_prob.set(i, fun_n[fs_]*fun_n[ss_]);
                        
                if(!allow_multiple) {
                    // All couples which "virtual edges" are equal to the just created "virtual edges" are set to 0 
                    if(ff_ == ff && ss_ == ss || ff_ == sf && ss_ == fs ||
                       sf_ == sf && fs_ == fs || sf_ == ff && fs_ == ss) {
                        nw_prob.set((*it), 0);                        
                    }

                    if(!directed)   // 4 more cases for undirected networks
                        if(ff_ == ss && ss_ == ff || ff_ == fs && ss_ == sf ||
                           sf_ == fs && fs_ == sf || sf_ == ss && fs_ == ff) {
                            nw_prob.set((*it), 0);
                        }
                }        
            }
        }

        edge_available.erase(i);
        edge_available.erase(j);
        for(size_t k=0; k<N; ++k) {
            edges_on_node[k].erase(i);
            edges_on_node[k].erase(j);
        }
    } 					 
}



/*
 * For a pan-sinha-network... (similar to above)
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 *
 * h          - number of binary levels
 * m          - Number of modules on level 0 that are put in one module on
 *              level 1 (on higher levels it is allways to / binary tree)
 * r          - Change linking probability by going one level up (linking
 *              probability between two nodes is choosen from the lowest 
 *              level on which they are in the same module)
 *        
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void couple_pan_sinha(std::vector< std::pair<size_t, size_t> > &couples,
                      size_t C, std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t h, size_t m, double r,
                      bool maintain_p=false, bool allow_multiple=false, 
                      bool self_loop=false, bool directed = false) {
    // Calculate N and M 
    size_t M=edges.size();
    size_t N=0;
    for(size_t i=0; i<M; ++i) {
        if(edges[i].first+1 > N)
            N = edges[i].first+1;

        if(edges[i].second+1 > N)
            N = edges[i].second+1;
    }

    std::vector<bool> edges_m; // edges matrix (necessary if you want to maintain (not) allow multiple) 
    if(maintain_p && !allow_multiple) {
        for(size_t i=0; i<N*N; ++i) 
            edges_m.push_back(false);

        for(size_t i=0; i<edges.size(); ++i)
            edges_m[edges[i].first*N+edges[i].second]=true;
    }

    size_t n_el_mod=pow(2,h)*m;                   // number of elementary modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in el. modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 to each until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }
    
    // Creating lookup table for asociating every node on every level with one module
    std::vector<size_t> module_no;
    size_t module_no_lvl=h+2;
    for(size_t i=0; i<module_no_lvl; ++i) 
        for(size_t j=0; j<N; ++j) 
	    module_no.push_back(bps_module_no(elem_c, m, i, j));
    

    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<M*M; ++count) {    
        // translate to edge
	size_t i=count/M;
	size_t j=count%M;
        
        if(j == 0 && i%1000 == 0) 
            std::cout << i << " / " << M << std::endl;

        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);  

        bool legal=true;  // is this pair legal with all boundary conditions	


        // both partners have to be different     
	if(i == j)
            legal=false;

        // If properties from network creation are to be maintained special
        // checks are necessary
        if(maintain_p) {                
            if(!self_loop && (ff == ss || fs == sf))
                legal=false;

            if(!allow_multiple) {
                if(edges_m[ff*N+ss] || edges_m[sf*N+fs]) 
                    legal=false;

                if(!directed)
                    if(edges_m[ss*N+ff] || edges_m[fs*N+sf]) 
                        legal=false;
            }
        }

        if(!legal)
            nw_prob.add(0);
        else
            nw_prob.add( pow(r, bps_nodes_level_comp_lu(module_no, N, ff, ss))*
                         pow(r, bps_nodes_level_comp_lu(module_no, N, sf, fs))*100000 );
    }


    for(size_t a=0; a<C; ++a) {  	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
        if(a%100 == 0) 
            std::cout << a << " / " << C << std::endl;

	// translate to edges and nodes
	size_t i=count/M;
	size_t j=count%M;
        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);   
			
        couples.push_back(std::pair<size_t, size_t> (i, j) );
            
        // UPDATE
        nw_prob.set(count, 0); // Don't want to see this again


        // If we want to maintain probabilities for substrate graphs we have to update
        // "virtual" functionalities und change some couples probabilities...     
        if(maintain_p) {
            for(size_t k=0; k<M*M; ++k) {
	        size_t e1=k/M;
	        size_t e2=k%M;
                size_t ff_(edges[e1].first), fs_(edges[e1].second), sf_(edges[e2].first), ss_(edges[e2].second);  


                // If we don't allow multiple we have to set all probabilities of couples
                // coinciding with the newly created couple to zero (may take a while)

                if(!allow_multiple) 
                    if(nw_prob.get(k)) {     // Only have to think about couples with non-zero probability
                        // All couples which "virtual edges" are equal to the just created "virtual edges" are set to 0 
                        if(ff_ == ff && ss_ == ss || ff_ == sf && ss_ == fs ||
                           sf_ == sf && fs_ == fs || sf_ == ff && fs_ == ss)
                            nw_prob.set(k, 0);

                        if(!directed)   // 4 more cases for undirected networks
                            if(ff_ == ss && ss_ == ff || ff_ == fs && ss_ == sf ||
                              sf_ == fs && fs_ == sf || sf_ == ss && fs_ == ff)
                                nw_prob.set(k, 0);
                    }
            }
        }

    } 				  
}



/*
 * For simple modular networks... (similar to above)
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 * m          - Number of modules 
 * r          - Change linking probability by going one level up (linking
 *              probability between two nodes is choosen from the lowest 
 *              level on which they are in the same module)
 *        
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 */

void couple_simple_modular(std::vector< std::pair<size_t, size_t> > &couples,
                      size_t C, std::vector< std::pair<size_t, size_t> >& edges, 
                      size_t m, double r,
                      bool maintain_p=false, bool allow_multiple=false, 
                      bool self_loop=false, bool directed = false) {
    // Calculate N and M 
    size_t M=edges.size();
    size_t N=0;
    for(size_t i=0; i<M; ++i) {
        if(edges[i].first+1 > N)
            N = edges[i].first+1;

        if(edges[i].second+1 > N)
            N = edges[i].second+1;
    }

    std::vector<bool> edges_m; // edges matrix (necessary if you want to maintain (not) allow multiple) 
    if(maintain_p && !allow_multiple) {
        for(size_t i=0; i<N*N; ++i) 
            edges_m.push_back(false);

        for(size_t i=0; i<edges.size(); ++i)
            edges_m[edges[i].first*N+edges[i].second]=true;
    }

    size_t n_el_mod=m;                            // number of modules 
    size_t n_floor = floor(double(N)/n_el_mod);   // min no nodes in modules
    
    std::vector<size_t> elem_c;      // calculate the size of elementary nodes starting by n_floor

    {                                // and adding 1 to each until N is reached
        size_t N_tmp=N;
   
        for(size_t i=0; i<n_el_mod; ++i) { 
            elem_c.push_back(n_floor);
	        N_tmp -= n_floor;
        }

        for(size_t i=0; i<N_tmp; ++i) 
            ++elem_c[i];
    }
    
    // The binary tree probability distribution object is created and initialized
    bt_draw nw_prob;
    for(size_t count=0; count<M*M; ++count) {    
        // translate to edge
	size_t i=count/M;
	size_t j=count%M;
        
        if(j == 0 && i%1000 == 0) 
            std::cout << i << " / " << M << std::endl;

        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);  

        bool legal=true;  // is this pair legal with all boundary conditions	


        // both partners have to be different     
	if(i == j)
            legal=false;

        // If properties from network creation are to be maintained special
        // checks are necessary
        if(maintain_p) {                
            if(!self_loop && (ff == ss || fs == sf))
                legal=false;

            if(!allow_multiple) {
                if(edges_m[ff*N+ss] || edges_m[sf*N+fs]) 
                    legal=false;

                if(!directed)
                    if(edges_m[ss*N+ff] || edges_m[fs*N+sf]) 
                        legal=false;
            }
        }

        if(!legal)
            nw_prob.add(0);
        else
            nw_prob.add( pow(r, size_t(!bps_same_module(elem_c, ff, ss)))*
                         pow(r, size_t(!bps_same_module(elem_c, sf, fs)))*100000 );
    }


    for(size_t a=0; a<C; ++a) {  	
	// draw random number and find asocciated entry
	size_t count=nw_prob.draw();
		
        if(a%100 == 0) 
            std::cout << a << " / " << C << std::endl;

	// translate to edges and nodes
	size_t i=count/M;
	size_t j=count%M;
        size_t ff(edges[i].first), fs(edges[i].second), sf(edges[j].first), ss(edges[j].second);   
			
        couples.push_back(std::pair<size_t, size_t> (i, j) );
            
        // UPDATE
        nw_prob.set(count, 0); // Don't want to see this again


        // If we want to maintain probabilities for substrate graphs we have to update
        // "virtual" functionalities und change some couples probabilities...     
        if(maintain_p) {
            for(size_t k=0; k<M*M; ++k) {
	        size_t e1=k/M;
	        size_t e2=k%M;
                size_t ff_(edges[e1].first), fs_(edges[e1].second), sf_(edges[e2].first), ss_(edges[e2].second);  


                // If we don't allow multiple we have to set all probabilities of couples
                // coinciding with the newly created couple to zero (may take a while)

                if(!allow_multiple) 
                    if(nw_prob.get(k)) {     // Only have to think about couples with non-zero probability
                        // All couples which "virtual edges" are equal to the just created "virtual edges" are set to 0 
                        if(ff_ == ff && ss_ == ss || ff_ == sf && ss_ == fs ||
                           sf_ == sf && fs_ == fs || sf_ == ff && fs_ == ss)
                            nw_prob.set(k, 0);

                        if(!directed)   // 4 more cases for undirected networks
                            if(ff_ == ss && ss_ == ff || ff_ == fs && ss_ == sf ||
                              sf_ == fs && fs_ == sf || sf_ == ss && fs_ == ff)
                                nw_prob.set(k, 0);
                    }
            }
        }

    } 				  
}





/*
 * Creates an Erdos Renyi network with an stochastical algorithm with N nodes
 * and on average M edges. From M and the maximal number of edges a global 
 * probability for creating edges is calculated.
 */

void create_erdos_renyi_stoch(std::vector< std::pair<size_t, size_t> > &edges, 
                              size_t N, size_t M, bool self_loop=false, 
                              bool directed = false) {
    // For each case, calculate first the probability for a arbitrary link
    // of existing d, so that the expectation value for number of links is M
    // and then iterate over all possible links

    if(!self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5);  
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=i+1; j<N; ++j)
	            if(double(rand())/RAND_MAX < d)
		            edges.push_back(std::pair<size_t, size_t> (i, j) );        
    } else if (self_loop && !directed) {
        double d=double(M)/(N*(N-1)*0.5+N);  
      
        for(size_t i=0; i<N; ++i) 
	        for(size_t j=i; j<N; ++j)
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



#endif // NETWORK_TOOLS_H
