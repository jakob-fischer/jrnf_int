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
#include <iml.h>
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




     void printout() {
       for(size_t i=0; i<=last; ++i) {
           std::cout << i << " : ";
           for(size_t j=0; j<probs[i].size(); ++j)
               std::cout << probs[i][j] << " / ";

           std::cout << std::endl;
       }

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
 *
 * TODO:      check!
 */

// Macro helper to determine "distance" 
size_t cws_links_cy_dist(size_t a, size_t b, size_t N) {
    size_t a_amb = std::abs(a-b);
    size_t a_amb2 = std::abs(N-a_amb);

    return std::min(a_amb, a_amb2);
}

 
void couple_watts_strogatz(std::vector< std::pair<size_t, size_t> > &couples,
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


    // Which edges are still available vor coupling?
    std::set<size_t> edge_available;
    
    for(size_t j=0; j<M; ++j)
	edge_available.insert(j);


    size_t close_v = ceil(double(M)/N); // needed to count which additional "connections" are far

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
    size_t far_links = 0;
    for(size_t i=0; i<edges.size(); ++i) 
        if(cws_links_cy_dist(edges[i].first, edges[i].second, N) > close_v)
            ++far_links;
    
    size_t prob_far_far=100000.0*double(far_links)/N*double(far_links)/N;
    size_t prob_clo_clo=100000.0*(1-double(far_links)/N)*(1-double(far_links)/N);
    size_t prob_far_clo=100000.0*double(far_links)/N*(1-double(far_links)/N);

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
 * 
 * TODO: check!
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
 * For a pan-sinha-network...
 *
 * couples    - reference to a vector of pairs to that the selected couples are
                written to 0-based 
 * C          - number of couples that are drawn
 * edges      - pairs vector containing the original networks edges
 *
 * h          -
 * m          -
 * r          -
 *        
 * maintain_p - maintain properties of original network in substrate graph
 *              (properties to maintain have to be given in following parameters)
 * allow_multiple  - the same link can occur multiple times
 * self_loop       - self loops are allowed
 * directed        - directed network is created
 * 
 * TODO: check!
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





/* bad written class for simple bidirectional network analysis. not 
 * developed and propperly commented at the moment because R and igraph 
 * can do better.
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
