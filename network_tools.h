#ifndef NETWORK_TOOLS_H
#define NETWORK_TOOLS_H

#include <vector>
#include <cstdlib>
#include <gmp.h>
#include <iml.h>
#include <iostream>



void create_barabasi_albert(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t m, size_t m0, bool allow_multiple=false) {
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
        //std::cout << "no_nodes=" << no_nodes << std::endl;
	size_t sum=0;
	for(size_t i=0; i<no_nodes; ++i) {
	    sum += degree[i];
	    available[i] = true;    
	}
	
	//std::cout << "sum=" << sum << std::endl;
	
	for(size_t j=0; j<m; ++j) {
	     //std::cout << "BLABLA j=" << j << std::endl;
	     double rnd=double(rand())/RAND_MAX;
	     size_t next=0;
	     
	     while(rnd > double(degree[next])/ sum || !available[next]) {
	         //std::cout << "next=" << next << "  rnd=" << rnd << "  sum=" << sum << std::endl;
		 //for(size_t i=0; i<N; ++i)
		 //    std::cout << i << "   " << degree[i] << "   " << available[i] << std::endl;
		
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
	     
	     //std::cout << "next=" << next << "  rnd=" << rnd << "  sum=" << sum << std::endl;
	     //for(size_t i=0; i<N; ++i)
             //   std::cout << i << "   " << degree[i] << "   " << available[i] << std::endl;
	}	      
    } 
	    
    if(no_nodes != N) {
        std::cout << "ERROR: Something went wrong!" << std::endl;  
    }
}



/*
 * Creates a simple erdos-renyi-network with n nodes and M edges. Edges are written
 * to the vector given as first parameter. Edges are indexed starting with zero. If 
 * allow_multiple is set true self loops and the multiple occurance of the same 
 * link is allowed.
 */

void create_erdos_renyi(std::vector< std::pair<size_t, size_t> > &edges, size_t n, size_t M, bool allow_multiple=false) {
    for(size_t j=edges.size(); j<M; ++j) {
        bool found=false;
	
	do {  
	    size_t first=rand()%n;  // Create node indices for link j randomly
            size_t second=rand()%n;
		
	    // If applicable check for self loops and multiple occorance
	    if(allow_multiple) {
	        found = true;	 
		edges.push_back(std::pair<size_t, size_t> (first, second) );
	    } else if(first != second) {
		if(find(edges.begin(), edges.end(), std::pair<size_t, size_t>(first, second)) == edges.end() &&
		   find(edges.begin(), edges.end(), std::pair<size_t, size_t>(second, first)) == edges.end()) {
		    edges.push_back(std::pair<size_t, size_t> (first, second) );
		    found = true;
		}
	    }
	} while (!found);
    }
}


void create_watts_strogatz(std::vector< std::pair<size_t, size_t> > &edges, size_t N, size_t K, double beta, bool allow_multiple=false) {
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