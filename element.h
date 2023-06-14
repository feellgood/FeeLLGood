#ifndef element_h
#define element_h

#include "node.h"

template <int N,int NPI>
class element
    {
    public: element(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */) : refNode(_p_node)
        {
        ind.resize(N);
        }
    
    /** indices to the nodes */
    std::vector<int> ind;
    
    double Kp[2*N][2*N];
    double Lp[2*N];
    
    inline void set_ind(const int i, const int idx)
        { ind[i] = idx-1; } // index convention Matlab/msh (one based) -> C++ (zero based)
    
    /** getter for N */
    inline constexpr int getN(void) const { return N; }

    /** getter for NPI */
    inline constexpr int getNPI(void) const { return NPI; }

    /** getter for node */
    inline const Nodes::Node &getNode(const int i)
        { return refNode[ind[i]]; }
    
    protected:
        /** vector of nodes */
        const std::vector<Nodes::Node> &refNode;
    
    };
    
#endif
