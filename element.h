#ifndef element_h
#define element_h

#include <execution>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "node.h"

/** \class element
template class, mother class for tetraedrons and facettes
*/

template <int N,int NPI>
class element
    {
    /** constructor */
    public:
    element(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */,
            const int _idx /**< index to params */,
            std::initializer_list<int> & _i /**< indices to the nodes */
            ) : idxPrm(_idx), refNode(_p_node)
        {
        ind.assign(_i);
        }

    /** indices to the nodes */
    std::vector<int> ind;
    
    /** index of the material parameters of the element */
    int idxPrm;
    
    /** matrix for integrales */
    double Kp[2*N][2*N];

    /** vector for integrales */
    double Lp[2*N];

    /** getter for N */
    inline constexpr int getN(void) const { return N; }

    /** getter for NPI */
    inline constexpr int getNPI(void) const { return NPI; }

    /** getter for node (only used by unit-tests ; not recommended) */
    inline const Nodes::Node &getNode(const int i) { return refNode[ind[i]]; }
    
    /** set all indices to zero */
    inline void indicesToZero(void)
        { std::fill(ind.begin(),ind.end(),0); }

    /** function to provide P matrix(2N,3N) coefficients, with respect to its block diagonal structure */
    double Pcoeff(const int i, const int j)
        {
        double val(0);
        int node_i = i % N;

        if (node_i == (j % N))
            {
            if (i < N)
                {
                val = refNode[ind[node_i]].ep(j / N);
                }
            else
                {
                val = refNode[ind[node_i]].eq(j / N);
                }
            }
        return val;
        }

/** make projection for tetra or facette. It computes Bp = P*B and stores result in inner vector Lp */
    void projection_vect(Pt::pt3D *B)
        {
        for (int i = 0; i < (2 * N); i++)
            {
            Lp[i] = 0;
            for (int k = 0; k < N; k++)
                {
                Lp[i] += Pt::pt3D(Pcoeff(i,k), Pcoeff(i,N + k), Pcoeff(i,2*N + k)).pScal(B[k]);
                }
            }
        }

/** make projection for tetra or facette. It computes Ap = (P*A)*trans(P) and stores result in inner matrix Kp */
    void projection_mat(double (&A)[3 * N][3 * N])
        {
        double PA[2 * N][3 * N];  // no need to initialize with zeros
        for (int i = 0; i < (2 * N); i++)
            {
            for (int k = 0; k < (3 * N); k++)
                {
                PA[i][k] = 0;
                for (int j = 0; j < (3 * N); j++)
                    {
                    PA[i][k] += Pcoeff(i, j) * A[j][k];
                    }
                }
            }

        for (int i = 0; i < (2 * N); i++)
            for (int k = 0; k < (2 * N); k++)
                {
                Kp[i][k] = 0;
                for (int j = 0; j < (3 * N); j++)
                    {
                    Kp[i][k] += PA[i][j] * Pcoeff(k, j);
                    }
                }
        }

    /** assemble the big sparse matrix K from tetra or facette inner matrix Kp */
    void assemblage_mat(const int NOD, std::vector<Eigen::Triplet<double>> &K) const
        {
        for (int i = 0; i < N; i++)
            {
            int i_ = ind[i];

            for (int j = 0; j < N; j++)
                {
                int j_ = ind[j];
                
                if(Kp[i][j] != 0)
                    { K.push_back(Eigen::Triplet<double>(NOD + i_, j_, Kp[i][j])); }
                
                if (Kp[i][N + j] != 0)
                    { K.push_back(Eigen::Triplet<double>(NOD + i_, NOD + j_, Kp[i][N + j])); }
                
                if (Kp[N + i][j] != 0)
                    { K.push_back(Eigen::Triplet<double>(i_, j_, Kp[N + i][j])); }
                
                if (Kp[N + i][N + j] != 0)
                    { K.push_back(Eigen::Triplet<double>(i_, NOD + j_, Kp[N + i][N + j])); }
                }
            }
        }

    /** assemble the big vector L from tetra or facette inner vector Lp */
    void assemblage_vect(const int NOD, Eigen::Ref<Eigen::VectorXd> L) const
        {
        for (int i = 0; i < N; i++)
            {
            const int i_ = ind[i];
            if(Lp[i] != 0)
                { L(NOD + i_) += Lp[i]; }
            if(Lp[N+i] != 0)
                { L(i_) += Lp[N + i]; }
            }
        }

    /** info: print node indices of the element and the vector index of the associated param */
    void infos() const
        {
        std::cout << "idxPrm: " << idxPrm << " ind: (";
        for(unsigned int i = 0; i < N-1; i++)
            { std::cout << ind[i] << ", "; }
        std::cout << ind[N-1] << ")\n";
        };
    
    protected:
        /** vector of nodes */
        const std::vector<Nodes::Node> &refNode;

        /** zeroBasing: index convention Matlab/msh (one based) -> C++ (zero based) */
        inline void zeroBasing(void)
            { for(int i=0;i<N;i++) {ind[i]--;} }
    
    private:
        /** a method to orientate the element must be provided in derived class */
        virtual void orientate() = 0;
    };
    
#endif
