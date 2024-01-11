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
            ) : idxPrm(_idx), Ms(0), refNode(_p_node)
        {
        ind.assign(_i);
        P.setZero();
        Lp.setZero();
        }

    /** indices to the nodes */
    std::vector<int> ind;
    
    /** index of the material parameters of the element */
    int idxPrm;
    
    /** magnetization at saturation of the element */
    double Ms;
    
    /** weights hat function of the element */
    Eigen::Vector<double,NPI> weight;

    /** matrix for integrales */
    Eigen::Matrix<double,2*N,2*N> Kp;

    /** vector for integrales */
    Eigen::Vector<double,2*N> Lp;

    /** block diagonal matrix for projections */
    Eigen::Matrix<double,2*N,3*N> P;

    /** getter for N */
    inline constexpr int getN(void) const { return N; }

    /** getter for NPI */
    inline constexpr int getNPI(void) const { return NPI; }

    /** set all indices to zero */
    inline void indicesToZero(void)
        { std::fill(ind.begin(),ind.end(),0); }

    /** build matrix P, must be refreshed when ep,eq are changing (setBasis) */
    void buildMatP(void)
        {
        for (int i = 0; i < N; i++)
            {
            const Eigen::Vector3d &ep = refNode[ind[i]].ep;
            P(i,i) = ep.x();
            P(i,N + i) = ep.y();
            P(i,2 * N + i) = ep.z();

            const Eigen::Vector3d &eq = refNode[ind[i]].eq;
            P(N + i,i) = eq.x();
            P(N + i,N + i) = eq.y();
            P(N + i,2 * N + i) = eq.z();
            }
        }

    /** assemble the big sparse matrix K from tetra or facette inner matrix Kp */
    void assemblage_mat(const int NOD /**< [in] nb nodes */,
                        std::vector<Eigen::Triplet<double>> &K /**< [out] COO matrix */ ) const
        {
        for (int i = 0; i < N; i++)
            {
            int i_ = ind[i];

            for (int j = 0; j < N; j++)
                {
                int j_ = ind[j];
                if(Kp(i,j) != 0) K.emplace_back( NOD + i_, j_, Kp(i,j) );
                if (Kp(i,N + j) != 0) K.emplace_back( NOD + i_, NOD + j_, Kp(i,N + j) );
                if (Kp(N + i,j) != 0) K.emplace_back( i_, j_, Kp(N + i,j) );
                if (Kp(N + i,N + j) != 0) K.emplace_back( i_, NOD + j_, Kp(N + i,N + j) );
                }
            }
        }

    /** assemble the big vector L from tetra or facette inner vector Lp */
    void assemblage_vect(const int NOD /**< [in] nb nodes */,
                        Eigen::Ref<Eigen::VectorXd> L /**< [out] vector */) const
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
    
    /** computes Gauss point of the element, return in result */
    virtual void getPtGauss(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result) const = 0;

    protected:
        /** vector of nodes */
        const std::vector<Nodes::Node> &refNode;

        /** zeroBasing: index convention Matlab/msh (one based) -> C++ (zero based) */
        inline void zeroBasing(void)
            { for(int i=0;i<N;i++) {ind[i]--;} }

    private:
        /** a method to orientate the element */
        virtual void orientate() = 0;
    };
    
#endif
