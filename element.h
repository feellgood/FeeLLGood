#ifndef element_h
#define element_h

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <execution>
#pragma GCC diagnostic pop

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "node.h"
#include "algebra/algebra.h"

/** \class element
\brief Template abstract class, mother class for tetraedrons and facettes.

template parameters are N number of sommits and NPI number of interpolation points. It contains a list of indices to the N nodes of the element, a reference to the full nodes vector, and index refering to the associated material parameters. All indices are zero based, derived class constructor should call zerobasing() if needed. It contains also the vector and matrix weight, Kp, Lp, P related to finite element computations. weight values are not initialized, they have to be set by derived class constructor.
Member function getPtGauss() returns Gauss points.
orientate() is a pure virtual function, it should manipulate indices to orientate positively the element.
*/

template <int N,int NPI>
class element
    {
    /** constructor */
    public:
    explicit element(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */,
            const int _idx /**< index to params */,
            std::initializer_list<int> & _i /**< indices to the nodes */
            ) : idxPrm(_idx), refNode(_p_node)
        {
        if(_i.size() == N)
            { ind.assign(_i); }
        else
            {
            std::cout<<"Warning: element constructor is given an init list with size() != N\n";
            }
        Kp.setZero();
        Lp.setZero();
        }

    /** indices to the nodes */
    std::vector<int> ind;
    
    /** index of the material parameters of the element */
    int idxPrm;
    
    /** weights hat function of the element */
    Eigen::Matrix<double,NPI,1> weight;

    /** matrix for integrales */
    Eigen::Matrix<double,2*N,2*N> Kp;

    /** vector for integrales */
    Eigen::Matrix<double,2*N,1> Lp;

    /** getter for N */
    inline constexpr int getN(void) const { return N; }

    /** getter for NPI */
    inline constexpr int getNPI(void) const { return NPI; }

/** build matrix P direcly from ep,eq in nodes
* P is block diagonal:
( Epx Epy Epz )
( Eqz Eqy Eqz )
with each block E(p|q)(x|y|z) a N*N diagonal matrix
* see here http://eigen.tuxfamily.org/dox-devel/TopicTemplateKeyword.html for the wierd template syntax
*/
    void buildMatP(Eigen::Ref<Eigen::Matrix<double,2*N,3*N>> P /**< [out] block diagonal matrix */)
        {
        P.setZero();
        Eigen::Matrix<double,Nodes::DIM,N,Eigen::RowMajor> tempo;
        for (int i = 0; i < N; i++) { tempo.col(i) = getNode(i).ep; }

        P.template block<N,N>(0,0).diagonal() = tempo.row(Nodes::IDX_X);
        P.template block<N,N>(0,N).diagonal() = tempo.row(Nodes::IDX_Y);
        P.template block<N,N>(0,2*N).diagonal() = tempo.row(Nodes::IDX_Z);

        for (int i = 0; i < N; i++) { tempo.col(i) = getNode(i).eq; }

        P.template block<N,N>(N,0).diagonal() = tempo.row(Nodes::IDX_X);
        P.template block<N,N>(N,N).diagonal() = tempo.row(Nodes::IDX_Y);
        P.template block<N,N>(N,2*N).diagonal() = tempo.row(Nodes::IDX_Z);
        }

    /** computes Gauss point of the element, return in result */
    virtual void getPtGauss(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result) const = 0;

    /** print node indices of the element and the vector index of the associated param */
    void infos(void)
        {
        std::cout << "idxPrm: " << idxPrm << " ind: {";
        for(unsigned int i = 0; i < N-1; i++)
            { std::cout << ind[i] << ": " << refNode[ind[i]].p << std::endl; }
        std::cout << ind[N-1] <<": " << refNode[ind[N-1]].p << "}\n";
        }

    protected:
        /** returns reference to node at ind[i] from mesh node vector */
        inline const Nodes::Node & getNode(const int i) const { return refNode[ind[i]]; }

        /** returns true if mesh node vector is not empty */
        inline bool existNodes(void) const { return (refNode.size() > 0); }

        /** zeroBasing: index convention Matlab/msh (one based) -> C++ (zero based) */
        inline void zeroBasing(void)
            { std::for_each(ind.begin(),ind.end(),[](int & _i){ _i--; } ); }

    private:
        /** vector of nodes */
        const std::vector<Nodes::Node> & refNode;

        /** a method to orientate the element */
        virtual void orientate() = 0;
    };
#endif
