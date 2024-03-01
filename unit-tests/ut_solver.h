template<int DIM,int NbRand>
void build_coeffs( std::vector<Eigen::Triplet<double>> &coeffs ,Eigen::Ref<Eigen::VectorXd> v)
    {
    for(int i=0;i<DIM;i++)
        {
        coeffs.push_back(Eigen::Triplet<double>(i,i,1.0) );
        v(i) = 1.0;
        }

    std::mt19937 gen(my_seed());
    std::uniform_int_distribution<> distrib(0, DIM-1);

    for (int nb=0; nb<NbRand; nb++)
        {
        int i = distrib(gen);
        int j = distrib(gen);
        coeffs.push_back(Eigen::Triplet<double>(i,j,1.0) );
        }
    }

template<int NOD,int NbRand>
        void loc_build_coeffs(std::vector<Eigen::Triplet<double>> &coeffs ,Eigen::Ref<Eigen::VectorXd> v,const double t)
            {
            build_coeffs<NOD,NbRand>(coeffs,v);
            
            std::mt19937 gen(my_seed());
            std::uniform_int_distribution<> distrib(0, NOD-1);
            
            std::for_each(coeffs.begin(),coeffs.end(),
                [&distrib,&gen]( Eigen::Triplet<double> &c)
                    {
                    if (c.col()==c.row()) 
                        {
                        double val = c.value() + distrib(gen)/(0.1+NOD);
                        c =  Eigen::Triplet(c.col(),c.row(),val); // a Triplet is a read only object, no possibility to alter its value once constructed
                        }
                    }
                );
            
            for(int i=0;i<NOD;i++)
                {
                double val = distrib(gen)/(0.1+NOD);
                v(i) += sin(val);
                }
         
            int i = distrib(gen);
            int j = distrib(gen);
            coeffs.push_back(Eigen::Triplet<double>(i,j,cos(t)));
            }

template <int NOD>
class DummyLinAlgebra
    {
    public:
        inline DummyLinAlgebra(const int _MAX_ITER,const double _TOL): MAX_ITER(_MAX_ITER), S_TOL(_TOL) {}
    
        int solve(double &t)
            {
            Eigen::VectorXd x(NOD), b(NOD);
            Eigen::SparseMatrix<double,Eigen::RowMajor> A(NOD,NOD);
            std::vector<Eigen::Triplet<double>> coeffs;

            loc_build_coeffs<NOD,30000>(coeffs,b,t);
            
            A.setFromTriplets(coeffs.begin(),coeffs.end());
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::IncompleteLUT<double> > _solver;
            _solver.setMaxIterations(MAX_ITER);
            _solver.setTolerance(S_TOL);
            _solver.compute(A);
            x = _solver.solveWithGuess(b, Eigen::VectorXd::Constant(NOD,2.0));
            std::cout << "#iterations:     " << _solver.iterations() << std::endl;
            std::cout << "estimated error: " << _solver.error()      << std::endl;
            BOOST_CHECK( _solver.iterations() >= 1);
            BOOST_CHECK( (x-b).norm() != 0.0 );
            double errorCheck = ((A*x)-b).norm();
            std::cout << "error = " << errorCheck << std::endl;
            BOOST_CHECK( errorCheck*errorCheck <= S_TOL );
            
            t += 0.01;
            return 0;
            }
    
    private:
        const int MAX_ITER;
        const double S_TOL;    
    }; // end class DummyLinAlgebra
