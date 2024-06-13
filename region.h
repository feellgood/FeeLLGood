#include <iterator>
#include <iostream>
#include <vector>
#include <string>

/** specifiy what is the problem to solve in the region */
const int LLG = 0x001;
const int DIFF_SPIN = 0x002;
const int ELECTROSTAT = 0x004;


/** \class region
region name, index and material constants
*/

template<class T,class PRM>
class region
    {
    public:
        /** constructor */
        inline region(const std::string regName,const int index,const int wts): name(regName), idx(index), whatToSolve(wts) {}

        /** iterator pointing to first element of the region */
        typename std::vector<T>::iterator first;

        /** iterator pointing to last element of the region */
        typename std::vector<T>::iterator last;

        /** print the region parameters */
        inline void infos() const
            {
            std::cout << "region name = " << name << "; what to solve: ";
            switch(whatToSolve)
                {
                case (LLG): std::cout << "\tLLG\n"; break;
                case (LLG | DIFF_SPIN): std::cout << "\tLLG + DIFF_SPIN\n"; break;
                case (DIFF_SPIN): std::cout << "\tDIFF_SPIN\n"; break;
                case (ELECTROSTAT): std::cout << "\tELECTROSTAT\n"; break;
                case (LLG | ELECTROSTAT): std::cout << "\tLLG + ELECTROSTAT\n"; break;
                default:
                    std::cout << "pb undefined: " << whatToSolve << std::endl;
                break;
                }
            };

        /** parameters associated to T */
        PRM p;

        /** getter: returns region name */
        inline std::string getName(void) { return name; }

        /** getter: returns region index */
        inline int getIndex(void) { return idx; }
        
    private:
        /** region name */
        const std::string name;

        /** index */
        const int idx;
        
        /** what problem(s) to solve in the region */
        const int whatToSolve;
    };
