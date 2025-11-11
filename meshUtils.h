#ifndef MESHUTILS_H
#define MESHUTILS_H

#include <gmsh.h>

/**
convenient template to check elements type while reading a gmsh mesh file
*/
template <int TYPE>
bool all_elems_are(std::vector<int> &container)
    { return std::all_of(container.begin(), container.end(), [](int &e_type){return e_type == TYPE;} ); }

/**
convenient template to check named regions(volume or surface) from feellgood settings exist in physical groups gmsh mesh inner structure built from file
*/
template <class T>
bool checkNamedObjects(std::vector<T> const &v_prm, const int dim_obj)
    {
    bool existRegions(true);
    std::vector<std::pair<int, int> > physGroups;
    gmsh::model::getPhysicalGroups(physGroups,dim_obj);
    std::vector<std::string> o_names;
    std::for_each(physGroups.begin(),physGroups.end(),[&o_names]( std::pair<int, int> &pGroup)
        {
        std::string physGroupName;
        gmsh::model::getPhysicalName(pGroup.first, pGroup.second, physGroupName);
        o_names.push_back(physGroupName);
        });

    std::for_each(v_prm.begin(),v_prm.end(),[&existRegions,&o_names](auto p)
        {
        if (p.regName != "__default__")
            {
            if (!std::any_of(o_names.begin(),o_names.end(),
                [&p] (std::string &elem) { return p.regName == elem; } ) )
                {
                std::cout << "Fatal Error: region " << p.regName << " not found.\n";
                existRegions = false;
                }
            } 
        } );
    return existRegions;
    }

/** order and suppress copies in the index input list */
template<typename T>
void suppress_copies(std::vector<T> &v_idx)
    {
    std::sort(v_idx.begin(), v_idx.end());
    v_idx.resize( std::distance(v_idx.begin(), std::unique(v_idx.begin(), v_idx.end())) );
    v_idx.shrink_to_fit();
    }

#endif
