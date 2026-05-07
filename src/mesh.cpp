#include "mesh.h"
#include "meshUtils.h"
#include <utility>

using namespace Mesh;

void mesh::infos(void) const
    {
    std::cout << "mesh:\n";
    std::cout << "  nodes:              " << getNbNodes() << '\n';
    std::cout << "  faces:              " << fac.size() << '\n';
    std::cout << "  tetraedrons:        " << tet.size() << '\n';
    std::cout << "  total volume:       " << vol << '\n';
    }

double mesh::thiele(const int region /** region index, or -1 for all magnetic regions */) const
    {
    double sum = std::transform_reduce(EXEC_POL,magTet.begin(),magTet.end(),0.0,std::plus<>(),
            [this, region](const int idxElem)
                {//the lambda computes sq(sum(grad(u))) on the tetrahedron indexed by idxElem
                double val(0);
                Tetra::Tet const &te = tet[idxElem];
                if((te.idxPrm == region) || (region == -1))
                    {
                    Eigen::Matrix<double,Nodes::DIM,Tetra::N> u_nod;
                    for (int i = 0; i< Tetra::N; i++)
                        u_nod.col(i) = getNode_u(te.ind[i]);
                
                    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> du_dz = u_nod * (te.da.col(Nodes::IDX_Z)).replicate(1,Tetra::NPI);
                
                    for(int npi=0;npi<Tetra::NPI;npi++)
                        val += du_dz.col(npi).dot(du_dz.col(npi)) * te.weight[npi];
                    }
                return val;
                });
    
    double volume_u(0.0);
    // to check: that sum seems to be the volume of the region, if so we should not recompute it 
    std::for_each(EXEC_POL,magTet.begin(),magTet.end(),
            [this,&volume_u,region](const int idxElem)
                {
                Tetra::Tet const &te = tet[idxElem];
                if((te.idxPrm == region) || (region == -1))
                    { volume_u += te.weight.sum(); }
                });
    double cross_section = volume_u/l.z();//average cross-section
    return 2.0*cross_section/sum;
    }


double mesh::avg(const std::function<double(Nodes::Node, Nodes::index)>& getter /**< [in] */,
                 Nodes::index d /**< [in] */,
                 int region /**< region index, or -1 for all magnetic regions */) const
    {
    double sum = std::transform_reduce(EXEC_POL, magTet.begin(), magTet.end(), 0.0,
                                       std::plus<>(),
                                       [this, getter, &d, region](const int &idxElem)
                                       {
                                       Tetra::Tet const &te = tet[idxElem];
                                       if (te.idxPrm != region && region != -1)
                                           return 0.0;
                                       Eigen::Matrix<double,Tetra::NPI,1> val;
                                       te.interpolation(getter, d, val);
                                       return te.weight.dot(val);
                                       });
    double volume = (region == -1) ? totalMagVol : paramTetra[region].volume;
    return sum / volume;
    }

double mesh::doOnNodes(const double init_val, const Nodes::index coord,
                     const std::function<bool(double, double)>& whatToDo) const
    {
    double result(init_val);
    std::for_each(node.begin(), node.end(),
                  [&result, coord, whatToDo](Nodes::Node const &n)
                  {
                  double val = n.p(coord);
                  if (whatToDo(val, result)) result = val;
                  });
    return result;
    }

void mesh::indexReorder()
    {
    std::set<Facette::Fac> sf;  // implicit use of operator< redefined in class Fac

    std::for_each(tet.begin(), tet.end(), [this,&sf](Tetra::Tet const &te)
                  {
                  const int ia = te.ind[0];
                  const int ib = te.ind[1];
                  const int ic = te.ind[2];
                  const int id = te.ind[3];

                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ic, ib} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ib, ic, id} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, id, ic} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ib, id} ));
                  });  // end for_each

    std::for_each(fac.begin(), fac.end(), [this, &sf](Facette::Fac &fa)
                  {
                  int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
                  auto it = sf.end();
                  for (int perm = 0; perm < 2; perm++)
                      {
                      for (int nrot = 0; nrot < 3; nrot++)
                          {
                          Facette::Fac fc(node, 0, 0, {0, 0, 0} );
                          fc.ind[(0 + nrot) % 3] = i0;
                          fc.ind[(1 + nrot) % 3] = i1;
                          fc.ind[(2 + nrot) % 3] = i2;
                          it = sf.find(fc);
                          if (it != sf.end()) break;
                          }

                      if (it != sf.end())
                          {  // found
                          Eigen::Vector3d p0p1 = node[it->ind[1]].p - node[it->ind[0]].p;
                          Eigen::Vector3d p0p2 = node[it->ind[2]].p - node[it->ind[0]].p;
                                
                          // fa.Ms will have the magnitude of first arg of copysign, with second arg
                          // sign
                          fa.dMs = std::copysign( paramTetra[it->idxPrm].Ms,
                               p0p1.dot(p0p2.cross(fa.calc_norm())) );  // carefull, calc_norm
                                                                        // computes the normal to
                                                                        // the face before idx swap
                          }
                      std::swap(i1, i2);  // it seems from ref archive we do not want to swap
                                        // inner fac indices but local i1 and i2
                      fa.n = fa.calc_norm(); // update normal vector
                      }                   // end perm
                  });  // end for_each
    }

void mesh::sortNodes(Nodes::index long_axis)
    {
    // Sort the nodes along this axis, indirectly through an array of indices.
    std::vector<int> permutation(node.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(),
              [this, long_axis](int a, int b)
              { return node[a].p(long_axis) < node[b].p(long_axis); });
    node_index.resize(node.size());
    for (size_t i = 0; i < node.size(); i++)
        node_index[permutation[i]] = i;

    // Actually sort the array of nodes.
    std::vector<Nodes::Node> node_copy(node);
    for (size_t i = 0; i < node.size(); i++)
        node[i] = node_copy[permutation[i]];

    // Update the indices stored in the elements.
    std::for_each(tet.begin(), tet.end(),
                  [this](Tetra::Tet &tetrahedron)
                  {
                      tetrahedron.ind[0] = node_index[tetrahedron.ind[0]];
                      tetrahedron.ind[1] = node_index[tetrahedron.ind[1]];
                      tetrahedron.ind[2] = node_index[tetrahedron.ind[2]];
                      tetrahedron.ind[3] = node_index[tetrahedron.ind[3]];
                  });
    std::for_each(fac.begin(), fac.end(),
                  [this](Facette::Fac &facette)
                  {
                      facette.ind[0] = node_index[facette.ind[0]];
                      facette.ind[1] = node_index[facette.ind[1]];
                      facette.ind[2] = node_index[facette.ind[2]];
                  });
    }

void mesh::checkFacettes() const 
    {
    std::vector<Facette::Fac> allFacCtnr;
    std::vector<std::pair<int, bool>> twoTetFacCtnr;    // Bool: if the facettes are 
                                                        // in the same region or not. 
    const int tetraN = 4;

    // Put all facettes into allFacCtnr
    std::for_each(tet.begin(), tet.end(),   // For each tetrahedron
            [this, &allFacCtnr, &twoTetFacCtnr](const Tetra::Tet &tetrahedron) 
            {
            for (int i = 0; i < tetraN; i++)  // For each 4 facettes
                {   
                Facette::Fac curFac(node, 0, tetrahedron.idxPrm, {tetrahedron.ind[i], 
                                 tetrahedron.ind[(i + 1) % tetraN], 
                                 tetrahedron.ind[(i + 2) % tetraN]});

                // Search for an equal facette, insert at the end if not found.
                int j = 0;
                while (j < allFacCtnr.size() && !(allFacCtnr[j] == curFac))
                    { j++; }
                if (j == allFacCtnr.size())
                    { allFacCtnr.push_back(curFac); } 
                else 
                    {   // If an equal facette is found
                    int k = 0;
                    while (k < twoTetFacCtnr.size() && twoTetFacCtnr[k].first != j)
                        { k++; }

                    if (k < twoTetFacCtnr.size())
                        {   // If one facette belongs to 3 tetrahedrons, throws an error
                        std::cout << "Error : wrong mesh generation";
                        exit(1);
                        }
                    else
                        {   // If only 2 tetrahedrons, updates twoTetFacCtnr
                        twoTetFacCtnr.push_back({j, (curFac.idxPrm == allFacCtnr[j].idxPrm)});
                        }
                    }
                }          
            });
    
    std::sort(twoTetFacCtnr.begin(), twoTetFacCtnr.end());
    
    int j = 0;
    for (int i = 0; i < allFacCtnr.size(); i++)
        {   // For each facette
        if (j < twoTetFacCtnr.size() && i == twoTetFacCtnr[j].first && twoTetFacCtnr[j].second)
            {   // If the current one belongs to two tetrahedrons of the same region
            int k = 0;
            while (k < fac.size() && !(fac[k] == allFacCtnr[i]))
                { k++; }
            if (k < fac.size())
                {
                std::cout << "Error : an internal facette belonging "
                                "to a surface region has been found";
                exit(1);
                }
            j++;
            }
        else if (j < twoTetFacCtnr.size() && i == twoTetFacCtnr[j].first && !twoTetFacCtnr[j].second)
            {   // If the current one belongs to two tetrahedrons of differents regions
            int k = 0;
            while (k < fac.size() && !(fac[k] == allFacCtnr[i]))
                { k++; }
            if (k < fac.size())
                {   // Looking for a second identical facette
                k++;
                while (k < fac.size() && !(fac[k] == allFacCtnr[i]))
                    { k++; }
                if (k < fac.size())
                    {
                    std::cout << "Error : an interface facette belonging "
                                    "to multiple surface regions has been found";
                    exit(1);
                    }
                }
            j++;
            }
        else
            {   // If the current one belongs to only one tetrahedron
            int k = 0;
            while (k < fac.size() && !(fac[k] == allFacCtnr[i]))
                { k++; }
            if (k < fac.size())
                {   // Looking for a second identical facette
                k++;
                while (k < fac.size() && !(fac[k] == allFacCtnr[i]))
                    { k++; }
                if (k < fac.size())
                    {
                    std::cout << "Error : an surface facette belonging "
                                    "to multiple surface regions has been found";
                    exit(1);
                    }
                }
            }
        }
    }

double mesh::surface(std::vector<int> &facIndices) const
    {
    double S(0);
    std::for_each(facIndices.begin(),facIndices.end(),[this,&S](int idx)
                  { S += fac[idx].calc_surf(); });
    return S;
    }

void mesh::setExtSpaceField(Settings &s /**< [in] */)
    { // see here for ref code
// /data/jc/st_feellgood_2024/src_Tube_scalfmm_zhang_ec_mu_oersted_spinHall_thiele_dyn20240320_dev
    extSpaceField.resize(magTet.size());
    int k(0);
    std::for_each(magTet.begin(), magTet.end(), [this,&s,&k](const int idx)
        {
        Tetra::Tet const &t = tet[idx];
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> pg;// gauss points
        t.getPtGauss(pg);
        for(int i=0;i<Tetra::NPI;i++)
            { extSpaceField[k].col(i) = s.getFieldSpace(pg.col(i)); }
        k++;
        });
    }

void mesh::init_distrib(Settings const &mySets)
{
    for (int nodeIdx = 0; nodeIdx < int(node.size()); ++nodeIdx)
        {
        Nodes::Node &n = node[nodeIdx];
        n.d[Nodes::NEXT].phi = 0.;
        n.d[Nodes::NEXT].phiv = 0.;
        // A non-magnetic node's magnetization is a vector of NAN.
        if (!magNode[nodeIdx])
            {
            n.d[Nodes::CURRENT].u = Eigen::Vector3d(NAN, NAN, NAN);
            n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
            continue;
            }
        // If the initial magnetization depends only on the node position, we do not need to
        // build the list of region names.
        if (mySets.getMagType() == POSITION_ONLY)
            {
            n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p);
            n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
            continue;
            }
        // Get the list of region indices this node belongs to.
        std::set<int> nodeRegions;
        // skip region 0 (__default__) which should be empty
        for (size_t regIdx = 1; regIdx < volumeRegions.size(); ++regIdx)
            {
            const std::vector<int> &regionTetras = volumeRegions[regIdx];
            bool node_in_region = std::any_of(EXEC_POL,
                regionTetras.begin(), regionTetras.end(),
                [this, nodeIdx](int tetIdx)
                {
                const Tetra::Tet &tetrahedron = tet[tetIdx];
                for (int i = 0; i < Tetra::N; ++i) // node of tetrahedron tetIdx
                    {
                    if (tetrahedron.ind[i] == nodeIdx)
                        { return true; }
                    }
                return false;
                });
            if (node_in_region)
                { nodeRegions.insert(regIdx); }
            }
        // Get the list of region names this node belongs to.
        std::vector<std::string> region_names;
        region_names.resize(nodeRegions.size());
        std::transform(nodeRegions.begin(), nodeRegions.end(), region_names.begin(),
            [this](int regIdx){ return paramTetra[regIdx].regName; });
        n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p, region_names);
        n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
        }
}