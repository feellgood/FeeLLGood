#include "feellgoodSettings.h"
#include "mesh.h"
#include "pt3D.h"

void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat)
    {
    if (f_in.fail())
        {
        std::cerr << strWhat << ": " << strerror(errno) << std::endl;
        SYSTEM_ERROR;
        }
    }

bool lookFor(const bool _b, std::ifstream &f_in, const std::string strWhat)
    {
    std::string symb = "";
    while ((f_in.peek() != EOF) && (symb != strWhat))
        {
        f_in >> symb;
        }

    if (_b) on_fail_msg_error(f_in, "could not find beacon " + strWhat);

    return !(f_in.fail());
    }

namespace Mesh
    {

void mesh::readMesh(Settings const &mySets)
    {
    if (mySets.verbose)
        {
        std::cout << "Reading mesh file " << mySets.getPbName() << ":\n";
        }
    std::string symb;
    std::ifstream msh(mySets.getPbName());

    on_fail_msg_error(msh, "cannot open file " + mySets.getPbName() );

    msh >> symb;
    if (symb == "$MeshFormat")
        {
        std::string mshFormat = "";
        msh >> mshFormat;
        if (mshFormat == "2.2")
            {
            if (mySets.verbose)
                {
                std::cout << "  file format: 2.2\n";
                }
            int tags, reg, TYP;

            bool beaconFound = lookFor(false, msh, "$PhysicalNames");

            int nbRegNames;
            if (beaconFound)
                {
                msh >> nbRegNames;

                while ((msh >> symb) && (symb != "$EndPhysicalNames") && (!msh.fail()))
                    {
                    std::string name;
                    msh >> tags >> name;

                    switch (stoi(symb))
                        {
                        case 2:
                            {
                            surfRegNames[tags] = name.substr(1, name.length() - 2);
                            break;
                            }
                        case 3:
                            {
                            volRegNames[tags] = name.substr(1, name.length() - 2);
                            break;
                            }
                        default:
                            std::cerr << "unknown type in mesh $PhysicalNames" << std::endl;
                            break;
                        }
                    }
                if (mySets.verbose)
                    {
                    std::cout << "  found " << nbRegNames << " regions:\n";
                    std::map<int, std::string>::iterator it;
                    for (it = surfRegNames.begin(); it != surfRegNames.end(); ++it)
                        {
                        std::cout << "    " << it->first << ": " << it->second << '\n';
                        }
                    for (it = volRegNames.begin(); it != volRegNames.end(); ++it)
                        {
                        std::cout << "    " << it->first << ": " << it->second << '\n';
                        }
                    }
                }
            else
                {
                if (mySets.verbose)
                    {
                    std::cout << "No $PhysicalNames defined.\n";
                    }

                msh.clear();
                msh.seekg(std::ios::beg);
                }

            beaconFound = lookFor(true, msh, "$Nodes");

            double scale = mySets.getScale();
            int nbNod;

            if (beaconFound)
                {
                msh >> nbNod;
                init_node(nbNod);

                for (int i = 0; i < nbNod; i++)
                    {
                    msh >> symb >> node[i].p;
                    node[i].p.rescale(scale);
                    }
                }

            on_fail_msg_error(msh, "error while reading nodes");

            for (auto it = surfRegNames.begin(); it != surfRegNames.end(); ++it)
                {
                s.push_back(Mesh::Surf(node, it->second));
                }

            beaconFound = lookFor(true, msh, "$Elements");

            int nbElem;

            if (beaconFound)
                {
                msh >> nbElem;
                if (mySets.verbose)
                    {
                    std::cout << "  element count: " << nbElem << '\n';
                    }
                while ((msh >> symb) && (symb != "$EndElements") && (!msh.fail()))
                    {
                    msh >> TYP >> tags >> reg;
                    for (int i = 1; i < tags; i++)
                        msh >> symb;
                    switch (TYP)
                        {
                        case 2:
                            {
                            int i0, i1, i2;
                            msh >> i0 >> i1 >> i2;
                            if (auto search = surfRegNames.find(reg); search != surfRegNames.end())
                                {  // found named surface
                                for (auto it = s.begin(); it != s.end(); ++it)
                                    {
                                    if (it->getName() == search->second)
                                        {
                                        it->push_back(Mesh::Triangle(node, i0, i1, i2));
                                        }
                                    }

                                int idx = mySets.findFacetteRegionIdx(search->second);
                                if (idx > -1)
                                    fac.push_back(
                                            Facette::Fac(node, nbNod, idx, i0, i1,
                                                         i2));  // we only want to store in fac
                                                                // vector the facettes for micromag
                                                                // problem, nothing for bc for stt
                                }
                            else
                                {
                                std::cout << "mesh reading error : unnamed surface region"
                                          << std::endl;
                                }  // unnamed surface

                            break;
                            }
                        case 4:
                            {
                            int i0, i1, i2, i3;
                            msh >> i0 >> i1 >> i2 >> i3;
                            if (auto search = volRegNames.find(reg); search != volRegNames.end())
                                {  // found named volume
                                int idx = mySets.findTetraRegionIdx(search->second);
                                if (idx > -1) tet.push_back(Tetra::Tet(node, idx, i0, i1, i2, i3));
                                }
                            else
                                {
                                std::cout << "mesh reading error : unnamed volume region"
                                          << std::endl;
                                }
                            break;
                            }
                        default: std::cerr << "unknown type in mesh $Elements" << std::endl; break;
                        }
                    }
                }
            on_fail_msg_error(msh, "error while reading elements");
            }
        else
            {
            std::cout << "mesh file format " << mshFormat << " not supported." << std::endl;
            SYSTEM_ERROR;
            }
        }

    if (!msh.fail())
        {
        if (mySets.verbose) std::cout << "  closing file\n";
        msh.close();
        }
    else
        {
        std::cout << "error before closing mesh." << std::endl;
        SYSTEM_ERROR;
        }

    for (unsigned int i = 0; i < tet.size(); i++)
        {
        tet[i].idx = i;
        }
    }

double mesh::readSol(bool VERBOSE, const std::string fileName)
    {
    double t(0);
    std::ifstream fin(fileName, std::ifstream::in);

    on_fail_msg_error(fin, "cannot open file: " + fileName);

    std::string str;
    getline(fin, str);
    while ((fin.peek() != EOF) && (str[0] == '#' && str[1] != '#'))
        {
        getline(fin, str);
        }  // to skip comments (if any)

    unsigned long idx;
    bool flag_t = false;
    while ((fin.peek() != EOF) && (str[0] == '#' && str[1] == '#'))
        {
        idx = str.find("## time:");
        if (std::string::npos != idx)
            {  // found beacon "## time:"
            t = stod(str.substr(idx + 7));
            flag_t = true;
            break;
            }
        getline(fin, str);
        }
    getline(fin, str);  // to skip "## columns:" line

    if (!flag_t)
        {
        std::cerr << "error: no ## time: beacon in input .sol file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    if (VERBOSE)
        {
        std::cout << ".sol file: " << fileName << " @ time t = " << t << std::endl;
        }

    for (unsigned int i = 0; i < node.size(); i++)
        {
        Nodes::Node &n = node[i];
        unsigned int i_;
        fin >> i_ >> n.u >> n.phi;
        if (i != i_)
            {
            std::cerr << "error: mesh node index mismatch between mesh and input .sol file"
                      << std::endl;
            fin.close();
            SYSTEM_ERROR;
            }
        }
    fin.close();

    return t;
    }

    }  // namespace Mesh
