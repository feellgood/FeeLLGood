#ifndef beacons_h
#define beacons_h

void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat)
    {
    if (f_in.fail())
        {
        std::cerr << strWhat << ": " << strerror(errno) << std::endl;
        SYSTEM_ERROR;
        }
    }

namespace beacons
    {
    namespace sol
        {
        const std::string time = "## time:";
        }
    
    namespace msh
        {
        const std::string format = "$MeshFormat";
        const std::string version = "2.2";
        const std::string begin_physical_names = "$PhysicalNames";
        const std::string end_physical_names = "$EndPhysicalNames";
        const std::string begin_nodes = "$Nodes";
        const std::string begin_elements = "$Elements";
        const std::string end_elements = "$EndElements";
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

    } // end namespace beacons
    
#endif
