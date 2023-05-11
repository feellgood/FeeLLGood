#ifndef tags_h
#define tags_h

void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat);

namespace tags
    {
    namespace sol
        {
        const std::string time = "## time:";
        const std::string rw_time = "## real-world time:";
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

        const int DIM_OBJ_2D = 2;
        const int DIM_OBJ_3D = 3;
        
        const int TYP_ELEM_TRIANGLE = 2;
        const int TYP_ELEM_TETRAEDRON = 4;
        }

    bool lookFor(const bool _b, std::ifstream &f_in, const std::string strWhat);

    } // end namespace tags
    
#endif
