#ifndef tags_h
#define tags_h

/**
error handler for input files
*/
void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat);

namespace tags
    {
    namespace sol
        {
        const std::string time = "## time:";
        const std::string rw_time = "## real-world time:";
        const std::string columns = "## columns:";
        const std::string defaultColumnsTitle = "idx\tmx\tmy\tmz\tphi";
        const std::string sColumnsTitle = "sx\tsy\tsz";
        }
    
    namespace evol
        {
        const std::string version = "## feeLLGood version:";
        const std::string hostname = "## hostname:";
        const std::string rw_time = "## real-world time:";
        const std::string settings_file = "## settings file:";
        const std::string columns = "## columns:";
        }

    namespace msh
        {
        const int DIM_OBJ_2D = 2;
        const int DIM_OBJ_3D = 3;
        
        const int TYP_ELEM_TRIANGLE = 2;
        const int TYP_ELEM_TETRAEDRON = 4;

        const int SIZE_TRIANGLE = 3;
        const int SIZE_TETRAEDRON = 4;
        }
    } // end namespace tags
    
#endif
