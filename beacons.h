#ifndef beacons_h
#define beacons_h

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

    } // end namespace beacons
    
#endif
