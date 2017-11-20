#include "fem.h"

/** convenient error handler */
class BadConversion : public std::runtime_error {
 public:
/**
do nothing constructor (to empty buffer for gmm?)
*/
   BadConversion(const std::string& s) : std::runtime_error(s) { }
 };

template <typename T> 
inline T str2num(const std::string& s)
 {
   std::istringstream iss(s);
   T x;
   if (!(iss >> x))
     throw BadConversion("convertToDouble(\"" + s + "\")");
   return x;
 }

void restoresol(Fem& fem, string *filename)  // filename may be NULL
{
string str("sol.in");

if (filename){
    str = *filename;
    }

ifstream fin(str, std::ifstream::in); //  *ct*
if (!fin){
    IF_VERBOSE(fem) cerr << "pb ouverture fichier: " << str << "en lecture" << endl;
    SYSTEM_ERROR;}

getline(fin, str); // 1eme ligne

unsigned long idx = str.find(":");
idx +=2;

fem.t = stod(str.substr(idx));

IF_VERBOSE(fem) cout << "fichier solution: " << str << " a l'instant : " << fem.t << endl;

const int NOD = fem.NOD;

for (int i=0; i<NOD; i++){
    Node &node = fem.node[i];
    Node node_;
    int i_;
    fin >> i_ >> node_.x >> node_.y >> node_.z;
    fin >> node.u[0] >> node.u[1] >> node.u[2] >> node.phi;

    node_.x*=fem.scale;
    node_.y*=fem.scale;
    node_.z*=fem.scale;

    double d2=sq(node.x-node_.x) + sq(node.y-node_.y) + sq(node.z-node_.z);

    if (d2 > sq(fem.diam * 1e-9)) // attention scaling écrit en dur ... 
	{
        IF_VERBOSE(fem){
        cerr << "WARNING difference dans la position des noeuds"<< endl;
        cerr << i  << "\t" << node.x  << "\t" << node.y  << "\t" << node.z << endl;
        cerr << i_ << "\t" << node_.x << "\t" << node_.y << "\t" << node_.z << endl;
        }
//        exit(1);
        }
    
    if (i!=i_){
        IF_VERBOSE(fem) cerr << "fichier solution incompatibilite de noeuds"<< endl;
        SYSTEM_ERROR;
        }
    }

fin.close();    			     
}
