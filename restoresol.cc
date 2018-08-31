#include "fem.h"

#include "pt3D.h"

using namespace std;

/** convenient error handler */
class BadConversion : public std::runtime_error {
 public:
/**
constructor
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

void Fem::restoresol(double scaling, string *filename)  // filename may be NULL
{
string str("sol.in");

if (filename){
    str = *filename;
    }

ifstream fin(str, std::ifstream::in); //  *ct*
if (!fin){
    IF_VERBOSE() cerr << "pb ouverture fichier: " << str << "en lecture" << endl;
    SYSTEM_ERROR;}

getline(fin, str); // 1eme ligne

unsigned long idx = str.find(":");
idx +=2;

t = stod(str.substr(idx));

IF_VERBOSE() cout << "fichier solution: " << str << " a l'instant t = " << t << endl;

//const int NOD = fem.NOD;

for (int i=0; i<NOD; i++){
    Node &n = node[i];
    Node node_;
    int i_;
    fin >> i_ >>  node_.p;// carefull! >> is overloaded for class pt3D
    fin >> n.u >> n.phi;

    //node_.x *= scaling;node_.y *= scaling;node_.z *= scaling;
node_.p *= scaling;

    double d2=sq(n.p.x()-node_.p.x()) + sq(n.p.y()-node_.p.y()) + sq(n.p.z()-node_.p.z());

    if (d2 > sq(diam * 1e-9)) // attention scaling écrit en dur ... 
	{
        IF_VERBOSE(){
        cerr << "WARNING difference dans la position des noeuds"<< endl;
        cerr << i  << "\t" << n.p << endl << i_ << "\t" << node_.p << endl;
        }
//        exit(1);
        }
    
    if (i!=i_){
        IF_VERBOSE() cerr << "fichier solution incompatibilite de noeuds"<< endl;
        SYSTEM_ERROR;
        }
    }

fin.close();    			     
}
