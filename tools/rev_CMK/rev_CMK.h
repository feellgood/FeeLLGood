#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

boost::char_separator<char> sep("=;:| \"\t");

struct Node { double x, y, z; };

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
boost::property<boost::vertex_color_t, boost::default_color_type, boost::property<boost::vertex_degree_t,int> > > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

typedef boost::graph_traits<Graph>::vertices_size_type size_type;

class BadConversion : public std::runtime_error {
 public:
   BadConversion(const std::string& s)
     : std::runtime_error(s)
     { }
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

template <int NBN>
void read_obj(tokenizer::iterator & tok_iter,Graph & G)
{
int ind[NBN];

for (int nbn=0; nbn<NBN; nbn++)
    {
    boost::advance(tok_iter, 1);
    ind[nbn]=str2num<int>(*tok_iter)-1;
    }

for (int i=0; i<NBN; i++)
    for (int j=0; j<NBN; j++)
        {
        int ori = ind[i];
        int ext = ind[j];
        if (ext > ori) add_edge(ori, ext, G);
        }    
}
 
void reverse_cmk(std::string nom, std::vector <int> &old2newlabel, std::vector <int> &new2oldlabel)
{
boost::property_map <Graph, boost::vertex_index_t>::type index_map;
std::vector <size_type> perm;
std::vector <Vertex>    inv_perm;

std::ifstream fin(nom.c_str());
if (!fin) {
   std::cerr << "Pb lecture fichier " << nom << std::endl;
   exit(1);
   }

std::string str;
for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    std::cout << str << std::endl;
    if (*tok_iter=="$Nodes") break; 
}

int NOD=0;
for (;;) {
    getline(fin, str);
    if (str.empty()) { continue; }
    
    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    boost::advance(tok_iter, 0);
    NOD=str2num<int>(*tok_iter);
    std::cout << NOD << std::endl;
    break;
    }

Graph G(NOD);
perm.resize(NOD);
inv_perm.resize(NOD);
old2newlabel.resize(NOD);
new2oldlabel.resize(NOD);

std::vector<Node> coord(NOD);
for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    if (*tok_iter=="$EndNodes") { break; }
    
    int i=str2num<int>(*tok_iter)-1; 
    boost::advance(tok_iter, 1);
    coord[i].x=str2num<double>(*tok_iter);
    boost::advance(tok_iter, 1);
    coord[i].y=str2num<double>(*tok_iter);
    boost::advance(tok_iter, 1);
    coord[i].z=str2num<double>(*tok_iter);
}



for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    std::cout << str << std::endl;
    if (*tok_iter=="$Elements") break; 
}

int NE=0;
for (;;) {
    getline(fin, str);
    if (str.empty()) { continue; }
    
    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    boost::advance(tok_iter, 0);
    NE=str2num<int>(*tok_iter);
    std::cout << NE << std::endl;
    break;
    }

int tags, TYP;

for (;;){
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    if (*tok_iter=="$EndElements") { break; } 

    boost::advance(tok_iter, 1);
    TYP=str2num<int>(*tok_iter);
    boost::advance(tok_iter, 1);
    tags=str2num<int>(*tok_iter);
    boost::advance(tok_iter, 1);
    boost::advance(tok_iter, tags-1);

    switch (TYP){
	case 15: break;
        case 1:{
            const int NBN=2;
            read_obj<NBN>(tok_iter, G);
            break;
	    }

        case 2:{
            const int NBN=3;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 3:{
            const int NBN=4;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 4:{
            const int NBN=4;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 5:{
            const int NBN=8;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 9:{
            const int NBN=6;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 16:{
            const int NBN=8;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
/* 3D */
        case 11:{
            const int NBN=10;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 17:{
            const int NBN=20;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        case 18:{
            const int NBN=15;
            read_obj<NBN>(tok_iter, G);
            break;
	    }
        default: {
	    std::cout << "unknown type : " << TYP << std::endl;
            exit(1);
	    }
	}
    }

index_map = get(boost::vertex_index, G);
std::cout << "\nOriginal Bandwidth: " << bandwidth(G) << std::endl;

// reverse cuthill_mckee_ordering
cuthill_mckee_ordering(G, inv_perm.rbegin(), get(boost::vertex_color, G), make_degree_map(G));

for (size_type c = 0; c != inv_perm.size(); ++c) perm[index_map[inv_perm[c]]] = c;

std::cout << "Final Bandwidth: " << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0])) << std::endl;

for (std::vector<Vertex>::const_iterator i=inv_perm.begin(); i!=inv_perm.end(); ++i){
    int oldlabel=inv_perm[index_map[*i]];
    int newlabel=index_map[*i];
    old2newlabel[oldlabel]=newlabel;
    new2oldlabel[newlabel]=oldlabel;	   
    }

for (;;){
    getline(fin, str); // 1eme ligne
    std::cout << str << std::endl;
    if (fin.eof()) break;
    }
fin.close();    
}



template <int NBN>
void write_obj(std::ofstream & fout, tokenizer::iterator & tok_iter, std::vector <int> &old2newlabel)
{
int ind[NBN];
    
for (int nbn=0; nbn<NBN; nbn++)
    {
    boost::advance(tok_iter, 1);
    ind[nbn]=str2num<int>(*tok_iter)-1;
    int i=old2newlabel[ind[nbn]]+1;
    fout << boost::format("%d ") % i;
    }
fout << std::endl;
}

void update_labelling(std::string nom, std::vector <int> &old2newlabel, std::vector <int> &new2oldlabel)
{
boost::char_separator<char> sep("=;:| \"\t");

std::string nom0=nom + ".orig";
int ier=rename(nom.c_str(), nom0.c_str());
if (ier) {
   perror( "Error renaming file" );
   exit(1);
   }

std::ifstream fin(nom0.c_str());
if (!fin) {
   std::cerr << "cannot open file " << nom << std::endl;
   exit(1);
   }

std::ofstream fout(nom.c_str());
if (!fout) {
   std::cerr << "cannot write file " << std::endl;
   exit(1);
   }

std::string str;
for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) {
       fout << std::endl;
       continue;
       }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    fout << str << std::endl;
    if (*tok_iter=="$Nodes") break; 
}

int NOD=0;
for (;;) {
    getline(fin, str);
    if (str.empty()) {
       fout << std::endl;
       continue;
       }
    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    boost::advance(tok_iter, 0);
    NOD=str2num<int>(*tok_iter);
    fout << NOD << std::endl;
    break;
    }

std::vector<Node> coord(NOD);
for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    if (*tok_iter=="$EndNodes") { break; }
    int i=str2num<int>(*tok_iter)-1;
    boost::advance(tok_iter, 1);
    coord[i].x=str2num<double>(*tok_iter);
    boost::advance(tok_iter, 1);
    coord[i].y=str2num<double>(*tok_iter);
    boost::advance(tok_iter, 1);
    coord[i].z=str2num<double>(*tok_iter);
}

for (int nod=1; nod<=NOD; nod++) {
    const int i = new2oldlabel[nod-1];
    fout << boost::format("%d %.18g %.18g %.18g") % nod % coord[i].x % coord[i].y % coord[i].z << std::endl;
    }
fout << "$EndNodes" << std::endl;

for (;;) {
    getline(fin, str); // 1eme ligne
    if (str.empty()) {
       fout << std::endl;
       continue;
       }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();

    if (*tok_iter=="$Elements") {
       fout << str << std::endl;
       break;
       } 
}

int NE=0;
for (;;) {
    getline(fin, str);
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    boost::advance(tok_iter, 0);
    NE=str2num<int>(*tok_iter);
    fout << NE << std::endl;
    break;
    }

int ne, tags, reg, TYP;

for (;;){
    getline(fin, str); // 1eme ligne
    if (str.empty()) { continue; }

    tokenizer tokens(str, sep);
    tokenizer::iterator tok_iter = tokens.begin();
    if (*tok_iter=="$EndElements") { break; }

    ne =str2num<int>(*tok_iter);
    boost::advance(tok_iter, 1);
    TYP=str2num<int>(*tok_iter);
    boost::advance(tok_iter, 1);
    tags=str2num<int>(*tok_iter);
    boost::advance(tok_iter, 1);
    reg=str2num<int>(*tok_iter);

    fout << boost::format("%d %d %d %d ") % ne % TYP % tags % reg;
    
    for (int ntag=1; ntag<tags; ntag++) {
        boost::advance(tok_iter, 1);
        int i=str2num<int>(*tok_iter);
        fout << boost::format("%d ") % i;
	    }
    
    switch (TYP){
        case 2:{
            const int NBN=3;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 3:{
            const int NBN=4;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 4:{
            const int NBN=4;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 5:{
            const int NBN=8;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 9:{
            const int NBN=6;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 16:{
            const int NBN=8;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
/* 3D */
        case 11:{
            const int NBN=10;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 17:{
            const int NBN=20;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }
        case 18:{
            const int NBN=15;
            write_obj<NBN>(fout, tok_iter, old2newlabel);
            break;
	        }

        default:
            exit(1);
	}
    }

fout << "$EndElements" << std::endl;

for (;;){
    getline(fin, str); // 1eme ligne
    fout << str << std::endl;
    if (fin.eof()) break;
    }
fin.close();
fout.close();    
}
