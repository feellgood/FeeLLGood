#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>
#include <unistd.h> // for getpid()

#include<boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "feellgoodSettings.h"

using namespace std;

void Settings::helloDialog()
{
cout << "\n\t ******************************\n";
cout <<   "\t *         feeLLGood          *\n";
cout <<   "\t *        version 2017        *\n";
cout <<   "\t *      cnrs Grenoble-INP     *\n";
cout <<   "\t ******************************\n";
cout << "\t process\t\t" << getpid() << endl;
}

void Settings::printToTerminal(std::vector<Seq> &seq)
{
cout << "\t name simul : "<< simName << endl;
cout << "\t name mesh file : " << pbName << endl;
cout << "\t temps final\t\t" << tf << endl;
cout << "\t pas de temps init\t" << dt << endl;
cout << endl;
cout << "\t sauve energies chaque " << n1 << " iterations" << endl;
cout << "\t photo config chaque " << n2 << " iterations" << endl;
cout << "\t restauration a partir d'un fichier sol.in: " << restore << endl;
cout << endl;
cout << "\t sequences de champ applique Bext" << endl;
for (vector<Seq>::iterator it = seq.begin(); it!=seq.end(); ++it) {
    cout << it->Bini << "\t" << it->Bfin << "\t" << it->dB << "\t" << it->a[0] << "\t" << it->a[1] << "\t" << it->a[2] << endl;
    }
}

void Settings::extract_comment(std::istream &flux)
{
string comment;
char c;
while (1)
      {
      while (flux.get(c))  // lecture des caracteres blancs (espaces, nl, tab, cr, np)
	    if (isspace(c)==0)
	       {
	       flux.putback(c);
	       break;
	       }

      flux.get(c);  
      if (c=='#') 
	 getline(flux, comment,'\n');  // lecture du commentaire
      else
	 {
	 flux.putback(c);
	 break;
	 }
      }
}

void Settings::dialog(std::vector<Seq> &seq)
{
int SEQ;
double theta; // removed unused variable tf  *ct*
pair <string,int>  p;


cout << "\t name simul ? "<<endl;
extract_comment(cin);   cin >> simName;
cout << "\t name simul\t\t"<< simName << endl;

cout << "\t file mesh ?"<< endl;
extract_comment(cin);   cin >> pbName;
cout << "\t name file pro\t" << pbName << endl;

ifstream fin(pbName);
if (!fin){
    cerr << "Cannot open problem file, fileName = " << pbName << endl;
    exit(1);}
fin.close();

extract_comment(cin);              cin >> tf;
      
extract_comment(cin);              cin >> dt;

extract_comment(cin);              cin >> theta;
p=make_pair("theta",-1);           param[p]=theta;  

extract_comment(cin);              cin >> n1;
p=make_pair("save_energies",-1);   param[p] = n1;

extract_comment(cin);              cin >> n2;
p=make_pair("take_photo",-1);      param[p] = n2;

extract_comment(cin);              cin >> restore;
p=make_pair("restore",-1);         param[p]=restore;

extract_comment(cin);              cin >> SEQ;

for (int nseq=0; nseq<SEQ; nseq++)
	{
    Seq s;
    extract_comment(cin);          cin >> s.Bini >> s.Bfin >> s.dB >> s.a[0] >> s.a[1] >> s.a[2];
    s.dB=fabs(s.dB)*(s.Bfin<s.Bini? -1.0: 1.0);
    normalize(s.a);
    seq.push_back(s);
    }

cout << "\t theta\t\t\t" << theta << endl;
printToTerminal(seq);
}

void Settings::read(std::vector<Seq> &seq)
{
	boost::property_tree::ptree root;
	boost::property_tree::read_json("settings.json",root);
	
	simName = root.get<string>("output_filename");
	pbName = root.get<string>("mesh_filename");
	tf = root.get<double>("final_time",0);
	dt = root.get<double>("initial_time_step",0);
	theta = root.get<double>("theta",0);
	n1 = root.get<int>("save_energies",0);
	n2 = root.get<int>("take_photo",0);
	restore = root.get<bool>("restore",0);
	
	int i=0;
	double trucs[6];
	for (boost::property_tree::ptree::value_type &s : root.get_child("field_sequence"))
		{
		int j=0;
		for(boost::property_tree::ptree::value_type &cell : s.second)
			{
			trucs[j] = cell.second.get_value<double>();
			j++;
			}
		Seq field;
		field.Bini = trucs[0];
		field.Bfin = trucs[1];
		field.dB = trucs[2];
		field.a[0] = trucs[3];
		field.a[1] = trucs[4];
		field.a[2] = trucs[5];
		
	seq.push_back(field);
		i++;	
		}
	
	cout << " reading parameters and settings from json file :\n" << endl;
}
