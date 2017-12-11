#include "fem.h"
#define DEBUG 0

void dialog(Fem &fem, vector<Seq> &seq)
{
cout << "\n\t ******************************\n";
cout <<   "\t *         feeLLGood          *\n";
cout <<   "\t *        version 2017        *\n";
cout <<   "\t *      cnrs Grenoble-INP     *\n";
cout <<   "\t ******************************\n";

int n1,n2,restore, SEQ;
double theta; // removed unused variable tf  *ct*
pair <string,int>  p;
map <pair<string,int>,double>  &param=fem.param;

cout << "\t process\t\t" << getpid() << endl;

cout << "\t nom simul ? "<<endl;
extract_comment(cin);   cin >> fem.simname;
cout << "\t nom simul\t\t"<< fem.simname << endl;

cout << "\t fichier mesh ?"<< endl;
extract_comment(cin);   cin >> fem.pbname;
cout << "\t nom fichier pro\t" << fem.pbname << endl;

ifstream fin(fem.pbname.c_str());
if (!fin){
    cerr << "Impossible d'ouvrir le fichier probleme " << fem.pbname << endl;
    exit(1);}
fin.close();

extract_comment(cin);              cin >> fem.tf;
      
extract_comment(cin);              cin >> fem.dt;

extract_comment(cin);              cin >> theta;
p=make_pair("theta",-1);           param[p]=theta;  

extract_comment(cin);              cin >> n1;
p=make_pair("save_energies",-1);   param[p] = n1;

extract_comment(cin);              cin >> n2;
p=make_pair("take_photo",-1);      param[p] = n2;

extract_comment(cin);              cin >> restore;
p=make_pair("restore",-1);         param[p]=restore;

extract_comment(cin);              cin >> SEQ;

for (int nseq=0; nseq<SEQ; nseq++) {
    Seq s;
    extract_comment(cin);          cin >> s.Bini >> s.Bfin >> s.dB >> s.a[0] >> s.a[1] >> s.a[2];
    s.dB=fabs(s.dB)*(s.Bfin<s.Bini? -1.0: 1.0);
    normalize(s.a);
    seq.push_back(s);
    }

cout << "\t temps final\t\t" << fem.tf << endl;
cout << "\t pas de temps init\t" << fem.dt << endl;
cout << "\t theta\t\t\t" << theta << endl;
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
