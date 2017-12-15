// #include "fem.h"

#include "linear_algebra.h"
#include "fmm_demag.h"

int main()
{
Settings mySettings = Settings();
LinAlgebra linAlg = LinAlgebra();
Fem fem;

OctreeClass *tree    = nullptr;
KernelClass *kernels = nullptr; 

#ifdef STAT
fem.stat.h = gsl_histogram_alloc (NCLASSES);
gsl_histogram_set_ranges_uniform (fem.stat.h, -HMAX, HMAX);
#endif

vector<Seq> seq;

mySettings.dialog(seq);
lecture(fem,mySettings, 0.0, nullptr);
femutil(fem,mySettings);
chapeaux(fem);
affichage(fem);

fem.t=0.;
pair <string,int> p; p= make_pair("restore",-1);
int restore = int(mySettings.param[p]);
if (restore)
    restoresol(fem,mySettings.getScale(), nullptr);
else
    init_distrib(fem);

/* direction de propagation de la paroi */
direction(fem);

fmm::init< CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem, tree, kernels);

double dt0= mySettings.dt;

int nseq=0;
for (vector<Seq>::iterator it = seq.begin(); it!=seq.end(); ++it) {
    double &Bini=it->Bini;
    double &Bfin=it->Bfin;
    double &dB=it->dB;
    triple &a=it->a;
    nseq++;
    fem.SEQ=nseq;

    cout << "dB : " << dB << endl;

    for (int loop=0; loop<=(int)((Bfin-Bini)/dB+0.5); loop++)
    	{
    	double Bext=Bini+dB*loop;
        fem.Bext=Bext;
        fem.Hext[0]=nu0*Bext*a[0];
        fem.Hext[1]=nu0*Bext*a[1];
        fem.Hext[2]=nu0*Bext*a[2]; 

	cout << "Bext : " << Bext*a[0] << "\t" << Bext*a[1] << "\t" << Bext*a[2] << endl;

        string str = mySettings.getSimName() +"_"+ to_string(fem.SEQ) + "_B" + to_string(fem.Bext) + ".evol";// +++ *ct*	
	/*	
	ostringstream ostr;
        ostr.str(""); ostr << fem.simname << boost::format("_%d_B%6f.evol") % fem.SEQ % fem.Bext; 
        str = ostr.str(); ofstream fout(str.c_str());
        */ // --- *ct*
	ofstream fout(str);
	if (!fout) {
            cerr << "erreur ouverture fichier" << endl; exit(1);
            }

        fmm::demag<0, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem,mySettings, tree, kernels);

#ifdef ORD2
        fmm::demag<1, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> (fem,mySettings, tree, kernels);
#endif
	    fem.DW_z  = 0.0;
        energy(fem,mySettings); 

	    evolution(fem);

	    int flag  = 0;
	    double dt = dt0;
	    double t= fem.t;
	    fem.vmax  = 0.0;
	    fem.DW_vz = fem.DW_vz0 = 0.0;
	    fem.DW_z  = 0.0;

	    int nt = 0;
	    saver(fem,mySettings,fout,nt);

        double start_cpu = cputime();

        while (t < mySettings.tf) {
            cout << "\n ------------------------------\n";
            if (flag) cout << "    t  : same (" << flag << ")" << endl;
		else cout << "nt = " << nt << ", t = " << t << endl; // *ct*
            //else cout << boost::format("nt = %d,  t = %2.8g") % nt % t << endl;
		            
		cout << "dt = " << dt << endl << endl;
            if (dt < DTMIN) {
                reset(fem);
                break;
                }

            int err = linAlg.vsolve(fem,mySettings, nt);  
            if (err) { cout << "err : " << err << endl;
		    	flag++; dt*= 0.5; mySettings.dt=dt; continue;}

            double dumax = dt*fem.vmax;
            //cout << boost::format("\t dumax = %2.2e,  vmax = %2.2e") % dumax % fem.vmax<< endl;// *ct*
		cout << "\t dumax = " << dumax << ",  vmax = "<< fem.vmax << endl;
            if (dumax < DUMIN) break; 
/*
            if (dumax > DUMAX) { 
			flag++; dt*= DUMAX/dumax; fem.dt=dt; continue;}
*/

            if (dumax > DUMAX) { 
			flag++; dt*= 0.5; mySettings.dt=dt; continue;}
          
        fmm::demag<0, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> 
		(fem,mySettings, tree, kernels); // Hd(u)
#ifdef ORD2
        fmm::demag<1, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass> 
		(fem,mySettings, tree, kernels); // Hd(v)
#endif

            energy(fem,mySettings);

         /*   cout << boost::format("\t energy %+2.3e") % fem.Etot << endl;
              cout << boost::format("\t   (dE/dt %+2.2e, av2 %+2.2e)")
                         % (fem.evol/dt) % fem.phy << endl;
         */
//            double dissip = 1. - fem.evol/dt/fem.phy;
//            cout << boost::format("\t dissipation %+3.1f") %(dissip*100) <<'%'<< endl;
            if (fem.evol > 0.0)
                {
		cout << "Warning energy increases! : " << fem.evol << endl;
//                flag++; dt*= 0.5; fem.dt=dt; continue;
                }

/* mise a jour de la vitesse du dernier referentiel et deplacement de paroi */
	    fem.DW_vz0 = fem.DW_vz; 
            fem.DW_z  += fem.DW_vz*dt;

            evolution(fem); t+=dt; fem.t=t; nt++; flag=0;
	    recentrage(fem, 0.1);

            saver(fem,mySettings,fout,nt);

            dt = min(1.1*dt, DTMAX); 
	    mySettings.dt=dt;
            }//endwhile

        if (dt < DTMIN) cout << " aborted:  dt < DTMIN";
        saver(fem,mySettings,fout,nt);
        savecfg_vtk(fem,mySettings.getSimName(),mySettings.getScale(),nt,nullptr);
        savesol(fem,mySettings.getSimName(),mySettings.getScale(),nt,nullptr);

        double end_cpu = cputime();
        double dureecpu = end_cpu - start_cpu;
        cout << "\n  * iterations: " << nt;
        cout << "\n  * duree cpu: " << dureecpu << " s"<< endl << endl;
        fout.close();

//      initialisation du temps pour l'iteration en champ suivante
        fem.t=0.;
        mySettings.dt=dt0;
        }//endfor Bext

    }//endfor it

cout << "--- the end ---" << endl;
delete tree;
delete kernels;

#ifdef STAT
gsl_histogram_free (fem.stat.h);
#endif
}
