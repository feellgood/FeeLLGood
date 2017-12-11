#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cctype>

using namespace std;

/* Extraction des commentaires commencant par un # en 1ere colonne */

void extract_comment(istream &flux)
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

