#include "fem.h"
#define DEBUG 0

void normalize(triple &a){
    double norme=sq(a[0])+sq(a[1])+sq(a[2]);
    norme=sqrt(norme);
    for (int d=0; d<3; d++)
        a[d]/= norme;
    }
