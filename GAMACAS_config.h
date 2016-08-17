#include "CAS_config.h"


/*When printing casjobs links, the inclusions of additional information can be turned on or off (for use with wget).
0=include info
1=SDSS CASJOBS
2=GZ4 Candels
*/
int wget=-1;

int matchSwitch=0;
float zslop=0.25;
float mslop=0.2;
int matchNum=5;

/*Turn on html code to display images as a webpage:
0 = off
1 = on
*/
int htmlSwitch=0;

int cason=1;
int cwacw=0;
float apMagMax=19.8;
//float zComp=1.0;
float zComp=0.977927;
//float zComp=0.954964; //redshift completeness - GAMA
//float zComp=0.878276722; //mass completeness - GAMA/
float z1=0.001;
float z2=0.2;


/*Angular spectroscopic correction bin step*/
float tStep=1.0;

/*Pair velocity difference limits*/
float minv=0.;
float maxv=500.;
float vStep=100;

/*Set z of galaxy to lowest in pair. For A and P_M contamination check*/
int contcheck=0;

/*Sample absolute magnitude limits*/
float m1=-10.0;
float m2=-30.0;
float lStep=1000000.0;

/*Sample mass limits*/
float mass1=6.5;
float mass2=12.5;
float mStep=0.5;

/*Use luminosity (magnitude) bins OR mass bins:
0 = mass 
1 = magnitude*/
int LorM=0;


/*Pair projected separation limits*/
float r1=0.0;
float r2=2.5;
float rStep=0.25;

/*Switch to bin by projected separation or relative separation 
0=physical projected separation in kpc
1=relative separation: rp/(r1+r2) */
int rpswitch = 1;

/*Amount to multiply combined pair radii by in order to get rid or overlap
The combined radii X radx is the min pair separation (in arcseconds)*/
float radx=0.0;

/*Maximum magnitude separation between pairs*/
float dMmin=0;
float dMmax=1000;

/*Global morphological class restrictions:
0=Off
201 => rtype[apos]>alimit(rmass) && rtype[apos]>rtype[spos]
211 => rtype[apos]>alimit(rmass) && rtype[spos]>rtype[apos]
*/
int grestrict=0;	// restrict all galaxies

/*Use GZ2 classification thersholds*/
/*0 = off
1=secondary
2=both
3=either*/
int prestrict=0;


/*Turn on or off Vmax weights (for flux limited sample)
0 = off - volume limited sample
1 = on - flux limited sample*/
int vmaxon=0;
int boundon=0;

/*Apply weight to account for interloper contamination (this correction is NOT mass dependent)
0 = off
1 = Kevin's method
2 = Steven's method*/
int removeint=0;

float a0=1.0;
float b0=20.1453319859;
float c0=0.0348624304216;


/*Switch to consider either the brighter of fainter galaxy only:
 0 = off
 1 = faint galaxy only
-1 = bright galaxy only*/
int fb=0;

/*Restrict bright member to a certain range*/
int brestrict=0;
float blim = 10.0;

/*Restrict heavy member to a certain range*/
int hrestrict=0;
float hlim = 5.0;

/*Mass ratios between galaxies*/
float ratiomin=1.0;
float ratiomax=100.0;

/*Switch to consider either heavy or light galaxy only:
 0 = off
 1 = heavy galaxy only
-1 = light galaxy only*/ 
int hl=0;

/* CAS Asymmetry limits*/
int Amass=0; //Mass dependent A_limit =1
float Adefault=0.35;

float pStep=0.1;

/*Switch to limit mergers using CAS by (1 OR 2) or (1 AND 2)
0= OFF
1 = OR
2 = AND
3 = PRIMARY ONLY*/
int ao=0;

/*Switch to select red, blue or red+blue pairs
0 = everything
1 = red only
2 = blue only
3 = ...*/
int rb=0;
