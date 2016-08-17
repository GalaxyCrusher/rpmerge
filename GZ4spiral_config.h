#include "GZ4_config.h"

/*When printing casjobs links, the inclusions of additional information can be turned on or off (for use with wget).
0=include info
1=SDSS CASJOBS
2=GZ4 Candels
*/
int wget=2;

int matchSwitch=1;
float zslop=0.25;
float mslop=0.2;
int matchNum=5;

char catCode[10]="UDS";
//char catCode[10]="GDS";

/*Turn on html code to display images as a webpage:
0 = off
1 = on
*/
int htmlSwitch=1;

int cason=0;
int cwacw=0;
float apMagMax=23.5;
float zComp=1.0;
//float zComp=0.846063;  //apMagMax=17.0  - GZ2
//float zComp=0.980763;   //apMagMax=17.77 - SDSS
//float zComp=0.743058;  //apMagMax=17.77 - GZ1
float z1=0.0;
float z2=4.0;


/*Angular spectroscopic correction bin step*/
float tStep=10.0;

/*Pair velocity difference limits*/
float minv=0.;
float maxv=1000000.;
float vStep=100;

/*Set z of galaxy to lowest in pair. For A and P_M contamination check*/
int contcheck=0;

/*Sample absolute magnitude limits*/
float m1=-10;
float m2=-30;
float lStep=1000000.0;

/*Sample mass limits*/
float mass1=8.0;
float mass2=11.5;
float mStep=1.0;

/*Use luminosity (magnitude) bins OR mass bins:
0 = mass 
1 = magnitude*/
int LorM=0;


/*Pair projected separation limits*/
float r1=0.0;
float r2=50;
float rStep=50;

/*Switch to bin by projected separation or relative separation 
0=physical projected separation in kpc
1=relative separation: rp/(r1+r2) */
int rpswitch = 0;

/*Amount to multiply combined pair radii by in order to get rid or overlap
The combined radii X radx is the min pair separation (in arcseconds)*/
float radx=0.0;

/*Maximum magnitude separation between pairs*/
float dMmin=0;
float dMmax=10;



/*GZ2 Class thresholds*/
float fesmin=0.1;
float ftmin=0.3;

float clmpmin=0.3;
float clmpmax=0.5;
float clmpspmin=0.2;

float egmin=0.5;

float spmin=0.3;
float spmax=0.8;

float twamax=0.4;
float twamin=0.8;

float mwamax=0.4;
float mwamin=0.8;

float lwamax=0.4;
float lwamin=0.6;

float awtight=0.33;
float awloose=0.66;

float oamin=0.0;	//N/A
float tamin=0.0;	//N/A
float oddmin=0.0;	//N/A

float pmin=0.2;
float tdmin=0.2;
float mrgtdmin=0.2;
float nomrgmin=0.2;

float irrmin=0.4;
float irrmax=0.1;

float dismin=0.4;
float dismax=0.1;

float barmin=0.2;
float barmax=0.2;

float noblgmin=0.5;		//No bulge

float jnblgmin=0.0;		//N/A

float obblgmin=0.5; 	//face-on obvious bulge

float doblgmin=0.5; 	//face-on dominant bulge



/*Global GZ morphological class restrictions:
0=Off

10=Pre pass				feturs>ftmin && edgeno>egmin && sparms>spmin && lwarms<lwamax

20=between passes		feturs>ftmin && edgeno>egmin && sparms>spmin && lwarms>lwamin
21=between passes v2	(feturs*edgeno*sparms)>fesmin  && lwarms>lwamin && (oarm>oamin || tarms>tamin)
22=simple pmerg cutoff	pmrg>pmin

30=merger proper		feturs>ftmin && edgeno>egmin && (irrgl>irrmin || dstrb>dismin) && pmrg>pmin && odd>oddmin

40=features and face on	feturs>ftmin && edgeno>egmin
41=All spirals

50=tight arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && twarms>twamin
51=tight arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && twarms>twamin && baryes>barmin

60=medium arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && mwarms>mwamin
61=medium arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && mwarms>mwamin && baryes>barmin

70=loose arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && lwarms>lwamin
71=loose arms   		feturs>ftmin && edgeno>egmin && sparms>spmin && lwarms>lwamin && baryes>barmin

80=No arms	   			feturs>ftmin && edgeno>egmin && sparms<spmax
81=No arms	   			feturs>ftmin && edgeno>egmin && sparms<spmax && baryes>barmin


110=No bulge			feturs>ftmin && edgeno>egmin && noblg>noblgmin
111=No bulge w bar 		feturs>ftmin && edgeno>egmin && noblg>noblgmin && baryes>barmin

120=JN bulge 			feturs>ftmin && edgeno>egmin && jnblg>jnblgmin
121=JN bulge w bar 		feturs>ftmin && edgeno>egmin && jnblg>jnblgmin && baryes>barmin

130=Ob bulge 			feturs>ftmin && edgeno>egmin && obblg>obblgmin
131=Ob bulge w bar 		feturs>ftmin && edgeno>egmin && obblg>obblgmin && baryes>barmin

140=Do bulge  			feturs>ftmin && edgeno>egmin && doblg>doblgmin
141=Do bulge w bar 		feturs>ftmin && edgeno>egmin && doblg>doblgmin && baryes>barmin


*/
int grestrict=41;	// morphological filter type

/*Filter the pair members instead of global filter:*/
/*0 = off (Global filter applies)
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
float ratiomax=100000.0;

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


