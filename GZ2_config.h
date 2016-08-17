/*Constants of nature*/
#define pi 3.1415927
#define cs 299792.458 
#define e 2.7188

/*Cosmological constants*/
#define omegaM 0.3
#define omegaL 0.7
#define omegaK 0.0
#define	hs 0.7

/*Hubble constant*/
float Ho = 100.*hs;
 
/*Hubble distance*/
float DH = cs/(100.*hs);

/*GZ2 switch*/
/*Catalogue length*/
#define ls 270000 //Sample length
#define ls2 270000 //Total catalogue length
#define ksize 1//Magnitude or Mass bins
#define lsize 1 //rp bins=
int Nvar = 44; //This number has to represent the exact number of GZ classification columns in the catalog
#define Nextra 5 //Number of calculated mophological index variables 
#define NvarArray 51 //This is the storage array index number and can be higher than Nvar

//pvalV[NvarArray][ksize][lsize][ls]

/*Columns for GZ2 morphological features of interest*/
#define smooth 0	//Smooth
#define feturs 1	//Features
#define edgeno 4	//Edge on - No
#define baryes 5	//Bar - Yes
#define sparms 7	//Spiral - Yes
#define noblg 9		//face-on no-bulge
#define	jnblg 10	//face-on just noticable bulge
#define obblg 11 	//face-on obvious bulge
#define doblg 12 	//face-on dominant bulge
#define odd 13		//Something Odd
#define dstrb 20	//Something Odd - Disturbed 
#define irrgl 21	//Something Odd - Irregular
#define pmrg 23		//Someting Odd - Merger
#define twarms 28	//Tight winding arms
#define mwarms 29	//Medium winding arms
#define lwarms 30	//Loose winding arms
#define oarm 31		//1 arm
#define tarms 32	//2 arms
#define threearms 33 //3 arms
#define fourarms 34	//4 arms
#define morearms 35	//4+ arms
#define unknarms 36	//? arms		


/*Columns for CW and ACW classifications*/
#define cw 38
#define acw 39

#define awidx 44    //arm winding tightness index
#define blgidx 45	//bulge domenance index
#define armnum 46	//most probable arm number
#define mclass 47	//most probable morpholgical class
#define spin 48

/*Columns for GZ4 Candels morphological features of interest*/
#define clmp 0		//N/A
#define clmpsp 0 	//N/A
#define tdbrs 0	//N/A
#define mrgtd 0	//N/A
#define nomrg 0	//N/A
#define sfr	0	//N/A

/*Columns for CAS measurments*/
#define cpos 0 		//N/A
#define apos 0 		//N/A
#define spos 0 		//N/A

/*Columns for G-M20 measurments*/
#define mpos 0 		//N/A
#define gpos 0 		//N/A


