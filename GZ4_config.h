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
#define ls 50000 //Sample length
#define ls2 50000 //Total catalogue length
#define ksize 10//Magnitude or Mass bins
#define lsize 10 //rp bins=
//int Nvar=50; //Without SFR column
int Nvar=51; //This number has to represent the exact number of GZ classification columns in the catalog
#define NvarArray 60 //This is the storage array index number and can be higher than Nvar

//pvalV[NvarArray][ksize][lsize][ls]

/*Columns for GZ4 Candels morphological features of interest*/
#define feturs 1	//Features

#define clmp 6	//Clumpy - Yes
#define clmpsp 17 //Clumpy spiral

#define edgeno 27	//Edge on - No
#define baryes 30	//Bar - Yes
#define sparms 32	//Spiral - Yes
#define twarms 34	//Tight winding arms
#define mwarms 35	//Medium winding arms
#define lwarms 36	//Loose winding arms
#define oarm 37		//1 arm
#define tarms 38	//2 arms
#define threearms 39 //3 arms
#define fourarms 40	//4 arms
#define morearms 41	//4+ arms
#define unknarms 42	//? arms
#define noblg 43		//face-on no-bulge
#define obblg 44 	//face-on obvious bulge
#define doblg 45 	//face-on dominant bulge

#define pmrg 46		//Merging
#define tdbrs 47		//Tidal debris
#define mrgtd 48		//Both
#define nomrg 49		//Neither

#define sfr	50

#define awidx 51    //arm winding tightness index
#define blgidx 52	//bulge domenance index

#define	jnblg 0		//N/A	
#define odd 0		//N/A
#define dstrb 0		//N/A 
#define irrgl 0		//N/A


/*Columns for CW and ACW classifications*/
#define cw 0		//N/A
#define acw 0		//N/A

/*Columns for CAS measurments*/
#define cpos 0		//N/A
#define apos 0		//N/A
#define spos 0		//N/A

/*Columns for G-M20 measurments*/
#define mpos 0		//N/A
#define gpos 0		//N/A



