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

/*GAMA+CAS switch*/
/*Catalogue length*/
#define ls 150000 //Sample length
#define ls2 150000 //Total catalogue length
#define ksize 15 //Magnitude or Mass bins
#define lsize 30 //rp bins
int Nvar=5;
#define Nextra 0 //Number of calculated mophological index variables 
#define NvarArray 6

//pvalV[NvarArray][ksize][lsize][ls]

/*Columns for CAS measurments*/
#define cpos 0
#define apos 1
#define spos 2

/*Columns for G-M20 measurments*/
#define mpos 3
#define gpos 4


/*Columns for GZ morphological features*/
#define smooth 0	//N/A
#define feturs 0	//N/A
#define clmp 0		//N/A
#define clmpsp 0 	//N/A
#define edgeno 0	//N/A
#define baryes 0	//N/A
#define sparms 0	//N/A
#define twarms 0	//N/A
#define mwarms 0	//N/A
#define lwarms 0	//N/A
#define oarm 0		//N/A
#define tarms 0		//N/A
#define threearms 0 //N/A
#define fourarms 0	//N/A
#define morearms 0	//N/A
#define unknarms 0	//N/A
#define noblg 0		//N/A
#define obblg 0 	//N/A
#define doblg 0 	//N/A
#define pmrg 0		//N/A
#define tdbrs 0		//N/A
#define mrgtd 0		//N/A
#define nomrg 0		//N/A
#define awidx 0		//N/A
#define blgidx 0	//N/A
#define	jnblg 0		//N/A	
#define odd 0		//N/A
#define dstrb 0		//N/A 
#define irrgl 0		//N/A
#define cw 0		//N/A
#define acw 0		//N/A
#define sfr	0		//N/A
#define awidx 	0		//N/A
#define blgidx 	0		//N/A
#define armnum 	0		//N/A
#define mclass 	0		//N/A
#define spin 	0		//N/A

char catCode[10]="";

/*GZ Class thresholds*/
float fesmin=0;
float ftmin=0;
float clmpmin=0;
float clmpmax=0;
float clmpspmin=0;
float egmin=0;
float spmin=0;
float spmax=0;
float twamax=0;
float twamin=0;
float mwamax=0;
float mwamin=0;
float lwamax=0;
float lwamin=0;
float oamin=0;
float tamin=0;
float oddmin=0;
float pmin=0;
float irrmin=0;
float irrmax=0;
float dismin=0;
float dismax=0;
float barmin=0;
float barmax=0;
float noblgmin=0;
float jnblgmin=0;
float obblgmin=0;
float doblgmin=0;
float tdmin=0;
float mrgtdmin=0;
float nomrgmin=0;
float awtight=0;
float awloose=0;
