/*rpmerge by Kevin Casteels
Version 1.0 modified Nov 21,2013
	-First post PhD version!
	-GalaxySelector subroutine

Pervious addditions:
	-Ability to bin by mass or magnitude
	-Fixed angular spec weight bug
	-Implimentation of GZ1_GAMA_CAS catalogue compatibility
	-Added GZ2 parameters such as "anything odd" values
	-Addes GZ1 parameters and CW and ACW pair determination
	-Placed all program variable definitions in config.h file
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"

/*Config file for the GAMA CAS paper*/
//#include "GAMACAS_config.h"

/*Config file for the GZ2 Merger Phases paper*/
#include "GZ2phases_config.h"

/*Config file for the GZ Candles Spiral Orgins paper*/
//#include "GZ4spiral_config.h"



/*Magnitude variable shared between main() and fuctions...*/
float Mi;

/*Functions included after main*/
float ran2(long *idum);
float trapzd(float (*func)(float), float a, float b, int n);
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
float qromb(float (*func)(float), float a, float b);
float e_func(float z);
float rtbis(float (*func)(float), float x1, float x2, float xacc);
float Mlim_func(float z);
float jerror(float total, int count, float values[]);
float xerror(float total, int count, float values[]);
float berror(float value, int count);
float serror(float total, int count, float values[]);
float edgewt(float coord,float coMax,float coMin,float maxRad,float minRad,float DA);
float alimit(float mass);

/*Main program starts:*/
int main(int argc, char **argv)
{
	printf("start\n");

	FILE *galcat;
	FILE *morphcat;
	FILE *matchcat;
	FILE *thetacorr;
	FILE *rpresults;
	FILE *pairlist;
	FILE *pairUrls;
	FILE *histout;
	FILE *casout;
	FILE *ncout;
	FILE *alllist;
	FILE *allUrls;
	FILE *asout;
	FILE *gm20out;

	char results[200],corrout[200],fileCode[200];
	float mlStep,ml1,ml2;

	char matchtmp[1000000]="";
	char urlHtml[10]="";
	char htmlCode1[1000]="";
	char htmlCode2[1000]="";
	char htmlCode3[1000]="";
	char htmlCode4[1000]="";

	int i,j,k,l,m,n,o,p,a,x,t,v,Nvar2;

	long long int rid;
	float rphoto_z;
	float rspec_z;
	float rvelDisp;
	float rap_g;
	float rap_r;
	float rra;
	float rdec;
	float rcol_gr;
	float rfracdev_r;
	float rpetro_g;
	float rpetro_r;
	float rpetroR50_r;
	float rpetroR90_r;
	float rkcorr_r;
	float rab_r;
	float rab_g;
	float rmass;

	long long int mid;
	float mphoto_z;
	float mspec_z;
	float mvelDisp;
	float map_g;
	float map_r;
	float mra;
	float mdec;
	float mcol_gr;
	float mfracdev_r;
	float mpetro_g;
	float mpetro_r;
	float mpetroR50_r;
	float mpetroR90_r;
	float mkcorr_r;
	float mab_r;
	float mab_g;
	float mmass;

	int rspi,rell,runc;

	//CAS readers
	float rc,rce,ras,rase,rs,rse,rcasflag;

	//Pair filter check switch
	int gpcheck[ls];

	long long int id[ls]; 
	float photo_z[ls];
	float spec_z[ls],spec_z2[ls2];
	float velDisp[ls];
	float ap_g[ls];
	float ap_r[ls];
	float ra[ls],ra2[ls2];
	float dec[ls],dec2[ls2];
	float col_gr[ls];
	float fracdev_r[ls];
	float radius[ls];
	float kcorr_r[ls];
	float ab_r[ls],ab_j;
	float ab_g[ls];
	float mass[ls],mass_2[ls2];
	float mlvar[ls];
	float lum[ls];
	float type[ls];
	float con[ls];
	float asy[ls];
	float smo[ls];
	float gin[ls];
	float m20[ls];
	


	//float rmerger,merger[ls],pmV[ksize][lsize][ls],pmT[ksize][lsize],pmB[ksize][lsize][ksize],pm_error;

	float rtype[NvarArray],mtype[NvarArray],pval[NvarArray][ls],pvalV[NvarArray][ksize][lsize][ls],pvalT[NvarArray][ksize][lsize],pvalB[NvarArray][ksize][lsize][ksize],p_error[NvarArray];
	float gvalV[NvarArray][ksize][ls],gvalT[NvarArray][ksize],gvalB[NvarArray][ksize][ksize],g_error[NvarArray];

	int galtot,galnum,galnum2;

	float DM,DA;
	float deltaz,z,deltav;
	float deltara,deltadec,avgdec;
	float theta,rp;
	int mcheck,bcheck,hcheck,pcheck,gcheck,gpcheckTemp,obcheck;
	float rbCut,coro;

	float tcorr[4][100000];

	float gbwtV[ksize][ls],gbradV[ksize][ls],gbmasV[ksize][ls],gblumV[ksize][ls] ={0};
	float wtV[ksize][lsize][ls],dvV[ksize][lsize][ls],rsep[ksize][lsize][ls] ={0};
	float radV[ksize][lsize][ls],radV2[ksize][lsize][ls],vdisV[ksize][lsize][ls] ={0};
	float assV[ksize][lsize][ls],masV[ksize][lsize][ls],lumV[ksize][lsize][ls],psep[ksize][lsize][ls] ={0};
	int gid1[ksize][lsize][ls],gid2[ksize][lsize][ls] ={0};

	float gbmagB[ksize];
	float gbmagBTotal[ksize];
	float gbwtT[ksize] = {0};
	float gbradT[ksize] = {0};
	float gbmasT[ksize] = {0};
	float gblumT[ksize] ={0};
	float radT[ksize][lsize] = {0};
	float radT2[ksize][lsize] = {0};
	float assT[ksize][lsize] = {0};
	float wtT[ksize][lsize] ={0};
	float rpB[ksize][lsize],magB[ksize][lsize];
	float rpT[ksize][lsize] = {0};
	float vdisT[ksize][lsize] ={0};
	float dvT[ksize][lsize] = {0};
	float masT[ksize][lsize] = {0};
	float lumT[ksize][lsize] ={0};

	int gbN[ksize] = {0};
	int gbNTotal[ksize] = {0};
	int gbID[ksize][ls];
	int rpN[ksize][lsize] ={0};
	int vdisN[ksize][lsize] ={0};
	
	int asyN[ksize] ={0};
	int asyID[ksize][ls];
	float asyT[ksize] = {0};
	float asyradT[ksize] = {0};
	float asymasT[ksize] = {0};
	float asylumT[ksize] = {0};
	float asyV[ksize][ls],asyradV[ksize][ls],asymasV[ksize][ls],asylumV[ksize][ls] = {0};

	int gm20N[ksize] ={0};
	int gm20ID[ksize][ls];
	float gm20T[ksize] = {0};
	float gm20radT[ksize] = {0};
	float gm20masT[ksize] = {0};
	float gm20lumT[ksize] = {0};
	float gm20V[ksize][ls],gm20radV[ksize][ls],gm20masV[ksize][ls],gm20lumV[ksize][ls] = {0};


	float tBin,rBin,mBin,vBin;
	int tCount=0;
	int rCount=0;
	int mCount=0;
	int maCount=0;
	int vCount=0;

	float minDist,maxDist,maxTheta;
	float Dmax,Dobj,Vmax,Vobj,Vwt,Zwt,DECwt,RAwt;

	int rawTal=0;
	int rawAss=0;
	int rawGm20=0;
	int rawPairs=0;
	float wtTal=0;
	float wtAss=0;
	float wtGm20=0;
	float wtPairs=0;

	float raMax=0;
	float decMax=0;
	float raMin=0;
	float decMin=0;

	float zsum=0;
	float colour=0;

	float dV_error,wt_error,a_error,r_error,r2_error,gbr_error,gbwt_error,gbm_error,m_error,asy_error,gbl_error,l_error;
	float nc_error;
	int argnum;
	long long int complete,total;

	argnum=argc;
	printf("argc = %d\n",argc);
	
	
	if(argnum==1)
	{	
		//printf("Please include the catalogue name as the first argument. The second argument is the angular spectroscopic completeness correction file.\n");
		printf("Usage: ./rpmerge_v* <input catalogue> <angular spectroscopic correction file>\n");
	}


	if(LorM==1)
	{
		mlStep=lStep;
		ml1=m1;
		ml2=m2;
	}
	else
	{
		mlStep=mStep;
		ml1=mass1;
		ml2=mass2;
	}

	printf("Creating output files\n");


	if(htmlSwitch==1){
		sprintf(urlHtml,"html");
	}
	else if(wget<0){
		sprintf(urlHtml,"wget");
		htmlSwitch=0;
	}
	else{
		sprintf(urlHtml,"urls");
	}

	sprintf(fileCode,"%.0f%.0f_%.1f-%.1f_%.1f:%.1f_%.1f_%.0f_%.2f-%.2f_ao%d_fb%d_hl%d_gr%d_pr%d_ck%d_ri%d_rx%.1f_rs%.1f_rp%d_%s",(fabs(ml1)),ml2,dMmin,dMmax,ratiomin,ratiomax,mlStep,maxv,z1,z2,ao,fb,hl,grestrict,prestrict,contcheck,removeint,radx,rStep,rpswitch,argv[1]);



	//strcpy(results,"1622_01_500_fb-1");
	sprintf(results,"%s.dat",fileCode);
	rpresults=fopen(results,"w");

	sprintf(results,"%s.nc",fileCode);
	ncout=fopen(results,"w");

	sprintf(results,"%s.pairs.list",fileCode);
	pairlist=fopen(results,"w");

	sprintf(results,"%s.all.list",fileCode);
	alllist=fopen(results,"w");

	sprintf(results,"%s.hist",fileCode);
	histout=fopen(results,"w");

	sprintf(results,"%s.cas",fileCode);
	casout=fopen(results,"w");

	sprintf(results,"%s.as",fileCode);
	asout=fopen(results,"w");

	sprintf(results,"%s.gm20",fileCode);
	gm20out=fopen(results,"w");

	sprintf(results,"%s.pairs.%s",fileCode,urlHtml);
	pairUrls=fopen(results,"w");

	sprintf(results,"%s.all.%s",fileCode,urlHtml);
	allUrls=fopen(results,"w");


	if(matchSwitch==1){
		sprintf(results,"%s.morph.cat",fileCode);
		morphcat=fopen(results,"w");
	}
	

	fprintf(rpresults,"\nPair velocity difference limits:\n");
	fprintf(rpresults,"minv=%f  maxv=%f\n\n",minv,maxv);

	fprintf(rpresults,"Sample magnitude limits:\n");
	fprintf(rpresults,"m1=%f  m2=%f  lStep=%f\n\n",m1,m2,lStep);

	fprintf(rpresults,"Redshift limits:\n");
	fprintf(rpresults,"z1=%f  z2=%f\n\n",z1,z2);

	fprintf(rpresults,"Pair projected separation limits:\n");
	fprintf(rpresults,"r1=%f  r2=%f  rStep=%f\n\n",r1,r2,rStep);

	fprintf(rpresults,"Maximum pair magnitude separation:\n");
	fprintf(rpresults,"dMmin=%f  dMmax=%f\n\n",dMmin,dMmax);

	fprintf(rpresults,"Maximum pair mass ratio:\n");
	fprintf(rpresults,"ratiomin=%f  ratiomax=%f\n\n",ratiomin,ratiomax);
	
	fprintf(rpresults,"fb=%d  hl=%d\n\n",fb,hl);

	fprintf(rpresults,"Adefault=%f\n",Adefault);

	fprintf(rpresults,"apMagMax=%f\n\n",apMagMax);
	printf("apMagMax=%f\n\n",apMagMax);


	if(wget>0)
	{
		fprintf(allUrls,"<h2>Apparent Mag Limit = %f</h2>",apMagMax);
		fprintf(pairUrls,"<h2>Apparent Mag Limit = %f</h2>",apMagMax);

		fprintf(allUrls,"<h2>Sample AB Mag Range: %f < log10(Mass) < %f</h2>",m1,m2);
		fprintf(pairUrls,"<h2>Sample AB Mag Range: %f < log10(Mass) < %f</h2>",m1,m2);

		fprintf(allUrls,"<h2>Sample Mass Range: %f < log10(Mass) < %f</h2>",mass1,mass2);
		fprintf(pairUrls,"<h2>Sample Mass Range: %f < log10(Mass) < %f</h2>",mass1,mass2);

		fprintf(allUrls,"<h2>Redshift Range: %f < z < %f</h2>",z1,z2);
		fprintf(pairUrls,"<h2>Redshift Range: %f < z < %f</h2>",z1,z2);
	}

	

	printf("%s\n",argv[1]);
	galcat=fopen(argv[1],"r");

	t=0;
	i=0;
	x=0;
	printf("Scanning file...\n");
	while (!feof(galcat))
	{
		fscanf(galcat,"%lld %f %f %f %f %f %f %f %f %f %f %f",&rid,&rra,&rdec,&rpetroR50_r,&rpetroR90_r,&rpetro_g,&rpetro_r,&rspec_z,&rab_g,&rab_r,&rvelDisp,&rmass);

		//printf("%lld %f %f\n",rid,rra,rdec);

		for(n=0;n<Nvar;n++)
		{
			fscanf(galcat," %f",&rtype[n]);
		}

		if(awidx !=0 && blgidx !=0 && armnum !=0){
			rtype[awidx] = 0.0*rtype[twarms] + 0.5*rtype[mwarms] + 1.0*rtype[lwarms];
			rtype[blgidx] = 0.0*rtype[noblg] + 0.5*rtype[obblg] + 1.0*rtype[doblg];


			/*if(rtype[oarm]>0.5){rtype[armnum] = 1.0;}
			else if(rtype[tarms]>0.5){rtype[armnum] = 2.0;}
			else if(rtype[threearms]>0.5){rtype[armnum] = 3.0;}
			else if(rtype[fourarms]>0.5){rtype[armnum] = 4.0;}
			else if(rtype[morearms]>0.5){rtype[armnum] = 10.0;}
			else {rtype[armnum] = 0.0;}
			*/

			if(rtype[oarm]>rtype[tarms] && rtype[oarm]>type[threearms] && rtype[oarm]>type[fourarms] && rtype[oarm]>type[morearms] && rtype[oarm]>type[unknarms]){rtype[armnum] = 1.0;}
			else if(rtype[tarms]>rtype[oarm] && rtype[tarms]>type[threearms] && rtype[tarms]>type[fourarms] && rtype[tarms]>type[morearms] && rtype[tarms]>type[unknarms]){rtype[armnum] = 2.0;}
			else if(rtype[threearms]>rtype[oarm] && rtype[threearms]>type[tarms] && rtype[threearms]>type[fourarms] && rtype[threearms]>type[morearms] && rtype[threearms]>type[unknarms]){rtype[armnum] = 3.0;}
			else if(rtype[fourarms]>rtype[oarm] && rtype[fourarms]>type[tarms] && rtype[fourarms]>type[threearms] && rtype[fourarms]>type[morearms] && rtype[fourarms]>type[unknarms]){rtype[armnum] = 4.0;}
			else if(rtype[morearms]>rtype[oarm] && rtype[morearms]>type[tarms] && rtype[morearms]>type[threearms] && rtype[morearms]>type[fourarms] && rtype[morearms]>type[unknarms]){rtype[armnum] = 5.0;}
			else if(rtype[unknarms]>rtype[oarm] && rtype[unknarms]>type[tarms] && rtype[unknarms]>type[threearms] && rtype[unknarms]>type[fourarms] && rtype[unknarms]>type[morearms]){rtype[armnum] = 10.0;}
			else {rtype[armnum] = 0.0;}


			if ((rtype[irrgl]>irrmin && rtype[odd]>oddmin)){rtype[mclass] = 7.0;} //irrgular
			else if(rtype[feturs]>ftmin && rtype[edgeno]>egmin && rtype[sparms]<spmax && rtype[irrgl]<irrmax && rtype[dstrb]<dismax){rtype[mclass] = 5.0;}
			else if(rtype[feturs]>0.2){rtype[mclass] = 3.0;}
			else if(rtype[smooth]>0.5){rtype[mclass] = 1.0;}
			else {rtype[mclass] = 0.0;}

			if(rtype[cw]>rtype[acw]){rtype[spin]=1.0;}
			else if(rtype[cw]<rtype[acw]){rtype[spin]=-1.0;}
			else{rtype[spin]=0.0;}
		}



		//printf("\n");
		//fscanf(galcat,"\n");

		//Global GZ2 Class Thresholdshttp://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=203.535095&dec=0.477853&scale=0.133479&width=424&height=424
		if(grestrict==0){gcheck=1;}

		else if(grestrict==10 && rtype[feturs]>ftmin && rtype[edgeno]>egmin && rtype[sparms]>spmin &&  rtype[lwarms]<lwamax){gcheck=1;}

		else if(grestrict==20 && rtype[feturs]>ftmin && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[lwarms]>lwamin){gcheck=1;}
		else if(grestrict==21 && ((rtype[feturs]*rtype[edgeno]*rtype[sparms])>fesmin)  && (rtype[lwarms]>lwamin) && (rtype[oarm]>oamin || rtype[tarms]>tamin) ){gcheck=1;}
		else if(prestrict==22 && (rtype[pmrg]>pmin)){gcheck=1;}

		else if(grestrict==30 && rtype[feturs]>ftmin && rtype[edgeno]>egmin && (rtype[irrgl]>irrmin || rtype[dstrb]>dismin) && rtype[pmrg]>pmin && rtype[odd]>oddmin){gcheck=1;}
		else if(grestrict==31 && (rtype[irrgl]>irrmin || rtype[dstrb]>dismin) && rtype[pmrg]>pmin && rtype[odd]>oddmin){gcheck=1;}
		else if(grestrict==32 && (rtype[irrgl]>irrmin || rtype[dstrb]>dismin) && rtype[pmrg]>pmin && rtype[odd]>oddmin && rtype[baryes]>barmin){gcheck=1;}

		else if(grestrict==35 && rtype[pmrg]>pmin){gcheck=1;}
		else if(grestrict==36 && rtype[tdbrs]>tdmin){gcheck=1;}
		else if(grestrict==37 && rtype[mrgtd]>mrgtdmin){gcheck=1;}

		else if(grestrict==38 && (rtype[tdbrs]>tdmin || rtype[mrgtd]>mrgtdmin)){gcheck=1;}


		else if(grestrict==40 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin){gcheck=1;}
		else if(grestrict==41 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin){gcheck=1;}
		else if(grestrict==42 && rtype[feturs]>ftmin && rtype[clmp]>clmpmin && rtype[clmpsp]>clmpspmin){gcheck=1;}
		else if(grestrict==43 && rtype[feturs]>ftmin && (rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin || (rtype[clmp]>clmpmin && rtype[clmpsp]>clmpspmin))) {gcheck=1;}

		else if(grestrict==44 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && (rtype[pmrg]>pmin || rtype[tdbrs]>tdmin || rtype[mrgtd]>mrgtdmin)){gcheck=1;}
		else if(grestrict==45 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && (rtype[pmrg]<pmin && rtype[tdbrs]<tdmin && rtype[mrgtd]<mrgtdmin)){gcheck=1;}



		else if(grestrict==50 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[twarms]>twamin){gcheck=1;}
		else if(grestrict==51 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[twarms]>twamin && rtype[baryes]>barmin){gcheck=1;}
		else if(grestrict==52 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]<awtight){gcheck=1;}
		else if(grestrict==53 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]<awtight && rtype[baryes]>barmin){gcheck=1;}


		else if(grestrict==60 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[twarms]<twamax && rtype[lwarms]<lwamax){gcheck=1;}
		else if(grestrict==61 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[twarms]<twamax && rtype[lwarms]<lwamax && rtype[baryes]>barmin){gcheck=1;}
		else if(grestrict==62 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]>awtight && rtype[awidx]<awloose){gcheck=1;}
		else if(grestrict==63 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]>awtight && rtype[awidx]<awloose && rtype[baryes]>barmin){gcheck=1;}


		else if(grestrict==70 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[lwarms]>lwamin){gcheck=1;}
		else if(grestrict==71 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[lwarms]>lwamin && rtype[baryes]>barmin){gcheck=1;}
		else if(grestrict==72 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]>awloose){gcheck=1;}
		else if(grestrict==73 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]>spmin && rtype[awidx]>awloose && rtype[baryes]>barmin){gcheck=1;}


		else if(grestrict==80 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]<spmax && rtype[irrgl]<irrmax && rtype[dstrb]<dismax){gcheck=1;}
		else if(grestrict==81 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[sparms]<spmax && rtype[irrgl]<irrmax && rtype[dstrb]<dismax && rtype[baryes]>barmin){gcheck=1;}

		else if(grestrict==110 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[noblg]>noblgmin){gcheck=1;}
		else if(grestrict==111 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[noblg]>noblgmin && rtype[baryes]>barmin){gcheck=1;}

		else if(grestrict==120 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[jnblg]>jnblgmin){gcheck=1;}
		else if(grestrict==121 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[jnblg]>jnblgmin && rtype[baryes]>barmin){gcheck=1;}

		else if(grestrict==130 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[obblg]>obblgmin){gcheck=1;}
		else if(grestrict==131 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[obblg]>obblgmin && rtype[baryes]>barmin){gcheck=1;}

		else if(grestrict==140 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[doblg]>doblgmin){gcheck=1;}
		else if(grestrict==141 && rtype[feturs]>ftmin && rtype[clmp]<clmpmax && rtype[edgeno]>egmin && rtype[doblg]>doblgmin && rtype[baryes]>barmin){gcheck=1;}



		// Limit to A>Alimit and A>S
		else if(grestrict==201 && rtype[apos]>alimit(rmass) && rtype[apos]>rtype[spos]){gcheck=1;}
		// Limit to A>Alimit and S>A
		else if(grestrict==211 && rtype[apos]>alimit(rmass) && rtype[spos]>rtype[apos]){gcheck=1;}
		
		else{gcheck=0;}


		//If we are only filtering paired galaxies:
		if(prestrict>0 && grestrict>0){
			gpcheckTemp = gcheck;
			gcheck=1;
		}

		//printf("%f  %f\n", rtype[feturs], rtype[edgeno]);

		//Set all redshifts to the z-range midpoint... useful for finding false pairs caused through over de-blending.
		//587722981736579107   171.99093627929688   -1.2140086889266968   8.175695       15.858483      16.691872    15.696057    9999         9999            9999            9999         9999         0.648         0.336         0.015         0.214         0.786         0.0           1.0           0.3           0.7           0.0           0.5           0.4           0.1           0.39          0.61          0.074         0.148         0.778         0.375         0.0           0.0           0.062         0.312         0.25          0.0           0.333         0.333         0.333         0.667         0.333         0.0           0.0           0.667         0.0           0.0           0.0           0.333         0.4      0.073    0.018     0.436      0.036    0.036    0.527   

		/*if(rspec_z==9999){
			rspec_z = 0.000099;
			rmass = 10.09999;
			rab_r = -20.09999;
			rab_g = -20.09999;
			rvelDisp = 100.09999;
		}*/
	
		//printf("%d %f %f %d\n",grestrict,rtype[irrgl],rtype[dstrb],gcheck);

		//Only select galaxies with low pmrg values
		//if (rtype[pmrg]<0.1 && rspec_z < z2){

		//printf("%f  %f  %f  %f  %d\n",rab_r,rspec_z, rpetro_r,rmass, gcheck );

		//Count the total number of galaxies in range, without considering morphology
		if((rab_r >= m2) && (rab_r <= m1) && (rspec_z >= z1) && (rspec_z <= z2) && 
		(rpetro_r<=apMagMax) && (rmass >= mass1) && (rmass <= mass2))
		{
			t = t+1;

			/*Absolute magnitude bin and array number*/
			if(LorM==1){
				k=floor((rab_r-ml2)/mlStep);
				mBin=(floor(rab_r/mlStep))*mlStep+mlStep;
			}
			else{
				k=floor((rmass-ml1)/mlStep);
				mBin=(floor(rmass/mlStep))*mlStep+mlStep;
			}
			if(k > maCount){maCount = k;}

			gbNTotal[k]++;                       //Add unweighted galaxy
			gbmagBTotal[k] = mBin - mlStep/2.;       //magnitude bin

		}


		obcheck=1;
		/*if((rra > (12.0*360/24.0)) && (rra < (24.0*360/24.0)) && (rdec < 15.0)){obcheck =1;}
		else{obcheck=0;}*/
		


		//Select the galaxies in the sample range and with correct morophology
		if((rab_r >= m2) && (rab_r <= m1) && (rspec_z >= z1) && (rspec_z <= z2) && 
		(rpetro_r<=apMagMax) && (rmass >= mass1) && (rmass <= mass2) && gcheck==1 && obcheck==1)
		{


			//Output a catalog of only the galaxies which meet the morphology, mass, luminosity and redshift requirements
			if(matchSwitch==1){
				fprintf(morphcat,"%lld %f %f %f %f %f %f %f %f %f %f %f",rid,rra,rdec,rpetroR50_r,rpetroR90_r,rpetro_g,rpetro_r,rspec_z,rab_g,rab_r,rvelDisp,rmass);

				for(n=0;n<(Nvar+Nextra);n++)
				{
					fprintf(morphcat," %f",rtype[n]);
				}

				fprintf(morphcat,"\n");
			}

			gpcheck[i]=gpcheckTemp;

			id[i]=rid;
			ra[i]=rra;
			dec[i]=rdec;
			spec_z[i]=rspec_z;
			//photo_z[i]=rphoto_z;
			ap_r[i]=rpetro_r;

			//radius[i]=rpetroR50_r;
			radius[i]=rpetroR90_r;
			ab_r[i]=rab_r;
			ab_g[i]=rab_g;

			//Use either luminosity or mass for binning:
			if(LorM==1)
			{
				mlvar[i]=rab_r;
			}
			else
			{
				mlvar[i]=rmass;
			}

			mass[i]=pow(10.0,rmass);
			lum[i]=pow(10.0,rpetro_r/2.5);
			velDisp[i]=rvelDisp; 

			col_gr[i]=rab_g-rab_r;

			fracdev_r[i]=1;

			for(n=0;n<Nvar;n++)
			{
				pval[n][i] = rtype[n];
			}

			if(awidx !=0 && blgidx !=0 && armnum !=0){
				pval[awidx][i] = rtype[awidx];
				pval[blgidx][i] = rtype[blgidx];
			}


			if(cason==1)
			{
				con[i]=pval[cpos][i];
				asy[i]=pval[apos][i];
				smo[i]=pval[spos][i];
				gin[i]=pval[gpos][i];
				m20[i]=pval[mpos][i];
			}
			else 
			{
				con[i]=0;
				asy[i]=0;
				smo[i]=0;
				gin[i]=0;
				m20[i]=0;
			}
			//printf("%lld  %f  %f\n",id[i],pval[cw][i],pval[acw][i]);

			/*Find RA and DEC boundary limits of sample*/
			if(rdec > decMax){decMax = rdec;}
			if(rdec < decMin){decMin = rdec;}
		
			if(rra > raMax){raMax = rra;}
			if(rra < raMin){raMin = rra;}

			zsum = zsum + rspec_z;
			colour = colour + rab_g - rab_r;
			
			i++;
		}

		/*Store RAs and DECs of all galaxies to calculate angular incompleteness*/
		if(rpetro_r<=apMagMax && argnum == 2 && zComp!=1.0)
		{
			ra2[x]=rra;
			dec2[x]=rdec;
			spec_z2[x]=rspec_z;
			mass_2[x]=rmass;
			x++;
		}



	}

	//Close the catalog file:
	fclose(galcat);

	//Record the galaxy counts
	galtot =t;
	galnum =i;
	galnum2=x;


	/**********************************************************************************************/
	/*Create a matched mass and redshift catalog a morphologically selected sample*/
	/**********************************************************************************************/

	if(matchSwitch==1)
	{

		galcat=fopen(argv[1],"r");

		fclose(morphcat);
		sprintf(results,"%s.morph.cat",fileCode);
		morphcat=fopen(results,"r");

		sprintf(results,"%s.match.cat",fileCode);
		matchcat=fopen(results,"w");

		i=0;
		printf("Creating matched random sample to morphlogical sample...\n");
		while (!feof(morphcat) && i<galnum)
		{
			fscanf(morphcat,"%lld %f %f %f %f %f %f %f %f %f %f %f",&mid,&mra,&mdec,&mpetroR50_r,&mpetroR90_r,&mpetro_g,&mpetro_r,&mspec_z,&mab_g,&mab_r,&mvelDisp,&mmass);

			for(n=0;n<Nvar;n++)
			{
				fscanf(morphcat," %f",&mtype[n]);
			}

			printf("Working on %lld  %f  %f\n",mid,mmass,mspec_z);

			int rescan=1;
			while(rescan>0){
				//sprintf(matchtmp,"");
				l=0;
				while (!feof(galcat) && l<matchNum)
				{
					fscanf(galcat,"%lld %f %f %f %f %f %f %f %f %f %f %f",&rid,&rra,&rdec,&rpetroR50_r,&rpetroR90_r,&rpetro_g,&rpetro_r,&rspec_z,&rab_g,&rab_r,&rvelDisp,&rmass);

					for(n=0;n<Nvar;n++)
					{
						fscanf(galcat," %f",&rtype[n]);
					}


					if(rid!=mid && (rspec_z < (mspec_z+zslop)) && (rspec_z > (mspec_z-zslop)) && (rmass<(mmass+mslop)) && (rmass>(mmass-mslop)) &&
						(rab_r >= m2) && (rab_r <= m1) && (rspec_z >= z1) && (rspec_z <= z2) && (rpetro_r<=apMagMax) && (rmass >= mass1) && (rmass <= mass2)){

						sprintf(matchtmp + strlen(matchtmp),"%lld %f %f %f %f %f %f %f %f %f %f %f",rid,rra,rdec,rpetroR50_r,rpetroR90_r,rpetro_g,rpetro_r,rspec_z,rab_g,rab_r,rvelDisp,rmass);

						for(n=0;n<Nvar;n++)
						{
							sprintf(matchtmp + strlen(matchtmp)," %f",rtype[n]);
						}
						sprintf(matchtmp + strlen(matchtmp),"\n");
						//printf("Match at %lld  %f  %f\n",rid,rmass,rspec_z);
						l++;
					}
				}
				if(l<matchNum){
					printf("%d matches found. Trying again\n",l);
					fclose(galcat);
					galcat=fopen(argv[1],"r");
					rescan=1;
				}
				else{
					printf("%d matches found\n",l);
					rescan=0;
				}

			}
			fprintf(matchcat,"%s",matchtmp);
			i++;			
		}
		fclose(morphcat);
		fclose(galcat);
		fclose(matchcat);
	}
	/**********************************************************************************************/



	//Add two array indicies if awidx and blgidx are active
	if(awidx !=0 && blgidx !=0){
		Nvar=Nvar+Nextra;
	}

	printf("mean z=%f\n",zsum/galnum);

	fprintf(rpresults,"Mean Redshift = %f\n",zsum/galnum);
	printf("galnum=%d\n",galnum);
	fprintf(rpresults,"Galaxy Count = %d\n",galnum);

	if(wget>0){
		fprintf(allUrls,"<h2>Mean Redshift = %f</h2>",zsum/galnum);
		fprintf(pairUrls,"<h2>Mean Redshift = %f</h2>",zsum/galnum);

		fprintf(allUrls,"<h2>Galaxy Count = %d</h2>",galnum);
		fprintf(pairUrls,"<h2>Galaxy Count = %d</h2>",galnum);

		fprintf(allUrls,"<h2>Total Galaxy Count = %d</h2>",galtot);
		fprintf(pairUrls,"<h2>Total Galaxy Count = %d</h2>",galtot);
	}


	/*Minimum and maximum distances*/
	minDist = DH*qromb(e_func,0,z1);
	maxDist = DH*qromb(e_func,0,z2);

	/*DA*/
	DA = minDist/(1+z1);

	/*Maximum angular separation*/
	maxTheta = (r2*360)/(2*pi*DA*1000);


	/**********************************************************************************************/
	/*Calculate angular spectroscopic incompleteness:*/
	/**********************************************************************************************/

	if(zComp==1.0)
	{
		printf("100 percent angular spectroscopic completeness\n");
		for(l=0;l<9999;l++){tcorr[3][l]=1.0;}
	}
	else if(argnum > 2)
	{
		thetacorr=fopen(argv[2],"r");

		printf("Scanning theta correction file...\n");
		l=0;
		while (!feof(thetacorr))
		{
			fscanf(thetacorr,"%f %f %f %f\n",&tcorr[0][l],&tcorr[1][l],&tcorr[2][l],&tcorr[3][l]);
			l++;
		}
		tCount=l-1;
	}	
	else
	{
		fprintf(rpresults,"galnum2=%d\n",galnum2);
	
		//zComp=galnum*1.0/galnum2;
	
		printf("galnum2=%d\n",galnum2);


		printf("Calculating angular separation corrections...\n");
		tCount=0;
		complete=0;
		total=0;
		for(i=0;i<galnum2;i++)
		{
			//printf("%d\n",i);
			for(j=0;j<galnum2;j++)
			{
				if(j==i)j++;
	
				/*To speed up search, only consider galaxies within certain coordinate range*/
				deltara = fabs(ra2[i]-ra2[j]);
				deltadec = fabs(dec2[i]-dec2[j]);
	
				if(deltara <= maxTheta && deltadec <=maxTheta)
				{

					avgdec = (dec2[j]+dec2[i])/2;
					theta = sqrt(pow(deltara*cos(avgdec*pi/180),2)+pow(deltadec,2))*3600;
		
					/*Calculate the angle bin and array number*/
					tBin = (floor(theta/tStep))*tStep+tStep;
					n = floor(theta/tStep);
	
					/*Place values in the array*/
					tcorr[0][n]=tBin;

					//if((mass_2[i]==9999 || mass_2[j]==9999) || (spec_z2[i]>10 || spec_z2[j]>10 || spec_z2[i]<0 || spec_z2[j]<0) ){printf("%f  %f  %f  %f\n",spec_z2[i],spec_z2[j],mass_2[i],mass_2[j]);}

					//if(spec_z2[i]!=9999 && spec_z2[j]!=9999){tcorr[1][n]++;}
					//if(spec_z2[i]<10 && spec_z2[j]<10 && spec_z2[i]>0 && spec_z2[j]>0){tcorr[1][n]++;}
					if(mass_2[i]!=9999 && mass_2[j]!=9999){tcorr[1][n]++;complete++;}
					tcorr[2][n]++;
					tcorr[3][n]=tcorr[1][n]/tcorr[2][n];
					total++;
				        //printf("%f %f %d %d\n",spec_z[j],tcorr[3][n],galnum,j);
					if(n > tCount){tCount = n;}
				}
			}	
		}
		zComp = sqrt((1.0*complete)/(1.0*total));

		sprintf(corrout,"tcorr_%s_%.0f.dat",argv[1],r2);

		thetacorr=fopen(corrout,"w");
		for(l=0;l<tCount;l++)
		{
			fprintf(thetacorr,"%f  %f  %f  %f\n",tcorr[0][l],tcorr[1][l],tcorr[2][l],tcorr[3][l]);
		}
			
	
	
	}

	printf("complete = %lld  total = %lld  zComp = %f\n",complete,total,zComp);
    	fprintf(rpresults,"zComp = %f\n",zComp);


	/**********************************************************************************************/

	printf("Finding galaxy pairs...\n");

	for(i=0;i<galnum;i++)
	{

		printf("\r%.0f %%",i*100.0/galnum);

		/*if(LorM==1){
			k=floor((mlvar[i]-ml2)/mlStep);
			mBin=(floor((mlvar[i])/mlStep))*mlStep + mlStep/2;
		}
		else{
			k=floor((mlvar[i]-ml1)/mlStep);
			mBin=(floor((mlvar[i])/mlStep))*mlStep + mlStep/2;
		}*/

		//printf("%d %f\n",k,mBin);

	        /*Calculate Vmax weights for current galaxy*/
		if(vmaxon==1)
		{		

			/*Vwt calculation*/		
			Mi=ab_r[i];		
			Dobj = rtbis(Mlim_func,z1,z2,0.0000001);			
			Dmax = DH*qromb(e_func,0,Dobj);
			Vmax = pow(maxDist,3)- pow(minDist,3);
			Vobj = pow(Dmax,3)- pow(minDist,3);
			Vwt = Vobj/Vmax;

			/*if the galaxy is observable throughout the sample its weight is 1*/
			if(Vwt > 1.) Vwt=1.;
		}
		else{Vwt=1.;}

		/*Combine Vmax and global incompleteness weight*/
		Vwt=Vwt*zComp;
		
		/*DA*/
		DA = DH*qromb(e_func,0,spec_z[i])/(1+spec_z[i]);

		/*Absolute magnitude bin and array number*/
		if(LorM==1){k=floor((mlvar[i]-ml2)/mlStep);}
		else{k=floor((mlvar[i]-ml1)/mlStep);}

		mBin=(floor(mlvar[i]/mlStep))*mlStep+mlStep;
		if(k > mCount){mCount = k;}


		gbN[k]++;                       //Add unweighted galaxy

		gbmagB[k] = mBin - mlStep/2.;       //magnitude bin

		//printf("%f  %f\n",gbwtT[k],1./Vwt);
		gbwtT[k] = gbwtT[k] + 1./Vwt;   //Add weighted galaxy
		//printf("%f  %f  %d  %d\n",gbwtT[k],1./Vwt,gbN[k],k);
		gbwtV[k][gbN[k]]=1./Vwt;	//All weights
		//printf("%f\n",gbwtV[k][gbN[k]]);

		gbradT[k] = gbradT[k] + (radius[i]*2*pi*DA*1000)/(360*3600);	//Total galaxy radii
		gbradV[k][gbN[k]] = (radius[i]*2*pi*DA*1000)/(360*3600);		//individual radii

		gbmasT[k] = gbmasT[k] + mass[i];	//Total galaxy mass
		gbmasV[k][gbN[k]] = mass[i];		//Individual masses

		gblumT[k] = gblumT[k] + lum[i];		//Total galaxy luminosities
		gblumV[k][gbN[k]] = lum[i];		//Individual luminosites

		//printf("\33[2K\r");
		//printf("%f  %f  %f  %f  %d  %d\n",gbwtT[k],gbwtV[k][gbN[k]],gbradT[k],gbradV[k][gbN[k]],k,i);
		gbID[k][gbN[k]] = i;


		//**************************************************************
		for(n=0;n<Nvar;n++)
		{
			gvalT[n][k] = gvalT[n][k] + pval[n][i];
			gvalV[n][k][gbN[k]] = pval[n][i];
			p=floor(pval[n][i]/pStep);
			gvalB[n][k][p]++;
		}
		//**************************************************************
	

		rawTal++;
		wtTal=wtTal + 1./Vwt;


		/*Count the assymetric galaxies:*/
		if(cason==1 && asy[i] > alimit(mass[i]) && asy[i]>smo[i])
		//if(cason==1 && (gin[i] > (-0.2)*m20[i] + 0.384) && asy[i] > alimit(mass[i]) && asy[i]>smo[i])
		{

			asyN[k]++;                       //Add unweighted galaxy
	
			asyT[k] = asyT[k] + 1./Vwt;   //Add weighted galaxy
			asyV[k][asyN[k]]=1./Vwt;	//All weights

			asyradT[k] = asyradT[k] + (radius[i]*2*pi*DA*1000)/(360*3600);	//Total galaxy radii
			asyradV[k][asyN[k]] = (radius[i]*2*pi*DA*1000)/(360*3600);		//individual radii
	
			asymasT[k] = asymasT[k] + mass[i];	//Total galaxy mass
			asymasV[k][asyN[k]] = mass[i];		//Individual masses

			asylumT[k] = asylumT[k] + lum[i];	//Total galaxy luminosities
			asylumV[k][asyN[k]] = lum[i];		//Individual luminosities

			asyID[k][asyN[k]] = i;

			rawAss++;
			wtAss = wtAss + 1./Vwt;
		}

		/*Count the G-M20 galaxies:*/
		if(cason==1 && (gin[i] > (-0.4)*m20[i] + 0.384))
		//if(cason==1 && (gin[i] > (-0.4)*asy[i] + 0.66))
		{

			gm20N[k]++;                       //Add unweighted galaxy
	
			gm20T[k] = gm20T[k] + 1./Vwt;   //Add weighted galaxy
			gm20V[k][gm20N[k]]=1./Vwt;	//All weights

			gm20radT[k] = gm20radT[k] + (radius[i]*2*pi*DA*1000)/(360*3600);	//Total galaxy radii
			gm20radV[k][gm20N[k]] = (radius[i]*2*pi*DA*1000)/(360*3600);		//individual radii
	
			gm20masT[k] = gm20masT[k] + mass[i];	//Total galaxy mass
			gm20masV[k][gm20N[k]] = mass[i];	//Individual masses

			gm20lumT[k] = gm20lumT[k] + lum[i];	//Total galaxy luminosities
			gm20lumV[k][gm20N[k]] = lum[i];		//Individual luminosities

			gm20ID[k][gm20N[k]] = i;

			rawGm20++;
			wtGm20 = wtGm20 + 1./Vwt;
		}


		/*Check if galaxy [i] is a pair of galaxy [j]*/
		for(j=0;j<galnum;j++)
		{

			if(j==i)j++;
			//printf("%lld %lld\n",id[i],id[j]);
			/*To speed up search, only consider galaxies within certain coordinate range*/
			deltara = fabs(ra[i]-ra[j]);
			deltadec = fabs(dec[i]-dec[j]);

			if(deltara <= maxTheta && deltadec <=maxTheta)
			{

				/*Now check velocity difference*/	
				deltaz = fabs(spec_z[j]-spec_z[i]);
				z = (spec_z[j]+spec_z[i])/2;
				deltav = cs*deltaz/(1+z);
	
				//printf("dV: %f\n",deltav);
	
				/*velocity condition*/
				if (deltav >= minv && deltav <= maxv)
				{
	


					if(cason!=1){mcheck=1;}
					else if((cason==1) && (ao==0)){mcheck=1;}
					else if((cason==1) && (ao==1) && ((asy[i]>=alimit(mass[i]) && asy[i]>smo[i]) || (asy[j]>=alimit(mass[j]) && asy[j]>smo[j]))){mcheck=1;}
					else if((cason==1) && (ao==2) && ((asy[i]>=alimit(mass[i]) && asy[i]>smo[i]) && (asy[j]>=alimit(mass[j]) && asy[j]>smo[j]))){mcheck=1;}
					else if((cason==1) && (ao==3) && (asy[i]>=alimit(mass[i]) && asy[i]>smo[i])){mcheck=1;}
					else {mcheck=0;}

		
					/*A simpler linear colour cut in (g-r) is used by Skibba et al. ():
					http://adsabs.harvard.edu/abs/2009MNRAS.399..966S
					(g-r) = 0.8 - 0.03 (Mr + 20)*/
			
					/*
					rbCut=0.8-0.03*(ab_r[i]+20);
			
			
					if((col_gr[i]-rbCut)>0)rb_i=1;
					else if((col_gr[i]-rbCut)<0)rb_i=2
					*/

		


					//*************************Control Sample*****************************
					//This is a check for asymm and pmerge contamination at small physical separations...

					if(contcheck==1)
					{
						//z is set to the lower z of the pair.
						//if(spec_z[i]>spec_z[j])z=spec_z[j];
						z=spec_z[i];
						ab_j=ab_r[j];
					}
					else if(contcheck==2)
					{
						//DM = DH*qromb(e_func,0,spec_z[j]);
						//ab_j = ap_r[j] -5*log10(DM*(1+spec_z[i])) -25;
						//printf("%f ",ab_j);
						z=spec_z[i];
						ab_j = ap_r[j] - ap_r[i] + ab_r[i];
					}
					else{ab_j=ab_r[j];}

					//*************************Control Sample*****************************

					/*Restrict bright galaxy to a absolute magnitude range.
					Faint galaxy can fall outside of this range*/
					if (brestrict!=1){bcheck=1;}
					else if(brestrict==1 && (ab_r[i] <= ab_j && ab_r[i] < (m2+blim)) ||
						(ab_j <= ab_r[i] && ab_j < (m2+blim))) {bcheck=1;}
					else{bcheck=0;}

					/*Restrict heavy galaxy to a mass range.
					Light galaxy can fall outside of this range*/
					if (hrestrict!=1){hcheck=1;}
					else if(hrestrict==1 && (mass[i] >= mass[j] && mass[i] > pow(10.0,mass2+hlim)) ||
						(mass[j] >= mass[i] && mass[j] > pow(10.0,mass2+hlim))) {hcheck=1;}
					else{hcheck=0;}

					//This applies the primary morphological filter:
					if(prestrict==0){pcheck=1;}
					else if(prestrict==1 && gpcheck[i]==1){pcheck=1;}						//Secondary
					else if(prestrict==2 && gpcheck[i]==1 && gpcheck[j]==1){pcheck=1;}		//Both
					else if(prestrict==3 && (gpcheck[i]==1 || gpcheck[j]==1)){pcheck=1;}	//Either
					else{pcheck=0;}
		

					/*Check if galaxies are within dM limits and mass ratio limits.
					A check is also done to see if only  considering bright/faint
					or heavy/light galaxy*/
					if((mcheck==1) && (bcheck==1) && (hcheck==1) && (pcheck==1) && (fb*ab_r[i] >= fb*ab_j) &&
						(((ab_j-ab_r[i]) <= dMmax && (ab_j-ab_r[i]) >= dMmin) || 
						((ab_r[i]-ab_j) <= dMmax && (ab_r[i]-ab_j) >= dMmin)) &&
						(hl*mass[i] >= hl*mass[j]) &&
						(((mass[i]/mass[j]>=ratiomin) && (mass[i]/mass[j]<=ratiomax)) ||
						((mass[j]/mass[i]>=ratiomin) && (mass[j]/mass[i]<=ratiomax))))
					{


						/*Calculate the projected distance between galaxies.*/
		
						/*Calculate angular separation*/
						//deltara = fabs(ra[i]-ra[j]);
						//deltadec = fabs(dec[i]-dec[j]);
						avgdec = (dec[j]+dec[i])/2;
						theta = sqrt(pow(deltara*cos(avgdec*pi/180),2)+pow(deltadec,2))*3600;
		
	
						/*Comoving distance*/
						DM = DH*qromb(e_func,0,z);
		
						/*DA*/
						DA = DM/(1+z);
		
						/*separation condition*/
						if(rpswitch==1)
						{
							//Relative separation of galaxies based angular separation divided by sum of radii
							rp= theta/(radius[i]+radius[j]);
						}
						else
						{
							//Physical projected spearation in kpc
							rp = (theta*2*pi*DA*1000)/(360*3600);
						}


						//Redshift Boundary Weights
						Zwt = 1.0;
						if(boundon==1)
						{
							if(z2*cs - cs*spec_z[i] < maxv)
							{
								if(spec_z[j]>spec_z[i]){rp = r2+1;}
								else{Zwt = 0.5;}
							}
							else if(cs*spec_z[i] - z1*cs < maxv)
							{
								if(spec_z[j]<spec_z[i]){rp = r2+1;}
								else{Zwt = 0.5;}
							}
							//printf("Zwt: %f\n",Zwt);
						}



						/*Check if pair is within search radius*/
						//if (rp < r2)
						if (rp < r2 && (theta > radx*(radius[i]+radius[j])))
						{
							//printf("%lld %lld %f %f\n",id[i],id[j],asy[i],asy[j]);

							/*Place galaxies in luminosity and separation bins:*/
	
							/*projected separation (rp) bin and array number*/
							l=floor(rp/rStep);
							rBin=(floor(rp/rStep))*rStep;
							if(l > rCount){rCount = l;}
		
							/*Absolute magnitude bin and array number*/
							if(LorM==1){k=floor((mlvar[i]-ml2)/mlStep);}
							else{k=floor((mlvar[i]-ml1)/mlStep);}

							mBin=(floor(mlvar[i]/mlStep))*mlStep+mlStep;
							if(k > mCount){mCount = k;}

							/*Velocity difference bin and array number*/
							//v=floor(deltav/vStep);
							//vBin=(floor(deltav/vStep))*vStep;
							//if(v > vCount){vCount = v;}
						
	
							rpN[k][l]++;						//Step up galaxy count
	
							magB[k][l] = mBin - mlStep/2.;		//magnitude bin
							rpB[k][l] = rBin + rStep/2.;		//radius bin


							//**************************************************************
							for(n=0;n<Nvar;n++)
							{
								pvalT[n][k][l] = pvalT[n][k][l] + pval[n][i];
								pvalV[n][k][l][rpN[k][l]] = pval[n][i];
								p=floor(pval[n][i]/pStep);
								pvalB[n][k][l][p]++;
							}
							//**************************************************************
	

							/*Check for co-rotating spirals in pairs*/
							if(cwacw==1)
							{
								if((pval[cw][i] > pval[acw][i] && pval[cw][j] > pval[acw][j]) ||
									(pval[cw][i] < pval[acw][i] && pval[cw][j] < pval[acw][j]))
									{coro=1.0;}
								else{coro=0.0;}
								//printf("%lld  %f  %f\n",id[i],pval[cw][i],pval[acw][i]);
								n=Nvar;

								pvalT[n][k][l] = pvalT[n][k][l] + coro;
								pvalV[n][k][l][rpN[k][l]] = coro;
								p=floor(coro/pStep);
								pvalB[n][k][l][p]++;
							}


							dvT[k][l] = dvT[k][l] + deltav;     //Total radial velocity difference
							dvV[k][l][rpN[k][l]] = deltav;     	//Individual radial velocity difference
	
							assT[k][l] = assT[k][l] + asy[i]; 	//Total GAMA CAS A measurement
							assV[k][l][rpN[k][l]] = asy[i];		//Individual GAMA CAS A measurement

							radT[k][l] = radT[k][l] + (radius[i]*2*pi*DA*1000)/(360*3600);	//Total Petro 90 radius measurements
							radV[k][l][rpN[k][l]] = (radius[i]*2*pi*DA*1000)/(360*3600);	//Individual Petro 90 radius measurements
	
							radT2[k][l] = radT2[k][l] + (radius[j]*2*pi*DA*1000)/(360*3600);	//Total Petro 90 radius of companion
							radV2[k][l][rpN[k][l]] = (radius[j]*2*pi*DA*1000)/(360*3600);	//Individual Petro 90 radius of companion
	
							masT[k][l] = masT[k][l] + mass[i];		//Total pairs masses
							masV[k][l][rpN[k][l]] = mass[i];		//Individual masses

							lumT[k][l] = lumT[k][l] + lum[i];		//Total pairs luminosities
							lumV[k][l][rpN[k][l]] = lum[i];		//Individual luminosities
	
							gid1[k][l][rpN[k][l]] = i;		//Record current galaxy id
							gid2[k][l][rpN[k][l]] = j;		//Record companion galaxy id
	
							psep[k][l][rpN[k][l]] = theta;
							rsep[k][l][rpN[k][l]] = rp;
	
	
	
							/*START****Find weights for luminosity function****************************/
							RAwt=1.0;
							DECwt=1.0;
							if(boundon==1)
							{
								RAwt=edgewt(ra[j],raMax,raMin,r2,r1,DA);
								DECwt=edgewt(dec[j],decMax,decMin,r2,r1,DA);
							}


							/*Find the Vmax weight*/
							if(vmaxon==1)
							{
								/*Use the Vmax weight of the faint galaxy for both*/
								if(ab_r[i]>ab_j){Mi=ab_r[i];}
								else{Mi=ab_j;}					
	
								/*Vwt calculation*/				
								Dobj = rtbis(Mlim_func,z1,z2,0.0000001);			
								Dmax = DH*qromb(e_func,0,Dobj);
								Vmax = pow(maxDist,3)- pow(minDist,3);
								Vobj = pow(Dmax,3)- pow(minDist,3);
								Vwt = Vobj/Vmax;
	
								/*if the galaxy is observable throughout the sample its weight is 1*/
								if(Vwt > 1.) Vwt=1.;
							}
							else{Vwt=1.;}

							/*Find the angular spectroscopic incompleteness weight bin*/
							n = floor(theta/tStep);
	
							/*Multiply by angular incompletness weight*/
							Vwt = Vwt*tcorr[3][n]*Zwt*RAwt*DECwt;

							if(removeint==1){Vwt = Vwt/(a0*pow(e,-rp/b0)/(a0+c0));}//Kevin's method
							if(removeint==2){Vwt = Vwt/(a0*pow(e,-rp/b0));}//Steven's method
							/*END*****Find weights for luminosity function****************************/
	
	
							/*Add the angular spectroscopic incompleteness and Vmax weights*/
							wtT[k][l] = wtT[k][l]+ 1./Vwt;		//Total combined weights
							wtV[k][l][rpN[k][l]] = 1./Vwt;		//Indiviual combined weights						
		
							rawPairs++; //Count up total raw pairs count
							wtPairs = wtPairs + 1./Vwt; //Count up weighted pairs count
		
							
							/*For galaxies with velocity dispersion measurements*/
							if(velDisp[i]>0)
							{
								vdisT[k][l] = vdisT[k][l] + velDisp[i]; //Total galaxy velocity dispersion
								vdisV[k][l][rpN[k][l]] = velDisp[i]; 	//Indiviual galaxy velocity dispersion
								vdisN[k][l]++;
							}

						}
					}
                }
			}
        }
	}

	printf("\nmaCount=%d  mCount=%d  rCount=%d\n",maCount,mCount,rCount);	

	fprintf(rpresults,"rawTal=%d  rawPairs=%d  rawNc=%f\n",rawTal,rawPairs,rawPairs*1.0/rawTal);
	fprintf(rpresults,"wtTal=%f  wtPairs=%f  wtNc=%f\n\n",wtTal,wtPairs,wtPairs/wtTal);	

	/**********************************************************************************************/
	/*Calculate average pmerge value for each bin.*/
	//fprintf(rpresults,"#M_r        r_p        N   wt         wt_error  pring     prg_err   parc      par_err   pdist     pdt_err   pirr      pir_err   pother    pot_err   pmerge    pm_error   dV          dV_error   asymm     a_error    pairsRad  r_error  compRad    r2_error  pairsMas  m_error   globalR   gbr_err   gbN   gbWt         gbWt_err  globalMas gbm_err\n");

	//fprintf(pairUrls,"#M_r        r_p       N       pmerge     asymm   id1     r90_1     id2    r90_2     rp\n");



	/*Add another variable if considering co-rotating galaxies*/
	if(cwacw==1){Nvar2=Nvar+1;}
	else{Nvar2=Nvar;}

	for(k=0;k<=mCount;k++)
	{			
		for(n=0;n<Nvar2;n++)
		{
			g_error[n]=serror(gvalT[n][k],gbN[k],gvalV[n][k]);
		}

		gbr_error=serror(gbradT[k],gbN[k],gbradV[k]);		//Global ave radius
		gbwt_error=jerror(gbwtT[k],gbN[k],gbwtV[k]);		//Global LF weighted
		gbm_error=serror(gbmasT[k],gbN[k],gbmasV[k]);		//Global ave mass
		gbm_error=serror(gblumT[k],gbN[k],gblumV[k]);		//Global ave luminosities
		asy_error=jerror(asyT[k],asyN[k],asyV[k]);			//Asymmetry fraction error


		//printf("%d\n",k);
		for(l=0;l<=rCount;l++)
		{

			
			for(n=0;n<Nvar2;n++)
			{
				p_error[n]=serror(pvalT[n][k][l],rpN[k][l],pvalV[n][k][l]);
			}

			dV_error=serror(dvT[k][l],rpN[k][l],vdisV[k][l]);	//radial velocity difference
			wt_error=jerror(wtT[k][l],rpN[k][l],wtV[k][l]);		//Pairs LF weighted
			a_error=serror(assT[k][l],rpN[k][l],assV[k][l]);	//GAMA CAS A measurement
			r_error=serror(radT[k][l],rpN[k][l],radV[k][l]);	//Primary ave radius
			r2_error=serror(radT2[k][l],rpN[k][l],radV2[k][l]);	//Companion ave radius
			m_error=serror(masT[k][l],rpN[k][l],masV[k][l]);	//Pairs ave mass
			l_error=serror(lumT[k][l],rpN[k][l],lumV[k][l]);	//Pairs ave luminosities



			//printf("s:%f  j:%f  %d\n",serror(asyT[k],asyN[k],asyV[k]),jerror(asyT[k],asyN[k],asyV[k]),asyN[k]);
			//printf("s:%f  j:%f  %d\n",serror(wtT[k][l],rpN[k][l],wtV[k][l]),jerror(wtT[k][l],rpN[k][l],wtV[k][l]),rpN[k][l]);

			nc_error=sqrt(pow(gbwt_error/gbwtT[k],2)+pow(wt_error/wtT[k][l],2))*(wtT[k][l]/gbwtT[k]); //Nc error

			/******************************rpresults**********************/
			fprintf(rpresults,"%f  %f  %f  %d  %f  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f",magB[k][l],gbradT[k]/gbN[k],gbr_error,gbN[k],gbwtT[k],gbwt_error,log10(gbmasT[k]/gbN[k]),log10(gbm_error),rpB[k][l],rpN[k][l],wtT[k][l],wt_error,dvT[k][l]/rpN[k][l],dV_error,assT[k][l]/rpN[k][l],a_error,radT[k][l]/rpN[k][l],r_error,radT2[k][l]/rpN[k][l],r2_error,log10(masT[k][l]/rpN[k][l]),log10(m_error));

			for(n=0;n<Nvar2;n++)
			{
				fprintf(rpresults,"  %f  %f",pvalT[n][k][l]/rpN[k][l],p_error[n]);
			}
			fprintf(rpresults,"\n");

			//fprintf(rpresults,"%f  %f\n",pvalT[baryes][k][l]/rpN[k][l],p_error[baryes]);

			/******************************rpresults**********************/


			/******************************ncout**********************/
			fprintf(ncout,"%f  %f  %f  %f\n",magB[k][l],rpB[k][l],(wtT[k][l]/gbwtT[k]),nc_error);
			/******************************ncout**********************/


			/******************************pairUrls**********************/	
			fprintf(pairUrls,"<br><br>Separation bin %f starts:<br><br>",rpB[k][l]);
			for(i=1;i<=rpN[k][l];i++)
			{
				

				//fprintf(pairlist,"%f  %lld  %lld  %f  %f  %f\n",rpB[k][l],id[gid1[k][l][i]],id[gid2[k][l][i]],psep[k][l][i],rsep[k][l][i],dvV[k][l][i]);
				fprintf(pairlist,"%lld  %f  %f\n",id[gid1[k][l][i]],ra[gid1[k][l][i]],dec[gid1[k][l][i]]);
				//printf("%f  %lld  %lld  %f  %f  %f\n",rpB[k][l],id[gid1[k][l][i]],id[gid2[k][l][i]],psep[k][l][i],rsep[k][l][i],dvV[k][l][i]);




				if(htmlSwitch==1){
					sprintf(htmlCode1,"<table border=\"0\"><tr><td><img src=\"");//http://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/GDS_61.jpg

					sprintf(htmlCode2,"\"></td><td><img src=\"");//http://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/GDS_61.jpg



					if(cason==1){
					sprintf(htmlCode3,"\"></td></tr><tr><td>Rp=%f, Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f<br>C=%f, A=%f, S=%f, M20=%f, Gini=%f</td>",rsep[k][l][i],gbmagB[k],id[gid1[k][l][i]],log10(mass[gid1[k][l][i]]),ap_r[gid1[k][l][i]],spec_z[gid1[k][l][i]],pval[cpos][gid1[k][l][i]],pval[apos][gid1[k][l][i]],pval[spos][gid1[k][l][i]],pval[mpos][gid1[k][l][i]],pval[gpos][gid1[k][l][i]]);


					sprintf(htmlCode4,"<td>Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f<br>C=%f, A=%f, S=%f, M20=%f, Gini=%f</td></tr></table>",gbmagB[k],id[gid2[k][l][i]],log10(mass[gid2[k][l][i]]),ap_r[gid2[k][l][i]],spec_z[gid2[k][l][i]],
						pval[cpos][gid2[k][l][i]],pval[apos][gid2[k][l][i]],pval[spos][gid2[k][l][i]],pval[mpos][gid2[k][l][i]],pval[gpos][gid2[k][l][i]]);
					}
					else{
					sprintf(htmlCode3,"\"></td></tr><tr><td>Rp=%f, Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f<br>Ftrs=%f, Clmpy=%f, FcOn=%f, Bar=%f, Sprl=%f, AWidx=%f, BLGidx=%f<br>1arm=%f, 2arm=%f, 3arm=%f, 4arm=%f, MRarm=%f, UNarm=%f</td>",rsep[k][l][i],gbmagB[k],id[gid1[k][l][i]],log10(mass[gid1[k][l][i]]),ap_r[gid1[k][l][i]],spec_z[gid1[k][l][i]],pval[feturs][gid1[k][l][i]],pval[clmp][gid1[k][l][i]],pval[edgeno][gid1[k][l][i]],pval[baryes][gid1[k][l][i]],pval[sparms][gid1[k][l][i]],pval[awidx][gid1[k][l][i]],pval[blgidx][gid1[k][l][i]],pval[oarm][gid1[k][l][i]],pval[tarms][gid1[k][l][i]],pval[threearms][gid1[k][l][i]],pval[fourarms][gid1[k][l][i]],pval[morearms][gid1[k][l][i]],pval[unknarms][gid1[k][l][i]]);


					sprintf(htmlCode4,"<td>Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f<br>Ftrs=%f, Clmpy=%f, FcOn=%f, Bar=%f, Sprl=%f, AWidx=%f, BLGidx=%f<br>1arm=%f, 2arm=%f, 3arm=%f, 4arm=%f, MRarm=%f, UNarm=%f</td></tr></table>",gbmagB[k],id[gid2[k][l][i]],log10(mass[gid2[k][l][i]]),ap_r[gid2[k][l][i]],spec_z[gid2[k][l][i]],pval[feturs][gid2[k][l][i]],pval[clmp][gid2[k][l][i]],pval[edgeno][gid2[k][l][i]],pval[baryes][gid2[k][l][i]],pval[sparms][gid2[k][l][i]],pval[awidx][gid2[k][l][i]],pval[blgidx][gid2[k][l][i]],pval[oarm][gid2[k][l][i]],pval[tarms][gid2[k][l][i]],pval[threearms][gid2[k][l][i]],pval[fourarms][gid2[k][l][i]],pval[morearms][gid2[k][l][i]],pval[unknarms][gid2[k][l][i]]);
					}

				}
				else{

					if(cason==1)
					{
						fprintf(pairUrls,"%f  %f_%f_%f_%f_%lld_%lld_%.2f-%.2f  %f  %f  %f  %f  %f  %f  %f  %d",magB[k][l],log10(mass[gid1[k][l][i]]),log10(mass[gid2[k][l][i]]),(mass[gid1[k][l][i]]/mass[gid2[k][l][i]]),rpB[k][l],id[gid1[k][l][i]],id[gid2[k][l][i]],asy[gid1[k][l][i]],asy[gid2[k][l][i]],radius[gid1[k][l][i]],radius[gid2[k][l][i]],psep[k][l][i],rsep[k][l][i],dvV[k][l][i],(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,rpN[k][l]);

					}
					else
					{
						fprintf(pairUrls,"%f  %f_%f_%f_%lld_%lld_%.2f-%.2f-%.2f  %f  %f  %f  %f  %f  %f  %f  %d",magB[k][l],log10(mass[gid1[k][l][i]]),log10(mass[gid2[k][l][i]]),rpB[k][l],id[gid1[k][l][i]],id[gid2[k][l][i]],(pval[feturs][gid1[k][l][i]]*pval[edgeno][gid1[k][l][i]]*pval[sparms][gid1[k][l][i]]),pval[lwarms][gid1[k][l][i]],pval[pmrg][gid1[k][l][i]],radius[gid1[k][l][i]],radius[gid2[k][l][i]],psep[k][l][i],rsep[k][l][i],dvV[k][l][i],(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,rpN[k][l]);
					}

					for(n=0;n<Nvar2;n++)
					{
						fprintf(pairUrls,"  %f",pvalV[n][k][l][i]);
					}
					fprintf(pairUrls,"\n");

				}


				
				if (wget==1){
					/*if(rpswitch==1)
					{
						DA = DH*qromb(e_func,0,spec_z[asyID[k][i]])/(1+spec_z[asyID[k][i]]);
						fprintf(pairUrls,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,40/DA);
					}
					else
					{*/
						//fprintf(pairUrls,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,(psep[k][l][i]/rsep[k][l][i])*0.34);


						DA = DH*qromb(e_func,0,spec_z[gid1[k][l][i]])/(1+spec_z[gid1[k][l][i]]);
						fprintf(pairUrls,"%shttp://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",htmlCode1,ra[gid1[k][l][i]],dec[gid1[k][l][i]],40/DA);

						DA = DH*qromb(e_func,0,spec_z[gid2[k][l][i]])/(1+spec_z[gid2[k][l][i]]);
						fprintf(pairUrls,"%shttp://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424%s%s\n",htmlCode2,ra[gid2[k][l][i]],dec[gid2[k][l][i]],40/DA,htmlCode3,htmlCode4);

					//}
				}

				//GZ4 Candels
				if (wget==2){

					fprintf(pairUrls,"%shttp://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/%s_%lld.jpg\n",
						htmlCode1,catCode,id[gid1[k][l][i]]);
					fprintf(pairUrls,"%shttp://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/%s_%lld.jpg%s%s\n",htmlCode2,catCode,id[gid2[k][l][i]],htmlCode3,htmlCode4);

				}	
			
			}
			/******************************pairUrls**********************/




			/******************************histout**********************/
			for(i=0;i<=10;i++)
			{
				fprintf(histout,"%f  %f  %f  %d",magB[k][l],rpB[k][l],(pStep*i+pStep/2.),rpN[k][l]);


				for(n=0;n<Nvar2;n++)
				{
					fprintf(histout,"  %f",pvalB[n][k][l][i]/rpN[k][l]);
				}
				fprintf(histout,"\n");

			}
			fprintf(histout,"\n");
			/******************************histout**********************/
		}
		fprintf(rpresults,"\n");

		/******************************casout**********************/

		fprintf(casout,"%f  %d  %f  %f  %d  %f  %f  %f  %f\n",magB[k][l-1],gbN[k],gbwtT[k],gbwt_error,asyN[k],asyT[k],asy_error,asyT[k]/gbwtT[k],sqrt(pow(asy_error/asyT[k],2)+pow(gbwt_error/gbwtT[k],2))*(asyT[k]/gbwtT[k]));

		//fprintf(casout,"%f  %d  %f  %f  %d  %f  %f  %f  %f\n",magB[k][l-1],gbN[k],gbwtT[k],gbwt_error,asyN[k],asyT[k],asy_error,asyT[k]/gbwtT[k],berror((asyT[k]/gbwtT[k]), gbN[k]) );

		/*for(i=1;i<=asyN[k];i++)
		{
			fprintf(casout,"\n",magB[k][1],
		*/


		/******************************casout**********************/




		/******************************allUrls, asout and gm20out**********************/

		a=k;
		if(mCount!=maCount){
			a=0;
			while(gbmagB[k]!=gbmagBTotal[a] || a==maCount){
				printf("%f  %f  %d  %d\n",gbmagB[k],gbmagBTotal[a],maCount,a);
				a++;
			}
		}


		if(wget==-1){
			fprintf(allUrls,"#!/bin/bash\n\n");
		}
		if(wget>0){
			fprintf(allUrls,"Number of Morph galaxies:%f:%d  %f:%d  %f</br>",gbmagB[k],gbN[k],gbmagBTotal[a],gbNTotal[a],(gbN[k]*1.0)/(gbNTotal[a]*1.0));

			fprintf(allUrls,"Mean Mass: %f +- %f</br>",log10(gbmasT[k]/gbN[k]),log10(gbm_error));

			fprintf(allUrls,"Mean Bar vote frac: %f +- %f</br>",gvalT[baryes][k]/gbN[k],g_error[baryes]);
			fprintf(allUrls,"Mean Arm Tightness Index: %f +- %f</br>",gvalT[awidx][k]/gbN[k],g_error[awidx]);

			fprintf(allUrls,"Mean Arm Number: %f %f   %f %f   %f %f   %f %f   %f %f  %f</br>",gvalT[oarm][k]/gbN[k],g_error[oarm],gvalT[tarms][k]/gbN[k],g_error[tarms],gvalT[threearms][k]/gbN[k],g_error[threearms],gvalT[fourarms][k]/gbN[k],g_error[fourarms],gvalT[unknarms][k]/gbN[k],g_error[unknarms],(gvalT[oarm][k]/gbN[k] + gvalT[tarms][k]/gbN[k] + gvalT[threearms][k]/gbN[k] + gvalT[fourarms][k]/gbN[k] + gvalT[unknarms][k]/gbN[k]));

			fprintf(allUrls,"Mean Bulge Dominance Index: %f +- %f</br>",gvalT[blgidx][k]/gbN[k],g_error[blgidx]);
			fprintf(allUrls,"Mean Star Formation Rate: %f +- %f</br>",gvalT[sfr][k]/gbN[k],g_error[sfr]);
			fprintf(allUrls,"Mean Colour: %f +- %f</br>",colour/gbN[k],(colour/gbN[k])*sqrt(gbN[k])/gbN[k]);

		}


		for(i=1;i<=gbN[k];i++)
		{
		

			fprintf(alllist,"%f  %lld  %f  %f\n",gbmagB[k],id[gbID[k][i]],ra[gbID[k][i]],dec[gbID[k][i]]);

		 	//fprintf(pairUrls,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,(psep[k][l][i]/rsep[k][l][i])*0.34);

			if(wget==-1){
				fprintf(allUrls,"wget -O \"%lld\" \"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\"\n",id[gbID[k][i]],ra[gbID[k][i]],dec[gbID[k][i]],40/DA);
			}
			else if(wget==0){
				fprintf(allUrls,"%f %lld %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",gbmagB[k],id[gbID[k][i]],log10(mass[gbID[k][i]]),ap_r[gbID[k][i]],ab_r[gbID[k][i]],spec_z[gbID[k][i]],pval[feturs][gbID[k][i]],pval[edgeno][gbID[k][i]],pval[baryes][gbID[k][i]],pval[sparms][gbID[k][i]],pval[lwarms][gbID[k][i]],pval[odd][gbID[k][i]],pval[dstrb][gbID[k][i]],pval[irrgl][gbID[k][i]],pval[pmrg][gbID[k][i]],pval[oarm][gbID[k][i]],pval[tarms][gbID[k][i]]);
			}
			else if(wget>0){

				if(htmlSwitch==1){
					sprintf(htmlCode2,"<table border=\"0\"><tr><td><img src=\"");//http://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/GDS_61.jpg
					//sprintf(htmlCode3,"\"></td></tr><tr><td>Bin=%f, ID=%lld, Mass=%f, z=%f</td></tr></table>",gbmagB[k],id[gbID[k][i]],log10(mass[gbID[k][i]]),spec_z[gbID[k][i]]);

			
					if(cason==1){
						sprintf(htmlCode3,"\"></td></tr><tr><td>Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f</td></tr></td></tr><tr><td>C=%f, A=%f, S=%f, M20=%f, Gini=%f</td></tr></table>",gbmagB[k],id[gbID[k][i]],log10(mass[gbID[k][i]]),ap_r[gbID[k][i]],spec_z[gbID[k][i]],pval[cpos][gbID[k][i]],pval[apos][gbID[k][i]],pval[spos][gbID[k][i]],pval[mpos][gbID[k][i]],pval[gpos][gbID[k][i]]);

					}
					else{

						sprintf(htmlCode3,"\"></td></tr><tr><td>Bin=%f, ID=%lld, Mass=%f, ApMag=%f, z=%f</td></tr></td></tr><tr><td>Ftrs=%f, Clmpy=%f, FcOn=%f, Bar=%f, Sprl=%f, AWidx=%f, BLGidx=%f</td></tr></td></tr><tr><td>1arm=%f, 2arm=%f, 3arm=%f, 4arm=%f, MRarm=%f, UNarm=%f</td></tr></td></tr><tr><td>irr=%f, dstb=%f, merg=%f, odd=%f</tr></table>",gbmagB[k],id[gbID[k][i]],log10(mass[gbID[k][i]]),ap_r[gbID[k][i]],spec_z[gbID[k][i]],pval[feturs][gbID[k][i]],pval[clmp][gbID[k][i]],pval[edgeno][gbID[k][i]],pval[baryes][gbID[k][i]],pval[sparms][gbID[k][i]],pval[awidx][gbID[k][i]],pval[blgidx][gbID[k][i]],pval[oarm][gbID[k][i]],pval[tarms][gbID[k][i]],pval[threearms][gbID[k][i]],pval[fourarms][gbID[k][i]],pval[morearms][gbID[k][i]],pval[unknarms][gbID[k][i]],pval[irrgl][gbID[k][i]],pval[dstrb][gbID[k][i]],pval[pmrg][gbID[k][i]],pval[odd][gbID[k][i]]);
					}


				}
				else{
				fprintf(allUrls,"%f_%lld_%f_%f_%f_%f %f %f %f %f %f %f %f %f %f %f %f\n",gbmagB[k],id[gbID[k][i]],log10(mass[gbID[k][i]]),ap_r[gbID[k][i]],ab_r[gbID[k][i]],spec_z[gbID[k][i]],pval[feturs][i],pval[edgeno][i],pval[baryes][i],pval[sparms][i],pval[lwarms][i],pval[odd][i],pval[dstrb][i],pval[irrgl][i],pval[pmrg][i],pval[oarm][i],pval[tarms][i]);
				}


				//SDSS CASJOBS:
				if (wget==1){
					DA = DH*qromb(e_func,0,spec_z[gbID[k][i]])/(1+spec_z[gbID[k][i]]);
					fprintf(allUrls,"%shttp://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424%s\n",htmlCode2,ra[gbID[k][i]],dec[gbID[k][i]],40/DA,htmlCode3);
				}

				//GZ4 Candels
				if (wget==2){
					fprintf(allUrls,"%shttp://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/%s_%lld.jpg%s\n",htmlCode2,catCode,id[gbID[k][i]],htmlCode3);
				}	

			}

		}
		fprintf(alllist,"\n");
		fprintf(allUrls,"\n");


		if(wget==1){fprintf(asout,"%f\n",gbmagB[k]);}

		for(i=1;i<=asyN[k];i++)
		{
		

		 //fprintf(pairUrls,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",(ra[gid1[k][l][i]]+ra[gid2[k][l][i]])/2.,(dec[gid1[k][l][i]]+dec[gid2[k][l][i]])/2.,(psep[k][l][i]/rsep[k][l][i])*0.34);
			if(wget==0){
				fprintf(asout,"%f_%lld_%f_%f_%f_%f_%f_%f\n",gbmagB[k],id[asyID[k][i]],log10(mass[asyID[k][i]]),ap_r[asyID[k][i]],ab_r[asyID[k][i]],spec_z[asyID[k][i]],asy[asyID[k][i]],smo[asyID[k][i]]);
			}

			DA = DH*qromb(e_func,0,spec_z[asyID[k][i]])/(1+spec_z[asyID[k][i]]);


			fprintf(asout,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",ra[asyID[k][i]],dec[asyID[k][i]],40/DA);
		
			//http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=221.998978&dec=1.981267&scale=0.001499&width=424&height=424
		}
		fprintf(asout,"\n");

		for(i=1;i<=gm20N[k];i++)
		{
		
			fprintf(gm20out,"%f_%lld_%f_%f_%f_%f_%f\n",gbmagB[k],id[gm20ID[k][i]],log10(mass[gm20ID[k][i]]),ap_r[gm20ID[k][i]],ab_r[gm20ID[k][i]],spec_z[gm20ID[k][i]],asy[gm20ID[k][i]]);


			DA = DH*qromb(e_func,0,spec_z[asyID[k][i]])/(1+spec_z[asyID[k][i]]);


			fprintf(gm20out,"http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=%f&dec=%f&scale=%f&width=424&height=424\n",ra[gm20ID[k][i]],dec[gm20ID[k][i]],40/DA);
		
			//http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra=221.998978&dec=1.981267&scale=0.001499&width=424&height=424
		}
		fprintf(gm20out,"\n");


	}





	fclose(rpresults);
	fclose(ncout);
	fclose(pairUrls);
	fclose(pairlist);
	fclose(alllist);
	fclose(histout);
	fclose(casout);
	fclose(allUrls);
	fclose(asout);
	fclose(gm20out);


	printf("All done!\n");

	return 0;
}



/*Random number generator*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



/*trapzd, polint and qromb are numerical recipies for doing integrations*/

#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define NRANSI


void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}
#undef NRANSI

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

float qromb(float (*func)(float), float a, float b)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	float ss,dss;
	float s[JMAXP],h[JMAXP+1];
	int j;

//	printf("qromb\n");

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K


float e_func(float z)
{

	float Ez;

	Ez = sqrt(omegaM*pow((1+z),3) + omegaK*pow((1+z),2) + omegaL);

	return 1./Ez;
}


#define JMAX 40

float rtbis(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(x1);
	fmid=(*func)(x2);
//	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
	        if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
#undef JMAX


float Mlim_func(float z)
{

	float mlim;
	
	//mlim = apMagMax-5*log10(DH*qromb(e_func,0,z)*(1+z)) - 25 - kmax(z) + 0.75*2.5*log10(1+z) - Mi;
	mlim = apMagMax-5*log10(DH*qromb(e_func,0,z)*(1+z)) - 25 - Mi;

	return mlim;
}




/*Jacknife errors*/
float jerror(float total,int count,float values[])
{
	int i;
	float deltaSum,error,jd;

	deltaSum =0;
	for(i=0;i<count;i++)
	{
		deltaSum = deltaSum + pow(total -(total - values[i]),2);
	}
	error= sqrt((count- 1)*deltaSum/count);

	//jd= sqrt((count- 1)*deltaSum/count);
	//error= jd/sqrt(count);

	return error;
}

/*Jacknife errors*/
/*float jerror(float total,int count,float values[])
{
	int i;
	float deltaSum,error,mean;

	mean = total/(count);
	
	for(i=0;i<count;i++)
	{
		deltaSum = deltaSum + pow((((total -values[i])/(count-1)-mean)),2);
	}

	error = sqrt((count-1)*deltaSum/count);

	return error;
}
*/



/*Poisson errors*/
float xerror(float total,int count,float values[])
{

	float error;

	error= total/sqrt(count-1);

	return error;
}


/*Poisson errors*/
float berror(float value,int count)
{

	// http://www.sigmazone.com/binomial_confidence_interval.htm
	float error;

	error= 1.0*sqrt((value*(1.0-value))/(count*1.0));

	return error;
}


/*Standard deviation of mean error*/
float serror(float total,int count,float values[])
{
	int i;
	float mean,total2,error,sd;

	mean = total/count;
	total2=0;
	for(i=0;i<count;i++)
	{
		total2 = total2 + pow(values[i]-mean,2);
	}
	sd = sqrt(total2/count);
	error= sd/sqrt(count);


	error = sqrt(total2/(count*(count-1)));
	
	return error;
}

//Geometric Boundary Weights

float edgewt(float coord,float coMax,float coMin,float maxRad,float minRad,float DA)
{
		float deltaCoord,edgeDist,area1,area2,area3,wt;



		if(coord > coMax - (coMax - coMin)/2)
		{
			deltaCoord = fabs(coord-coMax);
		}
		else
		{
			deltaCoord = fabs(coord-coMin);
		}
		
		edgeDist = deltaCoord*2*pi*DA*1000/360;
	
		if(edgeDist > maxRad){edgeDist = maxRad;}
		
		area1 = (pow(maxRad,2)/2)*(2*pi - acos(2*pow(edgeDist,2)/pow(maxRad,2) - 1)) + edgeDist*pow(pow(maxRad,2) - pow(edgeDist,2),0.5);
		
		if(edgeDist >= minRad)
		{
			area2 = pi*pow(minRad,2);
		}
		else
		{
			area2 = (pow(minRad,2)/2)*(2*pi - acos(2*pow(edgeDist,2)/pow(minRad,2) - 1)) + edgeDist*pow(pow(minRad,2) - pow(edgeDist,2),0.5);
		}

		area3 = pi*pow(maxRad,2) - pi*pow(minRad,2);
		
		wt=(area1-area2)/area3;	

		//printf("%f\n",wt);

		return wt;
}						

		
/*
//Redshift Boundary Weights
float redwt(float zi1,float zi2,float zmin,float zmax,float vmin,float vmax)
{
		float  Zwt = 1.0;
		if(zmax*cs - cs*zi1 < vmax)
		{
			if(zi2 > zi1){deltaRadius = maxRadius+1;}
			else{Zwt = 0.5;}
		}
		else if(cs*spec_z[i1] - minZ*cs < maxv)
		{
			if(zi2 < zi1){deltaRadius = maxRadius+1;}
			else{Zwt = 0.5;}
		}

		return Zwt;
}
*/

/*This is the mass dependent assymetry limit (A).
Masses below 10^10 are with a mx+b fit.*/
float alimit(float mass)
{
	float a= -0.0425091;
	float b= 0.424672;
	float Adefault = 0.35;
	float Alim;

	mass=log10(mass);
	if(Amass==1 && mass<9.99){Alim = a*mass + b + Adefault;}
	else{Alim=Adefault;}

	return Alim;
}





