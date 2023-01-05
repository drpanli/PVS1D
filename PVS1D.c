/*
A 1D Tissue Model of the Cardiac Purkinje-Ventricular System
by Mengya Yuan, Heqiang Lian, and Pan Li*, 
as in Spatiotemporal Patterns of Early Afterdepolarizations Underlying Abnormal T-wave Morphologies in A Tissue Model of the Purkinje-Ventricular System, PLoS One, Jan 2023
*Email: pan.li@nih.gov

To plot use matlab:
 [t,S,X,T,Si] = readLOOPmv('*.txt',201,1000,16);

*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define	 beats 5 
#define	 BCLMax   	1000 	
#define	 BCLMin   	200 	 
#define	 dBCL   	50 		

#define N 			50  /* number of P- cells */
#define ENDO		50	/* number of endo- cells */
#define MID			50	/* number of M- cells */
#define EPI			50	/* number of epi- cells */
#define M			ENDO+MID+EPI	/* total number of V- cells */

double ICKr, ICCaL, ICNa; // as scaling factor between 0 and 1 to GKr, GCaL, GNa; 
char DrugName[80]; 

const double Dp = 1.2; 		//gives 2m/s  
const double ER = 3.13;   	//gives 4.364ms P-V Delay
const double Dv = 0.1; 		//gives 0.5m/s 
int count1 = 0;
int run = 0;
double BCL,S2;

// GEOMETRY
double ageo,acap,vcell,vmyo,vnsr,vjsr,vsss,vcsr,vssl,vmito;

double radius;
double length;
double rcg;
const double pi = 3.14;
// length of the cell (cm)

double dx = 0.1;

// PHYSICAL CONSTANTS
const double frdy = 96485;
const double R = 8314;
const double temp = 310;

const double nao = 140;
const double cao = 1.8;
const double ko  = 5.4;
const double clo = 100;

const double zna = 1;
const double zk  = 1;
const double zcl = -1;
const double zca = 2;
const double ganai = 0.75;
const double ganao = 0.75;
const double gaki  = 0.75;
const double gako  = 0.75;
const double gacai = 1.0;
const double gacao = 0.341;

// REVERSAL POTENTIALS
double ena,ek,eca;

// TIME

double t,dt,tmax,T;
const double dtmin = 0.004;
const double dtmed = 0.004;
const double dtmax = 0.004;

// VOLTAGE
double dvdtclock;
const double dvdtthresh = 1;

double sum;
double R1, R2;
double pECG;

// STIMULUS CURRENT
const double stimdur = 2.5;
double istim = -300; 	
double tstim,stimtime;
int stimcount;


// Declare variables:

double itot, ma, mb, ha, hb, ja, jb, mtau, htau, jtau, mss, hss, ina, gna, mltau, ml3tau, mlss, jltau, jl3tau, jlss, inal2, gnal2, inal3, gnal3, inal;
double ibarca, ml3ss, hltau, hl3tau, hlss, hl3ss, jl3ss, jss, dss, dtau, fss, ftau, f2ss, f2tau, fca2ss;
double bss, gss, taub, taug, icat;
double atau, itau, i2tau, ass, i2ss, itos, itof, ito1;
double gkr, xrss, xrtau, rkr, ikr;
double eks, gks, xsss, xs1tau, xs2tau, iks;
double k1ss, gk1, ik1;
double yss, ytau, ifna, ifk, iftotal;
double allo, num, denommult, denomterm1, denomterm2, deltaE, inacass;
double fcass, fcatau, iss, fca2tau;
double fca_dtaucamkbar, camk0, alphacamk, betacamk, kmcam, kmcamk, trpnbar1, cmdnbar1, csqnbar1, trpnbar, kmtrpn, cmdnbar, kmcmdn, csqnbar, kmcsqn, bsrbar, kmbsr, bslbar, kmbsl, kmup, nsrbar, dkmplbbar,dqupcamkbar, IP3, k1, k1a, k0, k0a, k2, k2a, tauip3r, tautr1, tautr2, gaptau, sstau, ipcabar, kmpca, ibarnak, inacamax, kmcaact, kmnai1, kmnao, kmcai, kmcao, nu, ksat, gnab, pcab, pnab, prnak, gto, gtof, gcat, powtau, pca, gnal2, gnal3, gna, gtos, gtof;
double inaca, inak, ipca, icab, inab, inab, icatot, iktot, inatot, du, POip3, qip3, qdiff, REL, ireltau, irelss, qgap, dkmplb, dqupcamk, qup1, dkmplb, qtr1, qtr2, qdiffna, qgapna, dcasss, bsss, dcassl, trpn, cmdn, catotal, bmyo, cmyo, dmyo, dcajsr, qtr1, qtr2, csqn, bcsr, ccsr, dcansr, dnass, dnassl, dnai, dki, dcai, caavg, camkbound;
double ecl, alphaml, betaml, csqn1, bjsr, cjsr, dcacsr, dnasss, gnal, acttau, actss, icala, icalb, ivss, ivtau, icalx, icaly, istarca, istarvtau;
double gbarto1, gkp, kmto2, pcl, i2ftau, gbarclb, gnakbar, hill, kmnai2, kmko, ctkclbar, taudiffss, taudiff, ctnaclbar, qleakbar, qupbar, tautr, istarvss, icalxstar, icalystar, icaldelta, icaltheta, icalthetai, icaldeltai, icalca, icaloa, icalcia, icaloia, icalcstara, icalostara, icalcistara, icalcistara, icaloistara, dicalc, dicalo, dicalcstar, dicalostar, dicalci, dicaloi, dicalcistar, dicaloistar;
double alphaa, betaa, alphai, betai, alphai2, betai2, rto1;
double sa, sb, sr, sd, ikseta, ikstheta, iksomega, ikspsi, iksc1a, iksc2a, iksc3a, iksc4a, iksc5a, iksc6a, iksc7a, iksc8a, iksc9a, iksc10a, iksc11a, iksc12a, iksc13a, iksc14a, iksc15a, ikso1a, ikso2a, diksc1, diksc2, diksc3, diksc4, diksc5, diksc6, diksc7, diksc8, diksc9, diksc10, diksc11, diksc12, diksc13, diksc14, diksc15, dikso1, dikso2, oks;
double k1b, kp, ikp, ai2f, bi2f, i2fss, Kcaito2, ibarto2, ito2, iclb, nacax1, nacax2, nacax3, nacay1, nacay2, nacay3, inacasssr, inacai, ctkcl, ctnacl;
double ikp, icltot, ito2, iclb, relbtau, rela, qup, qleak, qtr, qdiffss, bsssr, bsscal, inacasssr, vsssr, vsscal, vsss, vss, cajsrb, cajsrc, dnasssr, qdiffcl, dcli, dclss;
double dvdt[N+M];
double v[N+M];
double m[N+M];
double h[N+M];
double j[N+M];
double d[N+M];
double f[N+M];
double f2[N+M];
double fca[N+M];
double fca2[N+M];
double xs1[N+M];
double xs2[N+M];
double xr[N+M];
double a[N+M];
double i[N+M];
double i2[N+M];
double ml[N+M];
double ml3[N+M];
double hl[N+M];
double hl3[N+M];
double jl[N+M];
double jl3[N+M];
double casss[N+M];
double cajsr[N+M];
double cacsr[N+M];
double cansr[N+M];
double cassl[N+M];
double nai[N+M];
double nassl[N+M];
double nasss[N+M];
double ki[N+M];
double cai[N+M];
double b[N+M];
double g[N+M];
double u[N+M];
double y[N+M];
double camktrap[N+M];
double cai[N+M];
double casssr[N+M];
double casscal[N+M];
double nasssr[N+M];
double ki[N+M];
double cli[N+M];
double clss[N+M];
double qrel[N+M];
double i2f[N+M];
double icalc[N+M];
double icalo[N+M];
double icalci[N+M];
double icaloi[N+M];
double icalcstar[N+M];
double icalostar[N+M];
double icalcistar[N+M];
double icaloistar[N+M];
double iksc1[N+M];
double iksc2[N+M];
double iksc3[N+M];
double iksc4[N+M];
double iksc5[N+M];
double iksc6[N+M];
double iksc7[N+M];
double iksc8[N+M];
double iksc9[N+M];
double iksc10[N+M];
double iksc11[N+M];
double iksc12[N+M];
double iksc13[N+M];
double iksc14[N+M];
double iksc15[N+M];
double ikso1[N+M];
double ikso2[N+M];
double camkactive[N+M];
double ical[N+M];
double qrel1[N+M];
double qrel2[N+M];
double qup2[N+M];
double S[36];
double Sv[48];

// ACTION POTENTIAL DURATIONS AND DIASTOLIC INTERVALS FOR EACH BEAT
double vmax1[beats+1];			// max voltage (mV)
double dvdtmax1[beats+1];		// max1 dV/dt (mV/ms)
double apd1[beats+1];			// Action Potential Duration (ms)
double toneapd1[beats+1];		// time of dv/dt max1 (ms)
double ttwoapd1[beats+1];		// time of 90% repolarization (ms)
double trep1[beats+1];			// time of full repolarization (ms)
double di1[beats+1];			// diastolic interval (ms)
double rmbp1[beats+1];			// resting membrane potential (mV)
double caimax1[beats+1];
double naimax1[beats+1];
int arrayloop;

double vmax79[beats+1];			// max79 voltage (mV)
double dvdtmax79[beats+1];		// max79 dV/dt (mV/ms)
double apd79[beats+1];			// Action Potential Duration (ms)
double toneapd79[beats+1];		// time of dv/dt max79 (ms)
double ttwoapd79[beats+1];		// time of 90% repolarization (ms)
double trep79[beats+1];			// time of full repolarization (ms)
double di79[beats+1];			// diastolic interval (ms)			// resting membrane potential (mV)
double caimax79[beats+1];
double naimax79[beats+1];
double rmbp79[beats+1];

double vmax80[beats+1];			// max80 voltage (mV)
double dvdtmax80[beats+1];		// max80 dV/dt (mV/ms)
double apd80[beats+1];			// Action Potential Duration (ms)
double toneapd80[beats+1];		// time of dv/dt max80 (ms)
double ttwoapd80[beats+1];		// time of 90% repolarization (ms)
double trep80[beats+1];			// time of full repolarization (ms)
double di80[beats+1];			// diastolic interval (ms)		// resting membrane potential (mV)
double caimax80[beats+1];
double naimax80[beats+1];
double rmbp80[beats+1];

double vmax159[beats+1];			// max voltage (mV)
double dvdtmax159[beats+1];		// max dV/dt (mV/ms)
double apd159[beats+1];			// Action Potential Duration (ms)
double toneapd159[beats+1];		// time of dv/dt max (ms)
double ttwoapd159[beats+1];		// time of 90% repolarization (ms)
double trep159[beats+1];			// time of full repolarization (ms)
double di159[beats+1];			// diastolic interval (ms)		// resting membrane potential (mV)
double caimax159[beats+1];
double naimax159[beats+1];
double rmbp159[beats+1];

// OUTPUT FILE
FILE *ap;
FILE *ecg;
FILE *apd;

FILE *VSState;
FILE *PSState;
int count;

// FUNCTIONS
void printSSP(char *infile);
void printSSV(char *infile);
void readSSP(char *readfile);
void readSSV(char *readfile);
void updateS();
void updateSv();
void pGeometry();
void pConstants();
void vGeometry();
void vConstants();
void delaytoscreen();
void record_variables();
void timestep ();
void Pcomp_revs (int n);
void Pcomp_ina (int n);
void Pcomp_ical (int n);
void Pcomp_icat (int n);
void Pcomp_ik1 (int n);
void Pcomp_iks (int n);
void Pcomp_ikp (int n);
void Pcomp_icab (int n);
void Pcomp_ikb (int n);
void Pcomp_inab (int n);
void Pcomp_ikr_block (int n, double block); // new: ikr = ikr*block;
void Pcomp_ikr (int n);
void Pcomp_inaca (int n);
void Pcomp_inak (int n);
void Pcomp_inal (int n);
void Pcomp_ipca (int n);
void Pcomp_ito1 (int n);
void Pcomp_ito1_block (int n, double block); // ito1 = ito1*block;
void Pcomp_istim (int n);
void Pcomp_itot (int n);
void Pcomp_qup1 (int n);
void Pcomp_qup2 (int n);
void Pcomp_qrel1 (int n);
void Pcomp_qrel2 (int n);
void Pcomp_qrel2DAD (int n); // for DAD simulations
void Pcomp_ip3 (int n);
void Pcomp_qtr1 (int n);
void Pcomp_qtr2 (int n);
void Pcomp_conc (int n);
void Pcomp_if (int n);


void Vcomp_revs (int n);			// function to calculate reversal potentials
void Vcomp_ina (int n);			// function to calculate the fast sodium current
void Vcomp_ical (int n);			// function to calculate the L-type calcium current
void Vcomp_ik1 (int n);			// function to calculate the time-independent potassium current
void Vcomp_iks (int n);			// function to calculate the slow delayed rectifier potassium current
void Vcomp_ikp (int n);			// function to calculate the plateau potassium current
void Vcomp_icab (int n);			// function to calculate the background calcium current
void Vcomp_iclb (int n);			// function to calculate the background chloride current
void Vcomp_ikr (int n);			// function to calculate the rapid delayedrectifier potassium current
void Vcomp_ikr_block (int n, double block);
void Vcomp_inaca (int n);			// function to calculate the sodium-calcium exchanger current
void Vcomp_inak (int n);			// function to calculate the sodium-potassium pump current
void Vcomp_inal (int n);			// function to calculate the slowly inactivating late sodium current
void Vcomp_ipca (int n);			// function to calculate the sarcolemmal calcium pump current
void Vcomp_ito1 (int n);
void Vcomp_ito1_block (int n, double block);			// function to calculate the 4AP-sensitive transient outward potassium current
void Vcomp_ito2 (int n);			// function to calculate the calcium-sensitive transient outward chloride current
void Vcomp_inab (int n);
void Vcomp_kcl (int n);			// function to calculate the potassium-chloride cotransporter
void Vcomp_nacl ();			// function to calculate the sodium-chloride cotransporter
void Vcomp_istim (int n);			// function to calculate the stimulus current
void Vcomp_itot (int n);			// function to calculate the total transmembrane current
void Vcomp_qleak (int n);			// function to calculate the SR leak flux
void Vcomp_qup (int n);			// function to calculate the SR uptake flux
void Vcomp_qrel (int n);			// function to calculate the SR release flux
void Vcomp_qrelDAD (int n);
void Vcomp_qtr (int n);			// function to calculate the SR transfer flux
void Vcomp_conc (int n);			// function to calculate concentration changes

void printtofile ();			// function for printing to file
void printAPD();

int n;
int nn;
int ii;
int jj;

double v1;   // for printing delay
double v2;



//WRITE SteadyState TO FILE? (values just before S2 beat) yes, write = 1. no, write = 0;
const int write = 0;
char filename1[255],filename2[255],filename3[255];

int main ()
{
	printf("Please enter£ºICKr ICCaL ICNa\n"); 
	scanf("%lf %lf %lf",&ICKr,&ICCaL,&ICNa);
    printf("They are£º%lf %lf %lf\n",ICKr,ICCaL,ICNa);
    
	// OPEN OUTPUT FILEs
	printf("Please enter£ºDrugName\n"); 
	scanf("%s",&DrugName);
	printf("They are£º%s\n",DrugName);
	sprintf(filename1, "%s-MV-adaptation.txt", DrugName);
	sprintf(filename2, "%s-ECG-adaptation.txt", DrugName);
	sprintf(filename3, "%s-APD-adaptation.txt", DrugName);

	ap =  fopen(filename1,"w");
	ecg = fopen(filename2,"w");
	apd = fopen(filename3,"w");
	
    for (run=0; run<(BCLMax-BCLMin)/dBCL; run++) {
	BCL=BCLMax-run*dBCL;
	S2=BCL;
    // initial values for cells (SS-files)
    //readSSP("SSfiles/PSState300msER321.txt");
    //readSSV("SSfiles/VSState300msER321.txt");

    //Values from Single Cell Steady State simulations
    
    for (n = 0 ; n < N; n++)
    {
        dvdt[n] = 0;
        
        //Single Cell SS 500ms P
        v		[n] = -84.697615;
        m		[n] = 0.000736;
        h		[n] = 0.996310;
        j		[n] = 0.999890;
        ical		[n] = -0.001113;
        d		[n] = 0.000015;
        f		[n] = 0.999244;
        f2		[n] = 0.884546;
        fca		[n] = 0.967652;
        fca2		[n] = 0.500641;
        xs1		[n] = 0.001064;
        xs2		[n] = 0.071040;
        xr		[n] = 0.014562;
        a		[n] = 0.000112;
        i		[n] = 0.958831;
        i2		[n] = 0.486621;
        ml		[n] = 0.000300;
        ml3		[n] = 0.042652;
        hl		[n] = 0.891296;
        hl3		[n] = 0.802605;
        jl		[n] = 0.000000;
        jl3     [n] = 0.000000;
        casss		[n] = 0.000147;
        cajsr		[n] = 2.484732;
        cacsr		[n] = 2.436541;
        cansr		[n] = 2.528173;
        cassl		[n] = 0.000147;
        nai		[n] = 14.272675;
        nassl		[n] = 14.272110;
        nasss		[n] = 14.272120;
        ki		[n] = 133.506492;
        cli		[n] = 20.000000;
        cai		[n] = 0.000063;
        b	[n] = 0.000399;
        g	[n] = 0.986255;
        u [n] = 0.557408;
        y [n] = 0.112337;
        camktrap [n] = 0.102607;
        camkactive[n] = 0;
        ical[n] = 0;
        qrel1[n] = 0;
        qrel2[n] = 0;
        qup2[n] = 0;
      
    }

    
    for (n = N; n < N+M; n++)
    {
        dvdt[n] = 0;
        //Single Cell SS 500ms V
        v		[n] = -87.266841;
        cai		[n] = 0.000126;
        casssr		[n] = 0.000173;
        casscal		[n] = 0.000185;
        cajsr		[n] = 1.136221;
        cansr		[n] = 1.185860;
        camktrap		[n] = 0.070306;
        nai		[n] = 8.863280;
        nasssr		[n] = 8.863319;
        ki		[n] = 143.954808;
        cli		[n] = 21.805698;
        clss		[n] = 21.805698;
        qrel		[n] = 0.001098;
        m		[n] = 0.001058;
        h		[n] = 0.990897;
        j		[n] = 0.996946;
        ml		[n] = 0.001058;
        hl		[n] = 0.239815;
        a		[n] = 0.000018;
        i		[n] = 0.997544;
        i2		[n] = 0.709993;
        i2f		[n] = 0.995940;
        xr		[n] = 0.004026;
        icalc		[n] = 0.999931;
        icalo		[n] = 0.000023;
        icalci		[n] = 0.000000;
        icaloi		[n] = 0.000000;
        icalcstar		[n] = 0.000044;
        icalostar		[n] = 0.000000;
        icalcistar		[n] = 0.000002;
        icaloistar		[n] = 0.000000;
        iksc1		[n] = 0.468140;
        iksc2		[n] = 0.251734;
        iksc3		[n] = 0.050763;
        iksc4		[n] = 0.004550;
        iksc5		[n] = 0.000153;
        iksc6		[n] = 0.139454;
        iksc7		[n] = 0.056250;
        iksc8		[n] = 0.007564;
        iksc9		[n] = 0.000339;
        iksc10		[n] = 0.015614;
        iksc11		[n] = 0.004203;
        iksc12		[n] = 0.000283;
        iksc13		[n] = 0.000791;
        iksc14		[n] = 0.000109;
        iksc15		[n] = 0.000023;
        ikso1		[n] = 0.000024;
        ikso2		[n] = 0.000006;
        camkactive[n] = 0;
        ical[n] = 0;
        qrel1[n] = 0;
        qrel2[n] = 0;
        qup2[n] = 0;
    }
    
    // TIME SETTINGS
	dvdtclock	= 1000;
	tmax		= BCL*(beats-1)+BCLMax+500;
	t		= 0;
	T       = 0;
	tstim		= 0;
	stimtime	= 1000;
	stimcount 	= -1;
	
	// MAIN LOOP
	while (t<=tmax)
	{
        for (n = 1; n < N; n++)
        {
            if (n == 1)
            {
                dt = dtmax;
                
                
                pGeometry();
                timestep();
                pConstants();
               
                
            }
        
            Pcomp_revs (n);
            Pcomp_ina (n);
            Pcomp_inal (n);
            Pcomp_inab (n);
            Pcomp_ical (n);
            Pcomp_icat (n);
            Pcomp_icab (n);
            Pcomp_ito1 (n);
            Pcomp_ikr (n);
            //Pcomp_ito1_block (n,0.2);
            //Pcomp_ikr_block (n,0.2);
            Pcomp_iks (n);
            Pcomp_ik1 (n);
            Pcomp_inaca (n);
            Pcomp_inak (n);
            Pcomp_ipca (n);
            Pcomp_if (n);
            Pcomp_istim (n);    // comment out if inserting stimulus in V
            Pcomp_itot (n);
            Pcomp_ip3 (n);
            Pcomp_qrel1 (n);
            Pcomp_qrel2 (n);
            //Pcomp_qrel2DAD(n);
            Pcomp_qup1 (n);
            Pcomp_qup2 (n);
            Pcomp_qtr1 (n);
            Pcomp_qtr2 (n);
            Pcomp_conc (n);
            
            if (n == N-1)
            {
                dvdt[n] = ((Dv*ER)*(v[n+1]-v[n]) + Dp*(v[n-1] - v[n]))/(dx*dx) - itot;
            }
            else
            {
                
            
                dvdt[n]	= Dp*(v[n-1] + v[n+1] - 2*v[n])/(dx*dx)-itot;
            }
        }
        
        for (n = N; n < N+M-1; n++)
        {
            if (n == N)
            {
                vGeometry();
                vConstants();
            }
            
            Vcomp_revs (n);			// function to calculate reversal potentials
            Vcomp_ina (n);			// function to calculate the fast sodium current
            Vcomp_ical (n);			// function to calculate the L-type calcium current
            Vcomp_ik1 (n);			// function to calculate the time-independent potassium current
            Vcomp_iks (n);			// function to calculate the slow delayed rectifier potassium current
            Vcomp_ikp (n);			// function to calculate the plateau potassium current
            Vcomp_icab (n);			// function to calculate the background calcium current
            Vcomp_iclb (n);			// function to calculate the background chloride current
            Vcomp_ikr (n);
            //Vcomp_ikr_block (n, 0.2);	// function to calculate the rapid delayed rectifier potassium current
            Vcomp_inaca (n);			// function to calculate the sodium-calcium exchanger current
            Vcomp_inak (n);			// function to calculate the sodium-potassium pump current
            Vcomp_inal (n);			// function to calculate the slowly inactivating late sodium current
            Vcomp_ipca (n);			// function to calculate the sarcolemmal calcium pump current
            Vcomp_ito1 (n);
            //Vcomp_ito1_block (n,0.2);			// function to calculate the 4AP-sensitive transient outward potassium current
            Vcomp_ito2 (n);			// function to calculate the calcium-sensitive transient outward chloride current
            Vcomp_inab (n);
            Vcomp_kcl (n);			// function to calculate the potassium-chloride cotransporter
            Vcomp_nacl (n);			// function to calculate the sodium-chloride cotransporter
            
            // remove comments if inserting stimulus in V
            //Vcomp_istim (n);			// function to calculate the stimulus current
            
            Vcomp_itot (n);			// function to calculate the total transmembrane current
            Vcomp_qleak (n);			// function to calculate the SR leak flux
            Vcomp_qup (n);			// function to calculate the SR uptake flux
            Vcomp_qrel (n);			// function to calculate the SR release flux
            //Vcomp_qrelDAD(n);
            Vcomp_qtr (n);			// function to calculate the SR transfer flux
            Vcomp_conc (n);
            
            if (n == N)
            {
                dvdt[n] = Dv*((1/ER)*(v[n-1]-v[n]) + (v[n+1] - v[n]))/(dx*dx) - itot;
            }
            else
            {
                dvdt[n]	= Dv*(v[n-1] + v[n+1] - 2*v[n])/(dx*dx)-itot;
            }
            
            
            
        }
        
        
        for (n = 1; n < N+M-1; n++)
        {
            v[n] += dvdt[n]*dt;
            
        }
        
        
        v[0] = v[1];
        v[N+M-1] = v[N+M-2];
        
        delaytoscreen();
        /*
        if ((v[79] >= -41) && rmbp79[stimcount] == 0)
        {
            rmbp79[stimcount] = v[79];
            //printf("79: %f \n", rmbp79[stimcount]);
        }
        if (v[80] >= -41 && rmbp80[stimcount] == 0)
        {
            rmbp80[stimcount] = v[80];
        }
        if (v[159] >= -41 && rmbp159[stimcount] == 0)
        {
            rmbp159[stimcount] = v[159];
        }
        */
		record_variables();
       	
		printtofile ();
		
        
		t	+= dt;
        
	}

    //printAPD();
    printf("numLines: %d \n",count);
	}
	fclose(ap);
	fclose(ecg);
	fclose(apd);
	return(1);
}


void printAPD()
{
    for (arrayloop=0; arrayloop<=beats; arrayloop++)
    {
		apd1[arrayloop] 	= ttwoapd1[arrayloop]-toneapd1[arrayloop];
        printf("APD90_1 : %f  %d \n", apd1[arrayloop], arrayloop);
        fprintf(apd, "APD90_1 : %f  %d \n", apd1[arrayloop], arrayloop);
    }
	for (arrayloop=1; arrayloop<=beats; arrayloop++)
    {
		di1[arrayloop] 	= toneapd1[arrayloop]-ttwoapd1[arrayloop-1];
        
    }
    for (arrayloop=0; arrayloop<=beats; arrayloop++){
		apd79[arrayloop] 	= ttwoapd79[arrayloop]-toneapd79[arrayloop];
        printf("APD90_80 : %f  %d \n", apd79[arrayloop], arrayloop);
        fprintf(apd, "APD90_80 : %f  %d \n", apd1[arrayloop], arrayloop);
    }
	for (arrayloop=1; arrayloop<=beats; arrayloop++)
		di79[arrayloop] 	= toneapd79[arrayloop]-ttwoapd79[arrayloop-1];
    
    
    for (arrayloop=0; arrayloop<=beats; arrayloop++)
    {
		apd80[arrayloop] 	= ttwoapd80[arrayloop]-toneapd80[arrayloop];
        printf("APD90_81 : %f  %d \n", apd80[arrayloop], arrayloop);
        fprintf(apd, "APD90_81 : %f  %d \n", apd1[arrayloop], arrayloop);
    }
	for (arrayloop=1; arrayloop<=beats; arrayloop++)
		di80[arrayloop] 	= toneapd80[arrayloop]-ttwoapd80[arrayloop-1];
    
    for (arrayloop=0; arrayloop<=beats; arrayloop++)
    {
		apd159[arrayloop] 	= ttwoapd159[arrayloop]-toneapd159[arrayloop];
        printf("APD90_160 : %f  %d \n", apd159[arrayloop], arrayloop);
        fprintf(apd, "APD90_160 : %f  %d \n", apd1[arrayloop], arrayloop);
    }
	for (arrayloop=1; arrayloop<=beats; arrayloop++)
		di159[arrayloop] 	= toneapd159[arrayloop]-ttwoapd159[arrayloop-1];

}

void timestep ()
{
    for (nn = 1; nn < M+N-1; nn++)
    {
        if (dt == dtmed)
        {
            if ((dvdt[nn]>dvdtthresh) || (t>(tstim-2)) || (stimtime<(stimdur+2)) || (dvdtclock<5))
            {
                dt = dtmin;
            }
        }
        
        else if (dt == dtmax)
        {
            if ((dvdt[nn]>dvdtthresh) || (t>(tstim-2)) || (stimtime<(stimdur+2)) || (dvdtclock<5))
            {
                dt = dtmin;
            }
            else if (dvdtclock>=5 && dvdtclock<20)
            {
                    dt = dtmed;
            }
        }
        
        if (dvdt[nn]>dvdtthresh) dvdtclock = 0;
    
    }
    
    dvdtclock += dt;
	
}


void Pcomp_revs (int n)
{
    eca	= (R*temp/(zca*frdy))*log(cao/cassl[n]);
	ena	= (R*temp/frdy)*log(nao/nassl[n]);
	ek	= (R*temp/frdy)*log(ko/ki[n]);
}

void Pcomp_ina (int n)
{
    ma	= 0.64*(v[n]+37.13)/(1-exp(-0.1*(v[n]+37.13)));
	mb	= 0.16*exp(-v[n]/11);
	if (v[n]<-40)
	{
		ha = 0.135*exp((70+v[n])/-6.8);
		hb = 3.56*exp(0.079*v[n])+310000*exp(0.35*v[n]);
		ja = (-127140*exp(0.2444*v[n])-0.003474*exp(-0.04391*v[n]))*(v[n]+37.78)/(1+exp(0.311*(v[n]+79.23)));
		jb = 0.1212*exp(-0.01052*v[n])/(1+exp(-0.1378*(v[n]+40.14)));
	}
	else
	{
		ha = 0.0;
		hb = 1/(0.13*(1+exp((v[n]+10.66)/-11.1)));
		ja = 0.0;
		jb = 0.3*exp(-0.0000002535*v[n])/(1+exp(-0.1*(v[n]+32)));
	}
	mtau	= 1/(ma+mb);
	htau	= 1/(ha+hb);
	jtau	= 1/(ja+jb);
	mss	= ma*mtau;
	hss	= ha*htau;
	jss	= 1*ja*jtau;
	m[n]	= mss-(mss-m[n])*exp(-dt/mtau);
	h[n]	= hss-(hss-h[n])*exp(-dt/htau);
	j[n]	= jss-(jss-j[n])*exp(-dt/jtau);
	ina	= ICNa*gna*pow(m[n],3)*h[n]*j[n]*(v[n]-ena);
}

void Pcomp_inal (int n)
{
	mltau	= 1/(0.64*(v[n]+37.13)/(1-exp(-0.1*(v[n]+37.13))) + 0.16*exp(-v[n]/11));
	ml3tau  = mltau;
	mlss	= 1/(1+exp(-(v[n]+28)/7));
	ml3ss   = 1/(1+exp(-(v[n]+63)/7));
	hltau   = 162+132/(1+exp(-(v[n]+28)/5.5));
	hl3tau  = 0.5*hltau;
	hlss	= 1/(1+exp((v[n]+28)/12));
	hl3ss	= 1/(1+exp((v[n]+63)/12));
	jltau   = 411;
	jl3tau  = 0.5*jltau;
	jlss	= hlss;
	jl3ss	= hl3ss;
	ml[n]	    = mlss-(mlss-ml[n])*exp(-dt/mltau);
	ml3[n]     = ml3ss-(ml3ss-ml3[n])*exp(-dt/ml3tau);
	hl[n]	    = hlss-(hlss-hl[n])*exp(-dt/hltau);
	hl3[n]     = hl3ss-(hl3ss-hl3[n])*exp(-dt/hl3tau);
	jl[n]	    = jlss-(jlss-jl[n])*exp(-dt/jltau);
	jl3[n]     = jl3ss-(jl3ss-jl3[n])*exp(-dt/jl3tau);
	inal2   = gnal2*ml[n]*hl[n]*jl[n]*(v[n]-ena);
	inal3   = gnal3*ml3[n]*hl3[n]*jl3[n]*(v[n]-ena);
	inal    = inal2 + inal3;
}

void Pcomp_ical (int n)
{
	ibarca		= pca*zca*zca*(((v[n]-15)*frdy*frdy)/(R*temp))*((gacai*casss[n]*exp((zca*(v[n]-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(v[n]-15)*frdy)/(R*temp))-1));
	dss		    = (1/(1.0+exp(-(v[n]-2.0)/7.8)));
	dtau		= (0.59+0.8*exp(0.052*(v[n]+13))/(1+exp(0.132*(v[n]+13))));
	fss	        = 1/(1.0 + exp((v[n]+16.5)/9.5));
	ftau        = 0.92/(0.125*exp(-(0.058*(v[n]-2.5))*(0.045*(v[n]-2.5)))+0.1);
	f2ss        = fss;
	f2tau       = 0.90/(0.02*exp(-(0.04*(v[n]-18.6))*(0.045*(v[n]-18.6)))+0.005);
	fcass		= 0.3/(1 - ical[n]/0.05) + 0.55/(1.0+casss[n]/0.003)+0.15;
	fcatau		= 10*camkactive[n]/(camkactive[n]+kmcam) + 0.5+1/(1.0+casss[n]/0.003);
	fca2ss		= 1.0/(1.0-ical[n]/0.01);
	fca2tau		= 1*(300.0/(1.0+exp((-ical[n]-0.175)/0.04))+125.0);
	d[n]		    = dss-(dss-d[n])*exp(-dt/dtau);
	f[n]		    = fss-(fss-f[n])*exp(-dt/ftau);
	f2[n]		    = f2ss-(f2ss-f2[n])*exp(-dt/f2tau);
	fca[n]		    = fcass-(fcass-fca[n])*exp(-dt/fcatau);
	fca2[n]		= fca2ss-(fca2ss-fca2[n])*exp(-dt/fca2tau);
	ical[n]		= ICCaL*d[n]*f[n]*f2[n]*fca[n]*fca2[n]*ibarca;
    
}

void Pcomp_icat (int n)
{
	bss	    = 1/(1+ exp (-(v[n]+30)/7));
	gss	    = 1/(1+exp((v[n]+61)/5));
	taub	= 1/(1.068*exp((v[n]+16.3)/30)+1.068*exp(-(v[n]+16.3)/30));
	taug    = 1/(0.015*exp(-(v[n]+71.7)/83.3)+0.015*exp((v[n]+71.7)/15.4));
	b[n]	    = bss-(bss-b[n])*exp(-dt/taub);
	g[n]	    = gss-(gss-g[n])*exp(-dt/taug);
	icat	= gcat*b[n]*g[n]*(v[n]-eca);
}

void Pcomp_ito1 (int n)
{
	atau	= 1/(25*exp((v[n]-82)/18)/(1+exp((v[n]-82)/18))+25*exp(-(v[n]+52)/18)/(1+exp(-(v[n]+52)/18)));
	itau	= 2.86+ 1/(exp(-(v[n]+125)/15)*0.1 + 0.1*exp((v[n]+2)/26.5));
	i2tau	= 21.5+ 1/(exp(-(v[n]+138.2)/52)*0.005 + 0.003*exp((v[n]+18)/12.5));
	ass	    = 1/(1+exp(-(v[n]-8.9)/10.3));
	iss	    = 1/(1+exp((v[n]+30)/11));
	i2ss	= iss;
	a[n]	    = ass-(ass-a[n])*exp(-dt/atau);
	i[n]	    = iss-(iss-i[n])*exp(-dt/itau);
	i2[n]	    = i2ss-(i2ss-i2[n])*exp(-dt/i2tau);
	itos    = gtos*a[n]*i[n]*i2[n]*(v[n]-ek);
	itof    = gtof*(v[n]-ek)/(1+exp(-(v[n]-3)/19.8));
	ito1	= itos + itof;
}

void Pcomp_ito1_block (int n, double block)
{
	atau	= 1/(25*exp((v[n]-82)/18)/(1+exp((v[n]-82)/18))+25*exp(-(v[n]+52)/18)/(1+exp(-(v[n]+52)/18)));
	itau	= 2.86+ 1/(exp(-(v[n]+125)/15)*0.1 + 0.1*exp((v[n]+2)/26.5));
	i2tau	= 21.5+ 1/(exp(-(v[n]+138.2)/52)*0.005 + 0.003*exp((v[n]+18)/12.5));
	ass	    = 1/(1+exp(-(v[n]-8.9)/10.3));
	iss	    = 1/(1+exp((v[n]+30)/11));
	i2ss	= iss;
	a[n]	    = ass-(ass-a[n])*exp(-dt/atau);
	i[n]	    = iss-(iss-i[n])*exp(-dt/itau);
	i2[n]	    = i2ss-(i2ss-i2[n])*exp(-dt/i2tau);
	itos    = gtos*a[n]*i[n]*i2[n]*(v[n]-ek);
	itof    = gtof*(v[n]-ek)/(1+exp(-(v[n]-3)/19.8));
	ito1	= itos + itof;
    ito1 = ito1*block;
}

void Pcomp_ikr (int n)
{
	gkr	    = 0.0326*sqrt(ko/5.4);
	xrss	= 1/(1+exp(-(v[n])/15));
	xrtau   = 400.0/(1.0+exp(v[n]/10.0)) + 100.0;
	rkr	    = 1/(1+exp((v[n])/35));
	xr[n]	    = xrss-(xrss-xr[n])*exp(-dt/xrtau);
	ikr	    = ICKr*gkr*xr[n]*rkr*(v[n]-ek); 
}
void Pcomp_ikr_block (int n, double block)
{
	gkr	    = 0.0326*sqrt(ko/5.4);
	xrss	= 1/(1+exp(-(v[n])/15));
	xrtau   = 400.0/(1.0+exp(v[n]/10.0)) + 100.0;
	rkr	    = 1/(1+exp((v[n])/35));
	xr[n]	    = xrss-(xrss-xr[n])*exp(-dt/xrtau);
	ikr	    = gkr*xr[n]*rkr*(v[n]-ek);
    ikr = ikr*block*ICKr;
}

void Pcomp_iks (int n)
{
	eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(ki[n]+prnak*nassl[n]));
	gks	    = 0.053*(1+0.6/(1+pow((0.000038/cassl[n]),1.4)));
	xsss	= 1/(1+exp(-(v[n]-9)/13.7));
	xs1tau	= 200/(exp(-(v[n]+10)/6) + exp((v[n]-62)/55));
	xs2tau	= 1500+ 350/(exp(-(v[n]+10)/4) + exp((v[n]-90)/58));
	xs1[n]	    = xsss-(xsss-xs1[n])*exp(-dt/xs1tau);
	xs2[n]	    = xsss-(xsss-xs2[n])*exp(-dt/xs2tau);
	iks	    = gks*xs1[n]*xs2[n]*(v[n]-eks);
}

void Pcomp_ik1 (int n)
{
    k1ss      = 1/(1+exp((v[n]+103-(2.9+ko*2.175))/10.15));
	gk1	      = 0.12*sqrt(ko);
	ik1	      = gk1*k1ss*(v[n]-ek);
}

void Pcomp_if (int n)
{
	yss       = 1/(1+exp((v[n]+87)/9.5));
	ytau      = 2000/(exp(-(v[n]+132)/10) + exp((v[n]+57)/60));
	y[n]         = yss - (yss-y[n])*exp(-dt/ytau);
	ifna	  = 0.012*y[n]*y[n]*(v[n]-ena);
	ifk       = 0.024*y[n]*y[n]*(v[n]-ek);
	iftotal   = ifna + ifk;
}

void Pcomp_inaca (int n)
{
	allo		= 1/(1+pow((kmcaact/(1.5*casss[n])),2));
	num		    = inacamax*(pow(nasss[n],3)*cao*exp(nu*v[n]*frdy/(R*temp))-pow(nao,3)*1.5*casss[n]*exp((nu-1)*v[n]*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v[n]*frdy/(R*temp));
	denomterm1	= kmcao*pow(nasss[n],3)+pow(kmnao,3)*1.5*casss[n]+pow(kmnai1,3)*cao*(1+1.5*casss[n]/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nasss[n]/kmnai1,3))+pow(nasss[n],3)*cao+pow(nao,3)*1.5*casss[n];
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inacass		= 0.2*allo*deltaE;
	
	allo		= 1/(1+pow((kmcaact/(1.5*cassl[n])),2));
	num		    = inacamax*(pow(nassl[n],3)*cao*exp(nu*v[n]*frdy/(R*temp))-pow(nao,3)*1.5*cassl[n]*exp((nu-1)*v[n]*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v[n]*frdy/(R*temp));
	denomterm1	= kmcao*pow(nassl[n],3)+pow(kmnao,3)*1.5*cassl[n]+pow(kmnai1,3)*cao*(1+1.5*cassl[n]/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nassl[n]/kmnai1,3))+pow(nassl[n],3)*cao+pow(nao,3)*1.5*cassl[n];
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inaca		= 0.8*allo*deltaE;
}

void Pcomp_inak (int n)
{
    inak	= ibarnak*(1/(1+exp(-1*(v[n]+92)*frdy/(R*temp))))*pow((nassl[n]/(nassl[n]+2.6)),3)*(ko/(ko+0.8));
}

void Pcomp_ipca (int n)
{
	ipca	= ipcabar/((kmpca/cassl[n])+1);
}

void Pcomp_icab (int n)
{
	icab	= pcab*zca*zca*((v[n]*frdy*frdy)/(R*temp))*((gacai*cassl[n]*exp((zca*v[n]*frdy)/(R*temp))-gacao*cao)/(exp((zca*v[n]*frdy)/(R*temp))-1));
}

void Pcomp_inab (int n)
{
    inab    = pnab*frdy*((frdy*v[n])/(R*temp))*(nassl[n]*exp((frdy*v[n])/(R*temp)) - nao)/(exp((frdy*v[n])/(R*temp))-1);
}

void Pcomp_istim (int n)
{
    if (n == 1) // Stimulus in cell no 1. (cell 0 is cell at boundary)
    {
        stimtime += dt;
    }
    
	if (t>=tstim && n == 1)
	{
		stimtime = 0.0;
		stimcount += 1;
		if (stimcount < beats-1)  tstim += BCL;
		else if (stimcount == beats-1) tstim += S2;
		else tstim = tmax+1;
		if (stimcount < beats) printf ("S1 Beat %d at time = %.2f ms in cell nr %d!\n", stimcount+1, t,n);
		else if (stimcount == beats)
        {
            printf ("S2 Beat at time = %.2f ms in cell nr %d!\n", t,n);
            if (write == 1)
            {
            printSSP("PSState340msER321.txt");
            printSSV("VSState340msER321.txt");
            printf("Steady State files written\n");
            }
        }
        rmbp1[stimcount] = v[1];
	}
}

void Pcomp_itot (int n)
{
    
    
    //////////////////////////////////////////////////
    
	if (stimtime>=0.0 && stimtime<stimdur && (n==1))
	{
		icatot	= ical[n]+icat+ipca+icab-2*inaca-2*inacass;
		iktot	= ikr+iks+ik1-2*inak+ito1+ifk+1*istim;
		inatot	= 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
		itot	= icatot+iktot+inatot;
	}
    
    ///////////////////////////////////////////////// comment out this part if inserting stimulus in V-Cell!
	else
	{
	
		icatot	= ical[n]+icat+ipca+icab-2*inaca-2*inacass;
		iktot	= ikr+iks+ik1-2*inak+ito1+ifk;
		inatot	= 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
		itot	= icatot+iktot+inatot;
	}
    
}

void Pcomp_ip3 (int n)
{
    du= dt*(casss[n]*k2*(1-u[n]) - k2a*u[n]);
    u[n] += du;
    POip3 = tauip3r*IP3*casss[n]*(1-u[n])/((1+IP3*k0/k0a)*(1+casss[n]*k1/k1a));
    qip3 = 10.920*(cajsr[n]-casss[n])*(POip3);
}

void Pcomp_qrel1 (int n)
{
	qdiff  = (casss[n]-cassl[n])/sstau;
    REL  = -((ical[n])*acap/(vsss*2.0*frdy) - (qrel1[n] + qip3)*vjsr/vsss + qdiff);
    ireltau = 2*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))/(1+(0.0123/cajsr[n]));
    if (REL > 0)
    {irelss  = 15*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))*REL/(1 + pow((1.0/cajsr[n]),8));}
    else {irelss = 0;}
    qrel1[n] += dt*((irelss-qrel1[n])/ireltau);
    
}

void Pcomp_qrel2 (int n)
{
	qgap  = (cassl[n]-cai[n])/gaptau;
    REL  = (-qup2[n]*vnsr/vmyo + qgap*vssl/vmyo+ (qrel2[n])*vcsr/vmyo);
    ireltau = 6*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))/(1+(0.0123/cacsr[n]));
    if (REL > 0)
    {irelss  = 91*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))*(REL)/(1 + pow((1/cacsr[n]),8));}
    else {irelss = 0;}
    qrel2[n] += dt*((irelss-qrel2[n])/ireltau);
}

void Pcomp_qrel2DAD(int n)
{
    qgap  = (cassl-cai)/gaptau;
    REL  = (-qup2[n]*vnsr/vmyo + qgap*vssl/vmyo+ (qrel2[n])*vcsr/vmyo);
    ireltau = 2*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))/(1+(0.0123/cacsr[n]));// DAD
    if (REL > 0)
    {irelss  = 91*(1+1*(1/(1+pow((0.28/camkactive[n]),8))))*(REL)/(1 + pow((0.25/cacsr[n]),8));}//DAD
    else {irelss = 0;}
    qrel2[n] += dt*((irelss-qrel2[n])/ireltau);
}

void Pcomp_qup1(int n)
{
    dkmplb		= dkmplbbar*camkactive[n]/(kmcamk+camkactive[n]);
	dqupcamk	= dqupcamkbar*camkactive[n]/(kmcamk+camkactive[n]);
	qup1		= 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cassl[n],1))-0.00105*cansr[n]/nsrbar;
}

void Pcomp_qup2 (int n)
{
    dkmplb		= dkmplbbar*camkactive[n]/(kmcamk+camkactive[n]);
	dqupcamk	= dqupcamkbar*camkactive[n]/(kmcamk+camkactive[n]);
	qup2[n]		= 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cai[n],1))-0.0042*cansr[n]/nsrbar;
}

void Pcomp_qtr1 (int n)
{
	qtr1		= (cansr[n]-cajsr[n])/tautr1;
}

void Pcomp_qtr2 (int n)
{
	qtr2		= (cansr[n]-cacsr[n])/tautr2;
}

void Pcomp_conc (int n)
{
	qdiff       = (casss[n]-cassl[n])/sstau;
	qgap        = (cassl[n]-cai[n])/gaptau;
    qdiffna     = (nasss[n]-nassl[n])/sstau;
    qgapna      = (nassl[n]-nai[n])/gaptau;
    
	dcasss		= dt*(-(ical[n]-2*inacass)*acap/(vsss*2.0*frdy)+(qrel1[n]+qip3)*vjsr/vsss-qdiff);
  
    
	bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+casss[n],2))+(bslbar*kmbsl/pow(kmbsl+casss[n],2)));
	casss[n]      += bsss*dcasss;
	
	dcassl		= dt*(-(qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(icat+ipca+icab-2*inaca)*acap/(vssl*2.0*frdy));
	trpn        = trpnbar1*(cassl[n]/(cassl[n]+kmtrpn));
	cmdn		= cmdnbar1*(cassl[n]/(cassl[n]+kmcmdn));
	catotal		= trpn+cmdn+dcassl+cassl[n];
	bmyo		= cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cassl[n]		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;
 	
	dcajsr		= dt*(qtr1-qrel1[n]-qip3);
	csqn1       = csqnbar1*(cajsr[n]/(cajsr[n]+kmcsqn));
	bjsr        = csqnbar1 - csqn1-cajsr[n]-dcajsr+kmcsqn;
	cjsr        = kmcsqn*(csqn1+cajsr[n]+dcajsr);
	cajsr[n]       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	
	dcacsr		= dt*(qtr2-qrel2[n]);
	csqn        = csqnbar*(cacsr[n]/(cacsr[n]+kmcsqn));
	bcsr        = csqnbar - csqn-cacsr[n]-dcacsr+kmcsqn;
	ccsr        = kmcsqn*(csqn+cacsr[n]+dcacsr);
	cacsr[n]       = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;
	
	dcansr	    = dt*(qup1+qup2[n]-qtr1*vjsr/vnsr-qtr2*vcsr/vnsr);
 	cansr[n]	   += dcansr;
 	
	dnasss	    = dt*((-(3*inacass)*acap)/((vsss)*zna*frdy)-qdiffna);
	nasss[n]      += dnasss;
	
	dnassl	    = dt*((-(3*inak+ina+inal+3*inaca+ifna+inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
	nassl[n]	   += dnassl;
	
	dnai        = dt*(qgapna*vssl/vmyo);
	nai[n]        += dnai;
	
	dki	        = dt*((-iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	ki[n]         += dki;
	
	dcai		= dt*(-(qup2[n])*vnsr/vmyo+qgap*vssl/vmyo+(qrel2[n])*vcsr/vmyo);
	trpn        = trpnbar*(cai[n]/(cai[n]+kmtrpn));
	cmdn		= cmdnbar*(cai[n]/(cai[n]+kmcmdn));
	catotal		= trpn+cmdn+dcai+cai[n];
	bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cai[n]		    = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;
	
	caavg       = (casss[n]*vsss+cassl[n]*vssl+cai[n]*vmyo)/(vsss+vmyo+vssl);
	
 	camkbound	= camk0*(1-camktrap[n])*1/(1+(kmcam/casss[n]));
	camktrap[n]	= dt*(alphacamk*camkbound*(camkbound+camktrap[n])-betacamk*camktrap[n]) + camktrap[n];
	camkactive[n]	= camkbound+camktrap[n];
	
}
    
    void Vcomp_revs (int n)
    {
        ena	= (R*temp/frdy)*log(nao/nai[n]);
        ek	= (R*temp/frdy)*log(ko/ki[n]);
        ecl	= -(R*temp/frdy)*log(clo/cli[n]);
    }
    
    // FAST SODIUM CURRENT
    void Vcomp_ina (int n)
    {
        ma	= 0.32*(v[n]+47.13)/(1-exp(-0.1*(v[n]+47.13)));
        mb	= 0.08*exp(-v[n]/11);
        if (v[n]<-40)
		{
            ha = 0.135*exp((80+v[n])/-6.8);
            hb = 3.56*exp(0.079*v[n])+310000*exp(0.35*v[n]);
            ja = (-127140*exp(0.2444*v[n])-0.00006948*exp(-0.04391*v[n]))*(v[n]+37.78)/(1+exp(0.311*(v[n]+79.23)));
            jb = 0.1212*exp(-0.01052*v[n])/(1+exp(-0.1378*(v[n]+40.14)));
		}
        else
		{
            ha = 0.0;
            hb = 1/(0.13*(1+exp((v[n]+10.66)/-11.1)));
            ja = 0.0;
            jb = 0.3*exp(-0.0000002535*v[n])/(1+exp(-0.1*(v[n]+32)));
		}
        mtau	= 1/(ma+mb);
        htau	= 1/(ha+hb);
        jtau	= 1/(ja+jb);
        mss	= ma*mtau;
        hss	= ha*htau;
        jss	= ja*jtau;
        m[n]	= mss-(mss-m[n])*exp(-dt/mtau);
        h[n]	= hss-(hss-h[n])*exp(-dt/htau);
        j[n]	= jss-(jss-j[n])*exp(-dt/jtau);
        ina	= ICNa*gna*pow(m[n],3)*h[n]*j[n]*(v[n]-ena);
    }
    
    // SLOWLY INACTIVATING LATE SODIUM CURRENT
    void Vcomp_inal (int n)
    {
        alphaml	= 0.32*(v[n]+47.13)/(1-exp(-0.1*(v[n]+47.13)));
        betaml	= 0.08*exp(-v[n]/11);
        mltau	= 1/(alphaml + betaml);
        mlss	= alphaml/(alphaml+betaml);
        hlss	= 1/(1+exp((v[n]+91.0)/6.1));
        ml[n]	= mlss-(mlss-ml[n])*exp(-dt/mltau);
        hl[n]	= hlss-(hlss-hl[n])*exp(-dt/hltau);
        //inal	= gnal*ml[n]*ml[n]*ml[n]*hl[n]*(v[n]-ena);
	
	/* Assigns parameter values to endo, M and epi cells */
	if (n <= N+ENDO)
		{
			inal	= 1.15*gnal*ml[n]*ml[n]*ml[n]*hl[n]*(v[n]-ena);	
		}
	else if (n <= N+ENDO+MID)
		{
			inal	= 1.7*gnal*ml[n]*ml[n]*ml[n]*hl[n]*(v[n]-ena);	
	    }
	else 
		{
			inal	= gnal*ml[n]*ml[n]*ml[n]*hl[n]*(v[n]-ena);	
		}
        
    }
    
    // L-TYPE CALCIUM CURRENT
    void Vcomp_ical (int n)
    {
        
        ibarca		= pca*zca*zca*((v[n]*frdy*frdy)/(R*temp))*((gacai*casscal[n]*exp((zca*v[n]*frdy)/(R*temp))-gacao*cao)/(exp((zca*v[n]*frdy)/(R*temp))-1));
        
        acttau = 0.59+0.8*((exp(0.052*(v[n]+13)))/(1+exp(0.132*(v[n]+13))));
        actss  = 1/(1.0+exp(-(v[n]-13.56)/9.45));
        icala  = actss/acttau;
        icalb  = (1-actss)/acttau;
        ivss   = (0.25+(1/(1+exp((v[n]+17.5)/3))))/1.25;
        ivtau  = 1/((1/(24.828*(1+exp((v[n]+49.10)/10.349))))+(1/(30.553*(1+exp(-(v[n]+0.213)/10.807)))));
        icalx  = ivss/ivtau;
        icaly  = (1-ivss)/ivtau;
        istarca = 25-(17.5/(1+pow((0.003/casscal[n]),4)));
        istarvtau = 1/((1/(24.828*(1+exp((v[n]+49.10)/10.349))))+(1/(istarca*(1+exp(-(v[n]+0.213)/10.807)))));
        istarvss = (0.0001+(1/(1+exp((v[n]+17.5)/3))))/1.0001;
        icalxstar = istarvss/istarvtau;
        icalystar = (1-istarvss)/istarvtau;
        icaldelta = 3/(1+pow((0.003/casscal[n]),4));
        icaltheta = 1.0;
        icalthetai = 0.000001;
        icaldeltai = icalthetai*(icalx*icalystar*icaldelta/(icaly*icalxstar*icaltheta));
        
        icalca = icalc[n];
        icaloa = icalo[n];
        icalcia = icalci[n];
        icaloia = icaloi[n];
        icalcstara = icalcstar[n];
        icalostara = icalostar[n];
        icalcistara = icalcistar[n];
        icaloistara = icaloistar[n];
        
        dicalc = dt*(icalb*icalo[n]+icaltheta*icalcstar[n]+icalx*icalci[n] - icalc[n]*(icala+icaly+icaldelta));
        dicalo = dt*(icala*icalc[n]+icaltheta*icalostar[n]+icalx*icaloi[n] - icalo[n]*(icalb+icaly+icaldelta));
        dicalcstar = dt*(icaldelta*icalc[n]+icalb*icalostar[n]+icalxstar*icalcistar[n] - icalcstar[n]*(icala+icaltheta+icalystar));
        dicalostar = dt*(icaldelta*icalo[n]+icala*icalcstar[n]+icalxstar*icaloistar[n] - icalostar[n]*(icalb+icaltheta+icalystar));
        dicalci = dt*(icaly*icalc[n]+icalthetai*icalcistar[n]+icalb*icaloi[n] - icalci[n]*(icala+icaldeltai+icalx));
        dicaloi = dt*(icaly*icalo[n]+icalthetai*icaloistar[n]+icala*icalci[n] - icaloi[n]*(icalb+icaldeltai+icalx));
        dicalcistar = dt*(icaldeltai*icalci[n]+icalystar*icalcstar[n]+icalb*icaloistar[n] - icalcistar[n]*(icala+icalthetai+icalxstar));
        dicaloistar = dt*(icaldeltai*icaloi[n]+icalystar*icalostar[n]+icala*icalcistar[n] - icaloistar[n]*(icalb+icalthetai+icalxstar));
        
        icalc[n] = icalca + dicalc;
        icalo[n] = icaloa + dicalo;
        icalci[n] = icalcia + dicalci;
        icaloi[n] = icaloia + dicaloi;
        icalcstar[n] = icalcstara + dicalcstar;
        icalostar[n] = icalostara + dicalostar;
        icalcistar[n] = icalcistara + dicalcistar;
        icaloistar[n] = icaloistara + dicaloistar;
        
        ical[n] = ICCaL*ibarca*(icalo[n]+icalostar[n]);
        
    }
    
    // BACKGROUND CALCIUM CURRENT
    void Vcomp_icab (int n)
    {
        icab	= pcab*zca*zca*((v[n]*frdy*frdy)/(R*temp))*((gacai*cai[n]*exp((zca*v[n]*frdy)/(R*temp))-gacao*cao)/(exp((zca*v[n]*frdy)/(R*temp))-1));
    }
    
    // 4AP-SENSITIVE TRANSIENT OUTWARD POTASSIUM CURRENT
    void Vcomp_ito1 (int n)
    {
        alphaa	= 1/(1.2089*(1+exp(-(v[n]-18.4099)/29.3814)));
        betaa	= 3.5/(1+exp((v[n]+100)/29.3814));
        alphai	= 0.025/(1+exp((v[n]+58)/5));
        betai	= 1/(9.7953*(1+exp(-(v[n]+19)/9)));
        //alphai2	= 1/(250*(1+exp((v[n]+60)/5)));
        
    if (n <= N+ENDO)
		{
			alphai2	= 0.0026/(1*(1+exp((v[n]+60)/5)));	
		}
	else if (n <= N+ENDO+MID)
		{
			alphai2	= 0.00225/(1*(1+exp((v[n]+60)/5)));	
	    }
	else 
		{
			alphai2	= 0.0039/(1*(1+exp((v[n]+60)/5)));	
		}
	
        betai2	= betai;
        atau	= 1/(alphaa+betaa);
        itau	= 1/(alphai+betai);
        i2tau	= 1/(alphai2+betai2);
        ass	= 1/(1+exp(-(v[n]+9.437)/7.133));
        iss	= alphai/(alphai+betai);
        i2ss	= alphai2/(alphai2+betai2);
        a[n]	= ass-(ass-a[n])*exp(-dt/atau);
        i[n]	= iss-(iss-i[n])*exp(-dt/itau);
        i2[n]	= i2ss-(i2ss-i2[n])*exp(-dt/i2tau);
        rto1	= exp(v[n]/550);
        //ito1	= gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);
        
        /* Assigns parameter values to endo, M and epi cells */
	if (n <= N+ENDO)
		{
			ito1	= 0.5*gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
		}
	else if (n <= N+ENDO+MID)
		{
			ito1	= 0.95*gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
	    }
	else 
		{
			ito1	= gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
		}
        
    }

void Vcomp_ito1_block (int n, double block)
{
    alphaa	= 1/(1.2089*(1+exp(-(v[n]-18.4099)/29.3814)));
        betaa	= 3.5/(1+exp((v[n]+100)/29.3814));
        alphai	= 0.025/(1+exp((v[n]+58)/5));
        betai	= 1/(9.7953*(1+exp(-(v[n]+19)/9)));
        //alphai2	= 1/(250*(1+exp((v[n]+60)/5)));
        
    if (n <= N+ENDO)
		{
			alphai2	= 0.0026/(1*(1+exp((v[n]+60)/5)));	
		}
	else if (n <= N+ENDO+MID)
		{
			alphai2	= 0.00225/(1*(1+exp((v[n]+60)/5)));	
	    }
	else 
		{
			alphai2	= 0.0039/(1*(1+exp((v[n]+60)/5)));	
		}
	
        betai2	= betai;
        atau	= 1/(alphaa+betaa);
        itau	= 1/(alphai+betai);
        i2tau	= 1/(alphai2+betai2);
        ass	= 1/(1+exp(-(v[n]+9.437)/7.133));
        iss	= alphai/(alphai+betai);
        i2ss	= alphai2/(alphai2+betai2);
        a[n]	= ass-(ass-a[n])*exp(-dt/atau);
        i[n]	= iss-(iss-i[n])*exp(-dt/itau);
        i2[n]	= i2ss-(i2ss-i2[n])*exp(-dt/i2tau);
        rto1	= exp(v[n]/550);
        //ito1	= gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);
        
        /* Assigns parameter values to endo, M and epi cells */
	if (n <= N+ENDO)
		{
			ito1	= block*0.5*gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
		}
	else if (n <= N+ENDO+MID)
		{
			ito1	= block*0.95*gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
	    }
	else 
		{
			ito1	= block*gbarto1*pow(a[n],3)*i[n]*i2[n]*rto1*(v[n]-ek);	
		}
}

    // RAPID DELAYED RECTIFIER POTASSIUM CURRENT
    void Vcomp_ikr (int n)
    {
        gkr	= 0.0138542*sqrt(ko/5.4);
        xrss	= 1/(1+exp(-(v[n]+10.0805)/4.25));
        //xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));
    if (n <= N+ENDO)
		{
			xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
		}
	else if (n <= N+ENDO+MID)
		{
			xrtau	= 2/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0006*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
	    }
	else 
		{
			xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
		}
	
		rkr	= 1/(1+exp((v[n]+10)/15.4));
        xr[n]	= xrss-(xrss-xr[n])*exp(-dt/xrtau);
        ikr	=ICKr* 1*gkr*xr[n]*rkr*(v[n]-ek);
    }

void Vcomp_ikr_block (int n, double block)
{
    gkr	= 0.0138542*sqrt(ko/5.4);
    xrss	= 1/(1+exp(-(v[n]+10.0805)/4.25));
    //xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));
    if (n <= N+ENDO)
		{
			xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
		}
	else if (n <= N+ENDO+MID)
		{
			xrtau	= 2/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0006*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
	    }
	else 
		{
			xrtau	= 1/(0.0006*(v[n] - 1.7384)/(1.0-exp(-0.136*(v[n] - 1.7384)))+ 0.0003*(v[n]+38.3608)/(exp(0.1522*(v[n] + 38.3608))-1.0));	
		}
	rkr	= 1/(1+exp((v[n]+10)/15.4));
    xr[n]	= xrss-(xrss-xr[n])*exp(-dt/xrtau);
    ikr	= 1*gkr*xr[n]*rkr*(v[n]-ek);
    ikr = ICKr*block*ikr;
}

    // SLOW DELAYED RECTIFIER POTASSIUM CURRENT
    void Vcomp_iks (int n)
    {
        
        eks	= (R*temp/frdy)*log((ko+prnak*nao)/(ki[n]+prnak*nai[n]));
        gks = 0.0826 * (1 + 0.6/(1 + pow((0.000038/cai[n]),1.4)));
        
        sa=1.4864e-2*exp(2.9877e-2*v[n]*(frdy/(R*temp)));
        sb=8.3986e-2*exp(-5.5461e-2*v[n]*(frdy/(R*temp)));
        sr=1.4601e-2*exp(2.4465e-1*v[n]*(frdy/(R*temp)));
        sd=3.1173e-3*exp(-4.2625e-1*v[n]*(frdy/(R*temp)));
        ikseta=7.732e-2*exp(-6.4726e-2*v[n]*(frdy/(R*temp)));
        ikstheta=8.9538e-2;
        iksomega=7.9405e-1*exp(-8.0174e-2*v[n]*(frdy/(R*temp)));
        ikspsi=5.8638e-1*exp(2.8206e-1*v[n]*(frdy/(R*temp)));
        
        iksc1a = iksc1[n];
        iksc2a = iksc2[n];
        iksc3a = iksc3[n];
        iksc4a = iksc4[n];
        iksc5a = iksc5[n];
        iksc6a = iksc6[n];
        iksc7a = iksc7[n];
        iksc8a = iksc8[n];
        iksc9a = iksc9[n];
        iksc10a = iksc10[n];
        iksc11a = iksc11[n];
        iksc12a = iksc12[n];
        iksc13a = iksc13[n];
        iksc14a = iksc14[n];
        iksc15a = iksc15[n];
        ikso1a = ikso1[n];
        ikso2a = ikso2[n];
        
        diksc1 = dt*(iksc2[n]*sb - iksc1[n]*4*sa);
        diksc2 = dt*(iksc1[n]*4*sa + iksc3[n]*2*sb + iksc6[n]*sd - iksc2[n]*(sb + 3*sa + sr));
        diksc3 = dt*(iksc2[n]*3*sa +iksc4[n]*3*sb + iksc7[n]*sd  - iksc3[n]*(2*sb + 2*sa + 2*sr));
        diksc4 = dt*(iksc3[n]*2*sa + iksc5[n]*4*sb + iksc8[n]*sd - iksc4[n]*(3*sb + sa + 3*sr));
        diksc5 = dt*(iksc4[n]*sa + iksc9[n]*sd - iksc5[n]*(4*sb + 4*sr));
        diksc6 = dt*(iksc2[n]*sr + iksc7[n]*sb - iksc6[n]*(sd + 3*sa));
        diksc7 = dt*(iksc6[n]*3*sa + iksc8[n]*2*sb + iksc3[n]*2*sr  + iksc10[n]*2*sd- iksc7[n]*(sb + 2*sa + sd + sr));
        diksc8 = dt*(iksc7[n]*2*sa + iksc9[n]*3*sb + iksc4[n]*3*sr + iksc11[n]*2*sd - iksc8[n]*(2*sb + sa + sd + 2*sr));
        diksc9 = dt*(iksc8[n]*sa + iksc5[n]*4*sr + iksc12[n]*2*sd - iksc9[n]*(3*sb + sd + 3*sr));
        diksc10 = dt*(iksc11[n]*sb + iksc7[n]*sr - iksc10[n]*(2*sa + 2*sd));
        diksc11 = dt*(iksc10[n]*2*sa + iksc12[n]*2*sb + iksc8[n]*2*sr + iksc13[n]*3*sd - iksc11[n]*(sb + sa + 2*sd + sr));
        diksc12 = dt*(iksc11[n]*sa + iksc9[n]*3*sr + iksc14[n]*3*sd - iksc12[n]*(2*sb + 2*sr + 2*sd));
        diksc13 = dt*(iksc14[n]*sb + iksc11[n]*sr - iksc13[n]*(sa + 3*sd));
        diksc14 = dt*(iksc13[n]*sa + iksc12[n]*2*sr + iksc15[n]*4*sd - iksc14[n]*(sb + 3*sd + sr));
        diksc15 = dt*(iksc14[n]*sr + ikso1[n]*ikseta - iksc15[n]*(4*sd + ikstheta));
        dikso1 = dt*(iksc15[n]*ikstheta + ikso2[n]*iksomega - ikso1[n]*(ikseta + ikspsi));
        dikso2 = dt*(ikso1[n]*ikspsi - ikso2[n]*iksomega);
        
        iksc1[n] = iksc1a + diksc1;
        iksc2[n] = iksc2a + diksc2;
        iksc3[n] = iksc3a + diksc3;
        iksc4[n] = iksc4a + diksc4;
        iksc5[n] = iksc5a + diksc5;
        iksc6[n] = iksc6a + diksc6;
        iksc7[n] = iksc7a + diksc7;
        iksc8[n] = iksc8a + diksc8;
        iksc9[n] = iksc9a + diksc9;
        iksc10[n] = iksc10a + diksc10;
        iksc11[n] = iksc11a + diksc11;
        iksc12[n] = iksc12a + diksc12;
        iksc13[n] = iksc13a + diksc13;
        iksc14[n] = iksc14a + diksc14;
        iksc15[n] = iksc15a + diksc15;
        ikso1[n] = ikso1a + dikso1;
        
        ikso2[n] = 1 - (iksc1[n] + iksc2[n] + iksc3[n] + iksc4[n] + iksc5[n] + iksc6[n] + iksc7[n] + iksc8[n] + iksc9[n] + iksc10[n] + iksc11[n] + iksc12[n] + iksc13[n] + iksc14[n] + iksc15[n] + ikso1[n]);
        
        oks = ikso1[n]+ ikso2[n];
        //iks = gks * oks * (v[n] - eks);
    if (n <= N+ENDO)
		{
			iks = 0.9*gks * oks * (v[n] - eks);
		}
	else if (n <= N+ENDO+MID)
		{
			iks = 0.3*gks * oks * (v[n] - eks);
	    }
	else 
		{
			iks = gks * oks * (v[n] - eks);
		}
    }
    
    // TIME-INDEPENDENT POTASSIUM CURRENT
    void Vcomp_ik1 (int n)
    {
        k1a	= 1.02/(1+exp(0.2385*(v[n]-ek-59.215)));
        k1b	= (0.49124*exp(0.08032*(v[n]-ek+5.476))+exp(0.06175*(v[n]-ek-594.31)))/(1+exp(-0.5143*(v[n]-ek+4.753)));
        k1ss	= k1a/(k1a+k1b);
        gk1	= 0.5*sqrt(ko/5.4);
        ik1	= gk1*k1ss*(v[n]-ek);
    }
    
    // PLATEAU POTASSIUM CURRENT
    void Vcomp_ikp (int n)
    {
        kp	= 1/(1+exp((7.488-v[n])/5.98));
        ikp	= gkp*kp*(v[n]-ek);
    }
    
    // CALCIUM-DEPENDENT TRANSIENT OUTWARD CHLORIDE CURRENT
    void Vcomp_ito2 (int n)
    {
        ai2f = 0.025/(1+exp((v[n]+58)/5));
        bi2f = 1/(5*(1+exp(-(v[n]+19)/9)));
        i2fss = ai2f/(ai2f+bi2f);
        i2f[n] = i2fss-(i2fss-i2f[n])*exp(-dt/i2ftau);
        Kcaito2 = 1-1/(1+pow((qrel[n]/kmto2),2));
        ibarto2	= pcl*zcl*zcl*((v[n]*frdy*frdy)/(R*temp))*(cli[n]-clo*exp(-zcl*v[n]*frdy/(R*temp)))/(1-exp(-zcl*v[n]*frdy/(R*temp)));
        ito2	= ibarto2*i2f[n]*Kcaito2;
    }
    
    // BACKGROUND CHLORIDE CURRENT
    void Vcomp_iclb (int n)
    {
        iclb	= gbarclb*(v[n]-ecl);
    }
    
    // SODIUM-CALCIUM EXCHANGER CURRENT
    void Vcomp_inaca (int n)
    {
        nacax1 = pow(nasssr[n],3)*cao*exp(0.27*v[n]*frdy/(R*temp)) - pow(nao,3)*casssr[n]*exp((0.27-1)*v[n]*frdy/(R*temp));
        nacax2 = (1+pow((0.000125/casssr[n]),2))*(1+0.32*exp((0.27-1)*v[n]*frdy/(R*temp)));
        nacax3 = 1.3*pow(nasssr[n],3) + pow(87.5,3)*casssr[n] + pow(12.3,3)*cao*(1+casssr[n]/0.0036)+0.0036*pow(nao,3)*(1+pow(nasssr[n],3)/pow(12.3,3))+pow(nasssr[n],3)*cao+pow(nao,3)*casssr[n];
        nacay1 = pow(nai[n],3)*cao*exp(0.27*v[n]*frdy/(R*temp)) - pow(nao,3)*cai[n]*exp((0.27-1)*v[n]*frdy/(R*temp));
        nacay2 = (1+pow((0.000125/cai[n]),2))*(1+0.32*exp((0.27-1)*v[n]*frdy/(R*temp)));
        nacay3 = 1.3*pow(nai[n],3) + pow(87.5,3)*cai[n] + pow(12.3,3)*cao*(1+cai[n]/0.0036)+0.0036*pow(nao,3)*(1+pow(nai[n],3)/pow(12.3,3))+pow(nai[n],3)*cao+pow(nao,3)*cai[n];
        
        inacasssr = 0.2*4.5*nacax1/(nacax2*nacax3);
        inacai    = 0.8*4.5*nacay1/(nacay2*nacay3);
        
        //inaca = inacasssr + inacai;
        
    if (n <= N+ENDO)
		{
			inaca = 0.9*inacasssr + 0.9*inacai;
		}
	else if (n <= N+ENDO+MID)
		{
			inaca = 1.3*inacasssr + 1.3*inacai;
	    }
	else 
		{
			inaca = inacasssr + inacai;
		}
        
    }
    
    // SODIUM-POTASSIUM PUMP CURRENT
    void Vcomp_inak (int n)
    {
        inak	= gnakbar*(1/(1+exp(-(v[n]+92)*frdy/(R*temp))))*pow((nai[n]/(nai[n]+kmnai2)),hill)*(ko/(ko+kmko));
    }
    
    // SARCOLEMMAL CALCIUM PUMP CURRENT
    void Vcomp_ipca (int n)
    {
        ipca	= ipcabar*cai[n]/(kmpca+cai[n]);
    }
    
    // BACKGROUND SODIUM CURRENT
    void Vcomp_inab (int n)
    {
        inab    = frdy*pnab*((frdy*v[n])/(R*temp))*(nai[n]*exp((frdy*v[n])/(R*temp)) - nao)/(exp((frdy*v[n])/(R*temp))-1);
    }
    
    // POTASSIUM-CHLORIDE COTRANSPORTER
    void Vcomp_kcl (int n)
    {
        ctkcl	= ctkclbar*(ek-ecl)/((ek-ecl)+87.8251);
    }
    
    // SODIUM-CHLORIDE COTRANSPORTER
    void Vcomp_nacl ()
    {
        ctnacl	= ctnaclbar*pow((ena-ecl),4)/(pow((ena-ecl),4)+pow(87.8251,4));
    }
    
    // STIMULUS CURRENT
    void Vcomp_istim (int n)
    {
        if (n == M+N - 2) // Stimulus inserted in Cell no. 158. (159 is cell at boundary)
        {
            stimtime += dt;
        }
        
        if (t>tstim+1.0 && n == M+N-2)
        {
            stimtime = 0.0;					// start applying stimulus
            stimcount += 1;					// increment the beat counter
            if (stimcount < beats-1)  tstim += BCL;		// time to apply next S1 stimulus
            else if (stimcount == beats-1) tstim += S2;	// time to apply next S1 stimulus
            else tstim = tmax+1;				// ensures only 1 S2 stimulus is applied
            
            if (stimcount < beats) printf ("S1 Stimulus %d applied at time = %.2f ms in cell no %d\n", stimcount+1, t,n);
            else if (stimcount == beats) printf ("S2 Stimulus applied at time = %.2f ms\n", t);
            
            rmbp1[stimcount] = v[n];
        }
		
    }
    
    // TOTAL TRANSMEMBRANE CURRENT
    void Vcomp_itot (int n)
    {
        
        // Don't comment out when inserting stimulus in V-Cell!
        /*
        if (stimtime>=0.0 && stimtime<stimdur && (n==M+N - 2))
		{
            icatot	= ical[n]+icab+ipca-2*inaca;
            iktot	= ikr+iks+ikp+ik1-2*inak+ito1+istim;
            inatot	= 3*inak+ina+3*inaca+inal+inab;
            icltot	= ito2+iclb;
            itot	= icatot+iktot+inatot+icltot;
		}
        
       else
       {
		*/
            icatot	= ical[n]+icab+ipca-2*inaca;
            iktot	= ikr+iks+ikp+ik1-2*inak+ito1;
            inatot	= 3*inak+ina+3*inaca+inal+inab;
            icltot	= ito2+iclb;
            itot	= icatot+iktot+inatot+icltot;
       //}
		
    }
    
    // EQUATIONS DESCRIBING SR FLUXES BEGIN HERE
    
    // SR RELEASE FLUX
    void Vcomp_qrel (int n)
    {
        relbtau = 4.75*(1+(1/(1+pow((0.28/camkactive[n]),10))));
        rela    = 0.1125*relbtau;
        ireltau = relbtau/(1+(0.0123/cajsr[n]));
        irelss  = rela*ical[n]/(1 + pow((1/cajsr[n]),8));
        qrel[n] += dt*(-(irelss+qrel[n])/ireltau);
    }

    void Vcomp_qrelDAD (int n)
    {
        relbtau = 4.75*(1+(1/(1+pow((0.28/camkactive[n]),10))));
        rela    = 0.1125*relbtau;
        ireltau = (1/3)*relbtau/(1+(0.0123/cajsr[n]));
        irelss  = rela*ical[n]/(1 + pow((0.25/cajsr[n]),8));
        qrel[n] += dt*(-(irelss+qrel[n])/ireltau);
    }

    // SR UPTAKE FLUX
    void Vcomp_qup (int n)
    {
        dkmplb		= dkmplbbar*camkactive[n]/(kmcamk+camkactive[n]);
        dqupcamk	= dqupcamkbar*camkactive[n]/(kmcamk+camkactive[n]);
        qup		= (dqupcamk+1)*qupbar*cai[n]/(cai[n]+kmup-dkmplb);
    }
    
    // SR LEAK FLUX
    void Vcomp_qleak (int n)
    {
        qleak	=  0.004375*cansr[n]/15;
    }
    
    // SR TRANSFER FLUX
    void Vcomp_qtr (int n)
    {
        qtr		= (cansr[n]-cajsr[n])/tautr;
    }
    
    // CONCENTRATION CHANGES
    void Vcomp_conc (int n)
    {
        // RESTRICTED SPACE CALCIUM CONCENTRATION
        qdiffss = (casssr[n] - casscal[n])/taudiffss;
        qdiff   = (casssr[n] - cai[n])/taudiff;
        
        bsssr = 1/ (1 + (0.047*0.00087/pow((0.00087+casssr[n]),2)) + (1.124*0.0087/pow((0.0087+casssr[n]),2)));
        bsscal = 1/ (1 + (0.047*0.00087/pow((0.00087+casscal[n]),2)) + (1.124*0.0087/pow((0.0087+casscal[n]),2)));
        
        casssr[n] = dt*bsssr*(2*inacasssr*acap/(zca*frdy*vsssr) + qrel[n]*vjsr/vsssr - qdiff - qdiffss) + casssr[n];
        casscal[n] = dt*bsscal*(-ical[n]*acap/(zca*frdy*vsscal) + qdiffss*vsssr/vsscal) + casscal[n];
        
        // JSR CALCIUM CONCENTRATION
        dcajsr		= dt*(qtr-qrel[n]);
        csqn        = csqnbar*(cajsr[n]/(cajsr[n]+kmcsqn));
        cajsrb      = csqnbar - csqn - cajsr[n] - dcajsr + kmcsqn;
        cajsrc      = kmcsqn*(csqn + cajsr[n] + dcajsr);
        cajsr[n]		= (sqrt(cajsrb*cajsrb+4*cajsrc)-cajsrb)/2;
        
        // NSR CALCIUM CONCENTRATION
        dcansr		= dt*(qup-qleak-qtr*vjsr/vnsr);
        cansr[n]		+= dcansr;
        
        // INTRACELLULAR CALCIUM CONCENTRATION
        dcai		= dt*(-1)*((icab+ipca-2*inacai)*acap/(zca*frdy*vmyo)+(qup-qleak)*vnsr/vmyo-qdiff*vsssr/vmyo);
        trpn        = trpnbar*(cai[n]/(cai[n]+kmtrpn));
        cmdn		= cmdnbar*(cai[n]/(cai[n]+kmcmdn));
        catotal		= trpn+cmdn+dcai+cai[n];
        bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
        cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
        dmyo		= -kmtrpn*kmcmdn*catotal;
        cai[n]		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;
        
        
        // INTRACELLULAR SODIUM, POTASSIUM AND CHLORIDE CONCENTRATIONS
        qdiffna = (nasssr[n] - nai[n])/taudiff;
        dnai		= dt*((-(3*inacai + 3*inak + ina + inal + inab)*acap)/(vmyo*zna*frdy)+ctnacl+qdiffna*vsssr/vmyo);
        nai[n]		+= dnai;
        //nai             = 10.0;
        dnasssr		= dt*((-3*inacasssr*acap)/(vsssr*zna*frdy)-qdiffna);
        nasssr[n]	+= dnasssr;
        //nasssr          = 10.0;
        dki		= dt*((-iktot*acap)/(vmyo*zk*frdy)+ctkcl);
        ki[n]		+= dki;
        
        qdiffcl     = (clss[n]-cli[n])/taudiff;
        dcli		= dt*((-iclb*acap)/(vmyo*zcl*frdy)+ctnacl+ctkcl+qdiffcl*vsssr/vmyo);     
        cli[n]		+= dcli;
        dclss		= dt*((-ito2*acap)/(vsssr*zcl*frdy)-qdiffcl);     
        clss[n]		+= dclss;
        
        
        // CALCIUM/CALMODULIN-DEPENDENT PROTEIN KINASE (CaMKII)
        camkbound	= camk0*(1-camktrap[n])*1/(1+(kmcam/casssr[n]));
        camktrap[n]	= dt*(alphacamk*camkbound*(camkbound+camktrap[n])-betacamk*camktrap[n]) + camktrap[n];
        camkactive[n]	= camkbound+camktrap[n];
    }

void printtofile ()
{
	count1   += 1;	
	
	//PRINT LAST 5 BEATS
	if (count1>=100)			
	{
	count1   = 0;
    printf("Run = %d Loading ... %f \n",run, 100*t/tmax);
    
      if (t>=(BCL*(beats-1))) {
		/* Determines pseudo ECG */
		sum = 0.0;
		for (n=1; n<N; n++)
			{  
			R1 = 1*(200+160-n+1);
			R2 = 1*(200+160-n-1);
			//sum -= 2*(v[n+1]-v[n-1])/(R1*R2); 
			sum -= (v[n+1]-v[n-1])*(1.0/R2-1.0/R1);
			}
		pECG = 100000000*0.00175*0.00175*sum/(64*dx*dx);
		sum = 0.0;
		for (n=N; n<N+M-1; n++)
			{  
			R1 = 1*(200+160-n+1);
			R2 = 1*(200+160-n-1);
			//sum -= 2*(v[n+1]-v[n-1])/(R1*R2); 
			sum -= (v[n+1]-v[n-1])*(1.0/R2-1.0/R1);
			}
		pECG += 100000000*0.0011*0.0011*sum/(64*dx*dx);
		
		fprintf(ecg,"%f \t", pECG);
    //print last 5 beats
    //if (t+2>=(BCL*(beats-4)))
    //{
        count += 1;

        fprintf(ap,"%f",run+1*t/(BCL*(beats-1)+BCLMax+500));
        for (n = 0; n < M+N; n++)
        {
            fprintf(ap," %f", v[n]);
        }
        fprintf(ap,"\n");
    //}
		}
	}
	
	if (t>=tmax-dtmax) {
		  fprintf(ecg,"\n", pECG);
	}
}

void delaytoscreen()
{
    if (v1 == 0 && v[N-1] > -41)
    {
        v1 = t;
    }
    if (v2 == 0 && v[N] > -41)
    {
        v2 = t;
        printf("delay: %f    ER: %f \n",v2-v1,ER);
    }
}

void printSSP (char *infile)
{
    FILE *gp = fopen(infile,"w");
    for (ii = 0; ii < N; ii++)
    {
        if (ii == 0)
        {
            updateS(1);
        }
        else
        {
            updateS(ii);
        }
        
        for(jj = 0; jj < 36; jj++)
        {
            fprintf(gp,"%f ",S[jj]);
        }
        fprintf(gp,"\n");
    }
    
    fclose(gp);
}

void printSSV(char *infile)
{
    FILE *gv = fopen(infile,"w");
    for (ii = N; ii < N+M; ii++)
    {
        if (ii == N+M-1)
        {
            updateSv(N+M-2);
        }
        else
        {
            updateSv(ii);
        }
        for (jj = 0; jj < 48; jj++)
        {
            fprintf(gv,"%f ",Sv[jj]);
        }
        fprintf(gv,"\n");
        
    }
    fclose(gv);
}

void readSSP(char *infile)
{
    FILE *fp = fopen(infile,"r");
    char line[1024];
    const char split[2] = " ";
    double value;
    nn = 0;
    
    while (fgets(line,1024,fp))
    {
        
        char *token;
        token = strtok(line,split);
        n = 0;
        while (token != NULL)
        {
            value = strtod(token,NULL);
            token = strtok(NULL,split);
            S[n] = value;
            n++;
        }
        
        dvdt[nn] = 0;
        
        v [nn] =S[0];
        m [nn] =S[1];
        h [nn] =S[2];
        j [nn] =S[3];
        d [nn] =S[4];
        f [nn] =S[5];
        f2 [nn] =S[6];
        fca [nn] =S[7];
        fca2 [nn] =S[8];
        xs1 [nn] =S[9];
        xs2 [nn] =S[10];
        xr [nn] =S[11];
        a [nn] =S[12];
        i [nn] =S[13];
        i2 [nn] =S[14];
        ml [nn] =S[15];
        ml3 [nn] =S[16];
        hl [nn] =S[17];
        hl3 [nn] =S[18];
        jl [nn] =S[19];
        jl3 [nn] =S[20];
        casss [nn] =S[21];
        cajsr [nn] =S[22];
        cacsr [nn] =S[23];
        cansr [nn] =S[24];
        cassl [nn] =S[25];
        nai [nn] =S[26];
        nassl [nn] =S[27];
        nasss [nn] =S[28];
        ki [nn] =S[29];
        cai [nn] =S[30];
        b [nn] =S[31];
        g [nn] =S[32];
        u [nn] =S[33];
        y [nn] =S[34];
        camktrap [nn] =S[35];
        
        
        camkactive[nn] = 0;
        ical[nn] = 0;
        qrel1[nn] = 0;
        qrel2[nn] = 0;
        qup2[nn] = 0;
        
        nn++;
        
        
        
    }
    fclose(fp);
}

void readSSV(char *infile)
{
    FILE *fp = fopen(infile,"r");
    char line[1024];
    const char split[2] = " ";
    double value;
    nn = 80;
    while (fgets(line,1024,fp))
    {
        
        char *token;
        token = strtok(line,split);
        n = 0;
        while (token != NULL)
        {
            value = strtod(token,NULL);
            token = strtok(NULL,split);
            Sv[n] = value;
            n++;
        }
        
        dvdt[nn] = 0;
        
        v[nn] = Sv[0];
        cai[nn] = Sv[1];
        casssr[nn] = Sv[2];
        casscal[nn] = Sv[3];
        cajsr[nn] = Sv[4];
        cansr[nn] = Sv[5];
        camktrap[nn] = Sv[6];
        nai[nn] = Sv[7];
        nasssr[nn] = Sv[8];
        ki[nn] = Sv[9];
        cli[nn] = Sv[10];
        clss[nn] = Sv[11];
        qrel[nn] = Sv[12];
        m[nn] = Sv[13];
        h[nn] = Sv[14];
        j[nn] = Sv[15];
        ml[nn] = Sv[16];
        hl[nn] = Sv[17];
        a[nn] = Sv[18];
        i[nn] = Sv[19];
        i2[nn] = Sv[20];
        i2f[nn] = Sv[21];
        xr[nn] = Sv[22];
        icalc[nn] = Sv[23];
        icalo[nn] = Sv[24];
        icalci[nn] = Sv[25];
        icaloi[nn] = Sv[26];
        icalcstar[nn] = Sv[27];
        icalostar[nn] = Sv[28];
        icalcistar[nn] = Sv[29];
        icaloistar[nn] = Sv[30];
        iksc1[nn] = Sv[31];
        iksc2[nn] = Sv[32];
        iksc3[nn] = Sv[33];
        iksc4[nn] = Sv[34];
        iksc5[nn] = Sv[35];
        iksc6[nn] = Sv[36];
        iksc7[nn] = Sv[37];
        iksc8[nn] = Sv[38];
        iksc9[nn] = Sv[39];
        iksc10[nn] = Sv[40];
        iksc11[nn] = Sv[41];
        iksc12[nn] = Sv[42];
        iksc13[nn] = Sv[43];
        iksc14[nn] = Sv[44];
        iksc15[nn] = Sv[45];
        ikso1[nn] = Sv[46];
        ikso2[nn] = Sv[47];
        
        
        camkactive[nn] = 0;
        ical[nn] = 0;
        qrel1[nn] = 0;
        qrel2[nn] = 0;
        qup2[nn] = 0;
        
        nn++;
    }
    fclose(fp);
    
}
void pGeometry()
{
    rcg = 1.54;
    radius = 0.00175;
    length = 0.0164;
    vcell = 1000*pi*radius*radius*length;
    ageo = 2*pi*radius*radius + 2*pi*radius*length;
    acap = rcg*ageo;
    vmyo = vcell * 0.60;
    vnsr = vcell * 0.04;
    vjsr	= vcell * 0.002;
    vcsr	= vcell * 0.008;
    vsss	= vcell * 0.02;
    vssl    = vcell * 0.15;
}

void vGeometry()
{
    radius = 0.0011;		// radius of the cell (cm)
    length = 0.01;		// length of the cell (cm)
    rcg = 2;
    vcell	= 38e-6;
    ageo	= 0.767e-4;
    acap	= 1.534e-4;
    vmyo	= 25.84e-6;
    vmito   = 9.88e-6;
    vnsr	= 2.098e-6;
    vjsr	= 0.182e-6;
    vsssr	= 0.76e-6;
    vsscal  = 0.076e-6;
    vss     = 0.836e-6;
}

void updateS(int ii)
{
        S[0] = v[ii];
        S[1] = m[ii];
        S[2] = h[ii];
        S[3] = j[ii];
        S[4] = d[ii];
        S[5] = f[ii];
        S[6] = f2[ii];
        S[7] = fca[ii];
        S[8] = fca2[ii];
        S[9] = xs1[ii];
        S[10] = xs2[ii];
        S[11] = xr[ii];
        S[12] = a[ii];
        S[13] = i[ii];
        S[14] = i2[ii];
        S[15] = ml[ii];
        S[16] = ml3[ii];
        S[17] = hl[ii];
        S[18] = hl3[ii];
        S[19] = jl[ii];
        S[20] = jl3[ii];
        S[21] = casss[ii];
        S[22] = cajsr[ii];
        S[23] = cacsr[ii];
        S[24] = cansr[ii];
        S[25] = cassl[ii];
        S[26] = nai[ii];
        S[27] = nassl[ii];
        S[28] = nasss[ii];
        S[29] = ki[ii];
        S[30] = cai[ii];
        S[31] = b[ii];
        S[32] = g[ii];
        S[33] = u[ii];
        S[34] = y[ii];
        S[35] = camktrap[ii];
}

void updateSv(int ii)
{
    Sv[0] = v[ii];
    Sv[1] = cai[ii];
    Sv[2] = casssr[ii];
    Sv[3] = casscal[ii];
    Sv[4] = cajsr[ii];
    Sv[5] = cansr[ii];
    Sv[6] = camktrap[ii];
    Sv[7] = nai[ii];
    Sv[8] = nasssr[ii];
    Sv[9] = ki[ii];
    Sv[10] = cli[ii];
    Sv[11] = clss[ii];
    Sv[12] = qrel[ii];
    Sv[13] = m[ii];
    Sv[14] = h[ii];
    Sv[15] = j[ii];
    Sv[16] = ml[ii];
    Sv[17] = hl[ii];
    Sv[18] = a[ii];
    Sv[19] = i[ii];
    Sv[20] = i2[ii];
    Sv[21] = i2f[ii];
    Sv[22] = xr[ii];
    Sv[23] = icalc[ii];
    Sv[24] = icalo[ii];
    Sv[25] = icalci[ii];
    Sv[26] = icaloi[ii];
    Sv[27] = icalcstar[ii];
    Sv[28] = icalostar[ii];
    Sv[29] = icalcistar[ii];
    Sv[30] = icaloistar[ii];
    Sv[31] = iksc1[ii];
    Sv[32] = iksc2[ii];
    Sv[33] = iksc3[ii];
    Sv[34] = iksc4[ii];
    Sv[35] = iksc5[ii];
    Sv[36] = iksc6[ii];
    Sv[37] = iksc7[ii];
    Sv[38] = iksc8[ii];
    Sv[39] = iksc9[ii];
    Sv[40] = iksc10[ii];
    Sv[41] = iksc11[ii];
    Sv[42] = iksc12[ii];
    Sv[43] = iksc13[ii];
    Sv[44] = iksc14[ii];
    Sv[45] = iksc15[ii];
    Sv[46] = ikso1[ii];
    Sv[47] = ikso2[ii];
}

void record_variables()
{
    if (stimcount>=0 && v[1]>vmax1[stimcount])
    {vmax1[stimcount] = v[1];}
    
    if (stimcount>=0 && cai[1]>caimax1[stimcount])
    {caimax1[stimcount] = cai[1];}
    
    if (stimcount>=0 && nai[1]>naimax1[stimcount])
    {naimax1[stimcount] = nai[1];}
    
    if (stimcount>=0 && dvdt[1]>dvdtmax1[stimcount])
    {
        dvdtmax1[stimcount] = dvdt[1];
        toneapd1[stimcount] = t;
    }
    if (stimcount>=0 && v[1]>=(vmax1[stimcount]-0.90*(vmax1[stimcount]-rmbp1[stimcount])))
        ttwoapd1[stimcount] = t;
    if (stimcount>=0 && v[1]>=(vmax1[stimcount]-0.98*(vmax1[stimcount]-rmbp1[stimcount])))
    {
        trep1[stimcount] = t;
    }
    
    
    
    if (stimcount>=0 && v[79] > vmax79[stimcount])
    {vmax79[stimcount] = v[79];}
    
    if (stimcount>=0 && cai[79] > caimax79[stimcount])
    {caimax79[stimcount] = cai[79];}
    
    if (stimcount>=0 && nai[79] > naimax79[stimcount])
    {naimax79[stimcount] = nai[79];}
    
    if (stimcount>=0 && dvdt[79] > dvdtmax79[stimcount])
    {
        dvdtmax79[stimcount] = dvdt[79];
        toneapd79[stimcount] = t;
    }
    if (stimcount>=0 && v[79] >=(vmax79[stimcount]-0.90*(vmax79[stimcount]-rmbp79[stimcount])))
        ttwoapd79[stimcount] = t;
    if (stimcount>=0 && v[79] >=(vmax79[stimcount]-0.98*(vmax79[stimcount]-rmbp79[stimcount])))
    {
        trep79[stimcount] = t;
    }
    
    
    
    if (stimcount>=0 && v[80]>vmax80[stimcount])
    {vmax80[stimcount] = v[80];}
    
    if (stimcount>=0 && cai[80]>caimax80[stimcount])
    {caimax80[stimcount] = cai[80];}
    
    if (stimcount>=0 && nai[80]>naimax80[stimcount])
    {naimax80[stimcount] = nai[80];}
    
    if (stimcount>=0 && dvdt[80]>dvdtmax80[stimcount])
    {
        dvdtmax80[stimcount] = dvdt[80];
        toneapd80[stimcount] = t;
    }
    if (stimcount>=0 && v[80]>=(vmax80[stimcount]-0.90*(vmax80[stimcount]-rmbp80[stimcount])))
    {
        ttwoapd80[stimcount] = t;
    }
    if (stimcount>=0 && v[80]>=(vmax80[stimcount]-0.98*(vmax80[stimcount]-rmbp80[stimcount])))
    {
        trep80[stimcount] = t;
    }
    
    
    
    if (stimcount>=0 && v[159]>vmax159[stimcount])
    {vmax159[stimcount] = v[159];}
    
    if (stimcount>=0 && cai[159]>caimax159[stimcount])
    {caimax159[stimcount] = cai[159];}
    
    if (stimcount>=0 && nai[159]>naimax159[stimcount])
    {naimax159[stimcount] = nai[159];}
    
    if (stimcount>=0 && dvdt[159]>dvdtmax159[stimcount])
    {
        dvdtmax159[stimcount] = dvdt[159];
        toneapd159[stimcount] = t;
    }
    if (stimcount>=0 && v[159]>=(vmax159[stimcount]-0.90*(vmax159[stimcount]-rmbp159[stimcount])))
        ttwoapd159[stimcount] = t;
    if (stimcount>=0 && v[159]>=(vmax159[stimcount]-0.98*(vmax159[stimcount]-rmbp159[stimcount])))
    {
        trep159[stimcount] = t;
    }
}

void pConstants()
{
    fca_dtaucamkbar = 10.0;
    camk0 = 0.05;
    alphacamk = 0.05;
    betacamk = 0.00068;
    kmcam = 0.0015;
    kmcamk = 0.15;
    trpnbar1 = 3.5e-3;
    cmdnbar1 = 1.25e-2;
    csqnbar1 = 1.2;
    trpnbar = 3.15e-2;
    kmtrpn = 0.5e-3;
    cmdnbar = 0.1125;
    kmcmdn = 2.38e-3;
    csqnbar = 2.88;
    kmcsqn = 0.8;
    bsrbar = 0.019975;
    kmbsr = 0.00087;
    bslbar = 0.4777;
    kmbsl = 0.0087;
    kmup   = 0.00028;
    nsrbar = 15.0;
    dkmplbbar = 0.00017;
    dqupcamkbar = 0.75;
    IP3 = 0.0001;
    k1 = 150000;
    k1a = 16.5;
    k0 = 96000;
    k0a = 9.6;
    k2 = 1800;
    k2a = 0.21;
    tauip3r = 3.7;
    tautr1 = 120;
    tautr2 = 120;
    gaptau = 12;
    sstau = 0.2;
    ipcabar = 0.0115;
    kmpca = 0.0005;
    ibarnak = 1.1004;
    inacamax = 2.52;
    kmcaact = 0.000125;
    kmnai1 = 12.3;
    kmnao = 87.5;
    kmcai = 0.0036;
    kmcao = 1.3;
    nu = 0.35;
    ksat = 0.27;
    gnab = 0.0025;
    pcab = 3.99e-8;
    pnab = 0.64e-8;
    prnak = 0.014;
    gtos = 0.1414;
    gtof = 0.042;
    gcat = 0.07875;
    powtau = 10;
    pca = 1.9926e-4;
    gnal2 = 0.052;
    gnal3 = 0.018;
    gna = 18;
}

void vConstants()
{
    gna = 9.075;
    gnal = 0.0065;
    hltau = 600;
    pca = 1.5552e-4;
    pcab = 1.995e-7;
    gbarto1 = 0.4975;
    prnak = 0.01833;
    gkp = 0.00276;
    kmto2 = 0.4;
    pcl = 9e-7;
    i2ftau = 8;
    pnab = 0.32e-8;
    gbarclb = 2.25e-4;
    inacamax = 4.5;		// maximal inaca (uA/uF)
    kmcaact = 0.000125;	// half-saturation concentration for cai activation (mM)
    kmnai1 = 12.3;		// half-saturation concentration for nai (mM)
    kmnao = 87.5;		// half-saturation concentration for nao  (mM)
    kmcai = 0.0036;		// half-saturation concentration for cai (mM)
    kmcao = 1.3;		// half-saturation concentration for cao (mM)
    nu = 0.27;			// position of energy barrier for inaca
    ksat = 0.32;
    gnakbar = 1.4;		// maximum inak (uA/uF)
    hill = 3;			// hill coefficient
    kmnai2 = 2.6;		// half-saturation concentration for nai (mM)
    kmko = 1.5;
    ipcabar = 0.0575;		// maximum ipca (uA/uF)
    kmpca = 0.0005;
    ctkclbar = 1.77e-5;
    taudiffss=2;
    taudiff=0.2;
    sstau=0.2;
    ctnaclbar = 2.46108e-5;
    qleakbar = 0.004375;	// maximum leak from nsr to myoplasm (mM/ms)
    nsrbar = 15.0;
    dqupcamkbar = 0.75;
    dkmplbbar = 0.00017;
    kmup = 0.00092;		// half-saturation concentration of qup (mmol)
    qupbar = 0.004375;
    tautr = 25.0;
    bsrbar = 0.047;		// maximum calcium binding by anionic binding sites in the subspace (mM)
    kmbsr = 0.00087;		// half saturation coefficient of anionic binding sites in the subspace (mM)
    bslbar = 1.124;		// maximum calcium binding by sarcolemmal binding sites in the subspace (mM)
    kmbsl = 0.0087;
    csqnbar = 10.0;		// maximum calcium buffered by calsequestrin (mM)
    kmcsqn = 0.8;
    cmdnbar = 50e-3;		// maximum calcium buffered by calmodulin (mM)
    kmcmdn = 2.38e-3;
    trpnbar = 70e-3;		// maximum calcium buffered by troponin (mM)
    kmtrpn = 0.5e-3;
    camk0 = 0.05;		// fraction of active CaMKII binding sites at equilibrium
    alphacamk = 0.05;		// phosphorylation rate of CaMKII (1/ms)
    betacamk = 0.00068;	// dephosphorylation rate of CaMKII (1/ms)
    kmcam = 0.0015;		// half saturation coefficient of CaM (mM)
    kmcamk = 0.15;
}
