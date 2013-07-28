/*
 *  main.c
 *  Assignment Four
 *
 *  Created by Derek Williams on 10-01-30.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *	This code numerically solves hyperbolic PDEs of the form: 
 *  
 *		Dt[u] + a Dx[u] + b Dy[u] + b Dxx[u] + u = F(t, x)
 *
 *		where Dt[], Dx[], Dy[], and Dxx[] are the differential operators for t, x, and y
 *
 *  The following explicit finite difference schemes are implemented:
 *		
 *		Forward-Time Forward-Space
 *		Forward-Time Backward-Space
 *		Forward-Time Central-Space
 *		Lax Fredrichs
 *		Leapfrog
 *		Equilibrium
 *      Lax Wendroff
 *
 *  The following implicit finite difference schemes are implemented:
 *
 *		Crank Nicolson
 *
 *  The solutions are animated in a window. Data is saved to text files when 's' is pressed.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <OpenGL/OpenGL.h>
#include <GLUT/GLUT.h>
#include <GLUI/GLUI.h>
#include <assert.h>

#include "main.h"

int lastFrameTime = 0;
int blowup = 0;
double dx = 0.1;
double dy = 0.1;
double dt = 0.0;
int Xm = 0;
int Ym = 0;
int Tn = 0;
float mytime = 0;
double **data = 0;
double *pdata;
double *qdata;
double *maxv = 0;
double maxmaxv = 0;
double *error = 0;
int maxVmN = 0;
int mypause = 1;			// start/pause animation
int dir = 1;				// direction of wave propagation
int udir = 1;				// direction of wave propagation
int zeroRightBoundary = 0;
int saving = 0;
int windMain = 0;
int drawSolution = 1;
int uTerm = 0;
int linespacing = 5;
int repeat = 1;
int first = 1;

int win_width;
int win_height;

enum cf_type {
	CFG_INT, 
	CFG_FLOAT, 
	CFG_STR
};

struct config {
	int	id;
	const char *name;
	cf_type type;
	int ival;
	double fval;
	char sval[BUFLEN];
};

config configs[] = {
	{CF_SCHEME,		"scheme",	CFG_INT,	DEFAULT_SCHEME, 0,				0},
	{CF_LAMBDA,		"lambda",	CFG_FLOAT,	0,				DEFAULT_LAMBDA, 0},
	{CF_MU,			"mu",		CFG_FLOAT,	0,				DEFAULT_MU,		0},
	{CF_TOP,		"top",		CFG_INT,	BOUNDS_TOP,		0,				0},
	{CF_BOTTOM,		"bottom",	CFG_INT,	BOUNDS_BOTTOM,	0,				0},
	{CF_LEFT,		"left",		CFG_INT,	BOUNDS_LEFT,	0,				0},
	{CF_RIGHT,		"right",	CFG_INT,	BOUNDS_RIGHT,	0,				0},
	{CF_BOUNDARY,	"boundary",	CFG_INT,	0,				0,				0},
	{CF_TIME,		"time",		CFG_FLOAT,	0,				TIME_END,		0},
	{CF_THREEDIM,	"threedim",	CFG_INT,	DEFAULT_THREEDIM, 0,			0},
	{CF_WIRE,		"wire",		CFG_FLOAT,	0,				0,				0},
	{CF_ACONST,		"aconst",	CFG_FLOAT,	1,				0,				0},
	{CF_BCONST,		"bconst",	CFG_FLOAT,	0,				0,				0},	
	{CF_A_FUNC,		"afunc",	CFG_INT,	0,				1,				0},	
	{CF_B_FUNC,		"bfunc",	CFG_INT,	0,				0,				0},	
	{CF_INIT_FUNC,	"initfunc",	CFG_INT,	0,				0,				0},	
	{CF_BOUNDS_FUNC,"boundsfunc",CFG_INT,	DEFAULT_BOUND,	0,				0},		
	{CF_H,			"h",		CFG_FLOAT,	0,				0.1,			0},		
	{CF_PBC,		"pbc",		CFG_INT,	0,				0,				0},		
	{CF_RHO,		"rho",		CFG_FLOAT,	0,				10,				0},		
	{CF_THETA,		"theta",	CFG_FLOAT,	0,				0,				0},		
	{CF_PHI,		"phi",		CFG_FLOAT,	0,				0,				0},			
	{CF_R,			"r",		CFG_FLOAT,	0,				0.1,			0},		
	{CF_SIGMA,		"sigma",	CFG_FLOAT,	0,				0.5,			0},		
	{CF_K,			"k",		CFG_FLOAT,	0,				100,			0},		
	{CF_DERIV,		"deriv",	CFG_INT,	0,				0,				0},		
	{-1,			NULL,		CFG_STR,	0,				0,				NULL},	
};

double ***data3d;
double **tmp3d;

double head;
double hmargin;
double vmargin;

//double scale = 1.0;

int xi = 8;

GLUI *glui;
GLUI_Spinner *spin;
GLUI_StaticText *pdetext;
int hListId;
float hTime;
float hLambda = 0;
float hMu = 0;
int hIC = 0;
float hAval = 0;
float hBval = 0;
float hSval = 0;
float hRval = 0;
float hKval = 0;
char configfile[BUFLEN] = "three";

double mouseinfo[4];

const char *schemeNames[NUM_SCHEMES] = {
	FORWARD_SPACE_NAME,
	BACKWARD_SPACE_NAME,
	CENTRAL_SPACE_NAME,
	LAX_FREDRICHS_NAME,
	LEAPFROG_NAME,
	EQUILIBRIUM_NAME,
	LAX_WENDROFF_NAME,
	CRANK_NICOLSON_NAME
};

const char *aFuncNames[2] = {
	A_CONST_NAME,
	A_SPACE_NAME
};

const char *boundNames[NUM_BOUND] = {
	BOUND_ONE_NAME,
	BOUND_TWO_NAME,
	BOUND_THREE_NAME,
	BOUND_FOUR_NAME,
	BOUND_FIVE_NAME
};

const char *initNames[NUM_INIT] = {
	INIT_ONE_NAME,
	INIT_TWO_NAME,
	INIT_THREE_NAME,
	INIT_FOUR_NAME,
	INIT_FIVE_NAME,
	INIT_SIX_NAME,
	INIT_SEVEN_NAME,
	INIT_EIGHT_NAME
};

const char *derivNames[NUM_DERIVS] = {
	DERIV_VANILLA_CALL_NAME,
	DERIV_VANILLA_PUT_NAME
};

int activeA = 0;

int initialize (void);
double applyScheme(int scheme, int n, int m);
void save (void);
void fdm (void);
void display(void);
void reshape(int width, int height);
void idle(void);
void buildMenu (void);
void Menu (int value);
void Mouse (int button, int state, int x, int y);
void Keyboard (unsigned char key, int x, int y);
int main(int argc, char** argv);
double worldx(int x);
double worldy(int y);
void drawString(char *s, void *font, int align, double x, double y);
void reshapeWorld();
void calcError(int n, int m, double val);

config *getConfigByName(char *name)
{
	int i;
	config *cp = configs; 
	
	while (cp->name != NULL)
	{
		if (strcmp(name, cp->name) == 0)
			return cp;
		cp++;
	}
	
	return NULL;
}

int writeConfig(const char *file)
{
	int i;
	config *cp = configs; 
	
	FILE *o;
	o = fopen(file, "r");
	
	if (o)
	{
		fprintf(stderr, "Config file already exists: %s\nOverwrite Y/N? ", file);

		char ch = fgetc(stdin);
		if (ch != 'Y')
			return 0;
	}

	o = fopen(file, "w");

	if (!o)
	{
		fprintf(stderr, "Error writing config file: %s\n", file);
		return 0;
	}
	
	printf("Writing config file: %s\n", file);
	
	while (cp->name != NULL)
	{
		switch (cp->type) {
			case CFG_INT:	
				fprintf(o, "%s\t%d\n", cp->name, cp->ival);
				break;
				
			case CFG_FLOAT:
				fprintf(o, "%s\t%f\n", cp->name, cp->fval);
				break;
				
			case CFG_STR:
				fprintf(o, "%s\t%s\n", cp->name, cp->sval);
				break;
				
			default: 
				break;
		}
		cp++;
	}
	
	fflush(o);
	fclose(o);
	
	return 1;
}

int loadConfig(const char *file)
{
	int i;
	char buf[BUFLEN], *pbuf, sbuf[BUFLEN];
	config *pconf;
	
	FILE *o;
	o = fopen(file, "r");
	
	if (!o)
	{
		fprintf(stderr, "Error opening config file: %s\n", file);
		return 0;
	}
	
	printf("Reading config file: %s\n", file);
	
	while (fgets(buf, BUFLEN, o))
	{
		pbuf = buf;
		while (*pbuf != '\t') 
			pbuf++;
		*pbuf = '\0';
		pbuf++;
		
		pconf = getConfigByName(buf);
		
		if (pconf == NULL)
			continue;
		
		switch (pconf->type) {
			case CFG_INT:	
				sscanf(pbuf, "%d", &pconf->ival);
				break;
				
			case CFG_FLOAT:
				sscanf(pbuf, "%lf", &pconf->fval);				
				break;
				
			case CFG_STR:
				strcpy(pconf->sval, pbuf);
				break;
				
			default:
				break;
		}
	}
	
	fclose(o);
		
	return 1;
}	
	

double FZero(double t, double x)
{
	return 0;
}

double FOne(double t, double x)
{
	return sin(M_PI*(x-t));
}

double FAOne(double t, double x)
{
	return t*sin(M_PI*(x-t));
}

double (*FActual)(double t, double x) = &FZero;
double (*FFunc)(double t, double x) = &FZero;

double BZero(double t, double *x)
{
	return 0;
}

double BCOne(double t, double *x)
{
	return -(1+t)*sin(M_PI*t);
}

double (*boundCondFunc)(double t, double *x) = DEFAULT_BOUNDARY;

/* functions to specify a(x, t) parameter */

double aConst(double t, double x)
{
	return configs[CF_ACONST].fval;
}

double aSpace(double t, double x)
{
	return 1 + 0.25 * (3 - x) * (1 + x);
}

double (*aFunc)(double t, double x) = DEFAULT_A;

double bConst(double t, double x)
{
	return configs[CF_BCONST].fval;
}

double (*bFunc)(double t, double x) = DEFAULT_B;

/* functions to specify initial conditions */

double initFunc(double t, double *x)
{
	double xabs = fabs(*x);

	switch (configs[CF_INIT_FUNC].ival) 
	{
		case INIT_ONE:
			if ( fabs(*x) <= 0.5)
				return pow(cos(M_PI * *x), 2);
			return 0;

		case INIT_TWO:
			if ( xabs <= 1.0 )
				return 1.0 - xabs;
			return 0;

		case INIT_THREE:
			return sin(2*M_PI**x);

		case INIT_FOUR:
			return sin(M_PI**x);

		case INIT_FIVE:
			if (fabs(*x) <= 1)
				return cos(xi * M_PI * *x)*pow(cos(0.5*M_PI**x), 2);
			else 
				return 0;

		case INIT_SIX:
			if ( xabs < 0.5 )
				return 1.0 - xabs;
			else if (xabs == 0.5)
				return 0.25;
			return 0;

		case INIT_SEVEN:
			return sin(1.2*(x[0] - x[1]))*cosh(x[0] + 2*x[1]);

		case INIT_EIGHT:
			return exp(-(x[0]*x[0] + x[1]*x[1]));
	}
}	
	
double n01(double t)
{
	return exp(-t*t/2.0)/sqrt(2.0*M_PI);
}

double N01_b[5] = {0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429};

double N01(double z)
{
	if (z > 6.0)
		return 1.0;
	else if (z < -6.0)
		return 0.0;
	
	double a = fabs(z);
	double t = 1.0/(1.0 + a * 0.2316419);
	double d = 0.3989423 * exp(-z*z/2.0);
	double n = N01_b[4];
	
	for (int i = 3; i >= 0; i--) 
		n = (n + N01_b[i])*t;
	n = 1.0 - d*n;
	if (z < 0.0)
		n = 1.0 - n;
	
	return n;
}

double solution (double t, double s)
{
	double res;
	double K = configs[CF_K].fval;
	double sigma = configs[CF_SIGMA].fval;
	double r = configs[CF_R].fval;
	double T = Tn*dt; // configs[CF_MU].fval;
	double Tsqrt;
	double d1;
	double d2;

	t = T - t;

	Tsqrt = sqrt(T - t);
	
	d1 = (log(s/K) + (r + sigma*sigma/2)*Tsqrt)/(sigma*Tsqrt);
	d2 = d1 - sigma*Tsqrt;
		
	switch (configs[CF_DERIV].ival) {
		case DERIV_VANILLA_CALL:
//			res = N01(s);
			res = s * N01(d1) - K*exp(-r*(T-t)) * N01(d2);
			break;

		case DERIV_VANILLA_PUT:
			res = K*exp(-r*(T-t)) * N01(-d2) - s * N01(-d1);
			break;
	}

	return res;
}

double solution3d (double t, double *x)
{
	assert(t <= configs[CF_TIME].fval);
		
	switch (configs[CF_INIT_FUNC].ival) 
	{
		case INIT_SEVEN:
			return exp(1.62*t) * sin(1.2*(x[0] - x[1])) * cosh(x[0] + 2*x[1]);
	
		case INIT_EIGHT:
			return exp(-(x[0]*x[0] + x[1]*x[1]));
	}
}

#define MALLOC(A, T, S)						\
	A = (T)malloc(S);						\
	if (A == NULL)							\
		perror("Error allocating memory.\n");

/* free allocated memory */

#define MFREE(X)	\
	if (X)			\
		free(X);	\
	X = NULL;

int cleanup (void)
{
	int i, j;
	
	MFREE(maxv);
	MFREE(error);	
	MFREE(pdata);
	MFREE(qdata);
	
	for (i = 0; i < Tn; i++) {
		if (data)
			MFREE(data[i]);
		
		if (data3d)
		{
			for (j = 0; j < Xm; j++)
				if (data3d[i])
					MFREE(data3d[i][j]);

			MFREE(data3d[i]);
		}
	}
	
	MFREE(data);
	MFREE(data3d);
	
	if (tmp3d)
	{
		for (i = 0; i < Xm; i++)
			MFREE(tmp3d[i]);
	}
	
	MFREE(tmp3d);
}

/* allocate memory */

int initialize (void) 
{
	int i, j, M;
	double x, y;
	
	if (!first) 
		cleanup();
	first = 0;
		
/*	if (hBval != 0 && hAval == 0)
		solution = &solution2;
	else
		solution = &solution1;	
*/
	
	dy = dx = configs[CF_H].fval;
	
	dt = dx*dx*dx*dx;
	
 	Xm = ceil((configs[CF_RIGHT].ival - configs[CF_LEFT].ival) / dx + 1);
	Ym = ceil((configs[CF_TOP].ival - configs[CF_BOTTOM].ival) / dy + 1);
	Tn = ceil(configs[CF_TIME].fval / dt + 1);
		
	MALLOC(data, double**, Tn * sizeof(double *))
	if (configs[CF_THREEDIM].ival) {
		MALLOC(data3d, double***, Tn * sizeof(double **))
	}
	MALLOC(maxv, double*, Tn * sizeof(double))
	MALLOC(error, double*, Tn * sizeof(double))
	
	for (i = 0; i < Tn; i++) {
		maxv[i] = TINY;
		error[i] = 0;
		
		MALLOC(data[i], double*, Xm * sizeof (double))
		
		for (j = 0; j < Xm; j++)
			data[i][j] = 0; 
		if (configs[CF_THREEDIM].ival) {
			MALLOC(data3d[i], double**, Xm * sizeof(double*))
		}
	}

	M = Xm;
	
	if (configs[CF_THREEDIM].ival) {
		for (i = 0; i < Tn; i++) {
			for (j = 0; j < Xm; j++) {
				MALLOC(data3d[i][j], double*, Ym * sizeof(double))
			}
		}

		if (Ym > M)
			M = Ym;
	}
		
	MALLOC(pdata, double*, M * sizeof(double))
	MALLOC(qdata, double*, M * sizeof(double))
	
	if (configs[CF_THREEDIM].ival)
	{	
		MALLOC(tmp3d, double**, Xm * sizeof(double*))
	
		for (i = 0; i < Xm; i++) {
			MALLOC(tmp3d[i], double*, Ym * sizeof(double))
		}	
	}
	
	for (i = 0; i < M; i++) {
		pdata[i] = 0;
		qdata[i] = 0;
	}	

	if (spin)
		spin->set_float_limits( 0, Tn/dt );
}

#undef MALLOC	

#define RES(X)			result += data[n-1][X];
#define RES0(X)			result += data[n-2][X];
#define RES1(X, A)		result += (A) * data[n-1][X];
#define RES2(X, FX, A)	result += (A) * ( data[n-1][X] + dt*(FFunc)(t-dt, FX) );
#define RES3(X, FX, A)	result += 2 * (A) * ( data[n-1][X] + dt*(FFunc)(t-dt, FX) );
#define RES4(X, FX, A)	result -= 2 * (A) * ( data[n-1][X] + dt*(FFunc)(t-dt, FX) );

void blackParams (int t, int mm, double alphas[3], double betas[3])
{
	double sigma2 = configs[CF_SIGMA].fval*configs[CF_SIGMA].fval;
	double r = configs[CF_R].fval;
	double m = mm;
	double m2 = m*m;
	
	alphas[0] = -(sigma2*m2/2.0 + r*m/2.0);
	alphas[1] = 2.0/dt + sigma2*m2 - r;
	alphas[2] = (-sigma2*m2/2.0 + r*m/2.0);
	
	betas[0] = -alphas[0];
	betas[1] = 2.0/dt - sigma2*m2 + r;
	betas[2] = -alphas[2];
}

double applyScheme(int scheme, int n, int m)
{
	double result = 0.0;	
	double lambda = configs[CF_LAMBDA].fval;
	double lambda2 = lambda*lambda;		
	double t = n * dt;
	double x = configs[CF_LEFT].ival + m * dx;

	double a = (aFunc)(t, x);
	double b = (bFunc)(t, x);
	
	double mu = configs[CF_MU].fval;
	
	int mm = m - 1;
	int mp = m + 1;

	if (configs[CF_PBC].ival)
	{
		mp = mp % Xm;
		if (mm < 0) 
			mm = Xm - 1;
	}
	else if (m == (Xm-1) && zeroRightBoundary) {
		return 0.0;
	}
	else if (m == (Xm-1) && configs[CF_SCHEME].ival != BACKWARD_SPACE) {
		mp = m;
	}
	
	if (scheme == LEAPFROG && n == 1)
		scheme = CENTRAL_SPACE;
	
	switch (scheme) {
		case FORWARD_SPACE:
			RES1(m, 1 + a * lambda)
			RES1(mp, - a * lambda)
			break;

		case BACKWARD_SPACE:
			RES1(mm, a * lambda)
			RES1(m, 1 - a * lambda)
			break;

		case CENTRAL_SPACE:
			RES(m)
			RES1(m, -2 * b * mu)
			RES1(mm, a * 0.5 * lambda)
			RES1(mp, - a * 0.5 * lambda)
			RES1(mm, b * mu)
			RES1(mp, b * mu)
			break;

		case LEAPFROG:
			RES0(m)
			RES1(mm, a * lambda)
			RES1(mp, - a * lambda)
			break;

		case LAX_FREDRICHS:
			RES1(mm, 0.5 * (1 + a * lambda))
			RES1(mp, 0.5 * (1 - a * lambda))
			break;			
					
		case EQUILIBRIUM:
			RES2(mm, x-dx,  0.5 * a * ( lambda + a*lambda2))
			RES2(m, x, 1 - a*lambda2)
			if (m != (Xm-1))
				RES2(mp, x+dx, 0.5 * a * (-lambda + a*lambda2))
			else {
				RES3(m, x, 0.5 * a * (-lambda + a*lambda2))
				RES4(mm, x-dx, 0.5 * a * (-lambda + a*lambda2))
			}
			break;
			
		case LAX_WENDROFF:
			RES1(m, 1 - a * a * lambda2)
			RES1(mm, 0.5 * a * lambda * ( 1 + a * lambda))
			RES1(mp, 0.5 * a * lambda * ( -1 + a * lambda))
			break;						

		default:
			break;
	}

	calcError(n, m, result);
	
	return result;
}

#undef RES(X)
#undef RES1(X, A)
#undef RES2(X, XF, A)
#undef RES3(X, XF, A)
#undef RES4(X, XF, A)

void calcError(const int n, const int m, const double val)
{	
	double t = n * dt;
	double x = configs[CF_LEFT].ival + m * dx;

	if (fabs(val) > BLOWUP_VAL && !blowup)
		blowup = n;
	
	if (val > maxv[n])
		maxv[n] = val;
	
	if (val > maxmaxv)
		maxmaxv = val;
	
	error[n] += pow(solution(t, x) - val, 2);	
}

void calcError3d(int n, int l, int m, double val)
{	
	double t = n * dt;
	double x[2];
	
	x[0] = configs[CF_LEFT].ival + l * dx;
	x[1] = configs[CF_BOTTOM].ival + m * dx;
	
	error[n] += pow(solution3d(t, x) - val, 2);
}

void save (void)
{
	char t[BUFLEN];
	int m, n, i = 0;
	
	FILE *o;

	saving = 1;
	idle();
		
	snprintf(t, BUFLEN, "V-%s.dat", derivNames[configs[CF_DERIV].ival]);
	o = fopen(t, "w");
	
	for (n = 0; n < Tn; n++) {		
		for (m = 1; m < Xm; m++)	
			fprintf(o, "%.4f ", data[n][m]);
		fprintf(o, "\n");
	}
	fflush(o);
	fclose(o);
	
	printf("Saved %s\n", t);
	
	snprintf(t, BUFLEN, "Max-%s.dat", derivNames[configs[CF_DERIV].ival]);
	o = fopen(t, "w");

	for (n = 0; n < Tn; n++) {		
		fprintf (o, "%.4f ", maxv[n]);
	}
	fprintf(o, "\n");

	fflush(o);
	fclose(o);

	printf("Saved %s\n", t);
	
	snprintf(t, BUFLEN, "Error-%d-%s.dat", (int)(1/dx), derivNames[configs[CF_DERIV].ival]);
	o = fopen(t, "w");
	
	for (n = 0; n < Tn; n++) {		
		if (configs[CF_THREEDIM].ival)
			fprintf (o, "%.4f ", sqrt(dy*dx*error[n]));
		else
			fprintf (o, "%.4f ", sqrt(dx*error[n]));
	}
	fprintf(o, "\n");
	
	fflush(o);
	fclose(o);
	
	printf("Saved %s\n", t);

	saving = 0;
}

void boundary (int n)
{
	double numer, denom;
	double b = bFunc(0, 0);   /* since b is constant */
	double lambda = configs[CF_LAMBDA].fval;
	double mu = configs[CF_MU].fval;
	int right = configs[CF_RIGHT].ival;
	
	switch (configs[CF_BOUNDARY].ival) {
		case BOUND_ONE:
			data[n][Xm-1] = data[n-1][Xm-2];
			break;
		
		case BOUND_TWO:
			denom = 1 + lambda * (1 - pdata[Xm-1]) - dt;
			numer = data[n-1][Xm-1] + lambda * qdata[Xm-1];
			data[n][Xm-1] = numer/denom;
			break;		
		
		case BOUND_THREE:
			denom = 1 - pdata[Xm-1]*(2 - pdata[Xm-2]);
			numer = qdata[Xm-1]*(2 - pdata[Xm-2]) - qdata[Xm-2];
			data[n][Xm-1] = numer/denom;
			break;
		
		/* Neumann Boundary */
		case BOUND_FOUR:
			denom = 1 + b*mu/2*(1 - pdata[Xm-1]);
			numer = data[n-1][Xm-1] * (1 - b*mu/2) 
				+ b*mu/2* ( data[n-1][Xm-2] + qdata[Xm-1] );
			data[n][Xm-1] = numer/denom;
			break;
		
		/* Exact Boundary */
		default:
			data[n][Xm-1] = solution(n*dt, right);
			break;
	}
}

inline double datumy (const int n, const int l, const int m)
{
	return tmp3d[l][m];
}

inline double datumx (const int n, const int l, const int m)
{
	return data3d[n][m][l];	
}

inline double gety (const int n, const int l, const int m)
{
	return data3d[n][l][m];	
}

inline void sety (double val, const int n, const int l, const int m)
{
	data3d[n][l][m] = val;
}

inline double getx (const int n, const int l, const int m)
{
	return tmp3d[m][l];	
}

inline void setx (double val, const int n, const int l, const int m)
{
	tmp3d[m][l] = val;
}

#define Y_CONST 1
#define X_CONST 2

void thomas3d (const int n, const int l, const int type, const double *aa, const double *bb)
{
	int m, M;
	double st, dd, denom, numer, val;
	double a = aFunc(0, 0);   /* since b is constant */
	double b = bFunc(0, 0);   /* since b is constant */
	double mu = configs[CF_MU].fval;
	void (*set)(double, const int, const int, const int);
	double (*get)(const int, const int, const int);		
	double (*datum)(const int, const int, const int);
	double b0, b1;
	
	int x0 = configs[CF_LEFT].ival;
	int x1 = configs[CF_RIGHT].ival;
	int y0 = configs[CF_BOTTOM].ival;
	int y1 = configs[CF_TOP].ival;
	
	double m0[3][2] = { {x0, y0}, {x0, y0},	{x0, y0} };
	double m1[3][2] = { {x1, y1}, {x1, y1}, {x1, y1} };
	
	if (type == Y_CONST)
	{
		M = Xm;
		get = getx;
		set = setx;
		datum = datumx;
		st = n*dt - dt/2.0;
	
		m0[0][1] += l*dy - dy;
		m0[1][1] += l*dy;
		m0[2][1] += l*dy + dy;
	
		m1[0][1] = m0[0][1];
		m1[1][1] = m0[1][1];
		m1[2][1] = m0[2][1];
	
		b0 =      b*mu/4 * solution3d(n*dt-dt, m0[0]);
		b0 += (1-b*mu)/2 * solution3d(n*dt-dt, m0[1]);
		b0 +=     b*mu/4 * solution3d(n*dt-dt, m0[2]);
		b0 +=    -b*mu/4 * solution3d(n*dt, m0[0]);
		b0 += (1+b*mu)/2 * solution3d(n*dt, m0[1]);
		b0 +=    -b*mu/4 * solution3d(n*dt, m0[2]);
		
		//	b0 = solution3d(st, m0[1]);

		b1 =      b*mu/4 * solution3d(n*dt-dt, m1[0]);
		b1 += (1-b*mu)/2 * solution3d(n*dt-dt, m1[1]);
		b1 +=     b*mu/4 * solution3d(n*dt-dt, m1[2]);
		b1 +=    -b*mu/4 * solution3d(n*dt, m1[0]);
		b1 += (1+b*mu)/2 * solution3d(n*dt, m1[1]);
		b1 +=    -b*mu/4 * solution3d(n*dt, m1[2]);

		//	b1 = solution3d(st, m1[1]);
	}
	else {
		M = Ym;
		get = gety;
		set = sety;
		datum = datumy;
		st = n*dt;
		
		m0[0][0] += l*dx - dx;
		m0[1][0] += l*dx;
		m0[2][0] += l*dx + dx;
		
		m1[0][0] = m0[0][0];
		m1[1][0] = m0[1][0];
		m1[2][0] = m0[2][0];
		
		b0 = solution3d(st, m0[1]);
		b1 = solution3d(st, m1[1]);
	}
		
	pdata[1] = 0;
	qdata[1] = b0;
	
	for (m = 1; m < M-1; m++) {
		dd =  bb[0] * datum(n-1, l-1, m);
		dd += bb[1] * datum(n-1, l, m);				
		dd += bb[2] * datum(n-1, l+1, m);
		
		denom = (aa[0] * pdata[m] + aa[1]);
		
		pdata[m+1] = -aa[2] / denom;
		qdata[m+1] = (dd - aa[0] * qdata[m]) / denom;
	}
	
	/* apply bounday condition */
			
	set(b1, n, l, M-1);
	
	for (m = M-2; m >= 0; m--) {
		val = pdata[m+1] * get(n, l, m+1) + qdata[m+1];
		set(val, n, l, m);
		if (type == X_CONST)
			calcError3d(n, l, m, val);
	}
}

void thomas (const int n)
{
	int m;
	double dd, denom, numer, val;
	double alphas[3], betas[3];
	int left = configs[CF_LEFT].ival;
	int right = configs[CF_RIGHT].ival;

	pdata[1] = 0;

	/* left boundary */
	qdata[1] = 0;

	for (m = 1; m < Xm-1; m++) {
		blackParams(n, m, alphas, betas);

		dd =  betas[0] * data[n-1][m+1];
		dd += betas[1] * data[n-1][m];
		dd += betas[2] * data[n-1][m-1];

		denom = (alphas[0] * pdata[m] + alphas[1]);

		pdata[m+1] = -alphas[2] / denom;
		qdata[m+1] = (dd - alphas[0] * qdata[m]) / denom;
	}

	/* right boundary */
	val = solution((double)n*dt, (double)right);
	data[n][Xm-1] = val;
	boundary(n);

	for (m = Xm-2; m >= 0; m--) {
		val = pdata[m+1]*data[n][m+1] + qdata[m+1];
		calcError(n, m, val);
		data[n][m] = val;
	}
}

void fdm (void)
{
	int m, n, i;
	double v0m, x;
	int x0 = configs[CF_LEFT].ival;
	int y0 = configs[CF_BOTTOM].ival;
	double k = configs[CF_K].fval;
	
	printf("Performing Finite Difference...\n");
	
	/* reset maxima and error */
	for (i = 0; i < Tn; i++) {
		maxv[i] = TINY;
		error[i] = 0;
	}
	maxmaxv = TINY;
	
	blowup = 0;	
	
	/* apply initial condition */
	
	for (m = 0; m < Xm; m++)
	{
		x = x0 + m*dx;
		v0m = ((x0 + m*dx) - k);
		data[0][m] = (v0m >= 0) ? v0m : 0;			
		calcError(0, m, (v0m >= 0) ? v0m : 0);
	}
	
	/* solve pde */
	
	for (n = 1; n < Tn; n++)
		thomas(n);
	
	reshapeWorld();
}

void drawString(char *s, void *font, int align, double x, double y)
{
	double width, wwidth;
	
	wwidth = (configs[CF_RIGHT].ival - configs[CF_LEFT].ival);
	width = glutBitmapLength(font, (unsigned char *)s);

	if (align == ALIGN_CENTER)
		glRasterPos2f(configs[CF_LEFT].ival + 0.5*wwidth*(1 - width/glutGet(GLUT_WINDOW_WIDTH)), y);
	else if (align == ALIGN_RIGHT)
		glRasterPos2f(x - wwidth*width/glutGet(GLUT_WINDOW_WIDTH), y);		
	else // (align == ALIGN_RIGHT)
		glRasterPos2f(x, y);
		
	while (*s != '\0')
		glutBitmapCharacter(font, *s++);	
}

void display3D (int t)
{
	int x, y, i;
	double soln, denom;
	char c[BUFLEN];
	int height = (configs[CF_TOP].ival-configs[CF_BOTTOM].ival);

	GLfloat v1[3], v2[3], v3[3], vn[3];
	GLfloat n1[3], n2[3];

	GLfloat ambient[] = {0.0, 0.0, 1.0, 0.0};
    GLfloat diffuse[] = {0.0, 0.0, 1.0, 0.0};
    GLfloat specular[] = {1.0, 1.0, 1.0, 0.0};
    GLfloat position[] = {-10.0, 10.0, -10.0, 1.0};

	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		
	if (!configs[CF_WIRE].ival)
	{
		glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		
		glEnable(GL_LIGHT0);
		glDepthFunc(GL_LESS);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);		
		
		glShadeModel(GL_SMOOTH);
		
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	}
	
	glPushMatrix();
	
//	glRotated(90, 1, 0, 0);
//	glRotated(90, 0, 0, 1);
	
	/* draw solution */
	double vx[2];
	for (x = 0; x < Xm-1; x++) {
		glBegin(GL_TRIANGLE_STRIP);

		vx[0] = v1[0] = configs[CF_LEFT].ival + x*dx;
		vx[1] = v1[2] = configs[CF_BOTTOM].ival;
		if (configs[CF_INIT_FUNC].ival > INIT_SIX)
			v1[1] = (drawSolution) ? solution3d(t*dt, vx) : data3d[t][x][0];
		else
			v1[1] = data[0][x];
		
		vx[0] = v2[0] = v1[0] + dx;
		vx[1] = v2[2] = v1[2];
		if (configs[CF_INIT_FUNC].ival > INIT_SIX)
			v2[1] = (drawSolution) ? solution3d(t*dt, vx) : data3d[t][x+1][0];
		else
			v2[1] = data[0][x+1];
		
		glVertex3fv(v1);
		glVertex3fv(v2);
		
		for (y = 1; y < Tn; y++)
		{			
			/* compute the normal */			
			vx[0] = v3[0] = configs[CF_LEFT].ival + x*dx;
			vx[1] = v3[2] = configs[CF_BOTTOM].ival + y*dy;			
			if (configs[CF_INIT_FUNC].ival > INIT_SIX)
				v3[1] = (drawSolution) ? solution3d(t*dt, vx) : data3d[t][x][y];
			else
				v3[1] = data[y][x];

			VSUB3(n2, v2, v1)
			VSUB3(n1, v3, v1)
			VCROSS3(vn, n1, n2)
			
			denom = vn[0]*vn[0] + vn[1]*vn[1] + vn[2]*vn[2];
			denom = sqrt(denom);
			if (!denom) denom = 1;
			
			glNormal3f(vn[0]/denom, vn[1]/denom, vn[2]/denom);

			VEQ3(v1, v3)
			
			vx[0] = v2[0] = v1[0] + dx;
			vx[1] = v2[2] = v1[2];						
			if (configs[CF_INIT_FUNC].ival > INIT_SIX)
				v2[1] = (drawSolution) ? solution3d(t*dt, vx) : data3d[t][x+1][y];
			else
				v2[1] = data[y][x+1];
			
			glVertex3fv(v1);
			glVertex3fv(v2);
		}

		glEnd();
		
	}
	
	if (!configs[CF_WIRE].ival)
		glDisable(GL_LIGHTING);		

	/* draw coordinate axes */
	
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	
	glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(4, 0, 0);	

		glColor3f(0, 1, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 4, 0);
	
		glColor3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 4);
	glEnd();	

	glPopMatrix();
	
	if (!configs[CF_WIRE].ival)
		glDisable(GL_DEPTH_TEST);
}

void display2D (int t)
{	
	int x, i;
	char c[BUFLEN];
	double left = configs[CF_LEFT].ival;
	int height = (configs[CF_TOP].ival-configs[CF_BOTTOM].ival);
	
	/* draw grid */
		
	glLineWidth (1);
	
	glBegin(GL_LINES);
	for (x = 0; x <= (configs[CF_RIGHT].ival - left)*linespacing; x++) {
		if (x % linespacing == 0)
			glColor3f(0.9f, 0.9f, 0.9f);
		else
			glColor3f(0.5f, 0.5f, 0.5f);
		glVertex2f(left + x/(double)linespacing, configs[CF_BOTTOM].ival);
		glVertex2f(left + x/(double)linespacing, configs[CF_TOP].ival);
	}
	glEnd();	
	
	for (x = 0; x <= (configs[CF_RIGHT].ival - left); x++) {
		snprintf(c, BUFLEN, "%d", (int)left + x);
		drawString(c, GLUT_BITMAP_HELVETICA_10, ALIGN_LEFT, left + x, configs[CF_BOTTOM].ival - vmargin/2);
		drawString(c, GLUT_BITMAP_HELVETICA_10, ALIGN_LEFT, left + x, configs[CF_TOP].ival + vmargin/4);		
	}
	
	double scaleh = 1;
	if (height > 10)
		scaleh = height/10.0;
	int lines = floor(height/scaleh);
	
	glBegin(GL_LINES);
	for (x = 0; x <= lines*linespacing; x++) {
		if (x % linespacing == 0)
			glColor3f(0.9f, 0.9f, 0.9f);
		else
			glColor3f(0.5f, 0.5f, 0.5f);
		glVertex2f(left, configs[CF_BOTTOM].ival + x*scaleh/(double)linespacing);
		glVertex2f(configs[CF_RIGHT].ival, configs[CF_BOTTOM].ival + x*scaleh/(double)linespacing);
	}
	glEnd();	
	
	for (x = 0; x <= lines; x++) {
		snprintf(c, BUFLEN, "%d", (int)floor(configs[CF_BOTTOM].ival + x*scaleh));
		drawString(c, GLUT_BITMAP_HELVETICA_10, ALIGN_LEFT, left - vmargin/2, configs[CF_BOTTOM].ival + x*scaleh);
		drawString(c, GLUT_BITMAP_HELVETICA_10, ALIGN_RIGHT, configs[CF_RIGHT].ival + vmargin/2, configs[CF_BOTTOM].ival + x*scaleh);		
	}
	
	/* draw error max(v_m(n)) */
	
	if (maxVmN)
	{
		glColor3f(0.0f, 0.0f, 1.0f);
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < Tn; i++) {
			glVertex2f(left + i * (configs[CF_RIGHT].ival - left)/(double)Tn, maxv[i]*3/maxmaxv - 1);
		}
		glEnd();	
	}
	else {
		
		/* draw actual solution */
		
		double soln = 0;
		if (drawSolution) {
			glPushAttrib(GL_LINE_BIT);
			glLineStipple (3, 0xAAAA);
			glColor3f(1.0f, 0.0f, 0.0f);
			glBegin(GL_LINE_STRIP);
			for (x = 0; x < Xm; x++) {
				soln = solution((t % Tn) * dt, left + x*dx);
				//				if (logy)
				//					soln = logf(soln)/t;
				glVertex2f(left + x*dx, soln);
			}
			glEnd();	
			glPopAttrib();
		}
		
		/* draw solution */
		
		glColor3f(0.0f, 1.0f, 0.0f);
		glBegin(GL_LINE_STRIP);
		for (x = 0; x < Xm; x++) {
			soln = data[t % Tn][x];
			//			if (logy)
			//				soln = logf(soln)/t;
			glVertex2f(left + x*dx, soln);
		}
		glEnd();		
	}
	
}

void getEyePos(double *pos)
{
	pos[0] = configs[CF_RHO].fval * cos(configs[CF_THETA].fval) * sin(configs[CF_PHI].fval);
	pos[1] = configs[CF_RHO].fval * sin(configs[CF_THETA].fval) * sin(configs[CF_PHI].fval);
	pos[2] = configs[CF_RHO].fval * cos(configs[CF_PHI].fval);
}

void pdestring(char *c)
{
	char aname[BUFLEN];
	char bname[BUFLEN];
	char tmp[BUFLEN];
	int	v;
	v = (configs[CF_THREEDIM].ival) ? -1 : 1;
			 
	bname[0] = aname[0] = '\0';

	if (activeA == 0 && configs[CF_ACONST].fval)
	{
		snprintf(tmp, BUFLEN, "D%s[u]", (configs[CF_THREEDIM].ival) ? "xx" : "x");
		
		if (configs[CF_ACONST].fval == v)
			snprintf(aname, BUFLEN, " + %s", tmp);
		else if (configs[CF_ACONST].fval == -v)
			snprintf(aname, BUFLEN, " - %s", tmp);
		else
			snprintf(aname, BUFLEN, " %s%.1f %s", (v*configs[CF_ACONST].fval > 0) ? "+ " : "", 
					 v*configs[CF_ACONST].fval, tmp);
	}
	else if (activeA != 0)
		snprintf(aname, BUFLEN, " + %s Dx[u]", aFuncNames[activeA]);

	snprintf(tmp, BUFLEN, "D%s[u]", (configs[CF_THREEDIM].ival) ? "yy" : "xx");	
			 
	if (configs[CF_BCONST].fval == v)
		snprintf(bname, BUFLEN, " + %s", tmp);
	else if (configs[CF_BCONST].fval == -v)
		snprintf(bname, BUFLEN, " - %s", tmp);
	else if (configs[CF_BCONST].fval) 
		snprintf(bname, BUFLEN, " %s%.1f %s", (-configs[CF_BCONST].fval > 0) ? "+ " : "", 
				 -configs[CF_BCONST].fval, tmp);
	
	if (uTerm)
		snprintf(c, BUFLEN, "Dt[u]%s%s%s = 0", aname, bname, (udir) ? " - u" : " + u");
	else
		snprintf(c, BUFLEN, "Dt[u]%s%s = 0", aname, bname);
	
}

void display(void)
{
	int x, t, i;
	char c[BUFLEN];
	int height = (configs[CF_TOP].ival-configs[CF_BOTTOM].ival);
	double v[3];
	
	if (lastFrameTime == 0)
        lastFrameTime = glutGet(GLUT_ELAPSED_TIME);
    
	int now = glutGet(GLUT_ELAPSED_TIME);
    int elapsedMilliseconds = now - lastFrameTime;
	float elapsedTime = elapsedMilliseconds / (1000.0f * dt);
    lastFrameTime = now;
	
	if (!mypause)
		mytime += elapsedTime;		
	
	t = mytime;

	if (t >= Tn) 
	{
		if (repeat) {
			t = 0;
			mytime = 0;
		}
		else {
			t = Tn - 1;
			mytime = TIME_END/dt;
		}

	}
	
	if (glui)
	{
		hTime = mytime*dt;
		glui->sync_live();
	}
		
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (saving)
		glutSetCursor(GLUT_CURSOR_WAIT);
	else
		glutSetCursor(GLUT_CURSOR_INHERIT);
	
	if (!configs[CF_THREEDIM].ival)
		display2D(t);
	else
	{
		reshapeWorld();	
	
		glPushMatrix();
		glLoadIdentity();
		
		getEyePos(v);
		gluLookAt(v[0], v[1], v[2], 0, 0, 0, 0, 1, 0);
		display3D(t);

		glPopMatrix();
	}
	
	/* draw title */
	
	glColor3f(1.0f, 1.0f, 1.0f);
	
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluOrtho2D(configs[CF_LEFT].ival, configs[CF_RIGHT].ival, configs[CF_BOTTOM].ival - vmargin, configs[CF_TOP].ival + vmargin + head);
    glMatrixMode(GL_MODELVIEW);

//	pdestring(c);
	
//	drawString(c, GLUT_BITMAP_HELVETICA_18, ALIGN_CENTER, 0, configs[CF_TOP].ival + vmargin + head*3/4);
	
//	pdetext->set_text(c);
	
//	snprintf(c, BUFLEN, "%s (Time %.2f)", schemeNames[configs[CF_SCHEME].ival], TIME_START + mytime*dt);
//	drawString(c, GLUT_BITMAP_HELVETICA_18, ALIGN_CENTER, 0, configs[CF_TOP].ival + vmargin + head*2/4);
	
//	snprintf(c, BUFLEN, "%.1f, %.1f (%.2f, %.2f)", mouseinfo[0], mouseinfo[1], mouseinfo[2], mouseinfo[3]);
//	drawString(c, GLUT_BITMAP_HELVETICA_18, ALIGN_CENTER, 0, configs[CF_TOP].ival + vmargin + head/4);

	if (configs[CF_THREEDIM].ival)
		snprintf(c, BUFLEN, "Error %.4f", sqrt(dy*dx*error[t]));
	else
		snprintf(c, BUFLEN, "Error %.4f", sqrt(dx*error[t]));
	
	drawString(c, GLUT_BITMAP_HELVETICA_18, ALIGN_CENTER, 0, configs[CF_TOP].ival + vmargin + head*3/4);
	
/*
	if (blowup > 0)
	{
		snprintf(c, BUFLEN, "Blow-up (Time %.2f)", TIME_START + (double)blowup*dt);
		drawString(c, GLUT_BITMAP_HELVETICA_18, ALIGN_CENTER, 0, configs[CF_TOP].ival + vmargin+ head*3/4);
	}
*/
	
	glutSwapBuffers();

}

void reshape(int width, int height)
{
	win_width = width;
	win_height = height;
	
    glViewport(0, 0, width, height);
    
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_LINE_STIPPLE);

	reshapeWorld();
}

void reshapeWorld()
{
	head = 0.3 * (configs[CF_TOP].ival - configs[CF_BOTTOM].ival);
	
	hmargin = 0.05 * (configs[CF_RIGHT].ival - configs[CF_LEFT].ival);
	vmargin = 0.05 * (configs[CF_TOP].ival - configs[CF_BOTTOM].ival);
	
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluPerspective(30, (GLfloat)win_width / (GLfloat)win_height, 1.0, 50.0);
    glMatrixMode(GL_MODELVIEW);		

}

void idle(void)
{
	if (glutGetWindow() != windMain)
		glutSetWindow(windMain);
    glutPostRedisplay();
}

void guiChange (int control)
{
	int dirty = 0;
	
	switch (control) {
		case SCHEME_CTRL:
			break;
			
		case H_CTRL:
			dirty = (dx!=1/(double)hListId); 
			configs[CF_H].fval = 1/(double)hListId;
			break;
	
		case A_CTRL:
			configs[CF_ACONST].fval = hAval;
/*			if (hAval == 0)
				solution = &solution2;
			else
				solution = &solution1;
*/			break;

		case B_CTRL:
			configs[CF_BCONST].fval = hBval;
/*			if (hBval == 0)
				solution = &solution1;
			else
				solution = &solution2;
*/			break;

		case TIME_CTRL:
			mytime = hTime/dt;
			break;

		case S_CTRL:
			configs[CF_SIGMA].fval = hSval;
			break;
			
		case R_CTRL:
			configs[CF_R].fval = hRval;
			break;

		case K_CTRL:
			configs[CF_K].fval = hKval;
			break;
			
		case MU_CTRL:
			dirty = (configs[CF_MU].fval != hMu);
			configs[CF_MU].fval = hMu;
			break;

		case LAMBDA_CTRL:
			dirty = (configs[CF_LAMBDA].ival != hLambda);
			configs[CF_LAMBDA].ival = hLambda;
			break;
			
		case XI_CTRL:
		case BC_CTRL:
			break;
			
		case TOP_CTRL:
		case BOTTOM_CTRL:
			reshapeWorld();
			break;

		case RIGHT_CTRL:
		case LEFT_CTRL:
			dirty = 1;
			reshapeWorld();
			break;			
			
		case SCALE_CTRL:
			reshapeWorld();
			break;
			
		case IC_CTRL:
			break;
			
		case SAVE_CTRL:
			writeConfig(configfile);
			break;

		case LOAD_CTRL:
			dirty = 1;
			loadConfig(configfile);
			break;

			
		default:
			break;
	}
	
	if (dirty)
	{
		initialize();
		glui->sync_live();			
	}
	
	fdm();
}

#define TEXT2PANL(N, P, T, Q, V, VV, C)											\
	GLUI_EditText *N = glui->add_edittext_to_panel(P, T, Q, &V, C, guiChange);	\

#define ITEXT2PANL(N, P, T, V, VV, C) TEXT2PANL(N, P, T, GLUI_EDITTEXT_INT, V, VV, C) \
	N->set_int_val(VV);

#define FTEXT2PANL(N, P, T, V, VV, C) TEXT2PANL(N, P, T, GLUI_EDITTEXT_FLOAT, V, VV, C) \
	N->set_float_val(VV);

#define SPIN2PANL(N, P, T, Q, V, C, MIN, MAX, S)								\
	GLUI_Spinner *N = glui->add_spinner_to_panel(P, T, Q, &V, C, guiChange);	\
	N->set_speed(S);															

#define ISPIN2PANL(N, P, T, V, C, MIN, MAX, S) SPIN2PANL(N, P, T, GLUI_SPINNER_INT, V, C, MIN, MAX, S)	\
	N->set_int_limits( MIN, MAX );																		\
	N->set_int_val(V);

#define FSPIN2PANL(N, P, T, V, C, MIN, MAX, S) SPIN2PANL(N, P, T, GLUI_SPINNER_FLOAT, V, C, MIN, MAX, S)	\
	N->set_float_limits( MIN, MAX );																		\
	N->set_float_val(V);

#define LIST2PANL(N, P, T, V, C, NUM, NMS, DEF)								\
	GLUI_Listbox *N = glui->add_listbox_to_panel(P, T, &V, C, guiChange);	\
	for (i = 0; i < NUM; i++)												\
		N->add_item(i, NMS[i]);												\
	N->set_int_val(DEF);	


void control (void)
{
	int i;
	char s[BUFLEN];
	
	glui = GLUI_Master.create_glui( "Controls", 0, 740, 0 );

	GLUI_Panel *panelparams = glui->add_panel ("Parameters");
	
	glui->add_edittext_to_panel(panelparams, "File: ", GLUI_EDITTEXT_TEXT, configfile);
	glui->add_button_to_panel(panelparams, "Save", SAVE_CTRL, guiChange);
	glui->add_button_to_panel(panelparams, "Load", LOAD_CTRL, guiChange);
	
	glui->add_column_to_panel(panelparams, false);
	
	glui->add_checkbox_to_panel(panelparams, "3D", &configs[CF_THREEDIM].ival);
//	glui->add_checkbox_to_panel(panelparams, "Wire Frame", &configs[CF_WIRE].ival );
	glui->add_checkbox_to_panel(panelparams, "Pause", &mypause );
	glui->add_checkbox_to_panel(panelparams, "Repeat", &repeat );
	glui->add_checkbox_to_panel(panelparams, "Solution", &drawSolution );
//	glui->add_checkbox_to_panel(panelparams, "Periodic Boundary", &configs[CF_PBC].ival );
//	glui->add_checkbox_to_panel(panelparams, "U(t, x)", &uTerm, UTERM_CTRL, guiChange );
	
	GLUI_Panel *panelpde = glui->add_panel ("Black Scholes");
	
//	pdetext = glui->add_statictext_to_panel(panelpde, "");

//	glui->add_separator_to_panel(panelpde);
	
//	LIST2PANL(listschemes, panelpde, "Scheme: ", configs[CF_SCHEME].ival, SCHEME_CTRL, NUM_SCHEMES, schemeNames, configs[CF_SCHEME].ival)
	LIST2PANL(derivs, panelpde, "Derivative Type: ", configs[CF_DERIV].ival, DERIV_CTRL, NUM_DERIVS, derivNames, configs[CF_DERIV].ival)

	spin = glui->add_spinner_to_panel(panelpde, "Time: ", GLUI_SPINNER_FLOAT, &hTime, TIME_CTRL, guiChange);
	spin->set_speed(0.1);
	spin->set_float_limits( 0, TIME_END );

	GLUI_Listbox *listbox2 = glui->add_listbox_to_panel(panelpde, "Div. Size:    ", &hListId, H_CTRL, guiChange);
	for (i = 1; i <= 10; i++)
	{
		snprintf(s, BUFLEN, "%d", i*10);
		listbox2->add_item(i*10, s);
	}
	listbox2->set_int_val(1.0/configs[CF_H].fval);

	glui->add_column_to_panel (panelpde, false);

	FTEXT2PANL(txtsig, panelpde, "Sigma: ", hSval, configs[CF_SIGMA].fval, S_CTRL)
	FTEXT2PANL(txtr, panelpde, "R: ", hRval, configs[CF_R].fval, R_CTRL)
	FTEXT2PANL(txtk, panelpde, "K: ", hKval, configs[CF_K].fval, K_CTRL)
		
//	FTEXT2PANL(txta, panelpde, "a(x, t): ", hAval, configs[CF_ACONST].fval, A_CTRL)
//	FTEXT2PANL(txtb, panelpde, "b(x, t): ", hBval, configs[CF_BCONST].fval, B_CTRL)
//	FTEXT2PANL(txtlam, panelpde, "Lambda: ", hLambda, configs[CF_LAMBDA].fval, LAMBDA_CTRL)

/*
	GLUI_Panel *panelinit = glui->add_panel ("Initial Conditions");

	LIST2PANL(listinit, panelinit, "U(0,x):     ", configs[CF_INIT_FUNC].ival, IC_CTRL, NUM_INIT, initNames, configs[CF_INIT_FUNC].ival)
	glui->add_column_to_panel (panelinit, false);
	ISPIN2PANL(spinxi, panelinit, "Xi:  ", xi, XI_CTRL, 0, 20, 0.1)	
*/
 
	GLUI_Panel *panelbounds = glui->add_panel ("Boundary");

	LIST2PANL(listboxbound, panelbounds, "Condition: ", configs[CF_BOUNDARY].ival, BC_CTRL, NUM_BOUND, boundNames, configs[CF_BOUNDARY].ival)
	
//	GLUI_Panel *panelbounds2 = glui->add_panel_to_panel(panelbounds, "Limits");

	ISPIN2PANL(spinleft, panelbounds, "Grid Left:  ", configs[CF_LEFT].ival, LEFT_CTRL, -100, 0, 0.1)
	ISPIN2PANL(spinbottom, panelbounds, "Grid Bottom:  ", configs[CF_BOTTOM].ival, BOTTOM_CTRL, -100, 0, 0.1)

	glui->add_column_to_panel (panelbounds, false);
	
	ISPIN2PANL(spinright, panelbounds, "Grid Right:  ", configs[CF_RIGHT].ival, RIGHT_CTRL, 1, 100, 0.1)	
	ISPIN2PANL(spintop, panelbounds, "Grid Top:  ", configs[CF_TOP].ival, TOP_CTRL, 1, 100, 0.1)

//	ISPIN2PANL(spindiv, panelbounds, "Grid Div.:  ", linespacing, DIV_CTRL, 1, 10, 0.1)
	
	glui->set_main_gfx_window( windMain );
	
	/* We register the idle callback with GLUI, *not* with GLUT */
	GLUI_Master.set_glutIdleFunc(idle);	
}


void Menu (int value) {
	int dirty = 0;
	
/*
	switch (value) {
		case FORWARD_SPACE:
		case BACKWARD_SPACE:
		case CENTRAL_SPACE:
		case LAX_FREDRICHS:
		case LEAPFROG:
		case EQUILIBRIUM:
		case LAX_WENDROFF:
		case CRANK_NICOLSON:
			configs[CF_SCHEME].ival = value;
			break;
			
		case H_10:
		case H_20:
		case H_40:
		case H_80:
			dirty = (dx!=1/(double)value); 
			dx = 1/(double)value;
			break;
			
		case LAMBDA_8:
			dirty = (lambda != 0.8);
			lambda = 0.8;
			break;

		case LAMBDA_10:
			dirty = (lambda != 1.0);
			lambda = 1.0;
			break;

		case LAMBDA_16:
			dirty = (lambda != 1.6);
			lambda = 1.6;
			break;
						
		case A_CONST:
			aFunc = &aConst;
			activeA = 0;
			break;		

		case A_SPACE:
			aFunc = &aSpace;
			activeA = 1;
			break;		

		case INIT_ONE:
			initFunc = &initOne;
			break;
			
		case INIT_TWO:
			initFunc = &initTwo;
			break;
			
		case INIT_THREE:
			initFunc = &initThree;
			break;		
			
		case INIT_FOUR:
			initFunc = &initFour;
			break;		

		case INIT_FIVE:
			initFunc = &initFive;
			break;		
			
		case F_ZERO:
			FFunc = &FZero;
			FActual = &FZero;
			break;
			
		case F_ONE:
			FFunc = &FOne;
			FActual = &FAOne;
			break;
			
		case MAX_VM_N:
			maxVmN = !maxVmN;
			break;
			
		case SAVE_ID:
			save();
			return;
			
		case SETTINGS_ID:
			control();
			return;
			
		default:
			return;
	}

	if (dirty) {
		initialize();
	}
	
	fdm();
 */

}

void buildMenu (void)
{
	int i, schemeMenuId, hMenuId, lambdaMenuId, aMenuId, initMenuId, fMenuId, menuId;
	char s[BUFLEN];
	
	schemeMenuId = glutCreateMenu(Menu);
	
	for (i = 0; i < NUM_SCHEMES; i++)
	{
		snprintf(s, BUFSIZ, "(%d) %s", i, schemeNames[i]);
		glutAddMenuEntry(s, i);
	}
		
	aMenuId = glutCreateMenu(Menu);
	
		glutAddMenuEntry(A_CONST_NAME, A_CONST);
		glutAddMenuEntry(A_SPACE_NAME, A_SPACE);

	initMenuId = glutCreateMenu(Menu);

		glutAddMenuEntry(INIT_ONE_NAME, INIT_ONE);
		glutAddMenuEntry(INIT_TWO_NAME, INIT_TWO);
		glutAddMenuEntry(INIT_THREE_NAME, INIT_THREE);
		glutAddMenuEntry(INIT_FOUR_NAME, INIT_FOUR);
		glutAddMenuEntry(INIT_FIVE_NAME, INIT_FIVE);

	fMenuId = glutCreateMenu(Menu);
	
		glutAddMenuEntry(F_ZERO_NAME, F_ZERO);
		glutAddMenuEntry(F_ONE_NAME, F_ONE);
	
	menuId = glutCreateMenu(Menu);
	
	glutAddSubMenu("Scheme", schemeMenuId);
	glutAddSubMenu("Grid Spacing (h)", hMenuId);
	glutAddSubMenu("λ (k/h)", lambdaMenuId);
	glutAddSubMenu("a(t, x)", aMenuId);
	glutAddSubMenu("Initial Condition (u0(x))v", initMenuId);
	glutAddSubMenu("F(t, x)", fMenuId);	
	glutAddMenuEntry("Save", SAVE_ID);
	
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	
}

void moveDown(void)
{
	configs[CF_THETA].fval -= M_PI/180;
}

void moveUp(void)
{
	configs[CF_THETA].fval += M_PI/180;
}

void moveLeft(void)
{
	configs[CF_PHI].fval += M_PI/180;
}

void moveRight(void)
{
	configs[CF_PHI].fval -= M_PI/180;
}

void moveForward(void)
{
	configs[CF_RHO].fval --;
}

void moveBackward(void)
{
	configs[CF_RHO].fval ++;
}

void Keyboard (unsigned char key, int x, int y)
{
	char p[2];
	
	p[0] = key;
	p[1] = '\0';
	
	switch (key) {
		case 'q':
			moveUp();
			break;
	
		case 'a':
			moveDown();
			break;
			
		case ' ':
			mypause = !mypause;
			break;

		case 'w':
			dir = !dir;
			fdm();
			break;

		case 'b':
			udir = !udir;
			fdm();
			break;
			
		case '+':
			if (mypause)
			{
				mytime += 1;
				if (mytime > Tn-1)
					mytime = 0;
			}
			break;

		case '-':
			if (mypause)
				mytime -= 1;
			if (mytime < 0)
				mytime = Tn-1;
			break;
			
		case '0':
		case '1':	
		case '2':	
		case '3':	
		case '4':	
		case '5':	
		case '6':	
		case '7':
			configs[CF_SCHEME].ival = atol(p);
			fdm();
			break;
			
		case 's':
			save();
			break;
			
		case 'p':
			configs[CF_PBC].ival = !configs[CF_PBC].ival;
			fdm();
			break;

		case 'r':
			zeroRightBoundary = !zeroRightBoundary;
			fdm();
			break;			
			
		case 'x':
			boundCondFunc = &BCOne;
			fdm();
			break;
			
		default:
			break;
	}
}

double worldx(int x)
{
	return configs[CF_LEFT].ival + (double)x*(configs[CF_RIGHT].ival - configs[CF_LEFT].ival + 0.2)/glutGet(GLUT_WINDOW_WIDTH) - 0.1;
}

double worldy(int y)
{
	double ddy = glutGet(GLUT_WINDOW_HEIGHT);
	return configs[CF_BOTTOM].ival + (double)(ddy - y)*(configs[CF_TOP].ival - configs[CF_BOTTOM].ival + 0.2)/ddy - 0.2;
}

#define MINFO(A, B, C, D)		\
	mouseinfo[0] = A;		\
	mouseinfo[1] = B;		\
	mouseinfo[2] = worldx(C);	\
	mouseinfo[3] = worldy(D);	

void Mouse (int button, int state, int x, int y)
{
	MINFO(button, state, x, y)
}

void PassiveMotionFunc (int x, int y)
{
	MINFO(0, 0, x, y)
}

void MotionFunc (int x, int y)
{
	MINFO(-1, 0, x, y)
}

#undef MINFO

void specialFunc(int key, int x, int y)
{
	switch (key) {
		case GLUT_KEY_LEFT:
			moveLeft();		
			break;
			
		case GLUT_KEY_RIGHT:
			moveRight();
			break;
		
		case GLUT_KEY_UP:
			moveForward();
			break;

		case GLUT_KEY_DOWN:
			moveBackward();
			break;
			
		default:
			break;
	}
	
}

int main(int argc, char** argv)
{
	int t, x;

    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(640, 640);
    
    windMain = glutCreateWindow("Numeric PDE Solver");
    
	control();

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
	glutSpecialFunc(specialFunc);
	glutMouseFunc(Mouse);
	glutKeyboardFunc(Keyboard);
	glutMotionFunc(MotionFunc);
	glutPassiveMotionFunc(PassiveMotionFunc);
	
	initialize();
	fdm();
	buildMenu();
	
    glutMainLoop();
	
    return EXIT_SUCCESS;
}