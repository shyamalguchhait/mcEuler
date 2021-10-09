/*******************************************************************************************
 *  Copyright Jessica C. Ramella-Roman, Steve L. Jacques and Scott A. Prahl 2005
 *  
 *  Stok1.c- Euler rotations
 *	Main program for Monte Carlo simulation of photon
 *	travelling into scattering media keeping track of 
 *  its status of polarization. Slab geometry.
 *  
 *  by Jessica C. Ramella-Roman
 *
 *  A report in preparation illustrates use of the program:
 *   
 *  J. Ramella-Roman, S. A. Prahl, S. L. Jacques 
 *  Three Monte Carlo programs of polarized light transport into scattering media: part I,
 *  2003 Optics Express, submitted 2005.
 * 
 *	This program is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU General Public License
 *	as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 ****/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "array.h"
#include "complex.h"
#include "mie.h"
#include "nrutil.h"
#include <time.h>

#define ALIVE 1
#define DEAD 0
#define NN 100
#define THRESHOLD 0.01 /* used in roulette */
#define CHANCE 0.1	   /* used in roulette */

#define RandomNum (double)RandomGen(1, 0, NULL)
#define SIGN(x) ((x) >= 0 ? 1 : -1)
#define InitRandomGen (double)RandomGen(0, 1, NULL)

/* Declare Subroutines */
double RandomGen(char Type, long Seed, long *Status);
void multS(double *S, double theta, double *S2);
void rotSphi(double *S, double phi, double *S2);
void rotateYYZZ(double *YY, double *ZZ, double phi, double theta, double *YY2, double *ZZ2);

/********************************** MAIN ****************************************/

int main()
{
	double pi = 3.1415926535897932384;

	/* Mie theory stuff */

	double radius, lambda, A;
	long nangles, i, j;
	struct complex m;
	struct complex *s1 = NULL;
	struct complex *s2 = NULL;
	double *mu = NULL;
	double x, qext, qsca, qback, g, rho, vol;
	double nre_p, nim_p, nre_med, nim_med;
	FILE *target;
	double IR_1, QR_1, UR_1, VR_1;
	double IT_1, QT_1, UT_1, VT_1;

	/* E field stuff */

	double phi, theta, I, I0;
	long ithedeg;

	/* Propagation parameters */

	double y, z, dy, dx, hw; /* photon position.  x already declared. Also, incrementals & max range. */
	double s;				 /* step sizes. s = -log(RND)/mus [cm] */
	long i_photon;			 /* current photon */
	short photon_status;	 /*  = ALIVE=1 or DEAD=0 */
	double mua;				 /* absorption coefficient [cm^-1] */
	double mus;				 /* scattering coefficient [cm^-1] */
	double albedo;			 /* albedo of tissue */
	double rnd;				 /* assigned random value 0-1 */
	int ix, iy, MM;
	double absorb, W;

	/**** allocate matrices and arrays *******/

	double **IR, **QR, **UR, **VR;
	double *U, *U2;
	double *S;	/* */
	double *S0; /* */
	double *S2; /* */
	double *s11 = NULL;
	double *s12 = NULL;
	double *s33 = NULL;
	double *s43 = NULL;
	double slabsize;
	int jjj;
	int Nphotons;
	double *XX, *XX2, *YY, *YY2, *ZZ, *ZZ2; /* Vyx, Vyy, Vyz, Vxx, Vxy, Vxz = coordinates of original x and y axes */
	double *IQUV;							/* [I, Q, U, V] Stokes Vector */
	double start_time, finish_time;

	MM = NN - 1;

	U = new_darray(3);
	U2 = new_darray(3);
	S = new_darray(4);
	S0 = new_darray(4);
	S2 = new_darray(4); /* dummy S*/

	XX = new_darray(3);
	YY = new_darray(3);
	ZZ = new_darray(3);

	XX2 = new_darray(6);
	YY2 = new_darray(6);
	ZZ2 = new_darray(6);

	IQUV = new_darray(4);

	/**** end  allocate matrices and arrays *******/

	start_time = clock();

	/* CHOOSE MIE SCATTERING parameters */
	radius = 2.0 / 2; /* microns */
	lambda = 0.6328;  /* microns */
	rho = 1.152e-4;	  /*Dilution 1*/
	Nphotons = 1e6;
	mua = 0.0; /*�a  */

	/* ------------------------*/
	nre_p = 1.59;
	nim_p = 0;
	nre_med = 1.33;
	nim_med = 0.0;
	nangles = 1000;

	IR = dmatrix(0, MM, 0, MM); /* [0:MM] */
	QR = dmatrix(0, MM, 0, MM); /* [0:MM] */
	UR = dmatrix(0, MM, 0, MM); /* [0:MM] */
	VR = dmatrix(0, MM, 0, MM); /* [0:MM] */

	/* Setup MIE SCATTERING parameters */
	mu = new_darray(nangles);
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);

	m.re = nre_p / nre_med;
	m.im = 0.0;
	x = 2 * pi * radius / (lambda / nre_med);
	vol = 4.0 / 3 * pi * radius * radius * radius;

	A = pi * radius * radius;

	for (i = 0; i <= nangles; i++)
		mu[i] = cos(pi * i / nangles);
	s11 = new_darray(nangles);
	s12 = new_darray(nangles);
	s33 = new_darray(nangles);
	s43 = new_darray(nangles);
	s1 = new_carray(nangles);
	s2 = new_carray(nangles);

	Mie(x, m, mu, nangles, s1, s2, &qext, &qsca, &qback, &g); /* <---- Call Scott's Mie program ----- */

	mus = qsca * A * rho * 1e4;
	albedo = mus / (mus + mua);

	free_darray(mu);

	for (i = 0; i <= nangles; ++i)
	{
		s11[i] = 0.5 * cabbs(s2[i]) * cabbs(s2[i]) + 0.5 * cabbs(s1[i]) * cabbs(s1[i]);
		s12[i] = 0.5 * cabbs(s2[i]) * cabbs(s2[i]) - 0.5 * cabbs(s1[i]) * cabbs(s1[i]);
		s33[i] = (cmul(conj(s1[i]), s2[i])).re;
		s43[i] = (cmul(conj(s1[i]), s2[i])).im;
		/*	printf("%5.5f\t %5.5f\t %5.5f\t %5.5f\n",s11[i],s12[i],s33[i],s43[i]);
	*/
	}
	IR_1 = 0; /*W*/
	QR_1 = 0;
	UR_1 = 0;
	VR_1 = 0;
	IT_1 = 0; /*W*/
	QT_1 = 0;
	UT_1 = 0;
	VT_1 = 0;

	for (iy = 0; iy < NN; iy++)
		for (ix = 0; ix < NN; ix++)
		{

			IR[iy][ix] = 0.0;

			QR[iy][ix] = 0.0;

			UR[iy][ix] = 0.0;

			VR[iy][ix] = 0.0;
		}

	/* SET UP Monte Carlo */

	hw = 7 / mus; /* [cm] , maximum range in x and y for output. */
	dx = 2.0 * hw / NN;
	dy = 2.0 * hw / NN;

	slabsize = 4 / mus;

	InitRandomGen;

	/* LAUNCHNOW*/

	printf("dia=%5.5f;\n mus=%5.5f;\n g=%5.5f;\n photons=%d;\n rho=%5.5f;\n", radius * 2, mus, g, Nphotons, rho);

	printf("slabsize =%1.7f\n", slabsize);

	for (jjj = 1; jjj <= 4; jjj++)
	{

		if (jjj == 1)
		{
			S0[0] = 1;
			S0[1] = 1;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch H\n");
		}

		if (jjj == 2)
		{
			S0[0] = 1;
			S0[1] = -1;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch V\n");
		}

		if (jjj == 3)
		{
			S0[0] = 1;
			S0[1] = 0;
			S0[2] = 1;
			S0[3] = 0;
			printf("launch P\n");
		}

		if (jjj == 4)
		{
			S0[0] = 1;
			S0[1] = 0;
			S0[2] = 0;
			S0[3] = 1;
			printf("launch R\n");
		}

		/* LAUNCH photon */
		for (i_photon = 1; i_photon <= Nphotons; i_photon++)
		{

			/* photon position */
			x = 0.0;
			y = 0.0;
			z = 0.0;

			/* V[0:5] = test */

			XX[0] = 1.0; /* original X-axis coord */
			XX[1] = 0.0;
			XX[2] = 0.0;

			/* V[0:5] = x */
			YY[0] = 0.0; /* original Y-axis coord */
			YY[1] = 1.0;
			YY[2] = 0.0;

			ZZ[0] = 0.0; /* original Y-axis coord */
			ZZ[1] = 0.0;
			ZZ[2] = 1.0;

			for (i = 0; i < 4; i++)
				S[i] = S0[i]; /* set incident Stokes vector to S0 */
			for (i = 0; i < 4; i++)
				S2[i] = 0.0; /* set incident Stokes vector to S0 */

			photon_status = ALIVE;
			W = 1; /* photon weight */

			while (photon_status == ALIVE)
			{

				/**** HOP */

				rnd = 0;
				while (rnd == 0)
					rnd = RandomNum; /* choose a step size */

				s = -log(RandomNum) / (mus + mua);

				x += ZZ[0] * s;
				y += ZZ[1] * s;
				z += ZZ[2] * s;

				/**** ABSORB */

				absorb = W * (1 - albedo);
				W -= absorb;

				if (z <= 0)
				{

					XX[1] = -(ZZ[0] * YY[2] - ZZ[2] * YY[0]);
					XX[2] = (ZZ[0] * YY[1] - ZZ[1] * YY[0]);

					if (fabs(YY[2]) < 1e-8 && fabs(XX[2]) < 1e-8)
					{
					}
					else
					{
						phi = atan2(YY[2], -XX[2]);
						rotSphi(S, phi, S2);

						phi = atan2(ZZ[1], ZZ[0]);
						rotSphi(S2, phi, S);
					}
					IR_1 += S[0];
					QR_1 += S[1];
					UR_1 += S[2];
					VR_1 += S[3];

					if (x >= -hw)
						ix = (int)(fabs(x + hw) / dx);
					if (y >= -hw)
						iy = (int)(fabs(y + hw) / dy);
					if (ix > MM)
						ix = MM; /* last bin [MM] is overflow */
					if (iy > MM)
						iy = MM; /* last bin [MM] is overflow */

					IR[iy][ix] += S[0];
					QR[iy][ix] += S[1];
					UR[iy][ix] += S[2];
					VR[iy][ix] += S[3];
					photon_status = DEAD;
				}

				else if (z >= slabsize)
				{

					XX[1] = -(ZZ[0] * YY[2] - ZZ[2] * YY[0]);
					XX[2] = (ZZ[0] * YY[1] - ZZ[1] * YY[0]);

					if (fabs(YY[2]) < 1e-8 && fabs(XX[2]) < 1e-8)
					{
					}
					else
					{
						phi = atan2(YY[2], -XX[2]);

						rotSphi(S, phi, S2);

						phi = atan2(ZZ[1], -ZZ[0]);

						rotSphi(S2, phi, S);
					}
					IT_1 += S[0];
					QT_1 += S[1];
					UT_1 += S[2];
					VT_1 += S[3];
					photon_status = DEAD;

				} /*z>slab size*/

				/* REJECTION METHOD to choose azimuthal angle phi and deflection angle theta */

				do
				{
					theta = acos(2 * RandomNum - 1);

					phi = RandomNum * 2.0 * pi;

					I0 = s11[0] * S[0] + s12[0] * (S[1] * cos(2 * phi) + S[2] * sin(2 * phi));

					ithedeg = floor(theta / pi * nangles);

					I = s11[ithedeg] * S[0] + s12[ithedeg] * (S[1] * cos(2 * phi) + S[2] * sin(2 * phi));

				} while (RandomNum * I0 >= I);

				rotSphi(S, phi, S2);

				for (i = 0; i < 4; i++)
					S[i] = S2[i];

				S2[0] = s11[ithedeg] * S[0] + s12[ithedeg] * S[1];

				S2[1] = s12[ithedeg] * S[0] + s11[ithedeg] * S[1];

				S2[2] = s33[ithedeg] * S[2] + s43[ithedeg] * S[3];

				S2[3] = -s43[ithedeg] * S[2] + s33[ithedeg] * S[3];

				S[1] = S2[1] / S2[0];
				S[2] = S2[2] / S2[0];
				S[3] = S2[3] / S2[0];
				S[0] = 1.0;

				rotateYYZZ(YY, ZZ, phi, theta, YY2, ZZ2); /* rotate coord of original x and y axis in local photon frame */

				for (i = 0; i < 3; i++)
				{
					YY[i] = YY2[i];
					ZZ[i] = ZZ2[i];
				}
				/*ROULETTE*/
				rnd = 0;
				while (rnd == 0)
					rnd = RandomNum;

				if (W < THRESHOLD)
				{
					if (rnd <= CHANCE)
						W /= CHANCE;
					else
						photon_status = DEAD;
				}

			} /* while (photon_status == ALIVE) */

		} /* end of single photon launching */

		printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ", IR_1 / (Nphotons), QR_1 / (Nphotons), UR_1 / (Nphotons), VR_1 / (Nphotons));
		printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ", IT_1 / (Nphotons), QT_1 / (Nphotons), UT_1 / (Nphotons), VT_1 / (Nphotons));

		IR_1 = 0;
		QR_1 = 0;
		UR_1 = 0;
		VR_1 = 0;

		IT_1 = 0;
		QT_1 = 0;
		UT_1 = 0;
		VT_1 = 0;

		if (jjj == 1)
		{

			target = fopen("outHI.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", IR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", IR[i][j]);
				fprintf(target, "\n");
			}

			fclose(target);

			target = fopen("outHQ.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", QR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", QR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outHU.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", UR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", UR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outHV.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", VR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", VR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);
			for (iy = 0; iy < NN; iy++)
				for (ix = 0; ix < NN; ix++)
				{
					IR[iy][ix] = 0.0;
					QR[iy][ix] = 0.0;
					UR[iy][ix] = 0.0;
					VR[iy][ix] = 0.0;
				}
		}
		/*111111111111111111111111111111111111111111111111111111111111111*/

		if (jjj == 2)
		{

			/* save data to file */
			target = fopen("outVI.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", IR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", IR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outVQ.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", QR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", QR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outVU.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", UR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", UR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outVV.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", VR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", VR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);
			for (iy = 0; iy < NN; iy++)
				for (ix = 0; ix < NN; ix++)
				{
					IR[iy][ix] = 0.0;
					QR[iy][ix] = 0.0;
					UR[iy][ix] = 0.0;
					VR[iy][ix] = 0.0;
				}
		}
		/* 222222222222222222222222222222222222222222222222222222222222222*/

		if (jjj == 3)
		{

			/* save data to file */

			target = fopen("outPI.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", IR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", IR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outPQ.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", QR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", QR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outPU.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", UR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", UR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outPV.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", VR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", VR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);
			for (iy = 0; iy < NN; iy++)
				for (ix = 0; ix < NN; ix++)
				{
					IR[iy][ix] = 0.0;
					QR[iy][ix] = 0.0;
					UR[iy][ix] = 0.0;
					VR[iy][ix] = 0.0;
				}
		}
		/* 33333333333333333333333333333333333333333333333333333333333333333333*/

		if (jjj == 4)
		{
			/* save data to file */
			target = fopen("outRI.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", IR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", IR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outRQ.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", QR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", QR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outRU.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", UR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", UR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);

			target = fopen("outRV.dat", "w");
			for (i = 0; i < NN; i++)
			{
				fprintf(target, "%5.5f", VR[i][0]);
				for (j = 1; j < NN; j++)
					fprintf(target, "\t%5.5f", VR[i][j]);
				fprintf(target, "\n");
			}
			fclose(target);
		}
		/* 444444444444444444444444444444444444444444444444444444444444444444444*/
	} /* end of 4 photon launchings */
	finish_time = clock();
	printf("Elapsed Time              = %10.2f seconds\n", (double)(finish_time - start_time) / CLOCKS_PER_SEC);
	fflush(NULL);
	return 0;
} /*main*/

/*************** end MAIN ************************************/
/*************************************************************/

/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double
RandomGen(char Type, long Seed, long *Status)
{
	static long i1, i2, ma[56]; /* ma[0] is not used. */
	long mj, mk;
	short i, ii;

	if (Type == 0)
	{ /* set seed. */
		mj = MSEED - (Seed < 0 ? -Seed : Seed);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++)
		{
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ)
				mk += MBIG;
			mj = ma[ii];
		}
		for (ii = 1; ii <= 4; ii++)
			for (i = 1; i <= 55; i++)
			{
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ)
					ma[i] += MBIG;
			}
		i1 = 0;
		i2 = 31;
	}
	else if (Type == 1)
	{ /* get a number. */
		if (++i1 == 56)
			i1 = 1;
		if (++i2 == 56)
			i2 = 1;
		mj = ma[i1] - ma[i2];
		if (mj < MZ)
			mj += MBIG;
		ma[i1] = mj;
		return (mj * FAC);
	}
	else if (Type == 2)
	{ /* get status. */
		for (i = 0; i < 55; i++)
			Status[i] = ma[i + 1];
		Status[55] = i1;
		Status[56] = i2;
	}
	else if (Type == 3)
	{ /* restore status. */
		for (i = 0; i < 55; i++)
			ma[i + 1] = Status[i];
		i1 = Status[55];
		i2 = Status[56];
	}
	else
		puts("Wrong parameter to RandomGen().");
	return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/************************************************************************************
 *	rotSphi(S,phi,S)
 *		Rotate S by phi [radians] and return as S
 *      multiply S for the rotational matrix of Chandrasekar or Boheren and Hoffman
 *		Uses invtan()
 ****/
void rotSphi(double *S, double phi, double *S2)
{
	double cos2phi, sin2phi;

	cos2phi = cos(2 * phi);
	sin2phi = sin(2 * phi);

	S2[0] = S[0];
	S2[1] = S[1] * cos2phi + S[2] * sin2phi;
	S2[2] = -S[1] * sin2phi + S[2] * cos2phi;
	S2[3] = S[3];
}

/**************************************************************************
 *	rotateYYZZ(YY,ZZ,phi, theta, YY,ZZ)
 *		Rotates coordinates of original x and y coord as seen by local photon frame.
 ****/
void rotateYYZZ(double *YY, double *ZZ, double phi, double theta, double *YY2, double *ZZ2)
{
	double ct, st, vt, temp;

	/* phi rotation around local z axis */
	ct = cos(phi);
	st = sin(phi);
	vt = 1 - ct;

	YY2[0] = (ZZ[0] * ZZ[0] * vt + ct) * YY[0] + (ZZ[0] * ZZ[1] * vt - ZZ[2] * st) * YY[1] + (ZZ[0] * ZZ[2] * vt + ZZ[1] * st) * YY[2];
	YY2[1] = (ZZ[0] * ZZ[1] * vt + ZZ[2] * st) * YY[0] + (ZZ[1] * ZZ[1] * vt + ct) * YY[1] + (ZZ[1] * ZZ[2] * vt - ZZ[0] * st) * YY[2];
	YY2[2] = (ZZ[0] * ZZ[2] * vt - ZZ[1] * st) * YY[0] + (ZZ[1] * ZZ[2] * vt + ZZ[0] * st) * YY[1] + (ZZ[2] * ZZ[2] * vt + ct) * YY[2];

	temp = 1 / sqrt(YY2[0] * YY2[0] + YY2[1] * YY2[1] + YY2[2] * YY2[2]);

	YY2[0] = YY2[0] * temp;
	YY2[1] = YY2[1] * temp;
	YY2[2] = YY2[2] * temp;

	ct = cos(theta);
	st = sin(theta);
	vt = 1 - ct;

	ZZ2[0] = (YY2[0] * YY2[0] * vt + ct) * ZZ[0] + (YY2[0] * YY2[1] * vt - YY2[2] * st) * ZZ[1] + (YY2[0] * YY2[2] * vt + YY2[1] * st) * ZZ[2];
	ZZ2[1] = (YY2[0] * YY2[1] * vt + YY2[2] * st) * ZZ[0] + (YY2[1] * YY2[1] * vt + ct) * ZZ[1] + (YY2[1] * YY2[2] * vt - YY2[0] * st) * ZZ[2];
	ZZ2[2] = (YY2[0] * YY2[2] * vt - YY2[1] * st) * ZZ[0] + (YY2[1] * YY2[2] * vt + YY2[0] * st) * ZZ[1] + (YY2[2] * YY2[2] * vt + ct) * ZZ[2];

	temp = 1 / fabs(sqrt(ZZ2[0] * ZZ2[0] + ZZ2[1] * ZZ2[1] + ZZ2[2] * ZZ2[2]));

	ZZ2[0] = ZZ2[0] * temp;
	ZZ2[1] = ZZ2[1] * temp;
	ZZ2[2] = ZZ2[2] * temp;
}
