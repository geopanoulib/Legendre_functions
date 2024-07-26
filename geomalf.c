/**
 *
 * Copyright (c) 2024, Georgios Panou
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * Authors: Koci J., Iossifidis C. and Panou G. <geopanou@survey.ntua.gr>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define format "%-10d%-10d%24.16le\n"

int main(void)
{
	register int p, q, t;
	int N, n, m, n2, np1, nm1, nm2, nm3;
	double theta, pihalf, raddeg, cosnm1x, cosnm2x, cosnx, sinnm1x, sinnm2x, sinnx;
	double A, B, C, coeff, coefn, coss, sins, Pi0, Pi1, pnn, pn0, Tni;						
	double *root, *rooti, *cosn, *sinn, *Pn0, *Pn1, *Pn_even, *Pn_odd, *Pn;
	double anm, bnn, bnm, cnn, cnm, bnn_1, cnn_1, Pn_2_m, Pn_2_m_2, Pn_m_2, Pnm, LFnum;
	double sqr2, sqr3, z;
	
	N = 6; /* The maximum degree */
	theta = 30.0; /* The co-latitude in degrees */
	pihalf = 2.0l * atan(1.0l);
	raddeg = pihalf / 90.0l;
	theta *= raddeg;
	LFnum = (double)(N + 1) * (N + 2) / 2; /* The number of Legendre Functions (in double type to avoid overflow) */
	
	if ((root = malloc((2 * N + 5) * sizeof(double))) == NULL)
		exit(1);
	else if ((rooti = malloc((2 * N + 4) * sizeof(double))) == NULL)
		exit(2);
	else if ((cosn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(3);
	else if ((sinn = malloc((N + 1) * sizeof(double))) == NULL)
		exit(4);
	else if ((Pn0 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(5);
	else if ((Pn1 = malloc((N + 1) * sizeof(double))) == NULL)
		exit(6);
	else if ((Pn_even = malloc((N + 3) * sizeof(double))) == NULL)
		exit(7);
	else if ((Pn_odd = malloc((N + 3) * sizeof(double))) == NULL)
		exit(8);
	else if ((Pn = malloc((N + 3) * sizeof(double))) == NULL)
		exit(9);
	
	printf("\nN = %d\nCo-latitude = %.8le deg\nThe number of Legendre Functions = %.0f\n", N, theta / raddeg, LFnum);
	/* Lookup tables for the square roots */
	for (t = 0; t <= 2 * N + 4; t++)
		root[t] = sqrt(1.0l * t);
	p = sqrt(2.0l * N + 4) + 1;
	for (t = 0; t < p; t++)
	{
		q = t * t;
		root[q] = 1.0l * t;
	}
	sqr2 = root[2];
	sqr3 = root[3];
	coeff = 0.0l;
	for (t = 0; t <= 2 * N + 3; t++)
	{
		z = root[t + 1];
		rooti[t] = coeff * z;
		coeff = z;
	}
	
	/**
	 * Computes all the multiple angle cosines by using the Chebyshev's method,
	 * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method
	 **/
	 	
	cosnm1x = cos(theta);
	coeff = 2.0l * cosnm1x;
	cosnm2x = 1.0l;
	cosn[0] = 1.0l;
	cosn[1] = cosnm1x;
	for (n = 2; n <= N; n++)
	{
		cosnx = coeff * cosnm1x - cosnm2x;
		cosn[n] = cosnx;
		cosnm2x = cosnm1x;
		cosnm1x = cosnx;
	}
	
	/**
	 * Computes all the multiple angle sines by using the Chebyshev's method,
	 * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method
	 **/

	coeff = 2.0l * cosn[1];
	sinnm2x = 0.0l;
	sinnm1x = sin(theta);
	sinn[0] = 0.0l;
	sinn[1] = sinnm1x;
	for (n = 2; n <= N; n++)
	{
		sinnx = coeff * sinnm1x - sinnm2x;
		sinn[n] = sinnx;
		sinnm2x = sinnm1x;
		sinnm1x = sinnx;
	}
	
	/* Computes the Legendre Polynomials (m = 0) and Functions (m = 1) */
	Pn0[0] = 1.0l;
	Pn0[1] = sqr3 * cosn[1];
	Pn1[0] = 0.0l;
	Pn1[1] = sqr3 * sinn[1];
	
	/* Computes the odd Legendre Polynomials (m = 0) and Functions (m = 1) */
	pnn = 1.0l;
	coss = cosn[1];
	sins = sinn[1];
	for (n = 3; n <= N; n += 2)
	{
		p = n - 1;
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75l * B;
		pnn *= 1.0l - A * coeff;
		Pi0 = coss;
		Pi1 = sins;
		q = n + 2;
		t = 3;
		while (p > 0) 
		{
			/* Computes the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni = B * C;
			Pi0 *= Tni;
			Pi1 *= Tni;
			Pi0 += cosn[t];
			Pi1 += t * sinn[t];
			t += 2;
			p -= 2;
			q += 2;
		}
		coeff = root[2 * n + 1];
		Pn0[n] = Pi0 * pnn * coeff;
		Pn1[n] = Pi1 * pnn * sqr2 * coeff / rooti[n];
	}
	
	/* Computes the even Legendre Polynomials (m = 0) and Functions (m = 1) */
	pn0 = 1.0l;
	pnn = 2.0l;
	coss = cosn[2];
	sins = sinn[2];
	for (n = 2; n <= N; n += 2)
	{
		p = n - 1;
		A = 1.0l / p;
		B = 1.0l / n;
		coeff = 1.0l - 0.75l * B;
		pnn *= 1.0l - A * coeff;
		/* Computes the tn ratio */
		A = 1.0l - B;
		pn0 *= A * A;
		Pi0 = coss;
		Pi1 = 2.0l * sins;
		q = n + 3;
		t = 4;
		p--;
		while (p > 0) 
		{
			/* Computes the Tni ratio */
			B = 1.0l - 1.0l / p;
			C = 1.0l + 1.0l / q;
			Tni = B * C;
			Pi0 *= Tni;
			Pi1 *= Tni;
			Pi0 += cosn[t];
			Pi1 += t * sinn[t];
			t += 2;
			p -= 2;
			q += 2;
		}
		coeff = root[2 * n + 1];
		Pn0[n] = (Pi0 * pnn + pn0) * coeff;
		Pn1[n] = Pi1 * pnn * sqr2 * coeff / rooti[n];
	}
	
	/* Computation of the Legendre Functions for m > 1 */
	printf("\n%-10s%-12s%1s\n", "n", "m", "Pnm");
	n = 0;
	printf(format, n, 0, Pn0[0]);
	
	n = 1;
	Pn_odd[0] = Pn0[n];
	Pn_odd[1] = Pn1[n];
	printf(format, n, 0, Pn_odd[0]);
	printf(format, n, 1, Pn_odd[1]);
	
	n = 2;
	Pn_even[0] = Pn0[n];
	Pn_even[1] = Pn1[n];
	Pn_even[2] = sqrt(3.0l) * (sqrt(5.0l) - Pn0[2]) / 3.0l;
	printf(format, n, 0, Pn_even[0]);
	printf(format, n, 1, Pn_even[1]);
	printf(format, n, 2, Pn_even[2]);
	
	n = 3;
	Pn_odd[0] = Pn0[n];
	Pn_odd[1] = Pn1[n];
	Pn_odd[2] = sqrt(15.0l) * (sqrt(21.0l) * Pn0[1] / 3.0l - Pn0[3]) / 5.0l;
	Pn_odd[3] = sqrt(15.0l) * (sqrt(14.0l) * Pn1[1] - Pn1[3]) / 15.0l;
	printf(format, n, 0, Pn_odd[0]);
	printf(format, n, 1, Pn_odd[1]);
	printf(format, n, 2, Pn_odd[2]);
	printf(format, n, 3, Pn_odd[3]);
	
	for (n = 4; n <= N; n++)
	{
		n2 = 2 * n;
		np1 = n + 1;
		nm1 = n - 1;
		nm2 = n - 2;
		nm3 = n - 3;
		coefn = root[n2 + 1] / root[n2 - 3];
		
		if(n % 2) /*========== n is odd ==========*/
		{
			/*<========= m is odd =========>*/
			Pn[1] = Pn1[n];
			for (m = 3; m <= nm2; m += 2)
			{
				z = rooti[nm1 + m];
				anm = coefn * rooti[nm1 - m] / z;
				bnm = coefn * rooti[nm3 + m] / z;
				cnm = rooti[np1 - m] / z;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_odd[m - 2];
				Pn_2_m = Pn_odd[m];
				Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
				Pn[m] = Pnm;
			}
			
			Pn_2_m_2 = Pn_odd[nm2];
			Pn_m_2 = Pn[nm2]; 
			cnn = 1.0l / root[n] / root[n2 - 1]; /* cnn coefficient */
			bnn = root[n2 + 1] * root[nm1] * cnn; /* bnn coefficient */
			Pnm = bnn * Pn_2_m_2 - cnn * Pn_m_2; /* Computes the Pnn Legendre function */
			Pn[n] = Pnm;
			
			/*<========= m is even =========>*/
			
			Pn[0] = Pn0[n];
			m = 2;
			z = rooti[nm1 + m];
			anm = coefn * rooti[nm1 - m] / z;
			bnm = sqr2 * coefn * rooti[nm3 + m] / z;
			cnm = sqr2 * rooti[np1 - m] / z;
			Pn_m_2 = Pn[m - 2];
			Pn_2_m_2 = Pn_odd[m - 2];
			Pn_2_m = Pn_odd[m];
			Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
			Pn[m] = Pnm;
			for (m = 4; m <= nm2; m += 2)
			{
				z = rooti[nm1 + m];
				anm = coefn * rooti[nm1 - m] / z;
				bnm = coefn * rooti[nm3 + m] / z;
				cnm = rooti[np1 - m] / z;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_odd[m - 2];
				Pn_2_m = Pn_odd[m];
				Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
				Pn[m] = Pnm;
			}
			
			Pn_2_m_2 = Pn_odd[nm3];
			Pn_m_2 = Pn[nm3];
			cnn_1 = sqr3 / root[nm1] / root[n2 - 1]; /* cn,n-1 coefficient */
			bnn_1 = cnn_1 * root[n2 + 1] * root[nm2] / sqr3; /* bn,n-1 coefficient*/
			Pnm = bnn_1 * Pn_2_m_2 - cnn_1 * Pn_m_2; /* Computes the Pn,n-1 Legendre function */
			Pn[nm1] = Pnm;
			
			for (m = 0; m <= n; m++)
			{
				Pnm = Pn[m];
				Pn_odd[m] = Pnm;
				printf(format, n, m, Pnm); /* Printing the odd Legendre Functions */
			}
		}
		else /*========== n is even ==========*/
		{
			/*<========= m is even =========>*/
			
			Pn[0] = Pn0[n];
			m = 2;
			z = rooti[nm1 + m];
			anm = coefn * rooti[nm1 - m] / z;
			bnm = sqr2 * coefn * rooti[nm3 + m] / z;
			cnm = sqr2 * rooti[np1 - m] / z;
			Pn_m_2 = Pn[m - 2];
			Pn_2_m_2 = Pn_even[m - 2];
			Pn_2_m = Pn_even[m];
			Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
			Pn[m] = Pnm;
			
			for (m = 4; m <= nm2; m += 2)
			{
				z = rooti[nm1 + m];
				anm = coefn * rooti[nm1 - m] / z;
				bnm = coefn * rooti[nm3 + m] / z;
				cnm = rooti[np1 - m] / z;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_even[m - 2];
				Pn_2_m = Pn_even[m];
				Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
				Pn[m] = Pnm;
			}
			
			Pn_2_m_2 = Pn_even[nm2];
			Pn_m_2 = Pn[nm2]; 
			cnn = 1.0l / root[n] / root[n2 - 1]; /* cnn coefficient */
			bnn = root[n2 + 1] * root[nm1] * cnn; /* bnn coefficient */
			Pnm = bnn * Pn_2_m_2 - cnn * Pn_m_2; /* Computes the Pnn Legendre function */
			Pn[n] = Pnm;
			
			/*<========= m is odd =========>*/
			
			Pn[1] = Pn1[n]; 
			for (m = 3; m <= nm2; m += 2)
			{
				z = rooti[nm1 + m];
				anm = coefn * rooti[nm1 - m] / z;
				bnm = coefn * rooti[nm3 + m] / z;
				cnm = rooti[np1 - m] / z;
				Pn_m_2 = Pn[m - 2];
				Pn_2_m_2 = Pn_even[m - 2];
				Pn_2_m = Pn_even[m];
				Pnm = anm * Pn_2_m + bnm * Pn_2_m_2 - cnm * Pn_m_2;
				Pn[m] = Pnm;
			}
			
			Pn_2_m_2 = Pn_even[nm3];
			Pn_m_2 = Pn[nm3];
			cnn_1 = sqr3 / root[nm1] / root[n2 - 1]; /* cn,n-1 coefficient */
			bnn_1 = cnn_1 * root[n2 + 1] * root[nm2] / sqr3; /* bn,n-1 coefficient*/
			Pnm = bnn_1 * Pn_2_m_2 - cnn_1 * Pn_m_2; /* Computes the Pn,n-1 Legendre function */
			Pn[nm1] = Pnm;
			
			for (m = 0; m <= n; m++)
			{
				Pnm = Pn[m];
				Pn_even[m] = Pnm;
				printf(format, n, m, Pnm); /* Printing the even Legendre Functions */
			}
		}
	}
	free(root);
	free(rooti);
	free(cosn);
	free(sinn);
	free(Pn0);
	free(Pn1);
	free(Pn_odd);
	free(Pn_even);
	free(Pn);
	return 0;
}
