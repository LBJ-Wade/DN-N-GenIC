#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


double  f_baryon;	/* Baryon fraction */
double	f_bnu;		/* Baryon + Massive Neutrino fraction */
double	f_cb;		/* Baryon + CDM fraction */
double	f_cdm;		/* CDM fraction */
double	f_hdm;		/* Massive Neutrino fraction */
double	obhh;		/* OmegaBaryon * HubbleParam^2 */
double	omega_curv;	/* = 1 - Omega - OmegaLambda */
double	omhh;		/* Omega * HubbleParam^2 */
double	onhh;		/* OmegaNu * HubbleParam^2 */
double	p_c;		/* The correction to the exponent before drag epoch */
double	p_cb;		/* The correction to the exponent after drag epoch */
double	theta_cmb;	/* The temperature of the CMB, in units of 2.7 K */


double b_k(double k)
/* Corrective multiplicative factor for free streaming scale */
{
   	theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */
	omega_curv = 1.0-Omega-OmegaLambda;
	omhh = Omega*HubbleParam*HubbleParam;
	obhh = OmegaBaryon*HubbleParam*HubbleParam;
	onhh = OmegaNu*HubbleParam*HubbleParam;
	f_baryon = OmegaBaryon/Omega;
	f_hdm = OmegaNu/Omega;
	f_cdm = 1.0-f_baryon-f_hdm;
	f_cb = f_cdm+f_baryon;
	f_bnu = f_baryon+f_hdm;


	/* Set up for the free-streaming & infall growth function */
	p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
	p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

	
	k*=(3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */	
	
	double	max_fs_correction;  /* Correction near maximal free streaming */
	double	qq;		/* Wavenumber rescaled by \Gamma */
	double	qq_nu;		/* Wavenumber compared to maximal free streaming */
	
	qq = k/(omhh*theta_cmb*theta_cmb);
	
	qq_nu = 3.92*qq*sqrt(NoOfDegenNu/f_hdm);
	
    max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(NoOfDegenNu,0.3+0.6*f_hdm)/
		(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
		
    return max_fs_correction;
}

double growth_cbnu(double k, double a)
/* Time dependent growth function for CDM + Baryons + Neutrinos: D_cbnu(z) */
{

	k*=(3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */	

	double	qq;		/* Wavenumber rescaled by \Gamma */
	double y_fs; 	/* The epoch of free-streaming for a given scale */
	
	qq = k/(omhh*theta_cmb*theta_cmb);
	y_fs = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*pow(NoOfDegenNu*qq/f_hdm,2);
	
	return pow(pow(f_cb,0.7/p_cb) + pow(growth(a)/(1 + y_fs),0.7),p_cb/0.7)*pow(growth(a),1-p_cb);
	
}

double gf_cbnu(double k,double a)
/* Time dependent growth factor for CDM + Baryons + Neutrinos
   D_cbnu(k,z) / D(z)
*/
{
	return growth_cbnu(k,a)/growth(a);
}
	
	

    
    
    
	

