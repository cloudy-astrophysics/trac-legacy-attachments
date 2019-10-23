/* This file is part of Cloudy and is copyright (C)1978-2014 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolEvaluate main routine to call others, to evaluate total cooling */
#include "cddefines.h"
#include "taulines.h"
#include "wind.h"
#include "coolheavy.h"
#include "radius.h"
#include "conv.h"
#include "h2.h"
#include "rt.h"
#include "opacity.h"
#include "ionbal.h"
#include "trace.h"
#include "dynamics.h"
#include "grainvar.h"
#include "atmdat.h"
#include "atoms.h"
#include "called.h"
#include "hmi.h"
#include "numderiv.h"
#include "magnetic.h"
#include "phycon.h"
#include "hyperfine.h"
#include "iso.h"
#include "thermal.h"
#include "cooling.h"
#include "pressure.h"
#include "mole.h"
#include "rfield.h"
#include "doppvel.h"
#include "freebound.h"
#include "dense.h"
#include "species.h"
#include "atmdat_gaunt.h"

/*fndneg search cooling array to find negative values */
STATIC void fndneg(void);

/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, double dc);

STATIC void CoolHyperfine (void);

/* set true to debug derivative of heating and cooling */
static const bool PRT_DERIV = false;



void CoolEvaluate(double *tot)
{
	static long int nhit = 0, 
	  nzSave=0;

	double 
	  qn, 
	  rothi=-SMALLFLOAT, 
	  rotlow=-SMALLFLOAT, 
	  x;

	static double oltcool=0., 
	  oldtemp=0.;

	long int coolnum, coolcal;

	DEBUG_ENTRY( "CoolEvaluate()" );

	/* returns tot, the total cooling,
	 * and dc, the derivative of the cooling */

	if( trace.lgTrace )
		fprintf( ioQQQ, "   COOLR TE:%.4e zone %li %li Cool:%.4e Heat:%.4e eden:%.4e edenTrue:%.4e\n", 
		phycon.te, 
		nzone, conv.nPres2Ioniz ,
		thermal.ctot , thermal.htot,dense.eden,dense.EdenTrue );

	/* must call TempChange since ionization has changed, there are some
	 * terms that affect collision rates (H0 term in electron collision) */
	TempChange(phycon.te , false);

	/* now zero out the cooling stack */
	CoolZero();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  0 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);
	if( gv.lgGrainPhysicsOn )
	{
		/* grain heating and cooling */
		/* grain recombination cooling, evaluated elsewhere
		* can either heat or cool the gas, do cooling here */
		CoolAdd("dust",0,MAX2(0.,gv.GasCoolColl));

		/* grain cooling proportional to temperature ^3/2 */
		thermal.dCooldT += MAX2(0.,gv.GasCoolColl)*3./(2.*phycon.te);

		/* these are the various heat agents from grains */
		/* options to force gas heating or cooling by grains to zero - for tests only ! */
		if( gv.lgDustOn() && gv.lgDHetOn )
		{
			/* rate dust heats gas by photoelectric effect */
			thermal.setHeating(0,13, gv.GasHeatPhotoEl);

			/* if grains hotter than gas then collisions with gas act
			* to heat the gas, add this in here
			* a symmetric statement appears in COOLR, where cooling is added on */
			thermal.setHeating(0,14, MAX2(0.,-gv.GasCoolColl) );

			/* this is gas heating due to thermionic emissions */
			thermal.setHeating(0,25, gv.GasHeatTherm );
		}
		else
		{
			thermal.setHeating(0,13,0.);
			thermal.setHeating(0,14,0.);
			thermal.setHeating(0,25,0.);
		}
	}
	else if( gv.lgBakesPAH_heat )
	{
		/* >>chng 06 jul 21, option to include Bakes PAH hack with grain physics off,
		 * needed to test dynamics models */
		thermal.setHeating(0,13,gv.GasHeatPhotoEl);
	}

	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  1 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* molecular molecules molecule cooling */
	if( mole_global.lgNoMole )
	{
		/* this branch - do not include molecules */
		hmi.hmicol = 0.;
		CoolHeavy.brems_cool_hminus = 0.;
		/* line cooling within simple H2 molecule - zero when big used */
		CoolHeavy.h2line = 0.;
		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = 0.;
		CoolHeavy.HD = 0.;

		/* thermal.heating(0,8) is heating due to collisions within X of H2 */
		thermal.setHeating(0,8,0.);
		/* thermal.heating(0,15) is H minus heating*/
		thermal.setHeating(0,15,0.);
		/* thermal.heating(0,16) is H2+ heating */
		thermal.setHeating(0,16,0.);
		hmi.HeatH2Dish_used = 0.;
		hmi.HeatH2Dexc_used = 0.;
		hmi.deriv_HeatH2Dexc_used = 0.;
	}

	else
	{
		/* save various molecular heating/cooling agent */
		thermal.setHeating(0,15, hmi.hmihet);
		thermal.setHeating(0,16, hmi.h2plus_heat);
		/* now get heating from H2 molecule, either simple or from big one */
		if( h2.lgEnabled  && hmi.lgH2_Thermal_BigH2 )
		{
			if( h2.lgEvaluated )
			{
				/* these are explicitly from big H2 molecule,
				 * first is heating due to radiative pump of excited states, followed by
				 * radiative decay into continuum of X, followed by dissociation of molecule
				 * with kinetic energy, typically 0.25 - 0.5 eV per event */
				hmi.HeatH2Dish_used = h2.HeatDiss;
				hmi.HeatH2Dexc_used = h2.HeatDexc;
				if (0)
					fprintf(ioQQQ,"DEBUG big %.2f\t%.5e\t%.2e\t%.2e\t%.2e\n", 
							  fnzone , phycon.te, hmi.HeatH2Dexc_used,
							  hmi.H2_total, dense.gas_phase[ipHYDROGEN] );
				/* negative sign because right term is really deriv of heating,
				 * but will be used below as deriv of cooling */
				hmi.deriv_HeatH2Dexc_used = -h2.HeatDexc_deriv;
			}
			else
			{
				hmi.HeatH2Dish_used = 0;
				hmi.HeatH2Dexc_used = 0;
				hmi.deriv_HeatH2Dexc_used = 0;
			}
		}

		else if( hmi.chH2_small_model_type == 'T' )
		{
			/* TH85 dissociation heating */
			/* these come from approximations in TH85, see comments above */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_TH85;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_TH85;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_TH85;
		}
		else if( hmi.chH2_small_model_type == 'H' )
		{
			/* Burton et al. 1990 */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BHT90;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BHT90;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BHT90;
		}
		else if( hmi.chH2_small_model_type == 'B')
		{
			/* Bertoldi & Draine */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BD96;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BD96;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BD96;
		}
		else if(hmi.chH2_small_model_type == 'E')
		{
			/* this is the default when small H2 used */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_ELWERT;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_ELWERT;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_ELWERT;
		}
		else
			TotalInsanity();

		/* heating due to photodissociation heating */
		thermal.setHeating(0,17,hmi.HeatH2Dish_used);

		/* heating due to continuum photodissociation */
		thermal.setHeating(0,28,0.);
		{
			double heat = 0.;
			for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
			{
				if( (*diatom)->lgEnabled && mole_global.lgStancil )
				{
					heat += (*diatom)->Cont_Diss_Heat_Rate();
				}
			}
			thermal.setHeating(0,28,MAX2( 0., heat ));
			CoolAdd("H2cD",0,MAX2(0.,-heat));
		}

		/* heating (usually cooling in big H2) due to collisions within X */
		/* add to heating is net heating is positive */
		thermal.setHeating(0,8, MAX2(0.,hmi.HeatH2Dexc_used) +
				// HD heating, cooling is CoolHeavy.HD
				MAX2(0.,hd.HeatDexc) + MAX2(0.,hd.HeatDiss));

		/* add to cooling if net heating is negative */
		CoolAdd("H2cX",0,MAX2(0.,-hmi.HeatH2Dexc_used));
		/*fprintf(ioQQQ,"DEBUG coolh2\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
			fnzone, phycon.te, dense.eden, hmi.H2_total, thermal.ctot, -hmi.HeatH2Dexc_used );*/
		/* add to net derivative */
		/*thermal.dCooldT += MAX2(0.,-hmi.HeatH2Dexc_used)* ( 30172. * thermal.tsq1 - thermal.halfte );*/
		/* >>chng 04 jan 25, check sign to prevent cooling from entering here,
		 * also enter neg sign since going into cooling stack (bug), in heatsum
		 * same term adds to deriv of heating */
		if( hmi.HeatH2Dexc_used < 0. )
			thermal.dCooldT -= hmi.deriv_HeatH2Dexc_used;

		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = (realnum)(MAX2((2.325*phycon.te-1875.)*1e-20,0.)*
		  dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHYDROGEN][1]*1.66e-11);

		if( h2.lgEnabled )
		{
			/* this is simplified approximation to H2 rotation cooling,
			 * big molecule goes this far better */
			CoolHeavy.h2line = 0.;
		}
		else
		{
			{
				// cooling from Glover & Abel, MNRAS, 388, 1627, Table 1-7
				GA1 low T
				-24.330855,4.4404496,-4.0460989,-1.1390725,9.8094223,8.6273872,0,0 ,
				-24.329086 , 4.6105087 , -3.9505350 , 12.363818 , -32.403165 , 48.853562 , -38.542008 , 12.066770,

				GAT2
				-24.216387 , 3.3237480 , -11.642384 , -35.553366 , -35.105689 , -10.922078 , 0 , 0
				-24.216387 , 4.2046488 , -1.3155285 , -1.6552763 , 4.1780102 , -0.56949697 , -3.3824407 , 1.09040270

				GAT3
				para H2
				-23.889798 , 1.8550774 , -0.55593388 , 0.28429361 , -0.20581113 , 0.13112378
				Ortho H2
				-23.748534 , 1.76676480 , -0.58634325 , 0.31074159 , -0.17455629 , 0.18530758

				GAT4
				Para-H2
				-24.126177 , 2.3258217 , -1.0082491 , 0.54823768 , -0.33679759 , 0.20771406
				Ortho-H2
				-24.020047 , 2.2687566 , -1.0200304 , 0.83561432 , -0.40772247 , 0.096025713

				GAT5
				Para-H2
				-23.489029 , 1.8210825 , -0.59110559 , 0.42280623 , -0.30171138 , 0.12872839
				Ortho-H2
				-23.7749 , 2.40654 , -1.23449 , 0.739874 , -0.258940 , 0.120573

				GAT6
				Para-H2
				-21.757160 , 1.3998367 -0.37209530 0.061554519 -0.37238286 0.23314157
				Ortho-H2
				-21.706641 1.3901283 -0.34993699 0.075402398 -0.23170723 0.068938876

				GAT7
				Para-H2T² 103K
				-22.817869 , 0.95653474 , 0.79283462 , 0.56811779 , 0.27895033 , 0.056049813
				Para-H2T > 103K
				-22.817869 , 0.66916141 , 7.1191428 , -11.176835 , 7.0467275 , -1.6471816
				Ortho-H2
				-21.703215 , 0.76059565 , 0.50644890 , 0.050371349 , -0.10372467 , -0.035709409

			}
			/*  rate for rotation lines from 
			*  >>refer	h2	cool	Lepp, S., & Shull, J.M. 1983, ApJ, 270, 578 */
			x = phycon.alogte - 4.;
			if( phycon.te > 1087. )
			{
				rothi = 3.90e-19*sexp(6118./phycon.te);
			}
			else
			{
				rothi = pow(10.,-19.24 + 0.474*x - 1.247*x*x);
			}

			/*  low density rotation cooling */
			/*&qn = pow(MAX2(findspecieslocal("H2")->den,1e-37),0.77) + 1.2*pow(MAX2(dense.xIonDense[ipHYDROGEN][0],1e-37),0.77);*/
			qn = pow(MAX2(hmi.H2_total,1e-37),0.77) + 1.2*pow(MAX2(dense.xIonDense[ipHYDROGEN][0],1e-37),0.77);
			/* these are equations 11 from LS83 */
			if( phycon.te > 4031. )
			{
				rotlow = 1.38e-22*sexp(9243./phycon.te)*qn;
			}
			else
			{
				rotlow = pow(10.,-22.90 - 0.553*x - 1.148*x*x)*qn;
			}

			CoolHeavy.h2line = 0.;
			if( rotlow > 0. )
				CoolHeavy.h2line += hmi.H2_total*rothi/(1. + rothi/rotlow);
			/* \todo 1 add this from LS83 or (better yet) update to another form.  See Galli & Palla 1998, A5-7. */
			//if( viblow > 0. )
			//	CoolHeavy.h2line += hmi.H2_total*vibhi/(1. + vibhi/viblow);
		}

		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && nzone>187&& iteration > 1/**/)
			{
				fprintf(ioQQQ,"h2coolbug\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					phycon.te, 
					CoolHeavy.h2line, 
					hmi.H2_total, 
					findspecieslocal("H-")->den, 
					hmi.HMinus_photo_rate,
					rothi,
					rotlow );
			}
		}

		if( hd.lgEnabled )
		{
			// heating was thermal.setHeating(0,8 above, with H2
			CoolHeavy.HD = -MIN2(0.,hd.HeatDexc) - MIN2(0.,hd.HeatDiss);
		}
		else
		{
			/* >>chng 22 aug 13, use Flower et al. (2000, MNRAS, 314, 753)
			 * HD cooling function
			 * http://ccp7.dur.ac.uk/cooling_by_HD/node3.html */
#if 0
			factor = 0.25*pow2( log10(double(dense.gas_phase[ipHYDROGEN])))+
					(0.283978*log10(double(dense.gas_phase[ipHYDROGEN]))-1.27333)*
					(sin(2.03657*log10(phycon.te)+4.63258))-2.08189*
					log10(double(dense.gas_phase[ipHYDROGEN]))+4.66288;

			CoolHeavy.HD = hmi.HD_total*pow(10., (0.5*log10(double(dense.gas_phase[ipHYDROGEN]))) +
					(-26.2982*pow(log10(phycon.te), -0.215807)) - pow(factor, 0.5));
#endif

			double aa=-26.2982, bb=-0.215807, omeg=2.03657, phi=4.63258,
					c1=0.283978, c2=-1.27333, d1=-2.08189, d2=4.66288;

			// the sum of hydrogen nuclei in all forms
			double y = log10(dense.gas_phase[ipHYDROGEN]);
			double x = phycon.alogte;

			double w = 0.5 * y + aa * pow(x,bb)
				- sqrt( 0.25 * y*y
				+ (c1*y+c2) * sin(omeg*x+phi) + (d1*y+d2) );

			CoolHeavy.HD = hmi.HD_total * pow(10.,w);

		}
	}

	fixit("test and enable chemical heating by default");
#if 0
	double chemical_heating = mole.chem_heat();	
	thermal.setHeating(0,29, MAX2(0.,chemical_heating) );
	/* add to cooling if net heating is negative */
	CoolAdd("Chem",0,MAX2(0.,-chemical_heating));
#endif

	/* cooling due to charge transfer ionization / recombination */
	CoolAdd("CT C" , 0. , thermal.char_tran_cool );

	/*  H- FB; H + e -> H- + hnu */
	/*  H- FF is in with H ff */
	CoolAdd("H-fb",0,hmi.hmicol);

	/* >>chng 96 nov 15, fac of 2 in deriv to help convergence in very dense
	 * models where H- is important, this takes change in eden into
	 * partial account */
	thermal.dCooldT += 2.*hmi.hmicol*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  2 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	CoolAdd("H2ln",0,CoolHeavy.h2line);
	/* >>chng 00 oct 21, added coef of 3.5, sign had been wrong */
	/*thermal.dCooldT += CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 mar 17, change 3.5 to 15 as per behavior in primal.in */
	/*thermal.dCooldT += 3.5*CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 may 18, from 15 to 30 as per behavior in primal.in - overshoots happen */
	/*thermal.dCooldT += 15.*CoolHeavy.h2line*phycon.teinv;*/
	/*>>chng 03 oct 03, from 30 to 3, no more overshoots in primalin */
	/*thermal.dCooldT += 30.*CoolHeavy.h2line*phycon.teinv;*/
	thermal.dCooldT += 3.0*CoolHeavy.h2line*phycon.teinv;

	{
		/* problems with H2 cooling */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /*&& nzone>300 && iteration > 1*/)
		{
			fprintf(ioQQQ,"CoolEvaluate debuggg\t%.2f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",
			fnzone, 
			phycon.te, 
			hmi.H2_total , 
			CoolHeavy.h2line,
			findspecieslocal("H-")->den , 
			dense.eden);
		}
	}

	CoolAdd("HDro",0,CoolHeavy.HD);
	thermal.dCooldT += CoolHeavy.HD*phycon.teinv;

	CoolAdd("H2+ ",0,CoolHeavy.H2PlsCool);
	thermal.dCooldT += CoolHeavy.H2PlsCool*phycon.teinv;

	/* heating due to three-body, will be incremented in iso_cool*/
	thermal.setHeating(0,3,0.);
	/* heating due to hydrogen lines */
	thermal.setHeating(0,23,0.);
	/* heating due to photoionization of all excited states of hydrogen species */
	thermal.setHeating(0,1,0.);

	/* isoelectronic species cooling, mainly lines, and level ionization */
	for( long int ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always call iso_cool since we must zero those variables
			* that would have been set had the species been present */
			iso_cool( ipISO , nelem );
		}
	}

	/* free-free free free brems emission for all ions */
	/* highest frequency where we have non-zero Boltzmann factors */
	long limit = min( rfield.ipMaxBolt, rfield.nflux );

	CoolHeavy.brems_cool_h = 0.;
	CoolHeavy.brems_cool_hminus = 0.;
	CoolHeavy.brems_cool_he = 0.;
	CoolHeavy.brems_cool_metals = 0.;
	CoolHeavy.brems_heat_total = 0.;

	if( CoolHeavy.lgFreeOn )
	{
		ASSERT(rfield.ipEnergyBremsThin < rfield.nflux_with_check);
		ASSERT(limit < rfield.nflux_with_check);

		t_brems_den sum;
		t_gaunt::Inst().brems_sum_ions(sum);
		double bfac = dense.eden*1.032e-11/phycon.sqrte*EN1RYD;

		/* do hydrogen first, before main loop since want to break out as separate
		 * coolant, and what to add on H- brems */
		CoolHeavy.brems_cool_h = sum.den_Hp*bfac*t_gaunt::Inst().brems_cool( 1, phycon.te );
		CoolHeavy.brems_cool_hminus = sum.den_Hm*bfac*t_gaunt::Inst().brems_cool( -1, phycon.te );

		/* now do helium, both He+ and He++ */
		CoolHeavy.brems_cool_he = sum.den_Hep*bfac*t_gaunt::Inst().brems_cool( 1, phycon.te ) +
			sum.den_Hepp*bfac*t_gaunt::Inst().brems_cool( 2, phycon.te );

		/* heavy elements */
		CoolHeavy.brems_cool_metals = 0.;
		for( long ion=1; ion < LIMELM+1; ++ion )
			if( sum.den_ion[ion] > 0. )
				CoolHeavy.brems_cool_metals +=
					sum.den_ion[ion]*bfac*t_gaunt::Inst().brems_cool( ion, phycon.te );

		/* ipEnergyBremsThin is index to energy where gas becomes optically thin to brems,
		 * so this loop is over optically thin frequencies 
		 * do not include optically thick part as net emission since self absorbed */
		CoolHeavy.brems_heat_total = 0.;
		for( long int i=rfield.ipEnergyBremsThin; i < limit; i++ )
		{
			/* the total heating due to bremsstrahlung */
			CoolHeavy.brems_heat_total +=
				opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i]*EN1RYD;
		}

		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>60 /*&& iteration > 1*/)
		{
			double sumfield = 0., sumtot=0., sum1=0., sum2=0.;
			for( long int i=rfield.ipEnergyBremsThin; i<limit;  i++ )
			{
				sumtot += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i];
				sumfield += rfield.flux[0][i]*rfield.anu[i];
				sum1 += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu[i];
				sum2 += opac.FreeFreeOpacity[i]*rfield.flux[0][i];
			}
			fprintf(ioQQQ,"DEBUG brems heat\t%.2f\t%.3e\t%.3e\t%.3e\t%e\t%.3e\t%.3e\n",
				fnzone,
				CoolHeavy.brems_heat_total,
				sumtot/SDIV(sumfield) ,
				sum1/SDIV(sum2),
				phycon.te , 
				t_gaunt::Inst().gauntff(1,phycon.te,rfield.anu[1218]),
				opac.FreeFreeOpacity[1218]);
		}
	}

	/* these two terms are both large, nearly canceling, near lte */
	CoolHeavy.brems_cool_net = 
		CoolHeavy.brems_cool_h + 
		CoolHeavy.brems_cool_he + 
		CoolHeavy.brems_cool_hminus + 
		CoolHeavy.brems_cool_metals - 
		CoolHeavy.brems_heat_total;
	/*fprintf(ioQQQ,"DEBUG brems\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
		fnzone,
		phycon.te,
		CoolHeavy.brems_cool_net,
		CoolHeavy.brems_cool_h ,
		CoolHeavy.brems_cool_he ,
		CoolHeavy.brems_cool_hminus,
		CoolHeavy.brems_cool_metals ,
		CoolHeavy.brems_heat_total);*/

	/* net free free brems cooling, count as cooling if positive */
	CoolAdd( "FF c" , 0, MAX2(0.,CoolHeavy.brems_cool_net) );

	/* now stuff into heating array if negative */
	thermal.setHeating(0,11, MAX2(0.,-CoolHeavy.brems_cool_net) );

	/* >>chng 96 oct 30, from HFFNet to just FreeFreeCool,
	 * since HeatSum picks up CoolHeavy.brems_heat_total */
	thermal.dCooldT += CoolHeavy.brems_cool_h*thermal.halfte;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  3 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* >>chng 02 jun 21, net cooling already includes this */
	/* end of brems cooling */

	/* heavy element recombination cooling, do not count hydrogenic since
	 * already done above, also helium singlets have been done */
	/* >>chng 02 jul 21, put in charge dependent rec term */
	CoolHeavy.heavfb = 0.;
	for( long int nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* note that detailed iso seq atoms are done in iso_cool */
			long limit_lo = MAX2( 1 , dense.IonLow[nelem] );
			long limit_hi = MIN2( nelem-NISO+1, dense.IonHigh[nelem] );
			for( long int ion=limit_lo; ion<=limit_hi; ++ion )
			{
				/* factor of 0.9 is roughly correct for nebular conditions, see
				 * >>refer	H	rec cooling	LaMothe, J., & Ferland, G.J., 2001, PASP, 113, 165 */
				/* note that ionbal.RR_rate_coef_used is the rate coef, cm3 s-1, needs eden */
				/* >>chng 02 nov 07, move rec arrays around, this now has ONLY rad rec,
				 * previously had included grain rec and three body */
				/* recombination cooling for iso-seq atoms are done in iso_cool */
				double one = dense.xIonDense[nelem][ion] * ionbal.RR_rate_coef_used[nelem][ion-1]*
					dense.eden * phycon.te * BOLTZMANN;
				/*fprintf(ioQQQ,"debugggfb\t%li\t%li\t%.3e\t%.3e\t%.3e\n", nelem, ion, one, 
					dense.xIonDense[nelem][ion] , ionbal.RR_rate_coef_used[nelem][ion]);*/
				CoolHeavy.heavfb += one;
				thermal.elementcool[nelem] += one;
			}
		}
	}

	/*fprintf(ioQQQ,"debuggg hvFB\t%i\t%.2f\t%.2e\t%.2e\n",iteration, fnzone,CoolHeavy.heavfb, dense.eden);*/

	CoolAdd("hvFB",0,CoolHeavy.heavfb);
	thermal.dCooldT += CoolHeavy.heavfb*.113*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  4 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* electron-electron brems from
	 * >>refer	ee	brems	Stepney and Guilbert, MNRAS 204, 1269 (1983)
	 * 2013 Mar Wang Ye fit to their table, roughly 3-10x larger at 1e9 - 1e10 K
	 * thank previous simple value */
	CoolHeavy.eebrm = ((SIGMA_THOMSON * FINE_STRUCTURE)/(ELECTRON_MASS * SPEEDLIGHT))*
			POW2(dense.eden*BOLTZMANN)*pow(phycon.te, 1.708064)*2633.23;

	/* >>chng 97 mar 12, added deriv */
	thermal.dCooldT += CoolHeavy.eebrm*thermal.halfte;
	CoolAdd("eeff",0,CoolHeavy.eebrm);

	/* add advective heating and cooling */
	/* this is cooling due to loss of matter from this region */
	CoolAdd("adve",0,dynamics.Cool() );
	/* >>chng02 dec 04, rm factor of 8 in front of dCooldT */
	thermal.dCooldT += dynamics.dCooldT();
	/* local heating due to matter moving into this location */
	thermal.setHeating(1,5, dynamics.Heat() );
	thermal.dHeatdT += dynamics.dHeatdT;

	/* total Compton cooling */
	CoolHeavy.tccool = rfield.cmcool*phycon.te;
	CoolAdd("Comp",0,CoolHeavy.tccool);
	thermal.dCooldT += rfield.cmcool;

	/* option for "extra" cooling, expressed as power-law in temperature, these
	 * are set with the CEXTRA command */
	if( thermal.lgCExtraOn )
	{
		CoolHeavy.cextxx = 
			(realnum)(thermal.CoolExtra*pow(phycon.te/1e4,(double)thermal.cextpw));
	}
	else
	{
		CoolHeavy.cextxx = 0.;
	}
	CoolAdd("Extr",0,CoolHeavy.cextxx);

	/* cooling due to wind expansion, only for winds expansion cooling */
	if( wind.lgBallistic() )
	{
		realnum dDensityDT = -(realnum)(wind.AccelTotalOutward/wind.windv + 2.*wind.windv/
		  radius.Radius);
		CoolHeavy.expans = -2.5*pressure.PresGasCurr*dDensityDT;
	}
	else if( dynamics.lgTimeDependentStatic && 
				iteration > dynamics.n_initial_relax)
	{
		realnum dens = scalingDensity();
		realnum dDensityDT = 
			(realnum)((dens-dynamics.Upstream_density)/
			(dynamics.timestep*0.5*(dens+dynamics.Upstream_density)));
		// pdV work term
		CoolHeavy.expans = -pressure.PresGasCurr*dDensityDT;
	}
	else
	{
		CoolHeavy.expans = 0.;
	}
	CoolAdd("Expn",0,CoolHeavy.expans);
	thermal.dCooldT += CoolHeavy.expans/phycon.te;

	/* cyclotron cooling */
	/* coef is 4/3 /8pi /c * vtherm(elec) */
	CoolHeavy.cyntrn = 4.5433e-25f*magnetic.pressure*PI8*dense.eden*phycon.te;
	CoolAdd("Cycl",0,CoolHeavy.cyntrn);
	thermal.dCooldT += CoolHeavy.cyntrn/phycon.te;

	/* heavy element collisional ionization
	 * derivative should be zero since increased coll ion rate
	 * decreases neutral fraction by proportional amount */
	CoolAdd("Hvin",0,CoolHeavy.colmet);


	double xIonDenseSave[LIMELM][LIMELM+1];
	if( atmdat.lgChiantiOn ||atmdat.lgStoutOn)
	{
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{		
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				xIonDenseSave[nelem][ion] = dense.xIonDense[nelem][ion];
				// zero abundance of species if we are using Chianti for this ion
				if( dense.lgIonChiantiOn[nelem][ion] || dense.lgIonStoutOn[nelem][ion] )
					dense.xIonDense[nelem][ion] = 0.;
			}
		}
	}

	(*TauDummy).Zero();
	(*(*TauDummy).Hi()).g() = 0.;
	(*(*TauDummy).Lo()).g() = 0.;
	(*(*TauDummy).Hi()).IonStg() = 0;
	(*(*TauDummy).Lo()).IonStg() = 0;
	(*(*TauDummy).Hi()).nelem() = 0;
	(*(*TauDummy).Lo()).nelem() = 0;
	(*TauDummy).Emis().Aul() = 0.;
	(*TauDummy).EnergyWN() = 0.;
	(*TauDummy).WLAng() = 0.;

	/* Iron cooling */
	coolnum = thermal.ncltot;
	CoolIron();
	for( coolcal = coolnum; coolcal < thermal.ncltot; coolcal++ )
		thermal.elementcool[ipIRON] += thermal.cooling[coolcal];
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 12 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	// reset abundances to original values, may have been set zero to protect against old cloudy lines
	if( atmdat.lgChiantiOn || atmdat.lgStoutOn)
	{
		// this clause, first reset abundances set to zero when Chianti included
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				dense.xIonDense[nelem][ion] = xIonDenseSave[nelem][ion];
			}
		}
	}

	/* opacity project lines Dima Verner added with g-bar approximation */
	coolnum = thermal.ncltot;
	CoolDima();
	for( int coolcal = coolnum; coolcal < thermal.ncltot; coolcal++ )
		thermal.dima += thermal.cooling[coolcal];

	/* do external database lines */
	dBaseTrim();
	dBaseUpdateCollCoeffs();
 	dBase_solve();

	/* Hyperfine line cooling */
	coolnum = thermal.ncltot;
	CoolHyperfine();
	for( coolcal = coolnum; coolcal < thermal.ncltot; coolcal++ )
		thermal.elementcool[ipHYDROGEN] += thermal.cooling[coolcal];
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  6 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Print number of levels for each species */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			static bool lgMustPrintHeader=true;
			if( lgMustPrintHeader )
			{
				lgMustPrintHeader = false;
				printf("DEBUG Levels\t%s",dBaseSpecies[0].chLabel );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf("\t%s",dBaseSpecies[ipSpecies].chLabel );
				}
				printf("\n" );
				printf("DEBUG Max\t%li" ,dBaseSpecies[0].numLevels_max );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf( "\t%li" ,dBaseSpecies[ipSpecies].numLevels_max );
				}
				printf("\n");
			}
			printf("DEBUG Local\t%li" ,dBaseSpecies[0].numLevels_local );
			for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
			{
				printf("\t%li" ,dBaseSpecies[ipSpecies].numLevels_local );
			}
			printf("\n");
		}
	}

	/* now add up all the coolants */
	CoolSum(tot);
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 14 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* negative cooling */
	if( *tot <= 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling is <=0, this is impossible.\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* bad derivative */
	if( thermal.dCooldT == 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling slope <=0, this is impossible.\n" );
		if( *tot > 0. && dense.gas_phase[ipHYDROGEN] < 1e-4 )
		{
			fprintf( ioQQQ, " Probably due to very low density.\n" );
		}
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fndstr(*tot,thermal.dCooldT);
	}

	/* lgTSetOn true for constant temperature model */
	if( (((((!thermal.lgTemperatureConstant) && *tot < 0.) && called.lgTalk) && 
	  !conv.lgSearch) && thermal.lgCNegChk) && nzone > 0 )
	{
		fprintf( ioQQQ, 
			" NOTE Negative cooling, zone %4ld, =%10.2e coola=%10.2e CHION=%10.2e Te=%10.2e\n", 
		  nzone, 
		  *tot, 
		  iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool, 
		  iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion, 
		  phycon.te );
		fndneg();
	}

	/* possibility of getting empirical cooling derivative
	 * normally false, set true with 'set numerical derivatives' command */
	if( NumDeriv.lgNumDeriv )
	{
		if( ((nzone > 2 && nzone == nzSave) && ! fp_equal( oldtemp, phycon.te )) && nhit > 4 )
		{
			/* hnit is number of tries on this zone - use to stop numerical problems
			 * do not evaluate numerical deriv until well into solution */
			thermal.dCooldT = (oltcool - *tot)/(oldtemp - phycon.te);
		}
		if( nzone != nzSave )
			nhit = 0;

		nzSave = nzone;
		nhit += 1;
		oltcool = *tot;
		oldtemp = phycon.te;
	}
	return;
}



STATIC void CoolHyperfine (void)
{
	static double	TeEvalCS = 0., TeEvalCS_21cm=0.;
	static double 	electron_rate_21cm,
		atomic_rate_21cm,
		proton_rate_21cm;


	for( int ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if ( ! dense.lgElmtOn[nelem] )
					continue;
			ionbal.ExcitationGround[nelem][nelem-ipISO] +=
					iso_sp[ipISO][nelem].fb[0].RateLevel2Cont;
		}
	}

	{
		enum {DEBUG_LOC=false};
		if ( DEBUG_LOC )
		{
			for (int nelem = 0; nelem < LIMELM; nelem++)
			{
				fprintf(ioQQQ, "ionbal.ExcitationGround: nelem = %d", nelem);
				for (int ion = 0; ion < nelem+1; ion++)
				{
					fprintf(ioQQQ,"\t%.6e", ionbal.ExcitationGround[nelem][ion]);
				}
				fprintf(ioQQQ, "\n");
			}
		}
	}

	/* evaluate H 21 cm spin changing collisions */
	if( !fp_equal(phycon.te,TeEvalCS_21cm) )
	{
		{
			/* this prints table of rates at points given in original data paper */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
#				define N21CM_TE	16
				int n;
				double teval[N21CM_TE]={2.,5.,10.,20.,50.,100.,200.,500.,1000.,
					2000.,3000.,5000.,7000.,10000.,15000.,20000.};
				for( n = 0; n<N21CM_TE; ++n )
				{
					fprintf(
						ioQQQ,"DEBUG 21 cm deex Te=\t%.2e\tH0=\t%.2e\tp=\t%.2e\te=\t%.2e\n",
						teval[n], 
						H21cm_H_atom( teval[n] ),
						H21cm_proton( teval[n] ),
						H21cm_electron( teval[n] ) );
				}
				cdEXIT(EXIT_FAILURE);
#				undef N21CM_TE
			}
		}
		/*only evaluate T dependent part when Te changes, but update
		 * density part below since densities may constantly change */
		atomic_rate_21cm = H21cm_H_atom( phycon.te );
		proton_rate_21cm = H21cm_proton( phycon.te );
		electron_rate_21cm = H21cm_electron( phycon.te );
		TeEvalCS_21cm = phycon.te;
	}
	/* H 21 cm emission/population,
	* cs will be sum of e cs and H cs converted from rate */
	double cs = (electron_rate_21cm * dense.eden + 
		atomic_rate_21cm * dense.xIonDense[ipHYDROGEN][0] +
		proton_rate_21cm * dense.xIonDense[ipHYDROGEN][1] ) *
		3./	dense.cdsqte;
	PutCS(  cs , HFLines[0] );

	/* do level pops for H 21 cm which is a special case since Lya pumping in included */
	RT_line_one_escape( HFLines[0], true,0.f, GetDopplerWidth(dense.AtomicWeight[(*HFLines[0].Hi()).nelem()-1]) );
	H21_cm_pops();

	hyperfine.cooling_total = HFLines[0].Coll().cool();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  5 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);


	if( !fp_equal(phycon.te,TeEvalCS) )
	{
		for( long int i=1; i < nHFLines; i++ )
		{
			cs = HyperfineCS( i );
			PutCS(  cs , HFLines[i] );
		}
		TeEvalCS = phycon.te;
	}


	/* now do level pops for all except 21 cm  */
	for( long int i=1; i < nHFLines; i++ )
	{
		/* remember current gas-phase abundance of this isotope */
		double save = dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1];

		/* bail if no abundance */
		if( save<=0. ) 
			continue;

		RT_line_one_escape( HFLines[i], true,0.f, GetDopplerWidth(dense.AtomicWeight[(*HFLines[i].Hi()).nelem()-1]) );

		/* set gas-phase abundance to total times isotope ratio */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] *= 
			hyperfine.HFLabundance[i];

		/* use the collision strength generated above and find pops and cooling */
		atom_level2( HFLines[i], true );

		/* put the correct gas-phase abundance back in the array */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] = save;

		/* find total cooling due to hyperfine structure lines */
		hyperfine.cooling_total += HFLines[i].Coll().cool();
	}

	return;
}



/*  */
#ifdef EPS
#	undef EPS
#endif
#define	EPS	0.01

/*fndneg search cooling array to find negative values */
STATIC void fndneg(void)
{
	long int i;
	double trig;

	DEBUG_ENTRY( "fndneg()" );

	trig = fabs(thermal.htot*EPS);
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( thermal.cooling[i] < 0. && fabs(thermal.cooling[i]) > trig )
		{
			fprintf( ioQQQ, " negative line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.cooling[i]/
			  thermal.htot );
		}

		if( thermal.heatnt[i] > trig )
		{
			fprintf( ioQQQ, " heating line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.heatnt[i]/
			  thermal.htot );
		}
	}
	return;
}

/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, 
  double dc)
{
	char chStrngLab[NCOLNT_LAB_LEN+1];
	long int i;
	realnum wl;
	double str, 
	  strong;

	DEBUG_ENTRY( "fndstr()" );

	strong = 0.;
	wl = -FLT_MAX;
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( fabs(thermal.cooling[i]) > strong )
		{
			/* this is the wavelength of the coolant, 0 for a continuum*/
			wl = thermal.collam[i];
			/* make sure labels are all valid*/
			/*>>chng 06 jun 06, bug fix, assert length was ==4, should be <=NCOLNT_LAB_LEN */
			ASSERT( strlen( thermal.chClntLab[i] ) <= NCOLNT_LAB_LEN );
			strcpy( chStrngLab, thermal.chClntLab[i] );
			strong = fabs(thermal.cooling[i]);
		}
	}

	str = strong;

	fprintf( ioQQQ, 
		"   fndstr cool: TE=%10.4e Ne=%10.4e C=%10.3e dC/dT=%10.2e ABS(%s %.1f)=%.2e nz=%ld\n", 
	  phycon.te, dense.eden, tot, dc, chStrngLab
	  , wl, str, nzone );

	/* option for extensive printout of lines */
	if( trace.lgCoolTr )
	{
		realnum ratio;

		/* flag all significant coolants, first zero out the array */
		coolpr(ioQQQ,thermal.chClntLab[0],1,0.,"ZERO");

		/* push all coolants onto the stack */
		for( i=0; i < thermal.ncltot; i++ )
		{
			/* usually positive, although can be neg for coolants that heats, 
			 * only do positive here */
			ratio = (realnum)(thermal.cooling[i]/thermal.ctot);
			if( ratio >= EPS )
			{
				/*>>chng 99 jan 27, only cal when ratio is significant */
				coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i], ratio,"DOIT");
			}
		}

		/* complete the printout for positive coolants */
		coolpr(ioQQQ,"DONE",1,0.,"DONE");

		/* now do cooling that was counted as a heat source if significant */
		if( thermal.heating(0,22)/thermal.ctot > 0.05 )
		{
			fprintf( ioQQQ, 
				"     All coolant heat greater than%6.2f%% of the total will be printed.\n", 
			  EPS*100. );

			coolpr(ioQQQ,"ZERO",1,0.,"ZERO");
			for( i=0; i < thermal.ncltot; i++ )
			{
				ratio = (realnum)(thermal.heatnt[i]/thermal.ctot);
				if( fabs(ratio) >=EPS )
				{
					coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i],
					  ratio,"DOIT");
				}
			}
			coolpr(ioQQQ,"DONE",1,0.,"DONE");
		}
	}
	return;
}

#undef	EPS
