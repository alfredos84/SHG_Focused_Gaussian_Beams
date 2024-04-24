/*---------------------------------------------------------------------------*/
// * This file contains all the parameters for the choosen nonlinear crystal.
// * By default the nonlinear crystal is MgO:sPPLT.
// * The properties such as crytal size, grid discretization, attenuation coefficients,
// * etc. are defined here as a global variables enabled to be used in device.
// * The funcion n(L, T) returns the MgO:sPPLT extraordinary refractive 
// * index from the Sellmeier Equation. 
// * Reference: Bruner et. al. "Temperature-dependent Sellmeier equation for 
// * the refractive index of stoichiometric lithium tantalate"
/*---------------------------------------------------------------------------*/


#ifndef _CRYSTAL
#define _CRYSTAL

#pragma once

// Define global constants
const real_t PI		= 3.14159265358979323846;	// pi
const real_t C		= 299792458*1E6/1E12;		// speed of ligth in vacuum [μm/ps]
const real_t EPS0	= 8.8541878128E-12*1E12/1E6;// vacuum pertivity [W.ps/V²μm]

std::string CrystalName = "MgO:sPPLT";

// Pump and signal Wavelengths [μm]
const real_t lp	= 1.064;
const real_t ls	= 0.5*lp;


// Crystal grid
const real_t LX	= 2000.0;		// Crystal width [μm]
const real_t LY	= 1000.0;		// Crystal heigth [μm]
const real_t LZ	= 30000.0;		// Crystal length [μm]
const real_t dx	= LX/(NX-1);		// x step [μm]
const real_t dy	= LY/(NY-1);		// y step [μm]
const real_t dz	= LZ/(NZ-1);		// z step [μm]



// Crystal properties
const real_t deff		= 11.00e-6;	// Eff. second-order susceptibility [um/V]
const real_t dQ			= (2.0/PI)*deff;	// Eff. second-order susceptibility [um/V]
const real_t Lambda	= 7.97;			// grating period for QPM  [μm] 
const real_t alpha_crp	= 0.17e-6;	// pump linear absorption [1/μm]
const real_t alpha_crs	= 1.57e-6;		// signal linear absorption [1/μm]
const real_t beta_crs	= 5e-5;		// signal 2-photons absorption [μm/W]
const real_t rhop		= 0;			// pump walk-off angle [rad] 
const real_t rhos		= 0;			// signal walk-off angle [rad] 


// Thermal properties
const real_t thermal_cond	= 8.4e-6;	// thermal conductivity [W/μm K]
const real_t heat_trf_coeff= 10e-12;	// Heat Transfer Coefficient in W/μm²K

real_t n(real_t L,real_t T)
{	// Sellmeier equation for the refractive index as a function of wavelength and temperature
	
	real_t A =  4.502483;
	real_t B =  0.007294;
	real_t C =  0.185087;
	real_t D =  -0.02357;
	real_t E =  0.073423;
	real_t F =  0.199595;
	real_t G =  0.001;
	real_t H =  7.99724;
	real_t b =  3.483933e-8 * pow(T + 273.15,2);
	real_t c =  1.607839e-8 * pow(T + 273.15,2);
	
	return sqrt( A + (B+b)/(pow(L,2)-pow((C+c),2)) + E/(pow(L,2)-pow(F,2)) + G/(pow(L,2)-pow(H,2))+ D*pow(L,2));
	
}


void GetCrystalProp( void )
{	// Print on screen all the parameters
	std::cout << "\n\nCrystal properties:\n" << std::endl;
	std::cout << "Crystal dimensions = ( " << LX*1e-3 << ", " << LY*1e-3 << ", " << LZ*1e-3 << " ) mm" << std::endl;
	std::cout << "Crystal discretization (dz, dy, dz) = ( " << dx << ", " << dy << ", " << dz << " ) μm " << std::endl;
	std::cout << "Nonlinear crystal = " << CrystalName << std::endl;
	std::cout << "Eff. second-order susc., deff = " << deff << " μm/V" << std::endl;
	std::cout << "Grating period for QPM, Λ = " << Lambda << " μm" << std::endl;
	std::cout << "Pump linear absorp., αp = " << alpha_crp << " 1/μm" << std::endl;
	std::cout << "Signal linear absorp., αs = " << alpha_crs << " 1/μm" << std::endl;
	std::cout << "Signal 2-photons absorp., βs = " << beta_crs << " 1/μm" << std::endl;
	std::cout << "Pump walk-off angle, ρ = " << rhop << " rad" << std::endl; 
	std::cout << "Signal walk-off angle, ρ = " << rhos << " rad" << std::endl;
	std::cout << "Thermal conductivity, k = " << thermal_cond << " W/μm K" << std::endl;
	std::cout << "Heat Transfer Coefficient, h = " << heat_trf_coeff << " W/μm² K\n\n" << std::endl;	
}

#endif // -> #ifdef _CRYSTAL
