/*---------------------------------------------------------------------------*/
// * This file contains the class Solver. 
// * The class contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution through the nonlinear crystal.
// * 
// * In particular, this file should be used when only two equation describes the 
// * problem, e.g., parametric down-convertion or second-harmonic generation.
// * Only two frequencies are involved in theses problems.
/*---------------------------------------------------------------------------*/


#ifndef _SOLVER
#define _SOLVER

#pragma once

void setZeros ( complex_t *A )
{	// Set 0
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			A[IDX(idx,idy,0)] = (0.0f, 0.0f);
		}
	}
	
	return ;
	
}

void complexProduct ( complex_t *C, complex_t *A, complex_t *B )
{	// Product of complex numbers in GPU
	
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			C[IDX(idx,idy,0)] = A[IDX(idx,idy,0)] * B[IDX(idx,idy,0)] ;
		}
	}
	
	return ;
	
}

#ifdef THERMAL
void dAdz( complex_t *dPump, complex_t *dSignal, complex_t *Pump, complex_t *Signal, real_t *Tfinal, real_t *DKint, real_t increment, uint s )
{	// Nonlinear term: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dPump/dSignal (x,y are different fields)
	
	real_t PI	= 3.14159265358979323846;     // pi
	complex_t Im (0.0f, 1.0f); // negative imaginary number -i
	
	uint idz = 0;
	
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			dPump[IDX(idx,idy,idz)] = Im*1.0*PI*dQ/(n(lp, Tfinal[IDX(idx,idy,s)])*lp) * Signal[IDX(idx,idy,idz)] * CpxConj(Pump[IDX(idx,idy,idz)]) * CpxExp(-(DKint[IDX(idx,idy,idz)]*(dz*increment+(s+1.0)*dz)/(s+1.0))) - 0.5*alpha_crp*Pump[IDX(idx,idy,idz)] ;
			// 		dPump[IDX(idx,idy,idz)] = (0.0f, 0.0f);
			dSignal[IDX(idx,idy,idz)] = Im*0.5*PI*dQ/(n(ls,Tfinal[IDX(idx,idy,s)])*ls)*Pump[IDX(idx,idy,idz)]*Pump[IDX(idx,idy,idz)] * CpxExp(+(DKint[IDX(idx,idy,idz)]*(dz*increment+(s+1.0)*dz)/(s+1.0))) - 0.5*(alpha_crs + beta_crs*0.5*EPS0*C*n(ls,Tfinal[IDX(idx,idy,s)])*CpxAbs2(Signal[IDX(idx,idy,idz)])) * Signal[IDX(idx,idy,idz)] ;
		}
	}
	
	return ;
	
}
#else
void dAdz( complex_t *dPump, complex_t *dSignal, complex_t *Pump, complex_t *Signal, real_t *DK, real_t Temp, real_t z, real_t increment )
{	// Nonlinear term: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dPump/dSignal (x,y are different fields)
	
	real_t PI	= 3.14159265358979323846;     // pi
	complex_t Im (0.0f, 1.0f); // negative imaginary number -i
	
	uint idz = 0;
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			dPump[IDX(idx,idy,idz)] = Im*1.0*PI*dQ/(n(lp, Temp)*lp) * Signal[IDX(idx,idy,idz)] * CpxConj(Pump[IDX(idx,idy,idz)]) * CpxExp(-Im*z*DK[IDX(idx,idy,idz)]) - 0.5*alpha_crp*Pump[IDX(idx,idy,idz)] ;
// 			dPump[IDX(idx,idy,idz)] = (0.0f, 0.0f);
			dSignal[IDX(idx,idy,idz)] = Im*0.5*PI*dQ/(n(ls, Temp)*ls) * Pump[IDX(idx,idy,idz)] * Pump[IDX(idx,idy,idz)] * CpxExp(+z*DK[IDX(idx,idy,idz)]) - 0.5*(alpha_crs + beta_crs*0.5*C*EPS0*n(ls, Temp) * CpxAbs2(Signal[IDX(idx,idy,idz)])) * Signal[IDX(idx,idy,idz)] ; 
		}
	}	
	
	return ;
	
}
#endif

void sumRK4( complex_t *Pump, complex_t *Signal, complex_t *k1p, complex_t *k1s, complex_t *k2p, complex_t *k2s,complex_t *k3p, complex_t *k3s,complex_t *k4p, complex_t *k4s )
{	// Final sum after appling the Rounge-Kutta algorithm
	
	uint idz = 0;
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			Pump[IDX(idx,idy,idz)] = Pump[IDX(idx,idy,idz)] + (k1p[IDX(idx,idy,idz)] + 2*k2p[IDX(idx,idy,idz)] + 2*k3p[IDX(idx,idy,idz)] + k4p[IDX(idx,idy,idz)]) * dz / 6;
			Signal[IDX(idx,idy,idz)] = Signal[IDX(idx,idy,idz)] + (k1s[IDX(idx,idy,idz)] + 2*k2s[IDX(idx,idy,idz)] + 2*k3s[IDX(idx,idy,idz)] + k4s[IDX(idx,idy,idz)]) * dz / 6;
		}
	}
	
	return ;
	
}


void linealCombination( complex_t *auxp, complex_t *auxs, complex_t *Pump, complex_t *Signal, complex_t *kp, complex_t *ks, real_t crk4 )
{	// Linear combination Ax + s.kx and saves the result in aux_x
	
	uint idz = 0;	
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			auxp[IDX(idx,idy,idz)] = Pump[IDX(idx,idy,idz)] + (kp[IDX(idx,idy,idz)] * crk4) ;
			auxs[IDX(idx,idy,idz)] = Signal[IDX(idx,idy,idz)] + (ks[IDX(idx,idy,idz)] * crk4) ;
		}
	}
	
	return ;
	
}


class Solver
{	// Difine the class Solver for modelling the fields propagation and heating
public:	
	complex_t *k1p, *k2p, *k3p, *k4p, *k1s, *k2s, *k3s, *k4s;
	complex_t *auxp, *auxs;	
	
	Solver()
	{	// Constructor
		//* RK4 (kx) and auxiliary (aux) GPU vectors 
		k1p = (complex_t*)malloc(nBytes2Dc);
		k2p = (complex_t*)malloc(nBytes2Dc);
		k3p = (complex_t*)malloc(nBytes2Dc);
		k4p = (complex_t*)malloc(nBytes2Dc);
		k1s = (complex_t*)malloc(nBytes2Dc);
		k2s = (complex_t*)malloc(nBytes2Dc);
		k3s = (complex_t*)malloc(nBytes2Dc);
		k4s = (complex_t*)malloc(nBytes2Dc);
		auxp = (complex_t*)malloc(nBytes2Dc);
		auxs = (complex_t*)malloc(nBytes2Dc);
	}
	
	~Solver()
	{	// Destructor
		free(k1p);	free(k2p);	free(k3p);		free(k4p);
		free(k1s);	free(k2s);	free(k3s);		free(k4s);	
		free(auxs);		free(auxp);		
	}
	
	void diffraction( Efields *A );
	void SSFM( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, real_t z, uint slice );
	void CWES( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, bool save_only_last );
	void solverRK4( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, real_t z, uint slice );
	void run( real_t Power, real_t waist, real_t focalpoint, real_t Temp, real_t T_inf, real_t Tpeltier1, real_t Tpeltier2, bool save_only_last, bool save_temperature );
};


void Solver::diffraction ( Efields *A )
{	// Applies the diffraction term to the electric fields
	
	
	// Set plan for FFT 2D//
	fftwf_plan plan2D; 
	
	plan2D = fftwf_plan_dft_2d(NX, NY, reinterpret_cast<fftwf_complex*>(A->Pump),
					   reinterpret_cast<fftwf_complex*>(A->PumpQ), FFTW_BACKWARD, FFTW_ESTIMATE); 
	fftwf_execute(plan2D);	
	
	plan2D = fftwf_plan_dft_2d(NX, NY, reinterpret_cast<fftwf_complex*>(A->Signal),
					   reinterpret_cast<fftwf_complex*>(A->SignalQ), FFTW_BACKWARD, FFTW_ESTIMATE); 
	fftwf_execute(plan2D);
	
	complexProduct (A->AuxQ, A->PropPump, A->PumpQ);
	copyMatrix<complex_t>( A->PumpQ, A->AuxQ );
	
	complexProduct (A->AuxQ, A->PropSignal, A->SignalQ);
	copyMatrix<complex_t>( A->SignalQ, A->AuxQ );
	
	plan2D = fftwf_plan_dft_2d(NX, NY, reinterpret_cast<fftwf_complex*>(A->PumpQ),
					   reinterpret_cast<fftwf_complex*>(A->Pump), FFTW_FORWARD, FFTW_ESTIMATE); 
	fftwf_execute(plan2D);
	FFTscale(A->Pump);
	
	plan2D = fftwf_plan_dft_2d(NX, NY, reinterpret_cast<fftwf_complex*>(A->SignalQ),
					   reinterpret_cast<fftwf_complex*>(A->Signal), FFTW_FORWARD, FFTW_ESTIMATE); 
	fftwf_execute(plan2D);
	FFTscale(A->Signal);
	
	fftwf_destroy_plan(plan2D);
	
	return ;
}


void Solver::solverRK4( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, real_t z, uint slice )
{	// Applies the Fourh-order Runge-Kutta Method with fixed step size dz
	// This function apply the fourth-order Runge-Kutta method	
	
	#ifdef THERMAL
	DK->setDKInt0();
	DK->IntegrateDK( T, slice );
	
	real_t increment = 0.0;
	dAdz( this->k1p, this->k1s, A->Pump, A->Signal, T->Tfinal, DK->DKint, 0.5*increment, slice );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k1p, this->k1s, 0.5 );
	
	increment = 0.5;
	dAdz( this->k2p, this->k2s, this->auxp, this->auxs, T->Tfinal, DK->DKint, 0.5*increment, slice );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k2p, this->k2s, 0.5 );
	dAdz( this->k3p, this->k3s, this->auxp, this->auxs, T->Tfinal, DK->DKint, 0.5*increment, slice );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k3p, this->k3s, 1.0 );
	
	increment = 1.0;
	dAdz( this->k4p, this->k4s, this->auxp, this->auxs, T->Tfinal, DK->DKint, 0.5*increment, slice );
	sumRK4( A->Pump, A->Signal, this->k1p, this->k1s, this->k2p, this->k2s, this->k3p, this->k3s,this->k4p, this->k4s );
	
	#else
	real_t increment = 0.0;
	dAdz( this->k1p, this->k1s, A->Pump, A->Signal, DK->DK, Temp, z, 0.5*increment );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k1p, this->k1s, 0.5 );

	increment = 0.5;
	dAdz( this->k2p, this->k2s, this->auxp, this->auxs, DK->DK, Temp, z, 0.5*increment );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k2p, this->k2s, 0.5 );
	dAdz( this->k3p, this->k3s, this->auxp, this->auxs, DK->DK, Temp, z, 0.5*increment );
	linealCombination( this->auxp, this->auxs, A->Pump, A->Signal, this->k3p, this->k3s, 1.0 );
	
	increment = 1.0;
	dAdz( this->k4p, this->k4s, this->auxp, this->auxs, DK->DK, Temp, z, 0.5*increment );
	sumRK4( A->Pump, A->Signal, this->k1p, this->k1s, this->k2p, this->k2s, this->k3p, this->k3s,this->k4p, this->k4s );
	#endif
	
	return ;
	
}


void Solver::SSFM( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, real_t z, uint slice )
{	// Applies the Split-Step Fourier Method using scheme N/2 - L - N/2
	solverRK4 ( A, T, DK, Temp, z, slice ); // Runge-Kutta 4 for dz/2
	diffraction ( A ) ; // diffraction for dz
	solverRK4 ( A, T, DK, Temp, z, slice ); // Runge-Kutta 4 for dz/2
	return ;
}


void Solver::CWES( Efields *A, Tfield *T, PhaseMatching *DK, real_t Temp, bool save_only_last )
{	// Solve the coupled-wave equations along the nonlinear crystal	
	real_t z = 0.0;
	uint s = 0;
	while( s < NZ) {
		if(save_only_last){// save fields in the last slice (save_only_last = true)
			if( s == NZ-1 ){
				// 				std::cout << "Saving only last slice" << std::endl;
				saveMatrixComplex ( A->Pump, "Pump_out" );
				saveMatrixComplex ( A->Signal, "Signal_out" );
			}
		}
		else{
			if( s <= NZ-1 ){// save fields in every slice (save_only_last = false)
				// 				std::cout << "Saving slice #" << s << std::endl;
				saveMatrixComplex ( A->Pump, "Pump_"+std::to_string(s) );
				saveMatrixComplex ( A->Signal, "Signal_"+std::to_string(s) );
			}				
		}
		// Split-Step Fourier Method (SSFM) for the z-th step
		SSFM( A, T,  DK, Temp, z, s );
		
		#ifdef THERMAL
		T->upDateQ( A,  s );
		#endif
		
		z+=dz; 
		s++;// next slice
	}
	
	return ;		
}


void Solver::run( real_t Power, real_t waist, real_t focalpoint, real_t Temp, real_t T_inf, real_t Tpeltier1, real_t Tpeltier2, bool save_only_last, bool save_temperature )
{	// Run the solver in the body of the main function.
	
	// 	Class instances
	Efields *A	= new Efields;
	Tfield *T		= new Tfield;
	PhaseMatching *DK	= new PhaseMatching;
	
	T->setTemperature( Temp );	// Set inicial conditions for temperature field
	
	#ifdef THERMAL
	// 	T->setOvenSurrounded( Tpeltier1, Tpeltier2 ); // set oven temperature
	T->setBottomOvens( Tpeltier1, Tpeltier2 ); // set oven temperature	
	// 	saveTensorRealGPU ( T->Tinic, "Tinic" );	
	T->setInitialQ();	// Set initial Q=0
	
	uint global_count = 0, num_of_iter = 5;	// accounts for the global iterations
	
	while (global_count < num_of_iter )
	{
		std::cout << "\n\nGlobal iteration #" << global_count << "\n\n" << std::endl;
		std::cout << "Solving 3D Heat Equation in the steady-state...\n" << std::endl;
		
		real_t tol = 5e-4/TSIZE; // tolerance for convergence
		std::cout << "\nTemperature tolerance = " << tol << "\n" << std::endl;
		
		real_t Reduced_sum = 1.0, aux_sum ;	// compares 2 consecutive steps in thermal calculations
		uint counter_T = 0; // accounts for the thermal iterations 
		std::cout.precision(20);
		
		while ( Reduced_sum >= tol )
		{
			T->upDate(T_inf); // update temperature
			uint chkconv = 5000; 
			if (counter_T % chkconv == 0)
			{	// Check temperature convergence every chkconv iterations
				if (counter_T > 0){aux_sum = Reduced_sum;}
				Reduced_sum = T->checkConvergence();
				std::cout << "\u03A3|Tf-Ti|²/N³ at #" << counter_T << " iteration is: " << Reduced_sum << std::endl;
				// check if reduced_sum reaches constant value
				if ( Reduced_sum==aux_sum ) Reduced_sum = 0.5*tol;
			}
			copyTensor<real_t>( T->Tinic, T->Tfinal);
			
			counter_T++;            
		}
		std::cout << counter_T << " iterations -> steady-state." << std::endl;
		
		// Set inicial conditions for electric fields
		DK->setDKFromTemperature( T );	// For simulations with thermal calculations
		A->setInputPump( Power, focalpoint, waist,  T->Tfinal );
		A->setNoisyField();
		A->setPropagators( T->Tfinal );
		
		std::cout << "\n\nSolving Coupled-Wave Equations (CWEs)...\n" << std::endl;
		this->CWES( A, T, DK, Temp, save_only_last );
		
		if (counter_T == 1 and global_count != 0){ global_count = num_of_iter; }
		global_count++;
	}
	
	if(save_temperature)	saveTensorReal ( T->Tfinal, "Tfinal" );
	#else
	// Set inicial conditions for electric fields
	A->setInputPump( Power, focalpoint, waist, T->Tfinal );	
	A->setNoisyField();
	A->setPropagators( T->Tfinal );
	std::cout << "\n\nEvolution in the crystal...\n" << std::endl;
	DK->setInicialDKConstant( Temp );	// for simulations without thermal calculations	
	this->CWES( A, T, DK, Temp, save_only_last );
	#endif
	
	std::cout << "Deallocating memory" << std::endl;
	delete A, T, DK;
	
	return ;
}

#endif // -> #ifdef _SOLVER
