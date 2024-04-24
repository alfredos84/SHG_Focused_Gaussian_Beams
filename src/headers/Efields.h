/*---------------------------------------------------------------------------*/
// * This file contains the class Efields to set all the electric fields involved in
// * the simulations.
/*---------------------------------------------------------------------------*/


#ifndef _EFIELDS
#define _EFIELDS

#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
/* Functions for FFT */

void FFTscale(complex_t *A)
{	// Scales a vector after Fourier transforms (FFT_INVERSE mode)	
	
	real_t size = NX*NY;
	uint idx;
	uint idy;
	uint idz = 0;
	
	for (idy = 0; idy <NY; idy++){
		for( idx = 0; idx < NX; idx++){
			A[IDX(idx,idy,idz)] = A[IDX(idx,idy,idz)] / size;
		}
	}
	return ;
	
}

void fftShift2DH( complex_t *Field, complex_t *aux)
{	 // Swap horizontally the values un a matrix
	
	uint c = (int) floor((real_t)NX/2);
	
	for (uint idy = 0; idy <NY; idy++){
		for( uint idx = 0; idx < c; idx++){ 		
			Field[IDX(idx+c,idy,0)]  =  aux[IDX(idx,idy,0)];
			Field[IDX(idx,idy,0)]  =  aux[IDX(idx+c,idy,0)];
		}
	}
	return ;
}


void fftShift2DV( complex_t *Field, complex_t *aux)
{	// Swap vertically the values un a matrix
	
	uint r = (int) floor((real_t)NY/2);
	
	for (uint idy = 0; idy <r; idy++)	{
		for( uint idx = 0; idx < NX; idx++){
			Field[IDX(idx,idy+r,0)]  =  aux[IDX(idx,idy,0)];
			Field[IDX(idx,idy,0)]  =  aux[IDX(idx,idy+r,0)];
		}  
	}
	
	return ;
}


void fftShift2D ( complex_t* in )
{	// Standard fftshift in 2D
	complex_t * aux = (complex_t*)malloc(nBytes2Dc);
	copyMatrix<complex_t>(aux, in); 	fftShift2DV(in, aux);
	copyMatrix<complex_t>(aux, in); 	fftShift2DH(in, aux);
	free(aux);
	
	return ;	
}


void setPump( complex_t *Pump, real_t Power, real_t focalpoint, real_t waistX, real_t *Tfinal )
{	// Set initial pump as a Gaussian beam
	
	real_t Temp = Tfinal[IDX(NX/2,NY/2,0)];
	
	real_t PI	  = 3.14159265358979323846;     // pi
	complex_t Im (0.0f, 1.0f); // imaginary number
	real_t uX     = 0.5*NX;
	real_t uY     = 0.5*NY;
	real_t waistY = waistX;
	real_t wX2    = waistX*waistX;
	real_t wY2    = waistY*waistY;
	real_t zRX    = PI*n(lp,Temp)*wX2/lp;
	real_t zRY    = PI*n(lp,Temp)*wY2/lp;
	real_t Ap0    = sqrtf(4*Power/(EPS0*C*PI*n(lp, Temp)*wX2));
	real_t etaX   = focalpoint/zRX;
	real_t etaY   = focalpoint/zRY;
	complex_t MX  = (1-Im*etaX);
	complex_t MY  = (1-Im*etaY);
	
	uint idz = 0;
	
	for( uint idx = 0; idx < NX; idx++){ 
		for (uint idy = 0; idy <NY; idy++){
			Pump[IDX(idx,idy,idz)] = (Ap0/MX) * CpxExp( (-powf((idx-uX)*dx,2) - powf((idy-uY)*dy,2) )/ (wX2*MX) ) ;
		}
	}
	
	return;
	
}


void beamPropagator ( complex_t *eiQz_pump, complex_t *eiQz_signal, real_t *Tfinal )
{	// This kernel set the beam propagator operators useful in diffraction calculations
	complex_t Im (0.0f,1.0f);
	real_t Temp = Tfinal[IDX(NX/2,NY/2,0)];
	real_t PI	= 3.14159265358979323846;     // pi
	real_t kp	= 2*PI*n(lp, Temp)/lp;
	real_t ks	= 2*PI*n(ls, Temp)/ls;
	real_t dfX= 1/dx/NX;	real_t dfY  = 1/dy/NY;
	real_t uX	= 0.5*NX;	real_t uY   = 0.5*NY;
	
	uint idz = 0;
	
	for (uint idy = 0; idy <NY; idy++){
		for( uint idx = 0; idx < NX; idx++){
			eiQz_pump[IDX(idx,idy,idz)] = CpxExp(-dz*(2*powf(PI,2)/kp * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2)) + 2*PI*dfX*tanf(rhop)*(idx-uX)) ); 
			eiQz_signal[IDX(idx,idy,idz)] = CpxExp(-dz*(2*powf(PI,2)/ks * ( dfX*dfX*powf(idx - uX,2) + dfY*dfY*powf(idy - uY,2))  + 2*PI*dfX*tanf(rhos)*(idx-uX)) ); 
		}
	}
	return ;	
}


void noiseGeneratorCPU ( complex_t *A )
{	// Noise generator for initial signal/idler vectors 
	uint seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<real_t> distribution(0.0,1.0e-20);
	
	real_t a1, a2;    
	for (int i=0; i<NX*NY; ++i) {
		a1=distribution(generator);
		a2=distribution(generator);
		complex_t a (a1,a2);
		A[i] = a;
	}
	
	return ;	
}


class Efields
{	
public:
	complex_t *Pump, *Signal;
	complex_t *PumpQ, *SignalQ;
	complex_t *PropPump, *PropSignal;		
	complex_t *AuxQ;
	
	Efields()
	{	// Constructor
		Pump = (complex_t*)malloc(nBytes2Dc);
		Signal = (complex_t*)malloc(nBytes2Dc);
		PumpQ = (complex_t*)malloc(nBytes2Dc);
		SignalQ = (complex_t*)malloc(nBytes2Dc);
		PropPump = (complex_t*)malloc(nBytes2Dc);
		PropSignal = (complex_t*)malloc(nBytes2Dc);
		AuxQ = (complex_t*)malloc(nBytes2Dc);
	}
	
	~Efields()
	{	// Destructor		
		free(Pump);		free(Signal);
		free(PumpQ);	free(SignalQ);
		free(PropPump);	free(PropSignal);
		free(AuxQ);
	}
	
	void setInputPump( real_t Power, real_t focalpoint, real_t waistX, real_t *Tfinal )
	{	// Set initial pump with a given power, beam waist and focal position
		setPump( this->Pump, Power, focalpoint, waistX, Tfinal );
		return ;
	}
	
	void setNoisyField()
	{	// Set signal vector as a noisy input
		noiseGeneratorCPU ( this->Signal );
		return ;
	}	
	
	void setPropagators( real_t *Tfinal )
	{	// Set vectors for beam propagators
		// 		* Diffraction operation ∇²_{xy}: Beam propagator for Pump and Signal
		beamPropagator ( this->PropPump, this->PropSignal, Tfinal );
		fftShift2D ( this->PropPump );	fftShift2D ( this->PropSignal );
		return ;
	}
	
};


#endif // -> #ifdef _EFIELDS
