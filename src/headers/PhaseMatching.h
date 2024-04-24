/*---------------------------------------------------------------------------*/
// * This file contains the class PhaseMatching.
/*---------------------------------------------------------------------------*/


#ifndef _PHASEMATCHING
#define _PHASEMATCHING

#pragma once


void initDK0( real_t *DKint )
{	// Set mismatch matrix filled of constant values
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			DKint[IDX(idx,idy,0)] = 0.0;
		}
	}
	
	return ;	
}


void integrateDK(real_t *DKint, real_t *DK,  real_t *Tfinal, uint s)
{	// This kernel integrates the DK matrix from z'=0 to z'=z
	
	for (uint idz = 0; idz <= s; idz++ ){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				DKint[IDX(idx,idy,0)] += DK[IDX(idx,idy,idz)];
			}
		}
	}
	
	return ;
	
}


void setInicialDKThermal( real_t *DK, real_t *Tfinal )
{	// Set mismatch matrix for thermal calculations
	
	for (uint idz = 0; idz < NZ; idz++ ){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				DK[IDX(idx,idy,idz)] = 2*PI*( 2*n(lp, Tfinal[IDX(idx,idy,idz)])/lp - n(ls, Tfinal[IDX(idx,idy,idz)])/ls + 1/(Lambda) );
			}
		}
	}
	
	return ;
	
}


void initDKConstant( real_t *DK, real_t Temp )
{	// Set mismatch matrix filled of constant values	
	
	for (uint idz = 0; idz < NZ; idz++ ){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				DK[IDX(idx,idy,idz)] = 2*PI*( 2*n(lp, Temp)/lp - n(ls, Temp)/ls + 1/(Lambda) );
			}
		}
	}
	
	return ;	
	
}


class PhaseMatching
{	// Difine the class PhaseMatching for calculations on the mismatch factor
public:
	real_t *DK, *DKint;
	
	PhaseMatching()
	{	// Constructor
		DK = (real_t*)malloc(nBytes3Dr);
		DKint = (real_t*)malloc(nBytes2Dr);
	}
	
	~PhaseMatching()
	{	// Destructor
		free(DK);	free(DKint);
	}
	
	void setDKInt0( void )
	{	// Set the DKint matrix = 0
		
		initDK0( this->DKint );
		
		return ;
	}
	
	void IntegrateDK( Tfield *T, uint s )
	{	// Sum in each slice the accumulated mismatch factor
		
		integrateDK( this->DKint, this->DK, T->Tfinal, s );
		
		return ;
	}
	
	void setDKFromTemperature( Tfield *T )
	{	// Set the DK matrix for thermal calculations
		
		setInicialDKThermal( this->DK, T->Tfinal );
		
		return ;
	}	
	
	void setInicialDKConstant( real_t Temp )
	{	// Set the DK matrix as a constant 
		
		initDKConstant( this->DK, Temp );
		
		return ;
	}
	
};

#endif // -> #ifdef _PHASEMATCHING
