/*---------------------------------------------------------------------------*/
// * This file contains the class Tfield to set the temperature in the same
// * grid than class Efields.
/*---------------------------------------------------------------------------*/


#ifndef _TFIELD
#define _TFIELD

#pragma once


real_t reducedSum( real_t *Tfinal, real_t *Tinic )
{	// Computes |Tfinal-Tinic|²/TSIZE.
	
	real_t sum = 0.0;
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){			
					sum += powf( Tfinal[IDX(idx,idy,idz)] - Tinic[IDX(idx,idy,idz)], 2);
			}
		}
	}
	
	return sum/TSIZE;
	
}


void initTemperature(real_t *Tinic, real_t Temp)
{	// Set the temperature of the nonlinear cristal at T = Temp that usually is the
	// phase-matching temperature. However, this is false in focused beams.
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				Tinic[IDX(idx,idy,idz)] = Temp;
			}
		}
	}
	
	return;
	
}


void initBottomOvensTemperature( real_t *Tfield, real_t TPeltier1, real_t TPeltier2 )
{	// Set the bottom face temperature/s. 
	// For a single oven set Tpeltier1 = Tpeltier2
	
	uint idy = 0;
	real_t nz = real_t(NZ);
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idx = 0; idx < NX; idx++){
			if ( idz < nz/2.0 ){
				Tfield[IDX(idx,idy,idz)] = TPeltier1;
			}
			if ( idz >= nz/2.0 ){
				Tfield[IDX(idx,idy,idz)] = TPeltier2;
			}
		}
	}
	
	return;
	
}


void initOvenSurrounded( real_t *Tfield, real_t TPeltier1, real_t TPeltier2 )
{	// Set the bottom, top and lateral faces temperature/s. 
	// For a single oven set Tpeltier1 = Tpeltier2
	
	real_t nz = real_t(NZ);
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				if ( (idx < NX) and (idy == 0) and (idz < nz/2.0) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier1;
				}
				if ( (idx < NX) and (idy == 0) and (idz >= nz/2.0) and (idz < NZ) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier2;
				}
				if ( (idx < NX) and (idy == NY-1) and (idz < nz/2.0) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier1;
				}
				if ( (idx < NX) and (idy == NY-1) and (idz >= nz/2.0) and (idz < NZ) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier2;
				}
				if ( (idy < NY) and (idx == 0) and (idz < nz/2.0) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier1;
				}
				if ( (idy < NY) and (idx == 0) and (idz >= nz/2.0) and (idz < NZ) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier2;
				}
				if ( (idy < NY) and (idx == NX-1) and (idz < nz/2.0) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier1;
				}
				if ( (idy < NY) and (idx == NX-1) and (idz >= nz/2.0) and (idz < NZ) ){
					Tfield[IDX(idx,idy,idz)] = TPeltier2;
				}
			}
		}
	}
	
	return;
	
}


void heatEquationTopOpen (real_t *Tf, real_t *Ti, real_t T_inf, real_t *Q )
{	// Change if other boundary conditions are required, Currently, the bottom face is in contact with
	// two ovens (Peltier cells) set in setPeltiers() class method.	
	
	// grid steps	
	real_t Bix = heat_trf_coeff*dx/thermal_cond, Biy = heat_trf_coeff*dy/thermal_cond, Biz = heat_trf_coeff*dz/thermal_cond;
	real_t ax = dx*dx/(2*dy*dy), ay = dy*dy/(2*dx*dx), az = dz*dz/(2*dx*dx); // for faces points
	real_t bx = dx*dx/(2*dz*dz), by = dy*dy/(2*dz*dz), bz = dz*dz/(2*dy*dy); // for faces points
	real_t av = dz/(2*dx), bv = dx/(2*dz), cv = dx*dz/(4*dy*dy), V = 1/(av+bv+2*cv+0.5*(Bix+Biz)); // for vertical sides
	real_t ah = dz/(2*dy), bh = dy/(2*dz), ch = dy*dz/(4*dx*dx), H = 1/(ah+bh+2*ch+0.5*(Biy+Biz)); // for horizontal sides
	real_t ad = dx/(2*dy), bd = dy/(2*dx), cd = dx*dy/(4*dz*dz), D = 1/(ad+bd+2*cd+0.5*(Bix+Biy)); // for depth sides
	real_t AV = dy*dz/dx, BV = dx*dz/dy, CV = dx*dy/dz, BiD = Bix*dy+Biy*dz+Biz*dx, GV = 1/(AV+BV+CV+BiD); //for vertices 
	real_t L2 = 1/(2/(dx*dx) + 2/(dy*dy) + 2/(dz*dz)); // for inner points
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				//* Inner points (x,y,z = 1 to N-2 ) */
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = L2 * ( (Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)])/(dx*dx) + (Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)])/(dy*dy) + (Ti[IDX(idx,idy,idz+1)] + Ti[IDX(idx,idy,idz-1)])/(dz*dz) + Q[IDX(idx,idy,idz)]/thermal_cond );
				}		
				
				// FACES /
				// Left  (x = 0)
				if ( (idx == 0) and (idy < NY-1) and (idy > 0) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx+1,idy,idz)] + ax*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + bx*(Ti[IDX(idx,idy,idz+1)] + Ti[IDX(idx,idy,idz-1)]) + Bix*T_inf ) / (1 + 2*ax + 2*bx + Bix) ;
				}
				// Right (x = Lx)
				if ( (idx == NX-1) and (idy < NY-1) and (idy > 0) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx-1,idy,idz)] + ax*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + bx*(Ti[IDX(idx,idy,idz+1)] + Ti[IDX(idx,idy,idz-1)]) + Bix*T_inf ) / (1 + 2*ax + 2*bx + Bix) ;
				}	
				// Top   (y = Ly)
				if ( (idx < NX-1) and (idx > 0) and (idy == NY-1) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx,idy-1,idz)] + ay*(Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)]) + by*(Ti[IDX(idx,idy,idz+1)] + Ti[IDX(idx,idy,idz-1)]) + Biy*T_inf ) / (1 + 2*ay + 2*by + Biy) ;
				}	
				// Back   (z = 0)
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx,idy,idz+1)] + az*(Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)]) + bz*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + Biz*T_inf ) / (1 + 2*az + 2*bz + Biz) ;
				}
				// Front  (z = Lcr)
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx,idy,idz-1)] + az*(Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)]) + bz*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + Biz*T_inf ) / (1 + 2*az + 2*bz + Biz) ;
				}
				
				
				// SIDES /	
				// 	VERTICAL (x=0, x=Lx, z=0, z=Lcr, y libre)
				// (x=0, z=0)
				if ( (idx == 0) and (idy < NY-1) and (idy > 0) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] =  (av*Ti[IDX(idx+1,idy,idz)] + bv*Ti[IDX(idx,idy,idz+1)] + cv*(Ti[IDX(idx,idy+1,idz)]+Ti[IDX(idx,idy-1,idz)]) + 0.5*T_inf*(Bix+Biz))*V ;
				}
				// (x=Lx, z=0)
				if ( (idx == NX-1) and (idy < NY-1) and (idy > 0) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] =  (av*Ti[IDX(idx-1,idy,idz)] + bv*Ti[IDX(idx,idy,idz+1)] + cv*(Ti[IDX(idx,idy+1,idz)]+Ti[IDX(idx,idy-1,idz)]) + 0.5*T_inf*(Bix+Biz))*V ;
				}
				// (x=Lx, z=Lcr)
				if ( (idx == NX-1) and (idy < NY-1) and (idy > 0) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] =  (av*Ti[IDX(idx-1,idy,idz)] + bv*Ti[IDX(idx,idy,idz-1)] + cv*(Ti[IDX(idx,idy+1,idz)]+Ti[IDX(idx,idy-1,idz)]) + 0.5*T_inf*(Bix+Biz))*V ;
				}
				// (x=0, z=Lcr)
				if ( (idx == 0) and (idy < NY-1) and (idy > 0) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] =  (av*Ti[IDX(idx+1,idy,idz)] + bv*Ti[IDX(idx,idy,idz-1)] + cv*(Ti[IDX(idx,idy+1,idz)]+Ti[IDX(idx,idy-1,idz)]) + 0.5*T_inf*(Bix+Biz))*V ;
				}
				
				// HORIZONTAL (y=0, y=Ly, z=0, z=Lcr)
				// (y=Ly, z=0)
				if ( (idy == NY-1) and (idx < NX-1) and (idx > 0) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] = (ah*Ti[IDX(idx,idy-1,idz)] + bh*Ti[IDX(idx,idy,idz+1)] + ch*(Ti[IDX(idx-1,idy,idz)]+Ti[IDX(idx+1,idy,idz)]) + 0.5*T_inf*(Biy+Biz))*H ;
				}
				// (y=Ly, z=Lcr)
				if ( (idy == NY-1) and (idx < NX-1) and (idx > 0) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] = (ah*Ti[IDX(idx,idy-1,idz)] + bh*Ti[IDX(idx,idy,idz-1)] + ch*(Ti[IDX(idx-1,idy,idz)]+Ti[IDX(idx+1,idy,idz)]) + 0.5*T_inf*(Biy+Biz))*H ;
				}
				
				// DEPTH (x=0, x=Lx, y=0, y=Ly)
				// (x=0, y=Ly)
				if ( (idx == 0) and (idy == NY-1) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = (ad*Ti[IDX(idx,idy-1,idz)] + bd*Ti[IDX(idx+1,idy,idz)] + cd*(Ti[IDX(idx,idy,idz+1)]+Ti[IDX(idx,idy,idz-1)]) + 0.5*T_inf*(Bix+Biy))*D ;
				}
				// (x=Lx, y=Lx)
				if ( (idx == NX-1) and (idy == NY-1) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = (ad*Ti[IDX(idx,idy-1,idz)] + bd*Ti[IDX(idx-1,idy,idz)] + cd*(Ti[IDX(idx,idy,idz+1)]+Ti[IDX(idx,idy,idz-1)]) + 0.5*T_inf*(Bix+Biy))*D ;
				}
				
				// VERTICES /
				// (x = 0, y = Ly, z = 0)
				if ( (idx == 0) and (idy == NY-1) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] =  ( AV*Ti[IDX(idx+1,idy,idz)] + BV*Ti[IDX(idx,idy-1,idz)] + CV*Ti[IDX(idx,idy,idz+1)] + T_inf*BiD )*GV ;
				}
				// (x = 0, y = Ly, z = Lcr)
				if ( (idx == 0) and (idy == NY-1) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] =  ( AV*Ti[IDX(idx+1,idy,idz)] + BV*Ti[IDX(idx,idy-1,idz)] + CV*Ti[IDX(idx,idy,idz-1)] + T_inf*BiD )*GV ;
				}
				
				// (x = Lx, y = Ly, z = 0)
				if ( (idx == NX-1) and (idy == NY-1) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] =  ( AV*Ti[IDX(idx-1,idy,idz)] + BV*Ti[IDX(idx,idy-1,idz)] + CV*Ti[IDX(idx,idy,idz+1)] + T_inf*BiD )*GV ;
				}
				// 	// (x = Lx, y = Ly, z = Lcr)
				if ( (idx == NX-1) and (idy == NY-1) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] =  ( AV*Ti[IDX(idx-1,idy,idz)] + BV*Ti[IDX(idx,idy-1,idz)] + CV*Ti[IDX(idx,idy,idz-1)] + T_inf*BiD )*GV ;
				}
				// 
			}
		}
	}
	
	return ;
}


void kernelHeatEquationOvenSurrounded (real_t *Tf, real_t *Ti, real_t T_inf, real_t *Q )
{	// Change if other boundary conditions are required, Currently, the bottom face is in contact with
	// two ovens (Peltier cells) set in setPeltiers() class method.	
	
	// grid steps	
	real_t Biz = heat_trf_coeff*dz/thermal_cond;
	real_t az = dz*dz/(2*dx*dx); // for faces points
	real_t bz = dz*dz/(2*dy*dy); // for faces points
	
	real_t L2 = 1/(2/(dx*dx) + 2/(dy*dy) + 2/(dz*dz)); // for inner points
	
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				//* Inner points (x,y,z = 1 to N-2 ) */
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz < NZ-1) and (idz > 0) ){
					Tf[IDX(idx,idy,idz)] = L2 * ( (Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)])/(dx*dx) + (Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)])/(dy*dy) + (Ti[IDX(idx,idy,idz+1)] + Ti[IDX(idx,idy,idz-1)])/(dz*dz) + Q[IDX(idx,idy,idz)]/thermal_cond );
				}		
				
				// FACES				
				// Back   (z = 0)
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz == 0) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx,idy,idz+1)] + az*(Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)]) + bz*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + Biz*T_inf ) / (1 + 2*az + 2*bz + Biz) ;
				}				
				// Front  (z = Lcr)
				if ( (idx < NX-1) and (idx > 0) and (idy < NY-1) and (idy > 0) and (idz == NZ-1) ){
					Tf[IDX(idx,idy,idz)] = ( Ti[IDX(idx,idy,idz-1)] + az*(Ti[IDX(idx+1,idy,idz)] + Ti[IDX(idx-1,idy,idz)]) + bz*(Ti[IDX(idx,idy+1,idz)] + Ti[IDX(idx,idy-1,idz)]) + Biz*T_inf ) / (1 + 2*az + 2*bz + Biz) ;
				}
			}
		}
	}
	
	return ;
}


void initQ0( real_t *Q )
{	// Set the initial internal heat source
	
	for (uint idz = 0; idz < NZ; idz++){
		for (uint idy = 0; idy < NY; idy++){
			for (uint idx = 0; idx < NX; idx++){
				Q[IDX(idx,idy,idz)] = 0.0;
			}
		}
	}
	
	return;
	
}


void setQ( real_t *Q, real_t *Tfinal, complex_t *Pump, complex_t *Signal, uint s )
{	// Set the internal heat source Q = αpIp +αsIs + βsIs²
	
	uint idz = 0;
	
	for (uint idy = 0; idy < NY; idy++){
		for (uint idx = 0; idx < NX; idx++){
			Q[IDX(idx,idy,s)] = alpha_crp*0.5*C*EPS0*n(lp, Tfinal[IDX(idx,idy,s)])*CpxAbs2(Pump[IDX(idx,idy,idz)]) + alpha_crs*0.5*C*EPS0*n(ls,Tfinal[IDX(idx,idy,s)])*CpxAbs2(Signal[IDX(idx,idy,idz)]) + beta_crs*powf(0.5*C*EPS0*n(ls,Tfinal[IDX(idx,idy,s)])*CpxAbs2(Signal[IDX(idx,idy,idz)]),2);
		}
	}
	
	return;
	
}


class Tfield
{	// Difine the class Tfield for calculation thermal profile
	
public:
	real_t *Tinic, *Tfinal, *Taux, *Q;
	
	Tfield()
	{	// Constructor
		Tinic = (real_t*)malloc(nBytes3Dr);
		Tfinal = (real_t*)malloc(nBytes3Dr);
		Taux = (real_t*)malloc(nBytes3Dr);
		Q = (real_t*)malloc(nBytes3Dr);
	}
	
	~Tfield()
	{	// Destructor
		free(Tinic);	free(Tfinal);
		free(Taux);	free(Q);
	}
	
	void setTemperature( real_t Temp )
	{	// Set inicial temperature
		initTemperature( this->Tinic, Temp );
		initTemperature( this->Tfinal, Temp );

		return ;
	}
	
	void setBottomOvens( real_t Tpeltier1, real_t Tpeltier2 )
	{	// Set ovens temperature in an open-top configuration
		initBottomOvensTemperature( this->Tinic, Tpeltier1, Tpeltier2 );
		initBottomOvensTemperature( this->Tfinal, Tpeltier1, Tpeltier2 );

		return ;
	}
	
	void setOvenSurrounded( real_t Tpeltier1, real_t Tpeltier2 )
	{	// Set ovens temperature in an surrounded-oven configuration
		initOvenSurrounded( this->Tinic, Tpeltier1, Tpeltier2 );
		initOvenSurrounded( this->Tfinal, Tpeltier1, Tpeltier2 );
		
		return ;
	}
	
	void upDate(real_t T_inf)
	{	// Update the temperature by solving the heat equation
		heatEquationTopOpen( this->Tfinal, this->Tinic, T_inf, this->Q );
		
		return ;
	}
	
	real_t checkConvergence(void)
	{	// Check whether the temperature changes during iterations
		return reducedSum(this->Tfinal, this->Tinic );
	}
	
	void setInitialQ(void)
	{	// Set the inicial  internal heat source Q=0
		initQ0( this->Q );
		
		return ;
	}
	
	void upDateQ( Efields *A, uint s )
	{	// Update the internal heat source Q
		setQ( this->Q, this->Tfinal, A->Pump, A->Signal, s );
		
		return ;
	}
	
};



#endif // -> #ifdef _TFIELD
