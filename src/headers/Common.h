/*---------------------------------------------------------------------------*/
// * This file contains a set of functions useful for timing the code
// * and indexing vectors and matrices
/*---------------------------------------------------------------------------*/


#ifndef _COMMON_H
#define _COMMON_H

inline double seconds()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

inline void TimingCode( double iElaps)
{
	if ( iElaps < 60. ){std::cout << "\n\n...time elapsed " <<  iElaps << " seconds\n\n " << std::endl;}
	else if ( iElaps >= 60. and iElaps < 3600. ){std::cout << "\n\n...time elapsed " <<  iElaps/60 << " minutes\n\n " << std::endl;}
	else{std::cout << "\n\n...time elapsed " <<  iElaps/3600 << " hours\n\n " << std::endl;}
	
	return ;
}

// Function for indexing matrices
inline int IDX( uint x,  uint y, uint z){
	return ((z*(NX*NY))+(y*NX)+x);
}

template<typename T>
inline void copyMatrix( T *out, T *in){
	
	uint idx, idy, idz = 0;
	
	for( idx = 0; idx < NX; idx++){ 
		for (idy = 0; idy <NY; idy++){
			out[IDX(idx,idy,idz)] = in[IDX(idx,idy,idz)];
		}
	}
	
	return ;
}

template<typename T>
inline void copyTensor( T *out, T *in){
	
	uint idx, idy, idz;
	
	for( idz = 0; idz < NZ; idz++){
		for( idx = 0; idx < NX; idx++){ 
			for (idy = 0; idy <NY; idy++){
				out[IDX(idx,idy,idz)] = in[IDX(idx,idy,idz)];
			}
		}
	}
	
	return ;
}

#endif // _COMMON_H
