/*---------------------------------------------------------------------------*/
// * This file contains two functions that save files in .dat extension
// * 1 - saveFile...Real()    : save real matrices or tensors
// * 2 - saveFile...Complex() : save complex matrices or tensors

// Inputs:
// - Vector   : matrix or tensor to save
// - Filename : name of the saved file

/*---------------------------------------------------------------------------*/


#ifndef _FILES
#define _FILES

#pragma once


void saveMatrixReal (real_t *Vector, std::string Filename)
{
	std::ofstream myfile;
	std::string extension = ".dat";
	myfile.open(Filename+extension);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Vector[IDX(ix,iy,0)] << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


void saveMatrixComplex (complex_t *Vector, std::string Filename)
{
	std::ofstream myfile;
	std::string filenamer = "_r.dat", filenamei = "_i.dat";
	myfile.open(Filename+filenamer);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Vector[IDX(ix,iy,0)].real() << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	myfile.open(Filename+filenamei);
	for (int iy = 0; iy < NY; iy++){
		for (int ix = 0; ix < NX; ix++)
			myfile << std::setprecision(20) << Vector[IDX(ix,iy,0)].imag() << "\t";
		myfile << "\n"; 
	}
	myfile.close();
	
	return;
	
}


void saveTensorReal ( real_t *Tensor, std::string Filename )
{
	
	std::ofstream myfile;
	std::string extension = ".dat";
	
	for (uint iz = 0; iz < NZ; iz++){
		myfile.open(Filename+"_"+std::to_string(iz)+extension);
		for (uint iy = 0; iy < NY; iy++){
			for (uint ix = 0; ix < NX; ix++){
				myfile << Tensor[IDX(ix,iy,iz)] << "\t";
			}
			myfile << "\n"; 
		}
		myfile.close();
	}
	
	return ;
	
}


#endif // -> #ifdef _FILES
