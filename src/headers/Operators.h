/*---------------------------------------------------------------------------*/
// * This file contains a set of overloaded operators to deal with complex numbers.
/*---------------------------------------------------------------------------*/


#ifndef _OPERATORS
#define _OPERATORS

#pragma once


/////////////////////////////////////     OPERATORS     ////////////////////////////////////////
inline complex_t  operator+(const real_t &a, const complex_t &b) {
	
	complex_t c (a   + b.real(), b.imag());
	return c;
}


inline complex_t  operator+(const complex_t &b, const real_t &a) {
	
	complex_t c(a   + b.real(), + b.imag());
	return c;
}


inline complex_t  operator+(const complex_t &a, const complex_t &b) {
	
	complex_t c (a.real() + b.real(), a.imag() + b.imag());
	return c;
}

inline complex_t  operator-(const complex_t &a) {
	
	complex_t c(-a.real(), -a.imag());
	return c;
}

inline complex_t  operator-(const real_t &a, const complex_t &b) {
	
	complex_t c (a   - b.real(), - b.imag());
	return c;
}


inline complex_t  operator-(const complex_t &b, const real_t &a) {
	
	complex_t c	(b.real() - a, b.imag()) ;
	return c;
}


inline complex_t  operator-(const complex_t &a, const complex_t &b) {
	
	complex_t c(a.real() - b.real(), a.imag() - b.imag());
	return c;
}


inline complex_t  operator*(const real_t &a, const complex_t &b) {
	
	complex_t c (a * b.real(), a * b.imag());	
	return c;
}


inline complex_t  operator*(const complex_t &b, const real_t &a) {
	
	complex_t c(a * b.real(), a * b.imag()) ;	
	return c;
}


inline complex_t  operator*(const complex_t &a, const complex_t &b) {
	
	complex_t c(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real()) ;
	return c;
}


inline complex_t  operator/(const complex_t &b, const real_t &a) {
	
	complex_t c (b.real()/a , b.imag()/a ) ;	
	return c;
}


inline complex_t operator/(const real_t& a, const complex_t& b) {
	
	real_t denominator = b.real() * b.real() + b.imag() * b.imag();
	complex_t c( (+a * b.real()) / denominator, (-a * b.imag()) / denominator);
	return c;
}


inline complex_t operator/(const complex_t& a, const complex_t& b) {
	
	real_t denominator = b.real() * b.real() + b.imag() * b.imag();
	complex_t c( (a.real() * b.real() + a.imag() * b.imag()) / denominator, (a.imag() * b.real() - a.real() * b.imag()) / denominator);
	return c;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//* Complex exponential e^(i*a) */
inline complex_t CpxExp (real_t a)
{
	complex_t b (cosf(a) , sinf(a)) ;
	return b;
}


//* Complex exponential e^(a+i*b) */
inline complex_t CpxExp (complex_t a)
{
	complex_t b( expf(a.real())*cosf(a.imag()), expf(a.real())*sinf(a.imag()) );
	return b;
}


//* Complex conjugate */
inline complex_t CpxConj (complex_t a)
{
	complex_t b (+a.real() , -a.imag()) ;	
	return b;
}


//* Complex absolute value  */
inline real_t CpxAbs (complex_t a)
{
	real_t b;
	b = sqrtf(a.real()*a.real() + a.imag()*a.imag());
	return b;
}


//* Complex square absolute value */
inline real_t CpxAbs2 (complex_t a)
{
	real_t b;
	b = a.real()*a.real() + a.imag()*a.imag();
	return b;
}


inline complex_t CpxSqrt(complex_t z)
{
	real_t magnitude = sqrtf(z.real() * z.real() + z.imag() * z.imag());
	real_t real = sqrtf(0.5f * (magnitude + z.real()));
	real_t imag = sqrtf(0.5f * (magnitude - z.real()));
	
	if (z.imag() < 0)	imag = -imag;
	
	complex_t b (real, imag);
	return b;
}

#endif // -> #ifdef _OPERATORS
