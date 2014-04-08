//*************************************
//	File: binary.cpp
//	Project: BEMRI in C++
//	Purpose: Implementation of binary class
//	Author: Eric Addison
//	Original date: 14Oct13
//*************************************

// This file holds the implementation of all the binary class methods

#include "binary.h"

namespace Astro
{

	binary::binary(ptMass *masses[])
	/// constructor with array of two pointers-to-ptMasses
	{
		m1 = masses[0];					// assign the ptMasses pointers
		m2 = masses[1];					// also for m2
		M = m1->mass + m2->mass;		// compute total mass M
		mu = (m1->mass * m2->mass) / M;	// reduced mass mu		
	}

	binary::binary(ptMass *a, ptMass *b)
	/// constructor with two individual ptMass pointers
	{
		m1 = a;							// assign the ptMasses pointers
		m2 = b;							// also for m2
		M = m1->mass + m2->mass;		// compute total mass M
		mu = (m1->mass * m2->mass) / M;	// reduced mass mu
	}

	double binary::calc_E()
	/// calculate internal binary energy, i.e. K + V
	{
		//Go into rest frame of mass 1 and calculate E	
		Vector3d muX = m2->pos - m1->pos;	// position of reduced mass relative to m1
		Vector3d muV = m2->vel - m1->vel;	// velocity of reduced mass relative to m1

		// assign new energy
		E = 0.5 * mu * muV.dot(muV) - G*mu*M/muX.norm();	// kinetic + potential
		return E;
	}

	// test function for binary class
	void testBinary()
	{
		ptMass m1, m2;
		m1.pos << 5,5,5;		
		m1.vel << -1, 2, 3;
		m2.vel << 1, 4, 1;
	
		binary(&m1, &m2);
		
	}

} // END namespace Astro
