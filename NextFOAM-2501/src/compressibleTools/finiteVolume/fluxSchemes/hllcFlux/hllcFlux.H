/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::compressible::godunovFluxes::hllcFlux

Description
    HLLC Arbitrary Lagrangian-Eulerian (ALE) flux class, computes numerical
    fluxes based on HLLC ALE scheme. Uses primitive variables as input and
    gives back conservative numerical fluxes.

    Luo, H.; Baum, Joseph D. and Löhner R. "On the computation of multi-
    material flows using ALE formulation."
    Journal of Computational Physics 194 (2004): 304–328.

Author
    Oliver Borm  All rights reserved.
    Modified JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    hllcFlux.C

\*---------------------------------------------------------------------------*/

#ifndef hllcFlux_H
#define hllcFlux_H

#include "godunovFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace godunovFluxes
{

/*---------------------------------------------------------------------------*\
                           Class hllcFlux Declaration
\*---------------------------------------------------------------------------*/

class hllcFlux
:
    public godunovFlux
{

protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("hllcFlux");

    // Constructors

        //- Construct from components
        hllcFlux
        (
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~hllcFlux()
    {}

    // Member Functions
    
        virtual void evaluate
        (
            scalar& rhoFlux,
            vector& rhoUFlux,
            scalar& rhoEFlux,
            const scalar pLeft,
            const scalar pRight,
            const vector ULeft,
            const vector URight,
            const scalar rhoLeft,
            const scalar rhoRight,
            const scalar kLeft,
            const scalar kRight,
            const scalar gammaLeft,
            const scalar gammaRight,
            const vector Sf,
            const scalar magSf,
            const vector dotX
        ) const;


        //- Update fluxes based on current state
        virtual void update();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace godunovFluxes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
