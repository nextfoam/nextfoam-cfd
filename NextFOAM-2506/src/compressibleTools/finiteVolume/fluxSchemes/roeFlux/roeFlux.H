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
    Foam::compressible::godunovFluxes::roeFlux

Description
    Implementation of inviscid numerical fluxes based on Roe & Pike FDS scheme
    with Harten's entropy fix, in an ALE formulation.

    ALE formulations:
    AIAA-98-2411, Darracq et.al
    AIAA-2006-2825, Gönc et.al

Author
    Hrvoje Jasak All rights reserved.
    Aleksandar Jemcov  All rights reserved.
    Oliver Borm  All rights reserved.
    Modified JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    roeFlux.C

\*---------------------------------------------------------------------------*/

#ifndef roeFlux_H
#define roeFlux_H

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
                           Class roeFlux Declaration
\*---------------------------------------------------------------------------*/

class roeFlux
:
    public godunovFlux
{

protected:

    // Protected data

        // Harten's Entropy-Fix parameter
        dimensionedScalar epsilon_;


public:

    //- Runtime type information
    TypeName("roeFlux");

    // Constructors

        //- Construct from components
        roeFlux
        (
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~roeFlux()
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
