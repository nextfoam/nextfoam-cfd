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
    Foam::compressible::godunovFluxes::rusanovFlux

Description
    Implementation of inviscid numerical fluxes based on Rusanov FDS scheme
    This is the first order scheme

Author
    Aleksandar Jemcov  All rights reserved.
    Oliver Borm  All rights reserved.
    Modified JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    rusanovFlux.C

\*---------------------------------------------------------------------------*/

#ifndef rusanovFlux_H
#define rusanovFlux_H

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
                           Class rusanovFlux Declaration
\*---------------------------------------------------------------------------*/

class rusanovFlux
:
    public godunovFlux
{

protected:

    // Protected data


public:

    //- Runtime type information
    TypeName("rusanovFlux");

    // Constructors

        //- Construct from components
        rusanovFlux
        (
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~rusanovFlux()
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
