/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::compressible::localTimeStep

Description
    The time-step is computed locally so that the local Courant number
    does not exceed the specified value.

    References:
    Blazek, Jiri. "Computational Fluid Dynamics: Principles and Applications",
    2001, Referex Engineering, Elsevier, Chapter 6.1.4 Method 2

Author
    Oliver Borm  All rights reserved.
    Modified JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    localTimeStep.C

\*---------------------------------------------------------------------------*/

#ifndef localTimeStep_H
#define localTimeStep_H

#include "fvCFD.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Class localTimeStep Declaration
\*---------------------------------------------------------------------------*/

class localTimeStep
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to Velocity field
        const  volVectorField& U_;

        //- CFL number limited local time step field
        volScalarField CoDeltaT_;

        tmp<volScalarField> trDeltaT_;

        //- projected area sum of faces for all cells
        volVectorField deltaS_;
       
        //- volScalarField of all cell volumes 
        volScalarField cellVolume_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        localTimeStep(const localTimeStep&);

        //- Disallow default bitwise assignment
        void operator=(const localTimeStep&);

        //- Update the characteristic length scale
        void updateDeltaS();

public:

    // Constructors

        //- Construct from mesh and Istream
        localTimeStep(const volVectorField& U);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the Courant-number limited time-step
        const volScalarField& CoDeltaT() const
        {
            return CoDeltaT_;
        }

        const tmp<scalarField> deltaT() const
        {
            return 1.0/trDeltaT_().primitiveField();
        }

        const tmp<scalarField> rDeltaT() const
        {
            return trDeltaT_().primitiveField();
        }

        //- Update the local deltaT field
        // and adjust it to the minimum value optionally
        void update
        (
            scalar maxCo, 
            const tmp<scalarField> Lambda, 
            Switch adjustLocalTimeStep = "no"
        );
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
