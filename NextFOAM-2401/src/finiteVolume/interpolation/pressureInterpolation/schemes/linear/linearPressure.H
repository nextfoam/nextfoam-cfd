/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::linearPressure

Group
    grpFvPressureInterpolationSchemes

Description
    Central-differencing interpolation scheme class for pressure

SourceFiles
    linearPressure.C

\*---------------------------------------------------------------------------*/

#ifndef linearPressure_H
#define linearPressure_H

#include "pressureInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class linearPressure Declaration
\*---------------------------------------------------------------------------*/

class linearPressure
:
    public pressureInterpolationScheme
{
    // Private Data

        mutable bool correct_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const linearPressure&) = delete;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from mesh
        linearPressure(const fvMesh& mesh)
        :
            pressureInterpolationScheme(mesh),
            correct_(false)
        {}

        //- Construct from Istream
        linearPressure(const fvMesh& mesh, Istream&)
        :
            pressureInterpolationScheme(mesh),
            correct_(false)
        {}

        //- Construct from faceFlux and Istream
        linearPressure
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            pressureInterpolationScheme(mesh),
            correct_(false)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const volScalarField&) const
        {
            if (hasActivePorousZone_)
            {
                //correct_ = true;

                return fvOptions_.linearWeights();
            }
            else
            {
                return this->mesh().surfaceInterpolation::weights();
            }
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return correct_;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<surfaceScalarField> correction
        (
            const volScalarField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
