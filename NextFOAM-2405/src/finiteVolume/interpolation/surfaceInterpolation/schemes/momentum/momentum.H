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
    Foam::momentum

Group
    grpFvSurfaceInterpolationSchemes

Description
    Momentum-differencing interpolation scheme class

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    momentum.C

\*---------------------------------------------------------------------------*/

#ifndef momentum_H
#define momentum_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class momentum Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class momentum
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const momentum&) = delete;


public:

    //- Runtime type information
    TypeName("momentum");


    // Constructors

        //- Construct from mesh
        momentum(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        momentum(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        momentum
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().surfaceInterpolation::AUWeights();
        }
};


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
momentumInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().surfaceInterpolation::AUWeights()
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
momentumInterpolate(const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp =
        momentumInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
