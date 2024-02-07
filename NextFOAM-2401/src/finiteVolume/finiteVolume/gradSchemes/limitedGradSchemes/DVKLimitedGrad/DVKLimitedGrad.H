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
    Foam::fv::DVKLimitedGrad

Group
    grpFvGradSchemes

Description
    DVKLimitedGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the maximum and minimum cell and cell neighbour values and is
    applied to all components of the gradient.

    This limiter function is based on Venkatakrishnan's limter.

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    DVKLimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef DVKLimitedGrad_H
#define DVKLimitedGrad_H

#include "gradScheme.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class DVKLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DVKLimitedGrad
:
    public fv::gradScheme<Type>
{
    // Private Data

        tmp<fv::gradScheme<Type>> basicGradScheme_;

        //- Limiter coefficient
        const scalar k_;

        //- Limiter epsilon coefficient
        scalar K_;


    // Private Member Functions

        void limitGradient
        (
            const Field<scalar>& limiter,
            Field<vector>& gIf
        ) const;

        void limitGradient
        (
            const Field<vector>& limiter,
            Field<tensor>& gIf
        ) const;

        //- No copy construct
        DVKLimitedGrad(const DVKLimitedGrad&) = delete;

        //- No copy assignment
        void operator=(const DVKLimitedGrad&) = delete;


public:

    //- RunTime type information
    TypeName("DVKLimited");


    // Constructors

        //- Construct from mesh and schemeData
        DVKLimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorInFunction(schemeData)
                    << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
            else
            {
                scalar A(10);
                scalar B(5);

                K_ = (1.0 - k_)*pow(A, B - B*k_);
            }
        }


    // Member Functions

        inline void limitFaceCmpt
        (
            scalar& limiter,
            const scalar eps,
            const scalar maxDelta,
            const scalar minDelta,
            const scalar extrapolate
        ) const;

        inline void limitFace
        (
            Type& limiter,
            const Type& eps,
            const Type& maxDelta,
            const Type& minDelta,
            const Type& extrapolate
        ) const;

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * Inline Member Function  * * * * * * * * * * * * * //

template<class Type>
inline void DVKLimitedGrad<Type>::limitFaceCmpt
(
    scalar& limiter,
    const scalar eps,
    const scalar maxDelta,
    const scalar minDelta,
    const scalar extrapolate
) const
{
    scalar sqrEpsilon = sqr(eps);

    if (extrapolate > 0.0)
    {
        limiter =
            min
            (
                limiter,
                (
                    sqr(maxDelta)
                  + sqrEpsilon
                  + 2.0*extrapolate*maxDelta
                )
               /stabilise
                (
                    sqr(maxDelta)
                  + 2.0*sqr(extrapolate)
                  + maxDelta*extrapolate + sqrEpsilon,
                  VSMALL
                )
            );
    }
    else if (extrapolate < 0.0)
    {
        limiter =
            min
            (
                limiter,
                (
                    sqr(minDelta)
                  + sqrEpsilon
                  + 2.0*extrapolate*minDelta
                )
               /stabilise
                (
                    sqr(minDelta)
                  + 2.0*sqr(extrapolate)
                  + minDelta*extrapolate + sqrEpsilon,
                  VSMALL
                )
            );
    }
}


template<class Type>
inline void DVKLimitedGrad<Type>::limitFace
(
    Type& limiter,
    const Type& eps,
    const Type& maxDelta,
    const Type& minDelta,
    const Type& extrapolate
) const
{
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
    {
        limitFaceCmpt
        (
            setComponent(limiter, cmpt),
            component(eps, cmpt),
            component(maxDelta, cmpt),
            component(minDelta, cmpt),
            component(extrapolate, cmpt)
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
