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
    Foam::solidBodyMotionFunctions::smoothStartLinearMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. 
    Smoothly starting constant velocity displacement.

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    smoothStartLinearMotion.C

\*---------------------------------------------------------------------------*/

#ifndef smoothStartLinearMotion_H
#define smoothStartLinearMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                  Class smoothStartLinearMotion Declaration
\*---------------------------------------------------------------------------*/

class smoothStartLinearMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Linear velocity
        vector velocity_;

        //- Ramping time scale
        scalar rampTime_;


    // Private Member Functions

        //- No copy construct
        smoothStartLinearMotion(const smoothStartLinearMotion&) = delete;

        //- No copy assignment
        void operator=(const smoothStartLinearMotion&) = delete;

        //- Velocity ramping factor resulting from rampTime_ value
        scalar rampFactor() const;


public:

    //- Runtime type information
    TypeName("smoothStartLinearMotion");


    // Constructors

        //- Construct from components
        smoothStartLinearMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new smoothStartLinearMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~smoothStartLinearMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
