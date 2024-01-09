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

\*---------------------------------------------------------------------------*/

#include "smoothStartLinearMotion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(smoothStartLinearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        smoothStartLinearMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::solidBodyMotionFunctions::smoothStartLinearMotion::rampFactor() const
{
    const scalar t = time_.value();

    if (t < rampTime_)
    {
        // Ramping region
        if (rampTime_ == 0.0)
        {
            return 1;
        }
        else
        {
            scalar t1 = rampTime_;
            scalar a = 30.0/pow(t1, 5);

            return a*(pow(t, 5)/5 - t1*pow(t, 4)/2 + sqr(t1)*pow(t, 3)/3);
        }
    }
    else
    {
        // Past ramping region
        return 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::smoothStartLinearMotion::smoothStartLinearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);

    if (rampTime_ < 0)
    {
        FatalIOErrorIn
        (
            "solidBodyMotionFunctions::smoothStartLinearMotion::smoothStartLinearMotion",
            SBMFCoeffs_
        )   << "Negative rampTime not allowed."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::smoothStartLinearMotion::transformation() const
{
    scalar t = time_.value();

    // Translation of centre of gravity with constant velocity
    vector displacement(0, 0, 0);

    if (t < rampTime_)
    {
        // Ramping region
        if (rampTime_ == 0.0)
        {
            displacement = velocity_*t;
        }
        else
        {
            scalar t1 = rampTime_;
            scalar a = 30.0/pow(t1, 5);

            displacement 
                = velocity_
                  *a*(pow(t, 6)/30 - t1*pow(t, 5)/10 + sqr(t1)*pow(t, 4)/12);
        }
    }
    else
    {
        // Past ramping region
        displacement 
            = velocity_*
              (
                // Displacement during the ramping region
                  0.5*rampTime_
                // Displacement during constant velocity after ramping region
                + (t - rampTime_)
              );
    }

    quaternion R(1);
    septernion TR(septernion(-displacement)*R);

    Info<< "solidBodyMotionFunctions::translation::transformation(): "
        << "Time = " << t << " velocity = " << rampFactor()*velocity_
        << " transformation = " << TR
        << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::smoothStartLinearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.readEntry("velocity", velocity_);
    SBMFCoeffs_.readEntry("rampTime", rampTime_);

    return true;
}


// ************************************************************************* //
