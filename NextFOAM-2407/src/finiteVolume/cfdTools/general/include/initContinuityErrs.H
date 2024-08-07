/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Global
    cumulativeContErr

Description
    Declare and initialise the cumulative continuity error.

\*---------------------------------------------------------------------------*/

#ifndef initContinuityErrs_H
#define initContinuityErrs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<uniformDimensionedScalarField> cumulativeContErrIO
    (
        fluidRegions.size()
    );

    List<scalar> cumulativeContErr(fluidRegions.size(), 0.0);
    
    forAll(fluidRegions, i)
    {
        cumulativeContErrIO.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                   "cumulativeContErr",
                    runTime.timeName(),
                    "uniform",
                    fluidRegions[i].thisDb(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                dimensionedScalar(word::null, dimless, Zero)
            )
        );

        cumulativeContErr[i] =  cumulativeContErrIO[i].value();
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    uniformDimensionedScalarField cumulativeContErrIO
    (
        IOobject
        (
            "cumulativeContErr",
            runTime.timeName(),
            "uniform",
            mesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(word::null, dimless, Zero)
    );
    scalar& cumulativeContErr = cumulativeContErrIO.value();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
