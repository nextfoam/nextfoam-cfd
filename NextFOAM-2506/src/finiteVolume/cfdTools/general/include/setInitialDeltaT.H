/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    setInitialDeltaT

Description
    Set the initial timestep corresponding to the timestep adjustment
    algorithm in setDeltaT but only if it would reduce the timestep.

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && ((CoNum > SMALL) || (DiNum > SMALL)))
    {
        if (CoNum < SMALL)
        {
            CoNum = SMALL;
        }

        if (DiNum < SMALL)
        {
            DiNum = SMALL;
        }

        runTime.setDeltaT
        (
            min
            (
                min(maxCo/CoNum, maxDi/DiNum)*runTime.deltaT().value(),
                min(runTime.deltaTValue(), maxDeltaT)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (CoNum > SMALL))
    {
        runTime.setDeltaT
        (
            min
            (
                maxCo*runTime.deltaTValue()/CoNum,
                min(runTime.deltaTValue(), maxDeltaT)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
