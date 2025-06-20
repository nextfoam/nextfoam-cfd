/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    readTimeControls

Description
    Read the control parameters used by setDeltaT

\*---------------------------------------------------------------------------*/

bool adjustTimeStep =
    runTime.controlDict().getOrDefault<bool>("adjustTimeStep", false);

bool adjustLocalTimeStep =
    runTime.controlDict().getOrDefault<bool>("adjustLocalTimeStep", false);

scalar maxCo =
    runTime.controlDict().getOrDefault<scalar>("maxCo", 1.0);

scalar scaleSpectralRadiiV = 
    runTime.controlDict().getOrDefault<scalar>("scaleSpectralRadiiV", 0.0);


// ************************************************************************* //
