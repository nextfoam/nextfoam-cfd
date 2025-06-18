/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::timeMarchingControl::nCorrPseudoTime() const
{
    return nCorrPseudoTime_;
}


inline bool Foam::timeMarchingControl::storeInitialResiduals() const
{
    // Start from second Pseudo-Time iteration
    //return corr_ == 2;
    // Store iniRes util the initialResidual is not an order of 1
    return corr_ < 5; 
}


inline bool Foam::timeMarchingControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::timeMarchingControl::finalIter() const
{
    return converged_ || (corr_ == nCorrPseudoTime_);
}


inline bool Foam::timeMarchingControl::finalInnerIter() const
{
    return finalIter();
}


// ************************************************************************* //
