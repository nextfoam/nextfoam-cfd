/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

inline Foam::label Foam::multiRegionPimpleControl::nCorrPIMPLE() const
{
    return nCorrPIMPLE_;
}


inline Foam::label Foam::multiRegionPimpleControl::nCorrPISO() const
{
    return nCorrPISO_;
}


inline Foam::label Foam::multiRegionPimpleControl::corrPISO() const
{
    return corrPISO_;
}


inline bool Foam::multiRegionPimpleControl::SIMPLErho() const
{
    return SIMPLErho_;
}


inline bool Foam::multiRegionPimpleControl::correct()
{
    setFirstIterFlag(mesh_);

    ++corrPISO_;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corrPISO = " << corrPISO_ << endl;
    }

    if (corrPISO_ <= nCorrPISO_)
    {
        return true;
    }

    corrPISO_ = 0;

    setFirstIterFlag(mesh_);

    return false;
}


inline bool Foam::multiRegionPimpleControl::storeInitialResiduals() const
{
    // Start from second PIMPLE iteration
    return (corr_ == 2) && (corrPISO_ == 0) && (corrNonOrtho_ == 0);
}


inline bool Foam::multiRegionPimpleControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::multiRegionPimpleControl::finalIter() const
{
    return converged_ || (corr_ == nCorrPIMPLE_);
}


inline bool Foam::multiRegionPimpleControl::finalInnerIter() const
{
    bool finalInnerIter =
        corrPISO_ == nCorrPISO_
     && corrNonOrtho_ == nNonOrthCorr_ + 1;

    if (finalOnLastPimpleIterOnly_)
    {
        return finalInnerIter && finalIter();
    }
    else
    {
        return finalInnerIter;
    }
}


inline bool Foam::multiRegionPimpleControl::turbCorr() 
{
    if (turbOnFinalIterOnly_)
    {
        if (finalIter())
        {
            setFirstIterFlag(mesh_, true, true);
        }
    }
    else
    {
        if (firstIter())
        {
            setFirstIterFlag(mesh_, true, true);
        }
        else
        {
            setFirstIterFlag(mesh_);
        }
    }

    return !turbOnFinalIterOnly_ || finalIter();
}


inline bool Foam::multiRegionPimpleControl::ddtCorr() const
{
    return ddtCorr_;
}


// ************************************************************************* //
