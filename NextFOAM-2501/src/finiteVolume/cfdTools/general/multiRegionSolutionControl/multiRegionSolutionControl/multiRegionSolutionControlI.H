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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::dictionary& Foam::multiRegionSolutionControl::dict
(
    const fvMesh* fvmesh
) const
{
    return fvmesh->solutionDict().subDict(algorithmName_);
}


inline Foam::label Foam::multiRegionSolutionControl::corr() const
{
    return corr_;
}


inline Foam::label Foam::multiRegionSolutionControl::corrNonOrtho() const
{
    return corrNonOrtho_;
}


inline Foam::label Foam::multiRegionSolutionControl::fluidIndex() const
{
    return fluidIndex_ - 1;
}


inline Foam::label Foam::multiRegionSolutionControl::solidIndex() const
{
    return solidIndex_ - 1;
}


inline bool Foam::multiRegionSolutionControl::solveFlow() const
{
    return solveFlow_;
}


inline bool Foam::multiRegionSolutionControl::solveEnergy() const
{
    return solveEnergy_;
}


inline bool Foam::multiRegionSolutionControl::solveSpecies() const
{
    return solveSpecies_;
}


inline Foam::label Foam::multiRegionSolutionControl::nFluidIndex() const
{
    return nFluidIndex_;
}


inline Foam::label Foam::multiRegionSolutionControl::nSolidIndex() const
{
    return nSolidIndex_;
}


inline Foam::label Foam::multiRegionSolutionControl::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline bool Foam::multiRegionSolutionControl::finalNonOrthogonalIter() const
{
    return corrNonOrtho_ == nNonOrthCorr_ + 1;
}


inline bool Foam::multiRegionSolutionControl::momentumPredictor() const
{
    return momentumPredictor_;
}


inline bool Foam::multiRegionSolutionControl::transonic() const
{
    return transonic_;
}


inline bool Foam::multiRegionSolutionControl::consistent() const
{
    return consistent_;
}


inline bool Foam::multiRegionSolutionControl::frozenFlow() const
{
    return frozenFlow_;
}


inline bool Foam::multiRegionSolutionControl::initialState() const
{
    bool initialState(false);

    IOobject phiHeader
    (
        "phi",
        Foam::Time::timeName(mesh_->time().startTime().value()),
        *mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if
    (
        !phiHeader.typeHeaderOk<surfaceScalarField>(true)
     && mesh_->time().timeIndex() == (mesh_->time().startTimeIndex() + 1)
    )
    {
        initialState = true;
    }

    return initialState;
}


inline bool Foam::multiRegionSolutionControl::correctNonOrthogonal()
{
    setFirstIterFlag(mesh_);

    ++corrNonOrtho_;

    if (debug)
    {
        Info<< algorithmName_ << " correctNonOrthogonal: corrNonOrtho = "
            << corrNonOrtho_ << endl;
    }

    if (corrNonOrtho_ <= nNonOrthCorr_ + 1)
    {
        return true;
    }

    corrNonOrtho_ = 0;

    setFirstIterFlag(mesh_);

    return false;
}


inline bool Foam::multiRegionSolutionControl::foamLog() const
{
    return foamLog_;
}


// ************************************************************************* //
