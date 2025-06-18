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

#include "timeMarchingControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeMarchingControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::timeMarchingControl::read()
{
    solutionControl::read(true);
    return true;
}


bool Foam::timeMarchingControl::readPseudoTime()
{
    solutionControl::read(false);

    const dictionary pseudoTimeDict(dict());

    nCorrPseudoTime_ = 
        pseudoTimeDict.getOrDefault<label>("nPseudoTimeIterations", 20);

    return true;
}


bool Foam::timeMarchingControl::criteriaSatisfied()
{
    if (residualControl_.empty())
    {
        return false;
    }

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.data().solverPerformanceDict();
    for (const entry& solverPerfDictEntry : solverDict)
    {
        const word& fieldName = solverPerfDictEntry.keyword();
        const label fieldi = applyToField(fieldName);

        if (fieldi != -1)
        {
            Pair<scalar> residuals = maxResidual(solverPerfDictEntry);
            checked = true;

            const bool absCheck =
                (residuals.first() < residualControl_[fieldi].absTol);

            achieved = achieved && absCheck;

            if (debug)
            {
                Info<< algorithmName_ << " solution statistics:" << endl;

                Info<< "    " << fieldName << ": tolerance = "
                    << residuals.first()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


bool Foam::timeMarchingControl::pseudoTimeCriteriaSatisfied()
{
    // no checks on first iteration - nothing has been calculated yet
    if ((corr_ == 1) || residualControl_.empty() || finalIter())
    {
        return false;
    }


    const bool storeIni = this->storeInitialResiduals();

    bool achieved = true;
    bool checked = false;    // safety that some checks were indeed performed

    const dictionary& solverDict = mesh_.data().solverPerformanceDict();
    for (const entry& solverPerfDictEntry : solverDict)
    {
        const word& fieldName = solverPerfDictEntry.keyword();
        const label fieldi = applyToField(fieldName);

        if (fieldi != -1)
        {
            Pair<scalar> residuals = maxResidual(solverPerfDictEntry);

            checked = true;

            scalar relative = 0.0;
            bool relCheck = false;

            const bool absCheck =
                (residuals.last() < residualControl_[fieldi].absTol);

            if (storeIni)
            {
                residualControl_[fieldi].initialResidual = residuals.first();
            }
            else
            {
                const scalar iniRes =
                    (residualControl_[fieldi].initialResidual + ROOTVSMALL);

                relative = residuals.last() / iniRes;
                relCheck = (relative < residualControl_[fieldi].relTol);
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< "Pseudo-Time" << " loop:" << endl;

                Info<< "    " << fieldName
                    << " Pseudo-Time iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldi].initialResidual
                    << ", abs tol = " << residuals.last()
                    << " (" << residualControl_[fieldi].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldi].relTol << ")"
                    << endl;
            }
        }
    }

    return checked && achieved;
}


void Foam::timeMarchingControl::setFirstIterFlag
(
    const bool check, 
    const bool force
)
{
    DebugInfo
        << "corr:" << corr_
        << endl;

    solutionControl::setFirstIterFlag(check, force);
}


void Foam::timeMarchingControl::storePrevIterFields() const
{
//    storePrevIter<label>();
    storePrevIter<scalar>();
    storePrevIter<vector>();
    storePrevIter<sphericalTensor>();
    storePrevIter<symmTensor>();
    storePrevIter<tensor>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeMarchingControl::timeMarchingControl
(
    fvMesh& mesh,
    const word& dictName,
    const bool pseudoTime,
    const bool verbose
)
:
    solutionControl(mesh, dictName),
    initialised_(false),
    nCorrPseudoTime_(0),
    converged_(false)
{
    if (pseudoTime)
    {
        readPseudoTime();

        if (verbose)
        {
            Info<< nl << "Pseudo-Time";

            if (residualControl_.empty())
            {
                Info<< ": no residual control data found. "
                    << "Calculations will employ " << nCorrPseudoTime_
                    << " pseudo-time loops" << nl;
            }
            else
            {
                Info<< ": max iterations = " << nCorrPseudoTime_ << nl;

                for (const fieldData& ctrl : residualControl_)
                {
                    Info<< "    field " << ctrl.name << token::TAB
                        << ": relTol " << ctrl.relTol
                        << ", tolerance " << ctrl.absTol
                        << nl;
                }
            }

            Info<< endl;
        }
    }
    else
    {
        read();

        if (verbose)
        {
            Info<< nl << algorithmName_;
    
            if (residualControl_.empty())
            {
                const scalar duration =
                    mesh_.time().endTime().value()
                - mesh_.time().startTime().value();
    
                Info<< ": no convergence criteria found. "
                    << "Calculations will run for " << duration << " steps."
                    << nl;
            }
            else
            {
                Info<< ": convergence criteria" << nl;
                for (const fieldData& ctrl : residualControl_)
                {
                    Info<< "    field " << ctrl.name << token::TAB
                        << " tolerance " << ctrl.absTol
                        << nl;
                }
            }
            Info<< endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::timeMarchingControl::loop()
{
    solutionControl::setFirstIterFlag(true, true);

    read();

    Time& runTime = const_cast<Time&>(mesh_.time());

    if (initialised_ && criteriaSatisfied())
    {
        Info<< nl
            << algorithmName_
            << " solution converged in "
            << runTime.timeName() << " iterations" << nl << endl;

        // Set to finalise calculation
        runTime.writeAndEnd();
    }
    else
    {
        initialised_ = true;
        storePrevIterFields();
    }

    return runTime.loop();
}


bool Foam::timeMarchingControl::pseudoTimeMarching()
{
    readPseudoTime();

    ++corr_;

    if (debug)
    {
        Info<< "Pseudo-Time" << " loop: corr = " << corr_ << endl;
    }

    setFirstIterFlag();

    if (corr_ == nCorrPseudoTime_ + 1)
    {
        if (!residualControl_.empty() && (nCorrPseudoTime_ != 1))
        {
            Info<< "Pseudo-Time" << ": not converged within "
                << nCorrPseudoTime_ << " iterations" << endl;
        }

        corr_ = 0;
        mesh_.data().setFinalIteration(false);
        return false;
    }

    bool completed = false;
    if (converged_ || pseudoTimeCriteriaSatisfied())
    {
        if (converged_)
        {
            Info<< "Pseudo-Time" << ": converged in " << corr_ - 1
                << " iterations" << endl;

            mesh_.data().setFinalIteration(false);
            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< "Pseudo-Time" << ": iteration " << corr_ << endl;
            storePrevIterFields();

            mesh_.data().setFinalIteration(true);
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            mesh_.data().setFinalIteration(true);
        }

        if (corr_ <= nCorrPseudoTime_)
        {
            Info<< "Pseudo-Time" << ": iteration " << corr_ << endl;
            storePrevIterFields();
            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
