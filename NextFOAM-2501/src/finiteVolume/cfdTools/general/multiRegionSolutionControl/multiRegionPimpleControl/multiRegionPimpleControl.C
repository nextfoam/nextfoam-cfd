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

#include "multiRegionPimpleControl.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiRegionPimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multiRegionPimpleControl::readControl()
{
    nCorrPIMPLE_ = 
    (
        globalSolutionDict_.subDict("PIMPLE").getOrDefault<label>
        (
            "nOuterCorrectors", 
            1
        )
    );

    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh = regions_[regioni];

        const dictionary& pimpleDict = this->dict(fvmesh);

        List<fieldData>& residualControl(residualControls_[regioni]);

        // Read residual information
        const dictionary residualDict
        (
            pimpleDict.subOrEmptyDict("residualControl")
        );

        DynamicList<fieldData> data(residualControl);

        for (const entry& dEntry : residualDict)
        {
            const word& fName = dEntry.keyword();
            const label fieldi = applyToField(residualControl, fName, false);
            if (fieldi == -1)
            {
                fieldData fd;
                fd.name = fName.c_str();

                if (dEntry.isDict())
                {
                    const dictionary& fieldDict = dEntry.dict();
                    fd.absTol = fieldDict.get<scalar>("tolerance");
                    fd.relTol = fieldDict.get<scalar>("relTol");
                    fd.initialResidual = 0.0;
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << dEntry.keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }

                data.append(fd);
            }
            else
            {
                fieldData& fd = data[fieldi];
                if (dEntry.isDict())
                {
                    const dictionary& fieldDict = dEntry.dict();
                    fd.absTol = fieldDict.get<scalar>("tolerance");
                    fd.relTol = fieldDict.get<scalar>("relTol");
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << dEntry.keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }
        }

        residualControl.transfer(data);
    }
}


bool Foam::multiRegionPimpleControl::criteriaSatisfied()
{
    bool checkedAndAchieved = true;

    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh(regions_[regioni]);
        List<fieldData>& residualControl(residualControls_[regioni]);

        // no checks on first iteration - nothing has been calculated yet
        if ((corr_ == 1) || residualControl.empty() || finalIter())
        {   
            return false;
        }

        bool storeIni = this->storeInitialResiduals();

        bool achieved = true;
        bool checked = false; // safety that some checks were indeed performed

        const dictionary& solverDict = fvmesh->data().solverPerformanceDict();
        for (const entry& solverPerfDictEntry : solverDict)
        {
            const word& fieldName = solverPerfDictEntry.keyword();
            const label fieldi = applyToField(residualControl, fieldName);

            if (fieldi != -1)
            {
                Pair<scalar> residuals = maxResidual
                (
                    fvmesh,
                    solverPerfDictEntry
                );

                checked = true;

                scalar relative = 0.0;
                bool relCheck = false;

                const bool absCheck = 
                    (residuals.last() < residualControl[fieldi].absTol);

                if (storeIni)
                {
                    residualControl[fieldi].initialResidual = 
                        residuals.first();
                }
                else
                {
                    const scalar iniRes =
                        residualControl[fieldi].initialResidual
                      + ROOTVSMALL;

                    relative = residuals.last() / iniRes;
                    relCheck = relative < residualControl[fieldi].relTol;
                }

                achieved = achieved && (absCheck || relCheck);

                if (debug)
                {
                    Info<< algorithmName_ << " loop:" << endl;

                    Info<< "    " << fieldName
                        << " PIMPLE iter " << corr_
                        << ": ini res = "
                        << residualControl[fieldi].initialResidual
                        << ", abs tol = " << residuals.last()
                        << " (" << residualControl[fieldi].absTol << ")"
                        << ", rel tol = " << relative
                        << " (" << residualControl[fieldi].relTol << ")"
                        << endl;
                }
            }
        }

        checkedAndAchieved = (checkedAndAchieved && checked && achieved);

    }

    return checkedAndAchieved;
}


void Foam::multiRegionPimpleControl::setFirstIterFlag
(
    fvMesh* fvMesh,
    const bool check, 
    const bool force
)
{
    DebugInfo
        << "corr:" << corr_
        << " corrPISO:" << corrPISO_
        << " corrNonOrtho:" << corrNonOrtho_
        << endl;

    multiRegionSolutionControl::setFirstIterFlag
    (
        fvMesh,
        check && corrPISO_ <= 1, 
        force
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionPimpleControl::multiRegionPimpleControl
(
    const Time& runTime, 
    const PtrList<fvMesh>& fluid, 
    const PtrList<fvMesh>& solid,
    bool foamLog
)
:
    multiRegionSolutionControl(runTime, "PIMPLE", fluid, solid, foamLog),
    nCorrPIMPLE_(0),
    corrPISO_(0),
    converged_(false)
{
	readControl();

    if (nCorrPIMPLE_ > 1)
    {
        Info<< nl;
        forAll(regions_, regioni)
        {
            const fvMesh* fvmesh = regions_[regioni];

            Info<< fvmesh->name() << " :" << endl;
            Info << algorithmName_;

            if (residualControls_[regioni].empty())
            {
                Info<< ": no residual control data found "
                    << "for region " << regions_[regioni]->name()
                    << " calculations will employ " << nCorrPIMPLE_
                    << " corrector loops" << nl << endl;
            }
            else
            {
                Info<< ": max iterations = " << nCorrPIMPLE_
                    << endl;
                Info<< regions_[regioni]->name() << " :" << endl;
                forAll(residualControls_[regioni], i)
                {
                    Info<< "    field " << residualControls_[regioni][i].name 
                        << token::TAB
                        << ": relTol " << residualControls_[regioni][i].relTol
                        << ", tolerance " 
                        << residualControls_[regioni][i].absTol
                        << nl;
                }
                Info<< endl;
            }
        }
    }
    else
    {
        Info<< nl << "PIMPLE" << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


Foam::multiRegionPimpleControl::multiRegionPimpleControl
(
    const Time& runTime, 
    const PtrList<dynamicFvMesh>& fluid, 
    const PtrList<fvMesh>& solid,
    bool foamLog
)
:
    multiRegionSolutionControl(runTime, "PIMPLE", fluid, solid, foamLog),
    nCorrPIMPLE_(0),
    corrPISO_(0),
    converged_(false)
{
	readControl();

    if (nCorrPIMPLE_ > 1)
    {
        Info<< nl;
        forAll(regions_, regioni)
        {
            const fvMesh* fvmesh = regions_[regioni];

            Info<< fvmesh->name() << " :" << endl;
            Info << algorithmName_;

            if (residualControls_[regioni].empty())
            {
                Info<< ": no residual control data found "
                    << "for region " << regions_[regioni]->name()
                    << " calculations will employ " << nCorrPIMPLE_
                    << " corrector loops" << nl << endl;
            }
            else
            {
                Info<< ": max iterations = " << nCorrPIMPLE_
                    << endl;
                Info<< regions_[regioni]->name() << " :" << endl;
                forAll(residualControls_[regioni], i)
                {
                    Info<< "    field " << residualControls_[regioni][i].name 
                        << token::TAB
                        << ": relTol " << residualControls_[regioni][i].relTol
                        << ", tolerance " 
                        << residualControls_[regioni][i].absTol
                        << nl;
                }
                Info<< endl;
            }
        }
    }
    else
    {
        Info<< nl << "PIMPLE" << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiRegionPimpleControl::loop()
{
    readControl();

    ++corr_;

    if (debug)
    {   
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    forAll(regions_, regioni)
    {
        fvMesh* fvmesh = const_cast<fvMesh*>(regions_[regioni]);

        setFirstIterFlag(fvmesh);
    }

    if (corr_ == nCorrPIMPLE_ + 1)
    {
		forAll(regions_, regioni)
        {
            fvMesh* fvmesh(const_cast<fvMesh*>(regions_[regioni]));

            fvmesh->data().setFinalIteration(false);
        }

        if (nCorrPIMPLE_ != 1)
        {
            Info<< nl << algorithmName_ << ": not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }

        corr_ = 0;
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< nl << algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;

            forAll(regions_, regioni)
            {
                fvMesh* fvmesh(const_cast<fvMesh*>(regions_[regioni]));

                fvmesh->data().setFinalIteration(false); 
            }

            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            if (foamLog_)
            {
                Info<< nl << "PIMPLE" << ": iteration " << corr_ << endl;
            }
            storePrevIterFields();

            forAll(regions_, regioni)
            {
                fvMesh* fvmesh(const_cast<fvMesh*>(regions_[regioni]));

                fvmesh->data().setFinalIteration(true); 
            }

            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            forAll(regions_, regioni)
            {
                fvMesh* fvmesh(const_cast<fvMesh*>(regions_[regioni]));

                fvmesh->data().setFinalIteration(true); 
            }
        }

        if (corr_ <= nCorrPIMPLE_)
        {
            if (foamLog_)
            {
                Info<< nl << "PIMPLE" << ": iteration " << corr_ << endl;
            }
            storePrevIterFields();
            completed = false;
        }
    }

    return !completed;
}


Foam::scalar Foam::multiRegionPimpleControl::compressibleCourantNo
(
    const fvMesh& fvmesh,
    const Time& runTime,
    const volScalarField& rho,
    const surfaceScalarField& phi
) const
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()/rho.primitiveField()
    );

    scalar CoNum = 0.5*gMax(sumPhi/fvmesh.V().field())*runTime.deltaTValue();

    scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(fvmesh.V().field()))*runTime.deltaTValue();

    if (foamLog_)
    {
        Info<< "Region: " << fvmesh.name() 
            << " Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;
    }

    return CoNum;
}


Foam::scalar Foam::multiRegionPimpleControl::solidRegionDiffusionNo
(
    const fvMesh& fvmesh,
    const Time& runTime,
    const volScalarField& Cprho,
    const volScalarField& kappa
) const
{
    surfaceScalarField kapparhoCpbyDelta
    (
        sqr(fvmesh.surfaceInterpolation::deltaCoeffs())
       *fvc::interpolate(kappa)
       /fvc::interpolate(Cprho)
    );

    const scalar DiNum = max(kapparhoCpbyDelta).value()*runTime.deltaTValue();
    const scalar meanDiNum =
        average(kapparhoCpbyDelta).value()*runTime.deltaTValue();

    if (foamLog_)
    {
        Info<< "Region: " << fvmesh.name() 
            << " Diffusion Number mean: " << meanDiNum
            << " max: " << DiNum << endl;
    }

    return DiNum;
}


// ************************************************************************* //
