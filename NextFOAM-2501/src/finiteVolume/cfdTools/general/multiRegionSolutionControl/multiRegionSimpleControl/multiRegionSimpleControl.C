/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "multiRegionSimpleControl.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiRegionSimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multiRegionSimpleControl::readControl()
{
    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh = regions_[regioni];

        const dictionary& simpleDict = this->dict(fvmesh);

        List<fieldData>& residualControl(residualControls_[regioni]);

        // Read residual information
        const dictionary residualDict
        (
            simpleDict.subOrEmptyDict("residualControl")
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

                fd.absTol = residualDict.get<scalar>(fName);
                fd.relTol = -1;
                fd.initialResidual = -1;

                data.append(fd);
            }
            else
            {
                fieldData& fd = data[fieldi];
                
                fd.absTol = residualDict.get<scalar>(fName);
            }
        }

        residualControl.transfer(data);
    }
}


bool Foam::multiRegionSimpleControl::criteriaSatisfied()
{
    bool checkedAndAchieved = true;

    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh(regions_[regioni]);
        List<fieldData>& residualControl(residualControls_[regioni]);

        if (residualControl.empty())
        {   
            return false;
        }

        bool achieved = true;
        bool checked = false; // safety that some checks were indeed performed

        const dictionary& solverDict = fvmesh->data().solverPerformanceDict();

        if (solverDict.empty())
        {
            checked = true;
            achieved = true;
        }

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

                const bool absCheck =
                    (residuals.first() < residualControl[fieldi].absTol);                
                achieved = achieved && absCheck;

                if (debug)
                {
                    Info<< algorithmName_ << " solution statistics:" << endl;

                    Info<< "    " << fieldName << ": tolerance = "
                        << residuals.first()
                        << " (" << residualControl[fieldi].absTol << ")"
                        << endl;
                }
            }
        }

        checkedAndAchieved = (checkedAndAchieved && checked && achieved);

    }

    return checkedAndAchieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionSimpleControl::multiRegionSimpleControl
(
    const Time& runTime, 
    const PtrList<fvMesh>& fluid, 
    const PtrList<fvMesh>& solid,
    bool foamLog
)
:
    multiRegionSolutionControl(runTime, "SIMPLE", fluid, solid, foamLog),
    initialised_(false)
{
    readControl();

    Info<< nl;

    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh = regions_[regioni];

        List<fieldData>& residualControl(residualControls_[regioni]);

        Info<< fvmesh->name() << " :" << endl;
        Info << algorithmName_;

        if (residualControl.empty())
        {
            const scalar duration =
                runTime_.endTime().value()
              - runTime_.startTime().value();

            Info<< ": no convergence criteria found. "
                << "Calculations will run for " << duration << " steps." 
                << nl;
        }
        else
        {
            Info<< ": convergence criteria" << nl;
            for (const fieldData& ctrl : residualControl)
            {
                Info<< "    field " << ctrl.name << token::TAB
                    << " tolerance " << ctrl.absTol
                    << nl;
            }
            Info<< endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiRegionSimpleControl::loop()
{
    forAll(regions_, regioni)
    {
        fvMesh* fvmesh = const_cast<fvMesh*>(regions_[regioni]);

        multiRegionSolutionControl::setFirstIterFlag(fvmesh, true, true);
    }

    readControl();

    Time& runTime = const_cast<Time&>(runTime_);

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


// ************************************************************************* //
