/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "multiRegionSolutionControl.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiRegionSolutionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::multiRegionSolutionControl::readFluidControl
(
    const fvMesh* fvmesh
)
{
    const dictionary& solutionDict = this->dict(fvmesh);

    if (globalSolutionDict_.found(algorithmName_))
    {
        const dictionary& algorithmDict
        (
            globalSolutionDict_.subDict(algorithmName_)
        );

        if (algorithmDict.found("solveFlow"))
        {
            const_cast<dictionary&>(solutionDict).add
            (
                "solveFlow",
                algorithmDict.get<bool>("solveFlow"),
                true
            );
        }

        if (algorithmDict.found("solveEnergy"))
        {
            const_cast<dictionary&>(solutionDict).add
            (
                "solveEnergy",
                algorithmDict.get<bool>("solveEnergy"),
                true
            );
        }
    }

    solveFlow_ = solutionDict.getOrDefault("solveFlow", true);
    solveEnergy_ = solutionDict.getOrDefault("solveEnergy", true);
    solveSpecies_ = solutionDict.getOrDefault("solveSpecies", true);
    nNonOrthCorr_ = 
        solutionDict.getOrDefault<label>("nNonOrthogonalCorrectors", 0);
    momentumPredictor_ = 
        solutionDict.getOrDefault("momentumPredictor", true);
    transonic_ = solutionDict.getOrDefault("transonic", false);
    consistent_ = solutionDict.getOrDefault("consistent", false);
    frozenFlow_ = solutionDict.getOrDefault("frozenFlow", false);

    if (algorithmName_ == "PIMPLE")
    {
        nCorrPISO_ = solutionDict.getOrDefault<label>("nCorrectors", 1);
        SIMPLErho_ = solutionDict.getOrDefault("SIMPLErho", false);
        turbOnFinalIterOnly_ =
            solutionDict.getOrDefault("turbOnFinalIterOnly", false);
        finalOnLastPimpleIterOnly_ =
            solutionDict.getOrDefault("finalOnLastPimpleIterOnly", false);
        ddtCorr_ = solutionDict.getOrDefault("ddtCorr", true);
    }

    solutionDict_ = solutionDict;
}


void Foam::multiRegionSolutionControl::readSolidControl
(
    const fvMesh* fvmesh
)
{
    const dictionary& solutionDict = this->dict(fvmesh);

    if (globalSolutionDict_.found(algorithmName_))
    {
        const dictionary& algorithmDict
        (
            globalSolutionDict_.subDict(algorithmName_)
        );

        if (algorithmDict.found("solveEnergy"))
        {
            const_cast<dictionary&>(solutionDict).add
            (
                "solveEnergy",
                algorithmDict.get<bool>("solveEnergy"),
                true
            );
        }
    }

    solveEnergy_ = solutionDict.getOrDefault("solveEnergy", true);
    nNonOrthCorr_ = 
        solutionDict.getOrDefault<label>("nNonOrthogonalCorrectors", 0);

    solutionDict_ = solutionDict;
}


Foam::label Foam::multiRegionSolutionControl::applyToField
(
    const List<fieldData>& residualControl,
    const word& fieldName,
    const bool useRegEx
) const
{
    forAll(residualControl, i)
    {
        if (residualControl[i].name.match(fieldName, !useRegEx))
        {
            return i;
        }
    }

    return -1;
}


void Foam::multiRegionSolutionControl::storePrevIterFields() const
{
//    storePrevIter<label>();
    storePrevIter<scalar>();
    storePrevIter<vector>();
    storePrevIter<sphericalTensor>();
    storePrevIter<symmTensor>();
    storePrevIter<tensor>();
}


bool Foam::multiRegionSolutionControl::fluidLoop()
{
    if (fluidIndex_ < nFluidIndex_)
    {
        mesh_ = const_cast<fvMesh*>(regions_[fluidIndex_]);

        readFluidControl(mesh_);

        fluidIndex_++;

        if (foamLog_)
        {
            Info<< "Solving for fluid region "
                << mesh_->name() << endl;
        }

        return true;
    }
    else
    {
        fluidIndex_ = 0;

        return false;
    }
}


bool Foam::multiRegionSolutionControl::solidLoop()
{
    if (solidIndex_ < nSolidIndex_)
    {
        mesh_ = const_cast<fvMesh*>(regions_[nFluidIndex_ + solidIndex_]);

        readSolidControl(mesh_);

        solidIndex_++;

        if (foamLog_)
        {
            Info<< "\nSolving for solid region "
                << mesh_->name() << endl;
        }

        return true;
    }
    else
    {
        solidIndex_ = 0;

        return false;
    }
}


void Foam::multiRegionSolutionControl::setFirstIterFlag
(
    fvMesh* fvmesh,
    const bool check,
    const bool force
)
{
    DebugInfo
        << "multiRegionSolutionControl: force:" << force
        << " check: " << check
        << " corr: " << corr_
        << " corrNonOrtho:" << corrNonOrtho_
        << endl;

    if (force || (check && corr_ <= 1 && corrNonOrtho_ == 0))
    {
        DebugInfo
            << "multiRegionSolutionControl: set firstIteration flag" << endl;

        fvmesh->data().setFirstIteration(true);
    }
    else
    {
        DebugInfo
            << "multiRegionSolutionControl: remove firstIteration flag" 
            << endl;

        fvmesh->data().setFirstIteration(false); 
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
bool Foam::multiRegionSolutionControl::maxTypeResidual
(   
    const fvMesh* fvmesh,
    const entry& solverPerfDictEntry,
    Pair<scalar>& residuals
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const word& fieldName = solverPerfDictEntry.keyword();

    if (fvmesh->foundObject<fieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(solverPerfDictEntry.stream());

        residuals.first() = cmptMax(sp.first().initialResidual());
        residuals.last()  = cmptMax(sp.last().initialResidual());

        return true;
    }
    
    return false;
}


Foam::Pair<Foam::scalar> Foam::multiRegionSolutionControl::maxResidual
(
    const fvMesh* fvmesh,
    const entry& solverPerfDictEntry
)
{
    Pair<scalar> residuals(0,0);

    // Check with builtin short-circuit
    const bool ok =
    (
        maxTypeResidual<scalar>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<vector>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<sphericalTensor>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<symmTensor>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<tensor>(fvmesh, solverPerfDictEntry, residuals)
    );

    if (!ok && multiRegionSolutionControl::debug)
    {
        Info<<"no residual for " << solverPerfDictEntry.keyword()
            << " on mesh " << fvmesh->name() << nl;
    }

    return residuals;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiRegionSolutionControl::multiRegionSolutionControl
(
    const Time& runTime,
    const word& algorithmName,
    const PtrList<fvMesh>& fluid,
    const PtrList<fvMesh>& solid,
    bool foamLog
)
:
    objectRegistry
    (
        IOobject
        (
            algorithmName,
            runTime.timeName(),
            runTime
        )
    ),
    runTime_(runTime),
    regions_(fluid.size() + solid.size()),
    residualControls_(fluid.size() + solid.size()),
    mesh_(nullptr),
    algorithmName_(algorithmName),
    solutionDict_(),
    globalSolutionDict_
    (
        IOobject
        (
            "fvSolution",
            runTime_.system(),
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true 
        )
    ),
    solveFlow_(true),
    solveEnergy_(true),
    solveSpecies_(true),
    nFluidIndex_(0),
    nSolidIndex_(0),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false),
    frozenFlow_(false),
    nCorrPISO_(0),
    SIMPLErho_(false),
    turbOnFinalIterOnly_(false),
    finalOnLastPimpleIterOnly_(false),
    ddtCorr_(true),
    fluidIndex_(0),
    solidIndex_(0),
    corr_(0),
    corrNonOrtho_(0),
    foamLog_(foamLog)
{
    if
    (
        IOobject
        (
            "fvSolution",
            runTime_.system(),
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ).typeHeaderOk<dictionary>(true)
    )
    {
        globalSolutionDict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        globalSolutionDict_.addWatch();
    }

    forAll(fluid, i)
    {
        regions_[i] = fluid.get(i);
        nFluidIndex_ ++;
    }

    forAll(solid, i)
    {
        regions_[i + nFluidIndex_] = solid.get(i);
        nSolidIndex_ ++;
    }
}


Foam::multiRegionSolutionControl::multiRegionSolutionControl
(
    const Time& runTime,
    const word& algorithmName,
    const PtrList<dynamicFvMesh>& fluid,
    const PtrList<fvMesh>& solid,
    bool foamLog
)
:
    objectRegistry
    (
        IOobject
        (
            algorithmName,
            runTime.timeName(),
            runTime
        )
    ),
    runTime_(runTime),
    regions_(fluid.size() + solid.size()),
    residualControls_(fluid.size() + solid.size()),
    mesh_(nullptr),
    algorithmName_(algorithmName),
    solutionDict_(),
    globalSolutionDict_
    (
        IOobject
        (
            "fvSolution",
            runTime_.system(),
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true 
        )
    ),
    solveFlow_(true),
    solveEnergy_(true),
    solveSpecies_(true),
    nFluidIndex_(0),
    nSolidIndex_(0),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false),
    frozenFlow_(false),
    nCorrPISO_(0),
    SIMPLErho_(false),
    turbOnFinalIterOnly_(false),
    finalOnLastPimpleIterOnly_(false),
    ddtCorr_(true),
    fluidIndex_(0),
    solidIndex_(0),
    corr_(0),
    corrNonOrtho_(0),
    foamLog_(foamLog)
{
    if
    (
        IOobject
        (
            "fvSolution",
            runTime_.system(),
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        ).typeHeaderOk<dictionary>(true)
    )
    {
        globalSolutionDict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        globalSolutionDict_.addWatch();
    }

    forAll(fluid, i)
    {
        regions_[i] = fluid.get(i);
        nFluidIndex_ ++;
    }

    forAll(solid, i)
    {
        regions_[i + nFluidIndex_] = solid.get(i);
        nSolidIndex_ ++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiRegionSolutionControl::~multiRegionSolutionControl()
{
    regions_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary Foam::multiRegionSolutionControl::dict() const
{
    return solutionDict_;
}


// ************************************************************************* //
