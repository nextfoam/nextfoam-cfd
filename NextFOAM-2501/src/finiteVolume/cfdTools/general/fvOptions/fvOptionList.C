/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "fvOptionList.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(optionList, 0);
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::dictionary& Foam::fv::optionList::optionsDict
(
    const dictionary& dict
)
{
    return dict.optionalSubDict("options", keyType::LITERAL);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::fv::optionList::readOptions(const dictionary& dict)
{
    checkTimeIndex_ = mesh_.time().timeIndex() + 2;

    bool allOk = true;
    for (fv::option& opt : *this)
    {
        bool ok = opt.read(dict.subDict(opt.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


void Foam::fv::optionList::checkApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        for (const fv::option& opt : *this)
        {
            opt.checkApplied();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::optionList::optionList(const fvMesh& mesh)
:
    PtrList<fv::option>(),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2),
    porousZoneCellMarker_(nullptr),
    porousInterfaceCells_(nullptr),
    activePorousCells_(nullptr),
    porousInterfaceFaces_(nullptr),
    porousInteriorFaces_(nullptr),
    porousFaces_(nullptr),
    linearWeights_(nullptr),
    porousGradpPtr_(nullptr),
    hasActivePorousZone_(false)
{}


Foam::fv::optionList::optionList(const fvMesh& mesh, const dictionary& dict)
:
    Foam::fv::optionList(mesh)
{
    reset(optionsDict(dict));

    hasActivePorousZone_ = findActivePorous();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::fv::optionList::~optionList()
{
    porousZoneCellMarker_.clear();
    porousInterfaceCells_.clear();
    activePorousCells_.clear();
    porousInterfaceFaces_.clear();
    porousInteriorFaces_.clear();
    porousFaces_.clear();
    linearWeights_.clear();
    porousGradpPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::optionList::reset(const dictionary& dict)
{
    // Count number of active fvOptions
    label count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            ++count;
        }
    }

    this->resize(count);

    count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            const word& name = dEntry.keyword();
            const dictionary& sourceDict = dEntry.dict();

            this->set
            (
                count++,
                option::New(name, sourceDict, mesh_)
            );
        }
    }
}


Foam::autoPtr<Foam::volScalarField> 
Foam::fv::optionList::getPorousZoneCellMarker() const
{
    if (debug)
    {
        Pout<< "optionsList::getPorousZoneCellMarker() : "
            << "Constructing porous media zone cell marker"
            << endl;
    }

    autoPtr<volScalarField> porousZoneCellMarkerPtr
    (
        new volScalarField
        (
            IOobject
            (
                "porousZoneCellMarker",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh_,
            dimensionedScalar("porousZoneCell", dimless, Zero)
        )
    );

    volScalarField& porousZoneCellMarker(porousZoneCellMarkerPtr());

    for (const fv::option& source : *this)
    {
        if (source.modelType() == "explicitPorositySource")
        {
            const word& cellZoneName = source.coeffs().get<word>("cellZone");

            const label zoneID = mesh_.cellZones().findZoneID(cellZoneName);

            if (zoneID != -1 && source.active())
            {
                const labelList& cellLabels = mesh_.cellZones()[zoneID];

                forAll(cellLabels, celli)
                {
                    porousZoneCellMarker.primitiveFieldRef()[cellLabels[celli]]
                        = 1.0;
                }
            }
        }
    }

    porousZoneCellMarker.correctBoundaryConditions();

    return porousZoneCellMarkerPtr;
}


const Foam::volScalarField& Foam::fv::optionList::porousZoneCellMarker() const
{
    if (!porousZoneCellMarker_)
    {
        porousZoneCellMarker_.reset(getPorousZoneCellMarker());
    }

    return porousZoneCellMarker_();
}


Foam::autoPtr<Foam::labelList> 
Foam::fv::optionList::getPorousInterfaceCells() const
{
    if (debug)
    {
        Pout<< "optionsList::getPorousInterfaceCells() : "
            << "Constructing porous media zone cell markers"
            << endl;
    }

    autoPtr<labelList> porousInterfaceCellsPtr(new labelList);

    labelList& porousInterfaceCells(porousInterfaceCellsPtr());
    
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    const surfaceScalarField& pBFs(porousInterfaceFaces());

    forAll(pBFs, facei)
    {
        if (pBFs[facei])
        {
            porousInterfaceCells.append(owner[facei]);
            porousInterfaceCells.append(neighbour[facei]);
        }
    }

    const surfaceScalarField::Boundary& pBFsBf = pBFs.boundaryField();

    forAll(pBFsBf, patchi)
    {
        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = p.faceCells();

        const scalarField& pPorousInterfaceFaces = pBFsBf[patchi];

        if (p.coupled())
        {
            forAll(pOwner, facei)
            {
                if (pPorousInterfaceFaces[facei])
                {
                    label own = pOwner[facei];
                    porousInterfaceCells.append(own);
                }
            }
        }
    }

    porousInterfaceCells.setSize(porousInterfaceCells.size());

    return porousInterfaceCellsPtr;
}


const Foam::labelList& Foam::fv::optionList::porousInterfaceCells() const
{
    if (!porousInterfaceCells_)
    {
        porousInterfaceCells_.reset(getPorousInterfaceCells());
    }

    return porousInterfaceCells_();
}


Foam::autoPtr<Foam::labelList> 
Foam::fv::optionList::getActivePorousCells() const
{
    if (debug)
    {
        Pout<< "optionsList::getActivePorousCells() : "
            << "Constructing porous media zone cell markers"
            << endl;
    }

    autoPtr<labelList> activePorousCellsPtr(new labelList);

    labelList& activePorousCells(activePorousCellsPtr());
   
    auto& zoneCellMarker = porousZoneCellMarker();

    const scalarField& zoneCells = zoneCellMarker.internalField();

    forAll(zoneCells, celli)
    {
        if (zoneCells[celli])
        {
            activePorousCells.append(celli);
        }
    }

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    const surfaceScalarField& pBFs(porousInterfaceFaces());

    forAll(pBFs, facei)
    {
        if (pBFs[facei])
        {
            if (zoneCells[owner[facei]])
            {
                activePorousCells.append(neighbour[facei]);
            }
            else
            {
                activePorousCells.append(owner[facei]);
            }
        }
    }

    const surfaceScalarField::Boundary& pBFsBf = pBFs.boundaryField();

    forAll(pBFsBf, patchi)
    {
        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = p.faceCells();

        const scalarField& pPorousInterfaceFaces = pBFsBf[patchi];

        if (p.coupled())
        {
            forAll(pOwner, facei)
            {
                if (pPorousInterfaceFaces[facei])
                {
                    label own = pOwner[facei];
                    activePorousCells.append(own);
                }
            }
        }
    }

    activePorousCells.setSize(activePorousCells.size());

    return activePorousCellsPtr;
}


const Foam::labelList& Foam::fv::optionList::activePorousCells() const
{
    if (!activePorousCells_)
    {
        activePorousCells_.reset(getActivePorousCells());
    }

    return activePorousCells_();
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::fv::optionList::getPorousInterfaceFaces() const
{
    if (debug)
    {
        Pout<< "optionsList::getPorousInterfaceFaces() : "
            << "Constructing porous media zone interface face markers"
            << endl;
    }

    auto& zoneCellMarker = porousZoneCellMarker();

    auto tPorousInterfaceFaces =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "porousInterfaceFaces",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh_,
            dimensionedScalar("porousFace", dimless, Zero)
        );

    auto& porousInterfaceFaces = tPorousInterfaceFaces.ref();
    porousInterfaceFaces.setOriented();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    scalarField& w = porousInterfaceFaces.primitiveFieldRef();
    const scalarField& porousZoneCells = zoneCellMarker.internalField();

    forAll(owner, facei)
    {
        const label& own = owner[facei];
        const label& nei = neighbour[facei];

        if (porousZoneCells[own] && !porousZoneCells[nei])
        {
            w[facei] = 1;
        }
        else if(!porousZoneCells[own] && porousZoneCells[nei])
        {
            w[facei] = 2;
        }
    }

    auto& wBf = porousInterfaceFaces.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        auto& pwBf = wBf[patchi];

        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();

        if (p.coupled())
        {
            const scalarField neighbourPorousZoneCells
            (
                zoneCellMarker.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pOwner, facei)
            {
                const label& own = pOwner[facei];

                if (porousZoneCells[own] && !neighbourPorousZoneCells[facei])
                {
                    pwBf[facei] = 1;
                }
                else if 
                (
                    !porousZoneCells[own] && neighbourPorousZoneCells[facei]
                )
                {
                    pwBf[facei] = 2;
                }
            }
        }
    }

    return tPorousInterfaceFaces;
}


const Foam::surfaceScalarField& 
Foam::fv::optionList::porousInterfaceFaces() const
{
    if (!porousInterfaceFaces_)
    {
        porousInterfaceFaces_.reset(getPorousInterfaceFaces().ptr());
    }

    return porousInterfaceFaces_();
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::fv::optionList::getPorousInteriorFaces() const
{
    if (debug)
    {
        Pout<< "optionsList::getPorousInteriorFaces() : "
            << "Constructing porous media zone interior face markers"
            << endl;
    }

    auto& zoneCellMarker = porousZoneCellMarker();

    auto tPorousInteriorFaces =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "porousInteriorFaces",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh_,
            dimensionedScalar("porousFace", dimless, Zero)
        );

    auto& porousInteriorFaces = tPorousInteriorFaces.ref();
    porousInteriorFaces.setOriented();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    scalarField& w = porousInteriorFaces.primitiveFieldRef();
    const scalarField& porousZoneCells = zoneCellMarker.internalField();

    forAll(owner, facei)
    {
        const label& own = owner[facei];
        const label& nei = neighbour[facei];

        if (porousZoneCells[own] && porousZoneCells[nei])
        {
            w[facei] = 1.0;
        }
    }

    auto& wBf = porousInteriorFaces.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        auto& pwBf = wBf[patchi];

        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();

        if (p.coupled())
        {
            const scalarField neighbourPorousZoneCells
            (
                zoneCellMarker.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pOwner, facei)
            {
                const label& own = pOwner[facei];

                if (porousZoneCells[own] && neighbourPorousZoneCells[facei])
                {
                    pwBf[facei] = 1.0;
                }
            }
        }
        else
        {
            forAll(pOwner, facei)
            {
                const label& own = pOwner[facei];

                if (porousZoneCells[own])
                {
                    pwBf[facei] = 1.0;
                }
            }
        }
    }

    return tPorousInteriorFaces;
}


const Foam::surfaceScalarField& 
Foam::fv::optionList::porousInteriorFaces() const
{
    if (!porousInteriorFaces_)
    {
        porousInteriorFaces_.reset(getPorousInteriorFaces().ptr());
    }

    return porousInteriorFaces_();
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::fv::optionList::getPorousFaces() const
{
    if (debug)
    {
        Pout<< "optionsList::getPorousFaces() : "
            << "Constructing porous media zone interior face markers"
            << endl;
    }

    auto& zoneCellMarker = porousZoneCellMarker();

    auto tPorousFaces =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "porousFaces",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh_,
            dimensionedScalar("porousFace", dimless, Zero)
        );

    auto& porousFaces = tPorousFaces.ref();
    porousFaces.setOriented();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    scalarField& w = porousFaces.primitiveFieldRef();
    const scalarField& porousZoneCells = zoneCellMarker.internalField();

    forAll(owner, facei)
    {
        const label& own = owner[facei];
        const label& nei = neighbour[facei];

        if (porousZoneCells[own] || porousZoneCells[nei])
        {
            w[facei] = 1.0;
        }
    }

    auto& wBf = porousFaces.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        auto& pwBf = wBf[patchi];

        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();

        if (p.coupled())
        {
            const scalarField neighbourPorousZoneCells
            (
                zoneCellMarker.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(pOwner, facei)
            {
                const label& own = pOwner[facei];

                if (porousZoneCells[own] || neighbourPorousZoneCells[facei])
                {
                    pwBf[facei] = 1.0;
                }
            }
        }
        else
        {
            forAll(pOwner, facei)
            {
                if (porousZoneCells[pOwner[facei]])
                {
                    pwBf[facei] = 1.0;
                }
            }
        }
    }

    return tPorousFaces;
}


const Foam::surfaceScalarField& Foam::fv::optionList::porousFaces() const
{
    if (!porousFaces_)
    {
        porousFaces_.reset(getPorousFaces().ptr());
    }

    return porousFaces_();
}


Foam::tmp<Foam::surfaceScalarField> 
Foam::fv::optionList::getLinearWeights() const
{
    auto tLinearWeights =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "linearWeights",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh_,
            dimensionedScalar("linearWeight", dimless, Zero)
        );

    auto& linearWeights = tLinearWeights.ref();
    linearWeights.setOriented();

    linearWeights = mesh_.surfaceInterpolation::weights();

    if (hasActivePorousZone())
    {
        const surfaceScalarField& pFs(porousFaces());
        const surfaceScalarField& pIFs(porousInterfaceFaces());

        forAll(pFs, facei)
        {
            if (pFs[facei])
            {
                if (pIFs[facei] == 1)
                {
                    linearWeights.primitiveFieldRef()[facei] = 0.0;
                }
                else if (pIFs[facei] == 2)
                {
                    linearWeights.primitiveFieldRef()[facei] = 1.0;
                }
                //else
                //{
                //    linearWeights.primitiveFieldRef()[facei] = 0.5;
                //}
            }
        }

        auto& lwBf = linearWeights.boundaryFieldRef();
        auto& pFsBf = pFs.boundaryField();

        forAll(mesh_.boundary(), patchi)
        {
            auto& pLwBf = lwBf[patchi];
            auto& pPorousFaces = pFsBf[patchi];
            auto& pPIFs = pIFs.boundaryField()[patchi];

            const polyPatch& p = mesh_.boundary()[patchi].patch();

            if (p.coupled())
            {
                forAll(pPorousFaces, facei)
                {
                    if(pPorousFaces[facei])
                    {
                        if (pPIFs[facei] == 1)
                        {
                            pLwBf[facei] = 0.0;
                        }
                        else if (pPIFs[facei] == 2)
                        {
                            pLwBf[facei] = 1.0;
                        }
                        //else
                        //{
                        //    pLwBf[facei] = 0.5;
                        //}
                    }
                }
            }
        }
    }

    return tLinearWeights;
}


const Foam::surfaceScalarField& Foam::fv::optionList::linearWeights() const
{
    if (!linearWeights_)
    {
        linearWeights_.reset(getLinearWeights().ptr());
    }

    return linearWeights_();
}


void Foam::fv::optionList::correctPorousFaceWeights
(
    surfaceScalarField& weights,
    bool flip
) const
{
    if (hasActivePorousZone())
    {
        const surfaceScalarField& pFs(porousFaces());
        const surfaceScalarField& pIFs(porousInterfaceFaces());

        forAll(pFs, facei)
        {
            if (pFs[facei])
            {
                if (pIFs[facei] == 1)
                {
                    weights.primitiveFieldRef()[facei] = 0.0;
                }
                else if (pIFs[facei] == 2)
                {
                    weights.primitiveFieldRef()[facei] = 1.0;
                }
                else if(flip)
                {
                    weights.primitiveFieldRef()[facei] = 
                        1.0 - weights.primitiveFieldRef()[facei];
                }
            }
        }

        auto& lwBf = weights.boundaryFieldRef();
        auto& pFsBf = pFs.boundaryField();

        forAll(mesh_.boundary(), patchi)
        {
            auto& pLwBf = lwBf[patchi];
            auto& pPorousFaces = pFsBf[patchi];
            auto& pPIFs = pIFs.boundaryField()[patchi];

            const polyPatch& p = mesh_.boundary()[patchi].patch();

            if (p.coupled())
            {
                forAll(pPorousFaces, facei)
                {
                    if (pPorousFaces[facei])
                    {
                        if (pPIFs[facei] == 1)
                        {
                            pLwBf[facei] = 0.0;
                        }
                        else if (pPIFs[facei] == 2)
                        {
                            pLwBf[facei] = 1.0;
                        }
                        else if (flip)
                        {
                            pLwBf[facei] = 1.0 - pLwBf[facei];
                        }
                    }
                }
            }
        }
    }
}


bool Foam::fv::optionList::appliesToField(const word& fieldName) const
{
    for (const fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            return true;
        }
    }

    return false;
}


bool Foam::fv::optionList::findActivePorous() const
{
    for (const fv::option& source : *this)
    {
        if (source.modelType() == "explicitPorositySource")
        {
            const word& cellZoneName = source.coeffs().get<word>("cellZone");

            const label zoneID = mesh_.cellZones().findZoneID(cellZoneName);

            if (zoneID != -1 && source.active())
            {
                if (!porousGradpPtr_ && mesh_.foundObject<volScalarField>("rUAp"))
                {
                    auto tPorousGradp
                    (
                        tmp<volVectorField>::New
                        (
                            IOobject
                            (
                                "porousGradp",
                                mesh_.time().timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE,
                                true
                            ),
                            mesh_,
                            dimensionedVector("porousGradp", dimless, Zero)
                        )
                    );
                    
                    porousGradpPtr_.reset(tPorousGradp.ptr());
                }

                return true;
            }
        }
    }

    return false;
}


void Foam::fv::optionList::adjustTransport
(
    volSymmTensorField& kappaEff,
    volSymmTensorField& alphaEff,
    const volScalarField& Cp
)
{
    const word fieldName("U");

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            const bool ok = source.isActive();

            if (ok)
            {
                source.adjustTransport(kappaEff, alphaEff, Cp);
            }
        }
    }
}


bool Foam::fv::optionList::read(const dictionary& dict)
{
    return readOptions(optionsDict(dict));
}


bool Foam::fv::optionList::writeData(Ostream& os) const
{
    // Write list contents
    for (const fv::option& opt : *this)
    {
        os  << nl;
        opt.writeHeader(os);
        opt.writeData(os);
        opt.writeFooter(os);
    }

    // Check state of IOstream
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const fv::optionList& options)
{
    options.writeData(os);
    return os;
}


// ************************************************************************* //
