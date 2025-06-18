/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "MRFZoneList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZoneList::MRFZoneList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<MRFZone>(),
    mesh_(mesh),
    coriolisForce_ // by Gill
    (
        IOobject
        (
            "MRFZoneList:coriolisForce",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimAcceleration, Zero)
    ),
    volumeAcceleration_ // by Gill
    (
        IOobject
        (
            "MRFZoneList:volumeAcceleration",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimArea*dimAcceleration, Zero)
    )
{
    reset(dict);

    active(true);

    volumeAcceleration_.setOriented(true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::MRFZoneList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "    No MRF zones active" << endl;
    }

    return a;
}


void Foam::MRFZoneList::reset(const dictionary& dict)
{
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
            const dictionary& modelDict = dEntry.dict();

            Info<< "    creating MRF zone: " << name << endl;

            this->set
            (
                count++,
                new MRFZone(name, mesh_, modelDict)
            );
        }
    }
}


const Foam::MRFZone& Foam::MRFZoneList::getFromName
(
    const word& name
) const
{
    DynamicList<word> names;
    for (const auto& mrf: *this)
    {
        if (mrf.name() == name)
        {
            return mrf;
        }

        names.append(mrf.name());
    }

    FatalErrorInFunction
        << "Unable to find MRFZone " << name
        << ". Available zones are: " << names
        << exit(FatalError);

    return first();
}


bool Foam::MRFZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    for (auto& mrf: *this)
    {
        bool ok = mrf.read(dict.subDict(mrf.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::MRFZoneList::writeData(Ostream& os) const
{
    for (const auto& mrf: *this)
    {
        os  << nl;
        mrf.writeData(os);
    }

    return os.good();
}


void Foam::MRFZoneList::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    for (const auto& mrf: *this)
    {
        mrf.addCoriolis(U, ddtU);
    }
}


void Foam::MRFZoneList::addAcceleration(fvVectorMatrix& UEqn) const
{
    for (const auto& mrf: *this)
    {
        mrf.addCoriolis(UEqn);
    }
}


void Foam::MRFZoneList::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    for (const auto& mrf: *this)
    {
        mrf.addCoriolis(rho, UEqn);
    }
}


Foam::tmp<Foam::volVectorField> Foam::MRFZoneList::DDt
(
    const volVectorField& U
) const
{
    auto tacceleration = volVectorField::New
    (
        IOobject::scopedName("MRFZoneList", "acceleration"),
        IOobject::NO_REGISTER,
        U.mesh(),
        dimensionedVector(U.dimensions()/dimTime, Zero)
    );
    auto& acceleration = tacceleration.ref();

    for (const auto& mrf: *this)
    {
        mrf.addCoriolis(U, acceleration);
    }

    return tacceleration;
}


Foam::tmp<Foam::volVectorField> Foam::MRFZoneList::DDt
(
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return rho*DDt(U);
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZoneList::phi() const
{
    auto tphi = surfaceScalarField::New
    (
        "phiMRF",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimVelocity*dimArea, Zero)
    );
    auto& phi = tphi.ref();

    for (const auto& mrf : *this)
    {
        mrf.makeAbsolute(phi);
    }

    return tphi;
}


Foam::tmp<Foam::volVectorField> Foam::MRFZoneList::coriolisForce
(
    const volVectorField& U
) const
{
    return coriolisForce_;
} // by Gill


Foam::tmp<Foam::volVectorField> Foam::MRFZoneList::coriolisForce
(
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return rho*coriolisForce_;
} // by Gill


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZoneList::volumeAcceleration() const
{
    return volumeAcceleration_;
} // by Gill


void Foam::MRFZoneList::makeRelative(volVectorField& U) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeRelative(U);
    }
}


void Foam::MRFZoneList::makeRelative(surfaceScalarField& phi) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeRelative(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZoneList::relative
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "relative(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeRelative(rphi.ref());

        tphi.clear();

        return rphi;
    }

    return tmp<surfaceScalarField>(tphi, true);
}


Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Foam::scalar>>
Foam::MRFZoneList::relative
(
    const tmp<FieldField<fvsPatchField, scalar>>& tphi
) const
{
    if (size())
    {
        tmp<FieldField<fvsPatchField, scalar>> rphi(New(tphi, true));

        for (const auto& mrf: *this)
        {
            mrf.makeRelative(rphi.ref());
        }

        tphi.clear();

        return rphi;
    }

    return tmp<FieldField<fvsPatchField, scalar>>(tphi, true);
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::MRFZoneList::relative
(
    const tmp<Field<scalar>>& tphi,
    const label patchi
) const
{
    if (size())
    {
        tmp<Field<scalar>> rphi(New(tphi, true));

        for (const auto& mrf: *this)
        {
            mrf.makeRelative(rphi.ref(), patchi);
        }

        tphi.clear();

        return rphi;
    }

    return tmp<Field<scalar>>(tphi, true);
}


void Foam::MRFZoneList::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeRelative(rho, phi);
    }
}


void Foam::MRFZoneList::addCoriolisFlux
(
    const surfaceScalarField& rAUf,
	surfaceScalarField& phi
) const
{
    for (const auto& mrf: *this)
    {
        mrf.addCoriolisFlux(rAUf, volumeAcceleration_, phi);
    }
} // by Gill


void Foam::MRFZoneList::addCoriolisFlux
(
    const surfaceScalarField& rhof,
    const surfaceScalarField& rhorAUf,
	surfaceScalarField& phi
) const
{
    for (const auto& mrf: *this)
    {
        mrf.addCoriolisFlux(rhof, rhorAUf, volumeAcceleration_, phi);
    }
} // by Gill


void Foam::MRFZoneList::makeAbsolute(volVectorField& U) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeAbsolute(U);
    }
}


void Foam::MRFZoneList::makeAbsolute(surfaceScalarField& phi) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeAbsolute(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFZoneList::absolute
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "absolute(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeAbsolute(rphi.ref());

        tphi.clear();

        return rphi;
    }

    return tmp<surfaceScalarField>(tphi, true);
}


void Foam::MRFZoneList::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeAbsolute(rho, phi);
    }
}


void Foam::MRFZoneList::correctBoundaryVelocity(volVectorField& U) const
{
    for (const auto& mrf: *this)
    {
        mrf.correctBoundaryVelocity(U);
    }
    
    if(mesh_.foundObject<volScalarField>("rUAp")) // by Gill
    {
        for (const auto& mrf: *this)
        {
            mrf.makeCoriolis(U, coriolisForce_);
        }
    }
}


void Foam::MRFZoneList::updateForce
(
    const surfaceScalarField& phi, 
    const volVectorField& U
) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeVolumeAcceleration(absolute(phi), U, volumeAcceleration_);
        mrf.makeCoriolis(U, coriolisForce_);
    }
} // by Gill


void Foam::MRFZoneList::updateForce
(
    const surfaceScalarField& rhof, 
    const surfaceScalarField& phi, 
    const volVectorField& U
) const
{
    for (const auto& mrf: *this)
    {
        mrf.makeVolumeAcceleration(absolute(phi/rhof), U, volumeAcceleration_);
        mrf.makeCoriolis(U, coriolisForce_);
    }
} // by Gill


void Foam::MRFZoneList::correctBoundaryFlux
(
    const volVectorField& U,
    surfaceScalarField& phi
) const
{
    FieldField<fvsPatchField, scalar> Uf
    (
        relative(mesh_.Sf().boundaryField() & U.boundaryField())
    );

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        if (isA<fixedValueFvsPatchScalarField>(phibf[patchi]))
        {
            phibf[patchi] == Uf[patchi];
        }
    }
}


void Foam::MRFZoneList::update()
{
    if (mesh_.topoChanging())
    {
        for (auto& mrf: *this)
        {
            mrf.update();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MRFZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
