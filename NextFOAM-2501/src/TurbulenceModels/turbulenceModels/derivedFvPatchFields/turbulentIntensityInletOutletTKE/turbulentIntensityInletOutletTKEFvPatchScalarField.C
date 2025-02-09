/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "turbulentIntensityInletOutletTKEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::
turbulentIntensityInletOutletTKEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    turbIntensity_(nullptr),
    UName_("U")
{}


Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::
turbulentIntensityInletOutletTKEFvPatchScalarField
(
    const turbulentIntensityInletOutletTKEFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    turbIntensity_(ptf.turbIntensity_.clone(p.patch())),
    UName_(ptf.UName_)
{}

Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::
turbulentIntensityInletOutletTKEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    turbIntensity_
    (
        PatchFunction1<scalar>::New(p.patch(), "turbIntensity", dict)
    ),
    UName_(dict.getOrDefault<word>("U", "U"))
{
    fvPatchFieldBase::readDict(dict);
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::
turbulentIntensityInletOutletTKEFvPatchScalarField
(
    const turbulentIntensityInletOutletTKEFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    turbIntensity_(ptf.turbIntensity_.clone(this->patch().patch())),
    UName_(ptf.UName_)
{}


Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::
turbulentIntensityInletOutletTKEFvPatchScalarField
(
    const turbulentIntensityInletOutletTKEFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    turbIntensity_(ptf.turbIntensity_.clone(this->patch().patch())),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);

    if (turbIntensity_)
    {
        turbIntensity_->autoMap(m);
    }
}


void Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);

    const turbulentIntensityInletOutletTKEFvPatchScalarField& tiptf =
        refCast<const turbulentIntensityInletOutletTKEFvPatchScalarField>(ptf);

    if (turbIntensity_)
    {
        turbIntensity_->rmap(tiptf.turbIntensity_(), addr);
    }
}


void Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    const auto& Up = patch().lookupPatchField<volVectorField>(UName_);

    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    this->refValue() = 1.5*sqr(turbIntensity_->value(t))*magSqr(Up);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentIntensityInletOutletTKEFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    if (turbIntensity_)
    {
        turbIntensity_->writeData(os);
    }

    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turbulentIntensityInletOutletTKEFvPatchScalarField
    );
}

// ************************************************************************* //
