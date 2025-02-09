/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "porousBafflePressureFvPatchField.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::autoPtr<Foam::scalarIOField> 
Foam::porousBafflePressureFvPatchField::UEqnDiag_(nullptr);

Foam::autoPtr<Foam::vectorIOField> 
Foam::porousBafflePressureFvPatchField::UEqnSource_(nullptr);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porousBafflePressureFvPatchField::calcUsource()
{
    if (!UEqnDiag_ && tangentialResistance_)
    {
        UEqnDiag_.reset
        (
            new scalarIOField
            (
                IOobject
                (
                    "UEqnDiag",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                scalarField(internalField().size(), Zero)
            )
        );

        UEqnSource_.reset
        (
            new vectorIOField
            (
                IOobject
                (
                    "UEqnSource",
                    db().time().timeName(),
                    db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                vectorField(internalField().size(), Zero)
            )
        );
    }

    if 
    (
        UEqnDiag_ && tangentialResistance_ 
     && db().foundObject<volVectorField>("U")
    )
    {
        const volVectorField& U(db().lookupObject<volVectorField>("U"));
        const fvMesh& mesh(U.mesh());
        const scalarField& V = mesh.V();

        const scalar t = mesh.time().timeOutputValue();

        const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

        const label patchi(patch().index());

        const vectorField& pUbf(U.boundaryField()[patchi]);
        const scalarField nu(turbModel.nu(patchi));

        tensor D(Zero);
        tensor F(Zero);

        if (mesh.nGeometricD() == 3)
        {
            D.xx() = D_->value(t);
            F.xx() = 0.5*I_->value(t);
        }

        D.yy() = D_->value(t);
        F.yy() = 0.5*I_->value(t);

        scalarField& UEqnDiag(UEqnDiag_());
        vectorField& UEqnSource(UEqnSource_());

        const labelUList& pOwner = patch().faceCells();
        vectorField axis(patch().nf());

        forAll(pOwner, facei)
        {
            coordSystem::cylindrical csys(Zero, axis[facei]);
            tensor d(csys.transform(D));
            tensor f(csys.transform(F));

            const label& celli = pOwner[facei];

            const tensor Cd = nu[facei]*d + mag(pUbf[facei])*f;

            const scalar isoCd = tr(Cd);

            UEqnDiag[celli] = V[celli]*isoCd;
            UEqnSource[celli] = V[celli]*(((Cd - I*isoCd) & pUbf[facei]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    D_(),
    I_(),
    length_(0),
    uniformJump_(false),
    tangentialResistance_(false),
    relaxFactor_(tangentialResistance_?-1:0.05)
{
    calcUsource();
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedJumpFvPatchField<scalar>(p, iF, dict, false),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    D_(Function1<scalar>::New("D", dict, &db())),
    I_(Function1<scalar>::New("I", dict, &db())),
    length_(dict.get<scalar>("length")),
    uniformJump_(dict.getOrDefault("uniformJump", false)),
    tangentialResistance_(dict.getOrDefault("tangentialResistance", false)),
    relaxFactor_
    (
        dict.getOrDefault<scalar>("relax", tangentialResistance_?-1:0.05)
    )
{
    if (valueRequired)
    {
        if (!this->readValueEntry(dict))
        {
            this->evaluate(Pstream::commsTypes::blocking);
        }
    }

    calcUsource();
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_),
    tangentialResistance_(ptf.tangentialResistance_),
    relaxFactor_(ptf.relaxFactor_)
{
    calcUsource();
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<scalar>(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_),
    tangentialResistance_(ptf.tangentialResistance_),
    relaxFactor_(ptf.relaxFactor_)
{
    calcUsource();
}


Foam::porousBafflePressureFvPatchField::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    D_(ptf.D_.clone()),
    I_(ptf.I_.clone()),
    length_(ptf.length_),
    uniformJump_(ptf.uniformJump_),
    tangentialResistance_(ptf.tangentialResistance_),
    relaxFactor_(ptf.relaxFactor_)
{
    calcUsource();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::porousBafflePressureFvPatchField::relaxFactor() const
{
    return relaxFactor_;
}


void Foam::porousBafflePressureFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (db().lookupObject<volScalarField>(internalField().name()).relaxed())
    {
        calcUsource();
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);

    scalarField Un(phip/patch().magSf());

    if (phip.internalField().dimensions() == dimMass/dimTime)
    {
        Un /= patch().lookupPatchField<volScalarField>(rhoName_);
    }

    if (uniformJump_)
    {
        Un = gAverage(Un);
    }
    scalarField magUn(mag(Un));

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalar t = db().time().timeOutputValue();
    const scalar D = D_->value(t);
    const scalar I = I_->value(t);

    setJump
    (
        -sign(Un)
        *(
            D*turbModel.nu(patch().index())
          + I*0.5*magUn
         )*magUn*length_
    );

    if (internalField().dimensions() == dimPressure)
    {
        setJump
        (
            jump()*patch().lookupPatchField<volScalarField>(rhoName_)
        );
    }

    this->relax();

    if (debug)
    {
        scalar avePressureJump = gAverage(jump());
        scalar aveVelocity = gAverage(Un);

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << " Average pressure drop :" << avePressureJump
            << " Average velocity :" << aveVelocity
            << endl;
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    D_->writeData(os);
    I_->writeData(os);
    os.writeEntry("length", length_);
    os.writeEntry("uniformJump", uniformJump_);
    os.writeEntry("tangentialResistance", tangentialResistance_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureFvPatchField
    );
}

// ************************************************************************* //
