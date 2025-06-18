if (thermo.includePressureWork())
{
    if (he.name() == "e")
    {
        tmp<volScalarField> tPOverRho
        (
            new volScalarField
            (
                IOobject
                (
                    "pOverRho",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar(dimPressure/dimDensity, Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        volScalarField& pOverRho = tPOverRho.ref();

        pOverRho.internalFieldRef() = p.internalField()/rho.internalField();
        pOverRho.correctBoundaryConditions();

        tmp<volScalarField> pressureWork
        (
            fvc::div(phi, tPOverRho, "div(phi," + he.name() + ")")
        );

        EEqn += pressureWork;
    }
}// Enthalpy formulation always includes pressure work

#ifdef PIMPLE_CONTROL
if (he.name() == "h" && thermo.dpdt())
{
    volScalarField dpdt(fvc::ddt(p));

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }

    EEqn -= dpdt;
}
#endif
