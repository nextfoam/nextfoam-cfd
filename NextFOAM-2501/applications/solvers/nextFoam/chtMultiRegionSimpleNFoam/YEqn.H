if (Y.size() && pvCoupling.solveSpecies())
{
    if (reactionEnabled)
    {
        reaction->correct();
        Qdot.reset(new volScalarField(reaction->Qdot()));
    }

    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::div(phi, Yi)
              - fvm::laplacian(rho*(Dm[i] + Dt()), Yi)
             ==
                fvOptions(rho, Yi)
            );

            #ifdef PIMPLE_CONTROL
            YiEqn += fvm::ddt(rho, Yi);
            #endif

            if (reactionEnabled)
            {
                YiEqn -= reaction->R(Yi);
            }

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
}
