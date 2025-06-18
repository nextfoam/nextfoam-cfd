while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix hEqn
    (
        (
            thermo.isotropic()
          ? -fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
          : -fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
        )
      ==
         fvOptions(rho, h)
    );

    #ifdef PIMPLE_CONTROL
    hEqn += fvm::ddt(betav*rho, h);
    #endif

    if (radiationEnabled)
    {
        hEqn -= radiation->Sh(thermo, h);
    }

    hEqn.relax();

    fvOptions.constrain(hEqn);

    hEqn.solve();

    fvOptions.correct(h);
}

thermo.correct();

if (radiationEnabled)
{
    radiation->correct();
}
