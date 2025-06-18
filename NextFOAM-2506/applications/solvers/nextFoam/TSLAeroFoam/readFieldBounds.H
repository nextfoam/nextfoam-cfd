    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p.dimensions(), 0);
    dimensionedScalar pMax("pMax", p.dimensions(), 0);

    fieldBounds.lookup(p.name()) >> pMin.value() >> pMax.value();

    // Temperature bounds
//     dimensionedScalar TMin("TMin", T.dimensions(), 0);
//     dimensionedScalar TMax("TMax", T.dimensions(), 0);
//
//     fieldBounds.lookup(T.name()) >> TMin.value() >> TMax.value();

    dimensionedScalar rhoEMin("rhoEMin", rhoE.dimensions(), 0);
    dimensionedScalar rhoEMax("rhoEMax", rhoE.dimensions(), 0);

    fieldBounds.lookup("rhoE") >> rhoEMin.value() >> rhoEMax.value();

    // rho static energy bounds
//     dimensionedScalar rhoeMin("rhoeMin", rhoE.dimensions(), 0);
//     dimensionedScalar rhoeMax("rhoeMax", rhoE.dimensions(), 0);
// 
//     fieldBounds.lookup("rhoe") >> rhoeMin.value() >> rhoeMax.value();

//     dimensionedScalar eMin("eMin", e.dimensions(), 0);
//     dimensionedScalar eMax("eMax", e.dimensions(), 0);
    
    dimensionedScalar eMin("eMin", he.dimensions(), 0);
    dimensionedScalar eMax("eMax", he.dimensions(), 0);

    fieldBounds.lookup("e") >> eMin.value() >> eMax.value();

    // static enthalpy bounds
    dimensionedScalar hMin("hMin", he.dimensions(), 0);
    dimensionedScalar hMax("hMax", he.dimensions(), 0);

    fieldBounds.lookup(he.name()) >> hMin.value() >> hMax.value();

    // rho bounds
    dimensionedScalar rhoMin("rhoMin", rho.dimensions(), 0);
    dimensionedScalar rhoMax("rhoMax", rho.dimensions(), 0);

    fieldBounds.lookup(rho.name()) >> rhoMin.value() >> rhoMax.value();

    // Velocity bound
    dimensionedScalar UMax("UMax", U.dimensions(), 0);

    fieldBounds.lookup(U.name()) >> UMax.value();
    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
