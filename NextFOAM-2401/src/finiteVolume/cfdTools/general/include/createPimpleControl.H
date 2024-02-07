#ifdef NFOAM
    pimpleControl pvCoupling(mesh);
#else
    pimpleControl pimple(mesh);
#endif
