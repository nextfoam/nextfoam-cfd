#ifdef NFOAM
    simpleControl pvCoupling(mesh);
#else  
    simpleControl simple(mesh);
#endif
