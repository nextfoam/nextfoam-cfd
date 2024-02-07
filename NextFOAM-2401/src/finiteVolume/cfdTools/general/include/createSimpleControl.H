#ifdef NFOAM
    simpleControl pvCoupling(mesh);
#else  
    #ifdef DBNFOAM
        simpleControl lusgs(mesh, "LU-SGS");
    #else
        simpleControl simple(mesh);
    #endif
#endif
