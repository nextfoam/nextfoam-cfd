/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Description
    Test file reading with broadcast

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"  // For fileSize()
#include "Fstream.H"
#include "Pstream.H"
#include "SpanStream.H"
#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

bool optUseSeek = false;
bool optVerbose = false;

// Get file contents. Usually master-only and broadcast
static List<char> slurpFile
(
    const fileName& pathname,
    const bool parallel = UPstream::parRun(),
    const bool masterOnly = true
)
{
    Info<< "slurp master-only:" << masterOnly
        << " broadcast:" << (masterOnly && parallel)
        << " seek:" << optUseSeek
        << " file: " << pathname << nl;

    if (optUseSeek)
    {
        Info<< "Rewinding gzstream does not work..." << nl;
    }

    // -------------------------

    List<char> buffer;

    ifstreamPointer ifp;

    if (UPstream::master() || !masterOnly)
    {
        ifp.open(pathname);
    }

    if (ifp && ifp->good())
    {
        Info<< "compressed:"
            << (IOstreamOption::COMPRESSED == ifp.whichCompression()) << nl;

        #if 0
        uint64_t inputSize = Foam::fileSize(pathname);

        if (IOstreamOption::COMPRESSED == ifp.whichCompression())
        {
            ifp->ignore(std::numeric_limits<std::streamsize>::max());

            const std::streamsize nread = ifp->gcount();

            if (nread == std::numeric_limits<std::streamsize>::max())
            {
                FatalErrorInFunction
                    << "Failed call to ignore()" << nl
                    << exit(FatalError);
            }
            inputSize = ifp->gcount();

            if (optUseSeek)
            {
                // Rewinding gzstream does not really work...
                ifp->rdbuf()->pubseekpos(0, std::ios_base::in);
            }
            else
            {
                // Open it again - gzstream rewinding is unreliable...
                ifp.open(pathname);
            }
        }

        buffer.resize(label(inputSize));
        ifp->read(buffer.data(), buffer.size_bytes());

        const std::streamsize nread = ifp->gcount();

        if (nread == std::numeric_limits<std::streamsize>::max())
        {
            FatalErrorInFunction
                << "Failed call to read()" << nl
                << exit(FatalError);
        }

        buffer.resize(label(nread));  // Extra safety (paranoid)

        #else

        if (IOstreamOption::COMPRESSED == ifp.whichCompression())
        {
            // For compressed files we do not have any idea how large
            // the result will be. So read chunk-wise.
            // Using the compressed size for the chunk size:
            // 50% compression = 2 iterations
            // 66% compression = 3 iterations
            // ...

            const auto inputSize = Foam::fileSize(pathname + ".gz");

            const uint64_t chunkSize =
            (
                (inputSize <= 1024)
              ? uint64_t(4096)
              : uint64_t(2*inputSize)
            );

            uint64_t beg = 0;

            bool normalExit = false;

            for (int iter = 1; iter < 100000; ++iter)
            {
                if (optVerbose)
                {
                    Info<< "iter " << iter << nl;
                    Info<< "chunk " << label(chunkSize) << nl;
                    Info<< "size " << label(iter * chunkSize) << nl;
                }

                buffer.resize(label(iter * chunkSize));
                ifp->read(buffer.data() + beg, chunkSize);

                const std::streamsize nread = ifp->gcount();

                if (optVerbose)
                {
                    Info<< "nread: " << nread << nl;
                }

                if
                (
                    nread < 0
                 || nread == std::numeric_limits<std::streamsize>::max()
                )
                {
                    if (iter == 0)
                    {
                        FatalErrorInFunction
                            << "Failed call to read()" << nl
                            << exit(FatalError);
                    }
                    break;
                }
                else
                {
                    beg += uint64_t(nread);
                    if (nread >= 0 && uint64_t(nread) < chunkSize)
                    {
                        normalExit = true;
                        if (optVerbose)
                        {
                            Info<< "stopped after "
                                << iter << " iterations" << nl;
                        }
                        buffer.resize(label(beg));
                        break;
                    }
                }
            }

            if (!normalExit)
            {
                FatalErrorInFunction
                    << "Abnormal exit" << nl
                    << exit(FatalError);
            }
        }
        else
        {
            const auto inputSize = Foam::fileSize(pathname);

            if (inputSize >= 0)
            {
                buffer.resize(label(inputSize));
                ifp->read(buffer.data(), buffer.size_bytes());

                const std::streamsize nread = ifp->gcount();

                if
                (
                    nread < 0
                 || nread == std::numeric_limits<std::streamsize>::max()
                )
                {
                    FatalErrorInFunction
                        << "Failed call to read()" << nl
                        << exit(FatalError);
                }

                buffer.resize(label(nread));  // Extra safety (paranoid)
            }
        }
        #endif
    }

    // Done with input file
    ifp.reset(nullptr);

    if (parallel && masterOnly)
    {
        // On the assumption of larger files,
        // prefer two broadcasts instead of serialization
        Pstream::broadcastList(buffer);
    }

    return buffer;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noCheckProcessorDirectories();
    argList::addBoolOption("seek", "seek with gzstream (fails!)");
    argList::addVerboseOption("addition information");
    argList::addBoolOption("seek", "seek with gzstream");
    argList::addBoolOption("no-broadcast", "suppress broadcast contents");

    argList::addNote("Test master-only reading (with broadcast)");

    argList::addArgument("srcFile");

    #include "setRootCase.H"

    const bool syncPar = (UPstream::parRun() && !args.found("no-broadcast"));
    optUseSeek = args.found("seek");
    optVerbose = args.verbose();

    auto srcName = args.get<fileName>(1);

    if (srcName.has_ext("gz"))
    {
        srcName.remove_ext();
        Info<< "stripping extraneous .gz ending" << endl;
    }

    ICharStream is;

    {
        List<char> buffer(slurpFile(srcName, syncPar));

        is.swap(buffer);
    }

    Pout<< "input:" << is.capacity() << endl;

    for (string line; is.getLine(line); /*nil*/)
    {
        Pout<< "L:" << is.lineNumber() << ": " << line.c_str() << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
