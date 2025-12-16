Clean: Delete DSEgenLib/build.

Configure:

Bash

mkdir build
cd build
cmake -A Win32 ..
Build:

Bash

cmake --build . --config Release
Result: You should find both DSEgenLib.dll and tdsEngLib.dll (copied automatically) in DSEgenLib/build/bin/Release/