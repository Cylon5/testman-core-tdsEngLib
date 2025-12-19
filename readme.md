Clean: Delete tdsEngLib/build.

Configure:

Bash

mkdir build
cd build
cmake -A Win32 ..
Build:

Bash

cmake --build . --config Release
cmake --build . --config Debug