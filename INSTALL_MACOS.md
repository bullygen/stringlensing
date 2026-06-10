# if you have macos, then proceed to build with this command

CXX=$(xcrun --find c++)
SDKROOT=$(xcrun --show-sdk-path)

"$CXX" -O3 -Wall -shared -std=c++17 \
    -undefined dynamic_lookup \
    -isysroot "$SDKROOT" \
    -I"$SDKROOT/usr/include/c++/v1" \
    $(python3 -m pybind11 --includes) \
    stringlensing.cpp \
    -o stringlensing$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")