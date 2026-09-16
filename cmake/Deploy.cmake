# Find all DLL dependencies to deploy on Windows.

# remove addendum release number from the ImageMagick version
string(REGEX MATCH "^[^-]+" ImageMagick_VERSION "${ImageMagick_VERSION}")

# list additional ImageMagick coders we want to deploy
set(ImageMagick_CODERS_PREFIX "$ENV{MSYSTEM_PREFIX}/lib/ImageMagick-${ImageMagick_VERSION}/modules-Q16HDRI/coders")
set(ImageMagick_CODERS
    "${ImageMagick_CODERS_PREFIX}/gif.dll"
    "${ImageMagick_CODERS_PREFIX}/heic.dll"
    "${ImageMagick_CODERS_PREFIX}/jpeg.dll"
    "${ImageMagick_CODERS_PREFIX}/jxl.dll"
    "${ImageMagick_CODERS_PREFIX}/png.dll"
    "${ImageMagick_CODERS_PREFIX}/webp.dll"
)
set(ImageMagick_DESCRIPTORS
    "${ImageMagick_CODERS_PREFIX}/gif.la"
    "${ImageMagick_CODERS_PREFIX}/heic.la"
    "${ImageMagick_CODERS_PREFIX}/jpeg.la"
    "${ImageMagick_CODERS_PREFIX}/jxl.la"
    "${ImageMagick_CODERS_PREFIX}/png.la"
    "${ImageMagick_CODERS_PREFIX}/webp.la"
)

# silence warnings about normalizing paths
cmake_policy(SET CMP0207 NEW)

file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES ${THERION} ${LOCH}
    LIBRARIES ${ImageMagick_CODERS}
    RESOLVED_DEPENDENCIES_VAR DLLS
    PRE_EXCLUDE_REGEXES "^api-ms-" "^ext-ms-"
    POST_EXCLUDE_REGEXES ".*system32/.*\\.dll"
    DIRECTORIES $ENV{PATH}
)

file(MAKE_DIRECTORY ${DLLS_DIR})

foreach(DEP ${DLLS} ${ImageMagick_CODERS} ${ImageMagick_DESCRIPTORS})
    message("Copying dependency: ${DEP}")
    file(COPY ${DEP} DESTINATION ${DLLS_DIR})
endforeach()
