#ifdef CATCH2_V3
#include <catch2/catch_test_macros.hpp>
#else
#include <catch2/catch.hpp>
#endif
#include "thexpmodel.h"

TEST_CASE("thexpmodel_parse_wallsrc")
{
    REQUIRE(thexpmodel_parse_wallsrc("maps") == TT_WSRC_MAPS);
    REQUIRE(thexpmodel_parse_wallsrc("centreline") == TT_WSRC_CENTERLINE);
    REQUIRE(thexpmodel_parse_wallsrc("centerline") == TT_WSRC_CENTERLINE);
    REQUIRE(thexpmodel_parse_wallsrc("all") == TT_WSRC_ALL);
    REQUIRE(thexpmodel_parse_wallsrc("maps,centreline") == (TT_WSRC_MAPS | TT_WSRC_CENTERLINE));
    REQUIRE(thexpmodel_parse_wallsrc("splay-tubes,splays,scans") == (TT_WSRC_SPLAYTUBES | TT_WSRC_SPLAYS | TT_WSRC_SCANS));
    REQUIRE(thexpmodel_parse_wallsrc("maps,") == TT_WSRC_MAPS);
    REQUIRE(thexpmodel_parse_wallsrc("maps,walls") == TT_WSRC_UNKNOWN);
    REQUIRE(thexpmodel_parse_wallsrc("") == TT_WSRC_UNKNOWN);
}
