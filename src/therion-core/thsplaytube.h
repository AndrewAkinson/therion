/**
 * @file thsplaytube.h
 * Passage tubes from cross-section splay shots.
 */

/* Copyright (C) 2026 Therion authors
 *
 * --------------------------------------------------------------------
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 * --------------------------------------------------------------------
 */

#ifndef thsplaytube_h
#define thsplaytube_h

#include "thdb3d.h"

/**
 * Builds passage tubes from splay shots.
 *
 * At each station with one leg in and one leg out (or a single leg at
 * a passage end) the cross-section plane bisects the legs. Splays within
 * a tolerance angle of that plane are projected onto it and form the
 * section polygon. Polygons of neighbouring stations are joined into
 * tubes along the legs; tubes taper to passage end stations without
 * splays and other open tube ends are closed. Polygon corners are at the
 * true splay ends. At junctions (and stations without a section) a body
 * is built from the station's splays and the neighbouring tubes are
 * extended into it. Where a station is drawn in a plan scrap and a drawn
 * wall is beyond the furthest point of the section on that side, a virtual horizontal
 * splay from the station to that wall is added to the section.
 *
 * @param walls Output tube surfaces.
 * @param section_splays Output splay shots used for the sections (lines).
 * @param polygons Output section polygons (closed lines).
 * @param junctions Output bodies at junctions.
 * @param virtual_splays Output virtual splays to the drawn walls (lines).
 */

void thsplaytube_build(class thdatabase * db, thdb3ddata & walls, thdb3ddata & section_splays,
  thdb3ddata & polygons, thdb3ddata & junctions, thdb3ddata & virtual_splays);

#endif
