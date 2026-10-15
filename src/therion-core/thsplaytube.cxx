/**
 * @file thsplaytube.cxx
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

#include "thsplaytube.h"
#include "thdatabase.h"
#include "thdb1d.h"
#include "thdataleg.h"
#include "thsurvey.h"
#include "thlog.h"
#include "thdb2dprj.h"
#include "thscrap.h"
#include "thscrapis.h"
#include "thpoint.h"
#include "thline.h"
#include "QuickHull.hpp"
#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

thsplaytube_options thsplaytube_opt;

namespace {

const double CANDIDATE_ANGLE = 45.0;  // max. angle of splays to the reference plane for plane fitting [deg]
const double SECTION_ANGLE = 25.0;  // max. angle of section splays to the fitted plane [deg]
const double TILTED_PLANE = 15.0;  // planes tilted more than this from vertical: splays judged in 3D only [deg]
const double LEVEL_SPLAY = 45.0;  // splays inclined less than this are judged by their horizontal angle [deg]
const double MAX_PLANE_TURN = 35.0;  // max. angle between fitted and reference plane [deg]
const double XSECTION_CLEAR = 0.1;  // overlap difference with the splay polygon that decides which way a drawn section faces
const double XSECTION_JOIN = 0.05;  // wall line ends of a drawn section this close are joined [m]
const double XSECTION_EDGE = 0.1;  // a drawn section's station may be this far outside the drawing (on the wall) [m]
const double JUNCTION_SITE = 0.5;  // spacing of the points along other legs that limit a junction footprint [m]
const double JUNCTION_STRAIGHT = 70.0;  // junction leg pairs bending more than this are not a through passage [deg]
const size_t JUNCTION_SPLAYS = 4;  // section splays needed for a cross section at a junction
const double STEEP = 0.7;  // |z| of the bisector above which the plane is not kept vertical
const double MAX_GAP = 150.0;  // larger angular gaps are closed at the station [deg]
const double VIRTUAL_MIN = 0.2;  // drawn wall this far beyond the section's furthest point gets a virtual splay [m]
const double OUTSIDE_REPORT = 0.5;  // splay ends this far outside the drawn walls are reported [m]
const double BRANCH_OPENING = 1.5;  // a branch stops where its drawn width grows by this factor (junction)
const double BRANCH_END = 1.0;  // and at the latest this far before the junction station [m]
const double BRANCH_OVERLAP = 0.5;  // a branch that cannot be joined to its junction runs this far into it [m]
const double ELEVATION_NEAR = 0.5;  // a leg outside the drawn elevation by up to this uses the nearest span [m]
const double ELEVATION_STEEP = 60.0;  // a steeper leg reads the extended elevation across it, not vertically [deg]
const double ELEVATION_LEVEL = 0.05;  // a leg shorter than this in plan or in the drawing gives no direction [m]
const int OVAL_POINTS = 24;  // points round an oval section of a leg without splay sections
const double OVAL_END = 1.0;  // an oval leg starts this far from a junction [m] (at most half the leg)
const double FILL_CELL = 0.1;  // grid for finding drawn plan not covered by the model [m]
const double FILL_LEVEL = 3.0;  // model pieces within this of a scrap's stations' heights count for it [m]
const int FILL_MIN_CELLS = 3;  // smaller pieces of uncovered drawn plan are ignored
const double FILL_PASSAGE_DEPTH = 0.5;  // passages without legs must be at least this deep [m]
const double WHOLE_STRETCH = 0.1;  // a side narrower than this times the other: the whole section is stretched
const double STRETCH_REPORT = 2.0;  // stretches by more than this factor are reported
const double PAIR_TOLERANCE = 20.0;  // max. position difference of corners joined directly (100 = quarter)
const double ANCHOR_TIE = 0.05;  // corners this close to an extreme are practically equal [m]

// Section polygon corner around the station: angle from u towards v and radius
// in the section plane, and offset along the plane normal (true splay end).
struct corner {
  double a, r, h;
  bool operator < (const corner & o) const { return a < o.a; }
};

struct section {
  bool ok = false;
  lxVec c, n, u, v;  // station, plane normal (along the passage), right and up directions in the plane
  double cx = 0.0, cy = 0.0;  // centre of the polygon in the plane (from the station along u, v)
  std::vector<lxVec> q;  // polygon corners (along u, v, n from the station), in polygon order
  std::vector<double> pos;  // corner positions 0-400: right 0, roof 100, left 200, floor 300 (see set_positions)
  std::vector<thdb3dvx *> vx;  // section polygon vertices
  std::string moved;  // anchors moved off a shared corner (for the log)
  int ntubes = 0;

  lxVec point(size_t k) const {
    return c + q[k].x * u + q[k].y * v + q[k].z * n;
  }
  lxVec centre() const {
    return c + cx * u + cy * v;
  }
  bool set_positions();
};

double deg(double rad);
double rad(double deg);

// Corner positions along the polygon, 0-400, anticlockwise looking along the
// passage. The rightmost (R), highest (T), leftmost (L) and lowest (B)
// corners are anchors at 0, 100, 200 and 300; between them a corner's
// position is the distance travelled along the polygon from the previous
// anchor, as a percentage of the quarter's length. Corners are rotated so
// that R is first.
// Returns false if the anchors are not in order around the polygon (the
// positions are then spread along the polygon by angle around its centre).
bool section::set_positions()
{
  size_t np = q.size(), k;
  double xmin = q[0].x, xmax = q[0].x, ymin = q[0].y, ymax = q[0].y;
  for (const lxVec & p : q) {
    xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
    ymin = std::min(ymin, p.y); ymax = std::max(ymax, p.y);
  }
  double tx = ANCHOR_TIE, ty = ANCHOR_TIE;
  // the extreme corner; among practically equal corners the one nearest the centre line
  auto anchor = [&](bool horiz, bool maxside) {
    size_t best = np;
    for (size_t i = 0; i < np; i++) {
      double val = horiz ? q[i].x : q[i].y, ext = horiz ? (maxside ? xmax : xmin) : (maxside ? ymax : ymin),
        tol = horiz ? tx : ty;
      if (fabs(val - ext) > tol) continue;
      double off = horiz ? fabs(q[i].y - cy) : fabs(q[i].x - cx),
        boff = (best < np) ? (horiz ? fabs(q[best].y - cy) : fabs(q[best].x - cx)) : 1e99;
      if (off < boff) best = i;
    }
    return best;
  };
  auto edge2 = [&](size_t k) {  // length of edge k in the section plane
    const lxVec & p0 = q[k % np], & p1 = q[(k + 1) % np];
    return std::hypot(p1.x - p0.x, p1.y - p0.y);
  };
  // anchor corners: right, roof, left, floor
  size_t a[4] = {anchor(true, true), anchor(false, true), anchor(true, false), anchor(false, false)};
  // Two anchors on one corner: in a section taller than wide the right or left
  // anchor moves, otherwise the roof or floor one. It moves along the stretch on
  // its side of the corner (towards the anchor before or after the pair) to where
  // the line at mid-height between roof and floor (mid-width between right and
  // left for roof or floor) crosses it, furthest out in its direction; if the
  // line misses the stretch, halfway along it. A point is inserted there.
  static const char * aname[4] = {"right", "roof", "left", "floor"};
  bool tall = (ymax - ymin) >= (xmax - xmin);
  for (int i = 0; i < 4; i++) {
    int j = (i + 1) % 4;
    if (a[i] != a[j]) continue;
    size_t X = a[i];
    int mv = (((i % 2) == 0) == tall) ? i : j;
    size_t s0 = (mv == i) ? a[(i + 3) % 4] : X, s1 = (mv == i) ? X : a[(i + 2) % 4];
    size_t nedges = (s1 + np - s0) % np;
    if (nedges == 0) continue;
    bool horizontal = ((mv % 2) == 0);
    double level = horizontal ? 0.5 * (q[a[1]].y + q[a[3]].y) : 0.5 * (q[a[0]].x + q[a[2]].x);
    double best = -1e99, bf = 0.0;
    size_t be = np;
    for (size_t e = 0; e < nedges; e++) {
      size_t k0 = (s0 + e) % np, k1 = (k0 + 1) % np;
      double c0 = horizontal ? q[k0].y : q[k0].x, c1 = horizontal ? q[k1].y : q[k1].x;
      if (((c0 - level) * (c1 - level) > 0.0) || (c0 == c1)) continue;
      double f = (level - c0) / (c1 - c0);
      lxVec pp = q[k0] + f * (q[k1] - q[k0]);
      double score = (mv == 0) ? pp.x : ((mv == 1) ? pp.y : ((mv == 2) ? -pp.x : -pp.y));
      if (score > best) {
        best = score;
        be = k0;
        bf = f;
      }
    }
    const char * how = "at mid-height";
    if (!horizontal) how = "at mid-width";
    if (be == np) {
      how = "halfway along";
      double len = 0.0, done = 0.0;
      for (size_t e = 0; e < nedges; e++) len += edge2(s0 + e);
      for (size_t e = 0; e < nedges; e++) {
        double el = edge2(s0 + e);
        if (done + el >= 0.5 * len) {
          be = (s0 + e) % np;
          bf = (el > 0.0) ? (0.5 * len - done) / el : 0.0;
          break;
        }
        done += el;
      }
    }
    size_t k1 = (be + 1) % np, na;
    if (bf <= 1e-9) na = be;
    else if (bf >= 1.0 - 1e-9) na = k1;
    else {
      na = be + 1;
      q.insert(q.begin() + na, q[be] + bf * (q[k1] - q[be]));
      np++;
      for (int m = 0; m < 4; m++)
        if (a[m] >= na) a[m]++;
    }
    a[mv] = na;
    moved += fmt::format(", {} anchor moved {}", aname[mv], how);
  }
  // corners from the right anchor
  size_t iR = a[0];
  std::rotate(q.begin(), q.begin() + iR, q.end());
  size_t ia[4];
  for (int m = 0; m < 4; m++) ia[m] = (a[m] + np - iR) % np;
  // an anchor that is also R (e.g. the lowest corner is the rightmost) ends the polygon
  for (int qq = 1; qq < 4; qq++)
    if ((ia[qq] == 0) && (ia[qq - 1] > 0)) ia[qq] = np;
  bool ordered = (ia[1] <= ia[2]) && (ia[2] <= ia[3]);
  pos.assign(np, 0.0);
  if (ordered) {
    // distance travelled along each quarter (in the section plane)
    auto edge = edge2;
    for (int qq = 0; qq < 4; qq++) {
      size_t k1 = (qq < 3) ? ia[qq + 1] : np;
      double len = 0.0, done = 0.0;
      for (k = ia[qq]; k < k1; k++) len += edge(k);
      for (k = ia[qq]; k < k1; k++) {
        if (k > ia[qq]) pos[k] = qq * 100.0 + ((len > 0.0) ? 100.0 * done / len : 0.0);
        done += edge(k);
      }
    }
    // every section has corners at 0, 100, 200 and 300: a corner that is
    // several anchors (e.g. the highest is also the leftmost) is listed once
    // for each, joined by zero-length edges
    std::vector<lxVec> nq;
    std::vector<double> npos;
    for (k = 0; k <= np; k++) {
      for (int qq = 0; qq < 4; qq++)
        if (ia[qq] == k) {
          nq.push_back(q[k % np]);
          npos.push_back(qq * 100.0);
        }
      bool isanchor = false;
      for (int qq = 0; qq < 4; qq++)
        if (ia[qq] == k) isanchor = true;
      if (!isanchor && (k < np)) {
        nq.push_back(q[k]);
        npos.push_back(pos[k]);
      }
    }
    q = nq;
    pos = npos;
    np = q.size();
  } else {
    double a0 = atan2(q[0].y - cy, q[0].x - cx);
    for (k = 1; k < np; k++) {
      double ang = atan2(q[k].y - cy, q[k].x - cx) - a0;
      while (ang < 0.0) ang += 2.0 * M_PI;
      pos[k] = ang / (2.0 * M_PI) * 400.0;
    }
  }
  for (k = 1; k < np; k++)
    pos[k] = std::max(pos[k], pos[k - 1]);
  return ordered;
}

double deg(double rad) {
  return rad / M_PI * 180.0;
}

double rad(double deg) {
  return deg / 180.0 * M_PI;
}

// Angle between a direction and a plane given by its unit normal [deg].
double plane_angle(const lxVec & w, const lxVec & n) {
  double wl = std::hypot(std::hypot(w.x, w.y), w.z);
  return (wl > 0.0) ? deg(asin(std::min(1.0, fabs(w * n) / wl))) : 90.0;
}

// Angle of a splay to a section plane as used for choosing section splays:
// splays within LEVEL_SPLAY of level by their direction seen from above (the
// angle between the splay and the plane's line in plan), steeper ones in 3D,
// as their direction from above means little; in 3D also for planes tilted
// more than TILTED_PLANE from vertical (steep passages).
double splay_angle(const lxVec & w, const lxVec & n) {
  double wh = std::hypot(w.x, w.y), nh = std::hypot(n.x, n.y);
  if ((wh < 1e-9) || (nh < 1e-9) || (deg(atan2(fabs(w.z), wh)) >= LEVEL_SPLAY)
      || (deg(asin(std::min(1.0, fabs(n.z)))) > TILTED_PLANE)) return plane_angle(w, n);
  return deg(asin(std::min(1.0, fabs(w.x * n.x + w.y * n.y) / (wh * nh))));
}

// Normal of the plane through the origin best fitting the unit directions
// (eigenvector of the smallest eigenvalue of their scatter matrix, Jacobi method).
lxVec fit_plane(const std::vector<lxVec> & dirs) {
  double a[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (const lxVec & w : dirs) {
    lxVec u = w;
    u.Normalize();
    double c[3] = {u.x, u.y, u.z};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        a[i][j] += c[i] * c[j];
  }
  for (int sweep = 0; sweep < 50; sweep++) {
    double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off < 1e-12) break;
    for (int p = 0; p < 2; p++) {
      for (int q = p + 1; q < 3; q++) {
        if (fabs(a[p][q]) < 1e-15) continue;
        double th = 0.5 * atan2(2.0 * a[p][q], a[q][q] - a[p][p]), c = cos(th), sn = sin(th);
        for (int k = 0; k < 3; k++) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - sn * akq;
          a[k][q] = sn * akp + c * akq;
        }
        for (int k = 0; k < 3; k++) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - sn * aqk;
          a[q][k] = sn * apk + c * aqk;
        }
        for (int k = 0; k < 3; k++) {
          double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - sn * vkq;
          v[k][q] = sn * vkp + c * vkq;
        }
      }
    }
  }
  int m = 0;
  for (int i = 1; i < 3; i++)
    if (a[i][i] < a[m][m]) m = i;
  lxVec n(v[0][m], v[1][m], v[2][m]);
  n.Normalize();
  return n;
}

double cross2(double ax, double ay, double bx, double by) {
  return ax * by - ay * bx;
}

// Adds triangle a, b, c to the face, ordered so that it faces along out
// (outwards), so that viewers light and show all faces correctly.
void tri(thdb3dfc * fc, thdb3dvx * a, thdb3dvx * b, thdb3dvx * c, const lxVec & out) {
  lxVec n = (b->get_vector() - a->get_vector()) ^ (c->get_vector() - a->get_vector());
  fc->insert_vertex(a);
  if (n * out >= 0.0) {
    fc->insert_vertex(b);
    fc->insert_vertex(c);
  } else {
    fc->insert_vertex(c);
    fc->insert_vertex(b);
  }
}

// Triangles of one piece of surface (built as one strip or fan, so all
// ordered the same way): they all get the order that most of them, weighted by
// size, need to face outwards.
struct piece {
  std::vector<std::array<thdb3dvx *, 3>> t;
  double vote = 0.0;
  void add(thdb3dvx * a, thdb3dvx * b, thdb3dvx * c, const lxVec & out) {
    lxVec n = (b->get_vector() - a->get_vector()) ^ (c->get_vector() - a->get_vector());
    vote += n * out;
    t.push_back({a, b, c});
  }
  void emit(thdb3dfc * fc) {
    bool flip = (vote < 0.0);
    for (const auto & x : t) {
      fc->insert_vertex(x[0]);
      fc->insert_vertex(flip ? x[2] : x[1]);
      fc->insert_vertex(flip ? x[1] : x[2]);
    }
  }
};

lxVec centroid(const std::vector<thdb3dvx *> & r) {
  lxVec c(0.0, 0.0, 0.0);
  for (thdb3dvx * v : r) c += v->get_vector();
  return (r.empty()) ? c : c / double(r.size());
}

// Joins two rings with corresponding vertices by triangles facing away from
// the line through their centres.
void loft_rings(thdb3dfc * fc, const std::vector<thdb3dvx *> & ra, const std::vector<thdb3dvx *> & rb) {
  lxVec ca = centroid(ra), cb = centroid(rb), axis = cb - ca;
  double al = axis.Length();
  if (al > 1e-9) axis = axis / al;
  size_t nv = std::min(ra.size(), rb.size());
  piece pc;
  for (size_t k = 0; k < nv; k++) {
    size_t k1 = (k + 1) % nv;
    lxVec m = 0.25 * (ra[k]->get_vector() + ra[k1]->get_vector() + rb[k1]->get_vector() + rb[k]->get_vector());
    lxVec out = m - (ca + ((m - ca) * axis) * axis);
    pc.add(ra[k], ra[k1], rb[k1], out);
    pc.add(ra[k], rb[k1], rb[k], out);
  }
  pc.emit(fc);
}

// Joins a ring to a point (a passage end) by triangles facing away from the line ring centre - point.
void taper_ring(thdb3dfc * fc, const std::vector<thdb3dvx *> & r, thdb3dvx * end) {
  lxVec c = centroid(r), axis = end->get_vector() - c;
  double al = axis.Length();
  if (al > 1e-9) axis = axis / al;
  piece pc;
  for (size_t k = 0, nv = r.size(); k < nv; k++) {
    lxVec m = (r[k]->get_vector() + r[(k + 1) % nv]->get_vector() + end->get_vector()) / 3.0;
    pc.add(r[k], r[(k + 1) % nv], end, m - (c + ((m - c) * axis) * axis));
  }
  pc.emit(fc);
}

// Closes a tube end: triangulates the section polygon (any shape) with the
// given vertices for its corners. Falls back to a fan from the station, which
// sees every corner of a polygon made of its splays.
bool cap(const section & s, const std::vector<thdb3dvx *> & ring, thdb3dvx * station, thdb3ddata & out,
  const lxVec & outward) {
  std::vector<std::unique_ptr<p2t::Point>> holder;
  std::vector<p2t::Point *> line;
  std::map<p2t::Point *, thdb3dvx *> vmap;
  for (size_t k = 0; k < s.q.size(); k++) {
    const lxVec & p = s.q[k];
    if (!line.empty() && (std::hypot(p.x - line.back()->x, p.y - line.back()->y) < 1e-6)) continue;
    if ((k + 1 == s.q.size()) && !line.empty() && (std::hypot(p.x - line[0]->x, p.y - line[0]->y) < 1e-6)) continue;
    holder.push_back(std::make_unique<p2t::Point>(p.x, p.y));
    line.push_back(holder.back().get());
    vmap[holder.back().get()] = ring[k];
  }
  thdb3dfc * fc = out.insert_face(THDB3DFC_TRIANGLES);
  bool done = false;
  if (line.size() >= 3) {
    try {
      p2t::CDT cdt(line);
      cdt.Triangulate();
      piece pc;
      for (p2t::Triangle * t : cdt.GetTriangles())
        pc.add(vmap.at(t->GetPoint(0)), vmap.at(t->GetPoint(1)), vmap.at(t->GetPoint(2)), outward);
      pc.emit(fc);
      done = true;
    } catch (...) {
      done = false;
    }
  }
  if (!done) {
    size_t nv = ring.size();
    piece pc;
    for (size_t k = 0; k < nv; k++)
      pc.add(station, ring[(k + 1) % nv], ring[k], outward);
    pc.emit(fc);
  }
  return done;
}

// Walls of the plan drawing: outline edges of each selected plan scrap in
// model coordinates (x, y), and the scraps in which each station is drawn.
typedef std::vector<std::array<double, 2>> ring2;

// point in a closed ring (even-odd)
bool in_ring(const ring2 & r, double x, double y) {
  bool in = false;
  for (size_t k = 0, n = r.size(); k < n; k++) {
    const auto & a = r[k], & b = r[(k + 1) % n];
    if (((a[1] > y) != (b[1] > y)) && (x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0])) in = !in;
  }
  return in;
}

// first crossing of the ray from (x, y) along (dx, dy) with a ring: distance and edge
bool ray_ring(const ring2 & r, double x, double y, double dx, double dy, double & dist, size_t & edge) {
  bool found = false;
  for (size_t k = 0, n = r.size(); k < n; k++) {
    const auto & a = r[k], & b = r[(k + 1) % n];
    double ex = b[0] - a[0], ey = b[1] - a[1], den = dx * ey - dy * ex;
    if (std::fabs(den) < 1e-12) continue;
    double t = ((a[0] - x) * ey - (a[1] - y) * ex) / den, u = ((a[0] - x) * dy - (a[1] - y) * dx) / den;
    if ((t > 1e-6) && (u >= 0.0) && (u <= 1.0) && (!found || (t < dist))) {
      dist = t;
      edge = k;
      found = true;
    }
  }
  return found;
}

struct drawn_walls {
  std::vector<std::vector<std::array<double, 4>>> scraps;
  std::vector<std::vector<ring2>> rings;  // outline rings of each scrap
  std::map<size_t, std::vector<size_t>> station_scraps;  // station index -> scraps
  std::vector<double> zlo, zhi;  // heights of the stations drawn in each scrap

  // a scrap at the level of height z (within FILL_LEVEL of its stations; any
  // scrap without stations)
  bool at_level(size_t si, double z) const {
    return (zlo[si] > zhi[si]) || ((z >= zlo[si] - FILL_LEVEL) && (z <= zhi[si] + FILL_LEVEL));
  }

  void collect(thdatabase * db) {
    thdb2dprjpr prjid = db->db2d.parse_projection("plan", false);
    if (prjid.newprj) return;
    db->db2d.process_projection(prjid.prj);
    for (thscrap * cs = prjid.prj->first_scrap; cs != NULL; cs = cs->proj_next_scrap) {
      if (!cs->fsptr->is_selected()) continue;
      thscrapis is(cs);
      is.outline_scan(cs->get_outline());
      std::vector<std::array<double, 4>> edges;
      std::vector<ring2> rr;
      double sx = cs->proj->shift_x, sy = cs->proj->shift_y;
      for (thscrapisolpt * ol = is.firstolseg; ol != NULL; ol = ol->next_segment) {
        ring2 r;
        for (thscrapisolpt * p = ol; p != NULL; p = p->next) {
          thscrapisolpt * q = (p->next != NULL) ? p->next : ol;
          if (q != p) edges.push_back({p->x + sx, p->y + sy, q->x + sx, q->y + sy});
          r.push_back({p->x + sx, p->y + sy});
        }
        if (r.size() >= 3) rr.push_back(r);
      }
      if (edges.empty()) continue;
      size_t si = scraps.size();
      scraps.push_back(edges);
      rings.push_back(rr);
      for (th2ddataobject * o2 = cs->fs2doptr; o2 != NULL; o2 = o2->nscrapoptr) {
        if (o2->get_class_id() != TT_POINT_CMD) continue;
        thpoint * pp = dynamic_cast<thpoint*>(o2);
        if (pp->station_name.id == 0) continue;
        size_t k = db->db1d.station_vec[pp->station_name.id - 1].uid - 1;
        std::vector<size_t> & v = station_scraps[k];
        if (std::find(v.begin(), v.end(), si) == v.end()) v.push_back(si);
      }
    }
    zlo.assign(scraps.size(), 1e99);
    zhi.assign(scraps.size(), -1e99);
    for (auto & ss : station_scraps)
      for (size_t si : ss.second) {
        zlo[si] = std::min(zlo[si], db->db1d.station_vec[ss.first].z);
        zhi[si] = std::max(zhi[si], db->db1d.station_vec[ss.first].z);
      }
  }

  bool inside(size_t si, double x, double y) const {
    bool in = false;
    for (const auto & e : scraps[si])
      if (((e[1] > y) != (e[3] > y)) && (x < (e[2] - e[0]) * (y - e[1]) / (e[3] - e[1]) + e[0])) in = !in;
    return in;
  }

  // nearest point on the drawn walls of a scrap
  void nearest(size_t si, double x, double y, double & nx, double & ny) const {
    double best = 1e99;
    for (const auto & e : scraps[si]) {
      double dx = e[2] - e[0], dy = e[3] - e[1], l2 = dx * dx + dy * dy;
      double t = (l2 > 0.0) ? std::max(0.0, std::min(1.0, ((x - e[0]) * dx + (y - e[1]) * dy) / l2)) : 0.0;
      double px = e[0] + t * dx, py = e[1] + t * dy, d = std::hypot(px - x, py - y);
      if (d < best) { best = d; nx = px; ny = py; }
    }
  }

  double wall_distance(size_t si, double x, double y) const {
    double best = 1e99;
    for (const auto & e : scraps[si]) {
      double dx = e[2] - e[0], dy = e[3] - e[1], l2 = dx * dx + dy * dy;
      double t = (l2 > 0.0) ? std::max(0.0, std::min(1.0, ((x - e[0]) * dx + (y - e[1]) * dy) / l2)) : 0.0;
      best = std::min(best, std::hypot(e[0] + t * dx - x, e[1] + t * dy - y));
    }
    return best;
  }

  // distance from (x, y) along (dx, dy) to the first wall of the scrap
  double ray(size_t si, double x, double y, double dx, double dy) const {
    double best = -1.0;
    for (const auto & e : scraps[si]) {
      double ex = e[2] - e[0], ey = e[3] - e[1], den = dx * ey - dy * ex;
      if (std::fabs(den) < 1e-12) continue;
      double t = ((e[0] - x) * ey - (e[1] - y) * ex) / den, u = ((e[0] - x) * dy - (e[1] - y) * dx) / den;
      if ((t > 1e-6) && (u >= 0.0) && (u <= 1.0) && ((best < 0.0) || (t < best))) best = t;
    }
    return best;
  }

  // drawn walls either side of a point at height z, along the horizontal
  // direction h: distances right and left (scraps at that level only)
  bool width(double x, double y, double z, const lxVec & h, double & right, double & left) const {
    right = left = -1.0;
    for (size_t si = 0; si < scraps.size(); si++) {
      if (!at_level(si, z) || !inside(si, x, y)) continue;
      double tr = ray(si, x, y, h.x, h.y), tl = ray(si, x, y, -h.x, -h.y);
      if ((tr > 0.0) && ((right < 0.0) || (tr < right))) right = tr;
      if ((tl > 0.0) && ((left < 0.0) || (tl < left))) left = tl;
    }
    return (right > 0.0) && (left > 0.0);
  }

  // drawn scraps containing the station: those where it is drawn, else any at
  // its level whose walls enclose it
  std::vector<size_t> around(size_t k, double x, double y, double z) const {
    std::vector<size_t> res;
    auto it = station_scraps.find(k);
    if (it != station_scraps.end())
      for (size_t si : it->second)
        if (inside(si, x, y)) res.push_back(si);
    if (res.empty())
      for (size_t si = 0; si < scraps.size(); si++)
        if (at_level(si, z) && inside(si, x, y)) res.push_back(si);
    return res;
  }
};

// Outlines of the extended elevation scraps, in extended coordinates (x along
// the extended elevation, y = height), and the scraps in which each station is drawn.
struct drawn_elevation {
  std::vector<std::vector<ring2>> rings;
  std::map<size_t, std::vector<size_t>> station_scraps;
  // the same scraps as drawn: scaled by their legs (-scale only without them)
  // and flipped only, not rotated or fitted to the stations (m); and the station
  // points drawn in each
  std::vector<std::vector<ring2>> raw;
  struct dpt { size_t si; double x, y; };
  std::map<size_t, std::vector<dpt>> drawn;

  void collect(thdatabase * db) {
    thdb2dprjpr prjid = db->db2d.parse_projection("extended", false);
    if (prjid.newprj) return;
    db->db2d.process_projection(prjid.prj);
    for (thscrap * cs = prjid.prj->first_scrap; cs != NULL; cs = cs->proj_next_scrap) {
      if (!cs->fsptr->is_selected()) continue;
      thscrapis is(cs);
      is.outline_scan(cs->get_outline());
      std::vector<ring2> rr;
      double sx = cs->proj->shift_x, sz = cs->proj->shift_z;
      for (thscrapisolpt * ol = is.firstolseg; ol != NULL; ol = ol->next_segment) {
        ring2 r;
        for (thscrapisolpt * p = ol; p != NULL; p = p->next) r.push_back({p->x + sx, p->y + sz});
        if (r.size() >= 3) rr.push_back(r);
      }
      if (rr.empty()) continue;
      size_t si = rings.size();
      rings.push_back(rr);
      // the outline as drawn, curves followed: the raw drawing with the scrap's
      // own -scale and -flip only (as Therion does before fitting it to the stations)
      double ms = cs->scale, mr = 0.0;
      if (cs->scale_p9) {
        double pd = std::hypot(cs->scale_p2x - cs->scale_p1x, cs->scale_p2y - cs->scale_p1y),
          rd = std::hypot(cs->scale_r2x - cs->scale_r1x, cs->scale_r2y - cs->scale_r1y);
        ms = (pd > 0.0) ? rd / pd : 1.0;
        mr = atan2(cs->scale_r2y - cs->scale_r1y, cs->scale_r2x - cs->scale_r1x) - atan2(cs->scale_p2y - cs->scale_p1y, cs->scale_p2x - cs->scale_p1x);
      }
      double ks = 1.0;  // the scale from the stations' legs
      auto drawn_xy = [&](const thdb2dpt * pt) {
        std::array<double, 2> q = {ms * (cos(mr) * pt->x - sin(mr) * pt->y), ms * (sin(mr) * pt->x + cos(mr) * pt->y)};
        if (cs->flip == TT_SCRAP_FLIP_HORIZ) q[0] = -q[0];
        if (cs->flip == TT_SCRAP_FLIP_VERT) q[1] = -q[1];
        return std::array<double, 2>{ks * q[0], ks * q[1]};
      };
      // the station points drawn in the scrap; then the scale from them, as Therion
      // scales a scrap by its stations (-scale only without them): the median of
      // its legs' true lengths in the elevation to their drawn lengths (the same
      // however the scrap is laid out; legs with a station drawn twice not counted)
      for (th2ddataobject * o2 = cs->fs2doptr; o2 != NULL; o2 = o2->nscrapoptr) {
        if (o2->get_class_id() != TT_POINT_CMD) continue;
        thpoint * pp = dynamic_cast<thpoint*>(o2);
        if (pp->station_name.id == 0) continue;
        size_t k = db->db1d.station_vec[pp->station_name.id - 1].uid - 1;
        std::vector<size_t> & v = station_scraps[k];
        if (std::find(v.begin(), v.end(), si) == v.end()) v.push_back(si);
        if (pp->type == TT_POINT_TYPE_STATION)
          drawn[k].push_back({si, drawn_xy(pp->point)[0], drawn_xy(pp->point)[1]});
      }
      auto once = [&](size_t k) -> const dpt * {
        const dpt * r = NULL;
        auto it = drawn.find(k);
        if (it != drawn.end())
          for (const dpt & p : it->second)
            if (p.si == si) {
              if (r != NULL) return NULL;
              r = &p;
            }
        return r;
      };
      std::vector<double> rs;
      for (thdb1dl & l : db->db1d.leg_vec) {
        if ((l.leg->flags & (TT_LEGFLAG_SPLAY | TT_LEGFLAG_SURFACE)) != 0) continue;
        size_t a = db->db1d.station_vec[l.leg->from.id - 1].uid - 1, b = db->db1d.station_vec[l.leg->to.id - 1].uid - 1;
        const dpt * pa = once(a), * pb = once(b);
        if ((pa == NULL) || (pb == NULL)) continue;
        const thdb1ds & sa = db->db1d.station_vec[a], & sb = db->db1d.station_vec[b];
        double d = std::hypot(pb->x - pa->x, pb->y - pa->y), t = std::hypot(std::hypot(sb.x - sa.x, sb.y - sa.y), sb.z - sa.z);
        if ((d > 0.0) && (t > 0.0)) rs.push_back(t / d);
      }
      if (!rs.empty()) {
        std::nth_element(rs.begin(), rs.begin() + rs.size() / 2, rs.end());
        ks = rs[rs.size() / 2];
        for (auto & d : drawn)
          for (dpt & p : d.second)
            if (p.si == si) { p.x *= ks; p.y *= ks; }
      }
      std::vector<ring2> rw;
      for (thscraplo * lo = cs->get_outline(); lo != NULL; lo = lo->next_outline) {
        ring2 r;
        for (thscraplo * l2 = lo; l2 != NULL; l2 = l2->next_line) {
          ring2 pts;
          thdb2dlp * prev = NULL;
          for (thdb2dlp * lp = l2->line->first_point; lp != NULL; lp = lp->nextlp) {
            if ((prev != NULL) && (lp->cp1 != NULL) && (lp->cp2 != NULL)) {
              auto p0 = drawn_xy(prev->point), c1 = drawn_xy(lp->cp1), c2 = drawn_xy(lp->cp2), p1 = drawn_xy(lp->point);
              for (int k = 1; k < 20; k++) {
                double t = k / 20.0, a = (1 - t) * (1 - t) * (1 - t), b = 3 * t * (1 - t) * (1 - t), c = 3 * t * t * (1 - t), d = t * t * t;
                pts.push_back({a * p0[0] + b * c1[0] + c * c2[0] + d * p1[0], a * p0[1] + b * c1[1] + c * c2[1] + d * p1[1]});
              }
            }
            pts.push_back(drawn_xy(lp->point));
            prev = lp;
          }
          if (l2->mode == TT_OUTLINE_REVERSED) std::reverse(pts.begin(), pts.end());
          for (auto & q : pts)
            if (r.empty() || (std::hypot(q[0] - r.back()[0], q[1] - r.back()[1]) > 1e-6)) r.push_back(q);
        }
        if (r.size() >= 3) rw.push_back(r);
      }
      raw.push_back(rw);
    }
  }

  // all drawn spans (floor, roof) at x
  std::vector<std::pair<double, double>> spans(const std::vector<size_t> & scraps, double x) const {
    std::vector<std::pair<double, double>> res;
    for (size_t si : scraps) {
      std::vector<double> ys;
      for (const ring2 & r : rings[si])
        for (size_t k = 0, n = r.size(); k < n; k++) {
          const auto & a = r[k], & b = r[(k + 1) % n];
          if (((a[0] - x) * (b[0] - x) > 0.0) || (a[0] == b[0]) || (b[0] == x)) continue;
          ys.push_back(a[1] + (x - a[0]) / (b[0] - a[0]) * (b[1] - a[1]));
        }
      std::sort(ys.begin(), ys.end());
      for (size_t k = 0; k + 1 < ys.size(); k += 2) res.push_back(std::make_pair(ys[k], ys[k + 1]));
    }
    return res;
  }

  // the leg i -> j as drawn: of all pairs of drawn points of its two stations
  // in one scrap (a station may be drawn twice: in two surveys' scraps where they
  // join or overlap, or at both ends of a loop), the pair that best fits the leg
  // (hl along, dz up), and the drawn position at fraction f along the leg
  // (along: how far j is drawn from i along the scrap, signed)
  bool drawn_leg(size_t i, size_t j, double f, double hl, double dz, size_t & si, double & x, double & y,
                 double * fit = NULL, double * along = NULL) const {
    auto ii = drawn.find(i), jj = drawn.find(j);
    if ((ii == drawn.end()) || (jj == drawn.end())) return false;
    double best = 1e99;
    const dpt * bi = NULL, * bj = NULL;
    for (const dpt & a : ii->second)
      for (const dpt & b : jj->second) {
        if ((a.si != b.si) || raw[a.si].empty()) continue;
        double sc = std::hypot(std::fabs(b.x - a.x) - hl, (b.y - a.y) - dz);
        if (sc < best) { best = sc; bi = &a; bj = &b; }
      }
    if (bi == NULL) return false;
    si = bi->si;
    x = bi->x + f * (bj->x - bi->x);
    y = bi->y + f * (bj->y - bi->y);
    if (fit != NULL) *fit = best;
    if (along != NULL) *along = bj->x - bi->x;
    return true;
  }

  // all spans of a scrap as drawn on the line through (x, y) in direction
  // (dx, dy) (a unit vector; (0, 1) = floor and roof above and below x), as
  // distances along it from (x, y)
  std::vector<std::pair<double, double>> raw_spans(size_t si, double x, double y, double dx, double dy) const {
    std::vector<std::pair<double, double>> res;
    std::vector<double> ts;
    for (const ring2 & r : raw[si])
      for (size_t k = 0, n = r.size(); k < n; k++) {
        const auto & a = r[k], & b = r[(k + 1) % n];
        double sa = dx * (a[1] - y) - dy * (a[0] - x), sb = dx * (b[1] - y) - dy * (b[0] - x);
        if ((sa * sb > 0.0) || (sa == sb) || (sb == 0.0)) continue;  // a vertex is counted once
        double u = sa / (sa - sb);
        ts.push_back((a[0] + u * (b[0] - a[0]) - x) * dx + (a[1] + u * (b[1] - a[1]) - y) * dy);
      }
    std::sort(ts.begin(), ts.end());
    for (size_t k = 0; k + 1 < ts.size(); k += 2) res.push_back(std::make_pair(ts[k], ts[k + 1]));
    return res;
  }

  // the drawn span (floor, roof) at x containing height z, or the nearest within ELEVATION_NEAR
  bool span(const std::vector<size_t> & scraps, double x, double z, double & floor, double & roof) const {
    double best = -1.0;
    for (size_t si : scraps) {
      std::vector<double> ys;
      for (const ring2 & r : rings[si])
        for (size_t k = 0, n = r.size(); k < n; k++) {
          const auto & a = r[k], & b = r[(k + 1) % n];
          if (((a[0] - x) * (b[0] - x) > 0.0) || (a[0] == b[0])) continue;
          if (b[0] == x) continue;  // count a vertex once
          ys.push_back(a[1] + (x - a[0]) / (b[0] - a[0]) * (b[1] - a[1]));
        }
      std::sort(ys.begin(), ys.end());
      for (size_t k = 0; k + 1 < ys.size(); k += 2) {
        double d = (z < ys[k]) ? ys[k] - z : ((z > ys[k + 1]) ? z - ys[k + 1] : 0.0);
        if ((d <= ELEVATION_NEAR) && ((best < 0.0) || (d < best))) {
          best = d;
          floor = ys[k];
          roof = ys[k + 1];
        }
      }
    }
    return best >= 0.0;
  }
};

// Body around a station from its splays: the hull of the splay directions
// with corners at the true splay ends (as therion's splay walls, but without legs).
bool splay_body(const lxVec & c, const std::vector<lxVec> & splays, thdb3ddata & out) {
  if (splays.size() < 3) return false;
  std::vector<quickhull::Vector3<double>> dirs;
  std::vector<lxVec> ends;
  dirs.push_back(quickhull::Vector3<double>(0.0, 0.0, 0.0));
  ends.push_back(c);
  for (const lxVec & w : splays) {
    lxVec d = w;
    d.Normalize();
    dirs.push_back(quickhull::Vector3<double>(d.x, d.y, d.z));
    ends.push_back(c + w);
  }
  quickhull::QuickHull<double> qh;
  auto hull = qh.getConvexHull(dirs, true, true);
  auto idx = hull.getIndexBuffer();
  if (idx.size() < 3) return false;
  std::vector<thdb3dvx *> vx(ends.size(), nullptr);
  thdb3dfc * fc = out.insert_face(THDB3DFC_TRIANGLES);
  for (size_t k = 0; k < idx.size(); k++)
    if (vx[idx[k]] == nullptr) vx[idx[k]] = out.insert_vertex(ends[idx[k]]);
  piece pc;  // the hull is ordered consistently; turn it outwards as a whole
  for (size_t k = 0; k + 2 < idx.size(); k += 3) {
    lxVec m = (ends[idx[k]] + ends[idx[k + 1]] + ends[idx[k + 2]]) / 3.0;
    pc.add(vx[idx[k]], vx[idx[k + 1]], vx[idx[k + 2]], m - c);
  }
  pc.emit(fc);
  return true;
}

}



namespace {

// Cross-sections drawn in the plan: scraps with -projection none, each placed
// with a `point section -scrap` and cut across the passage by a `line section`
// (the one just before the point in the plan scrap). The section's station is
// a station point in the section scrap if there is one, else the station point
// at the middle of the cut. Where the station is in the section drawing: that
// station point, else the `point section` if it lies within the drawing (the
// section drawn in the plan's own drawing units, as SexyTopo exports it).
struct drawn_xsection {
  std::string scrap, plan;  // section scrap, and the plan scrap placing it
  long station = -1;  // station index, -1 if not found
  std::string linked;  // how the station was found
  bool placed = false;  // where the station is in the drawing is known
  std::string at;  // how that was found
  lxVec cut;  // horizontal direction of the cut across the passage (model)
  std::vector<std::array<double, 2>> pts;  // outline (m) from the station: across (cut direction), up
};

// metres per drawing unit of a scrap (its -scale), 0 if not known
double scrap_metres(thscrap * cs) {
  if (cs->scale_p9) {
    double p = std::hypot(cs->scale_p2x - cs->scale_p1x, cs->scale_p2y - cs->scale_p1y),
      r = std::hypot(cs->scale_r2x - cs->scale_r1x, cs->scale_r2y - cs->scale_r1y);
    return (p > 0.0) ? r / p : 0.0;
  }
  return 0.0;
}

// points of a line (raw drawing units), curves followed
std::vector<std::array<double, 2>> line_points(thline * ln) {
  std::vector<std::array<double, 2>> res;
  thdb2dlp * prev = NULL;
  for (thdb2dlp * lp = ln->first_point; lp != NULL; lp = lp->nextlp) {
    if ((prev != NULL) && (lp->cp1 != NULL) && (lp->cp2 != NULL)) {
      double x0 = prev->point->x, y0 = prev->point->y;
      for (int k = 1; k < 6; k++) {
        double t = k / 6.0, a = (1 - t) * (1 - t) * (1 - t), b = 3 * t * (1 - t) * (1 - t), c = 3 * t * t * (1 - t), d = t * t * t;
        res.push_back({a * x0 + b * lp->cp1->x + c * lp->cp2->x + d * lp->point->x,
                       a * y0 + b * lp->cp1->y + c * lp->cp2->y + d * lp->point->y});
      }
    }
    res.push_back({lp->point->x, lp->point->y});
    prev = lp;
  }
  return res;
}

void collect_xsections(thdatabase * db, std::vector<drawn_xsection> & out) {
  thdb2dprjpr prjid = db->db2d.parse_projection("plan", false);
  if (prjid.newprj) return;
  db->db2d.process_projection(prjid.prj);
  auto station_of = [&](thpoint * pp) {
    return (long) db->db1d.station_vec[pp->station_name.id - 1].uid - 1;
  };
  for (thscrap * cs = prjid.prj->first_scrap; cs != NULL; cs = cs->proj_next_scrap) {
    if (!cs->fsptr->is_selected()) continue;
    double mp = scrap_metres(cs);
    thline * cut = NULL;
    for (th2ddataobject * o2 = cs->fs2doptr; o2 != NULL; o2 = o2->nscrapoptr) {
      if (o2->get_class_id() == TT_LINE_CMD) {
        thline * ln = dynamic_cast<thline*>(o2);
        if ((ln->type == TT_LINE_TYPE_SECTION) && (ln->first_point != NULL) && (ln->first_point != ln->last_point)) cut = ln;
        continue;
      }
      if (o2->get_class_id() != TT_POINT_CMD) continue;
      thpoint * sp = dynamic_cast<thpoint*>(o2);
      if ((sp->type != TT_POINT_TYPE_SECTION) || (sp->get_scrap() == nullptr)) continue;
      thscrap * xs = sp->get_scrap();
      drawn_xsection d;
      d.scrap = xs->name;
      d.plan = cs->name;
      double xm = scrap_metres(xs);
      // the station: a station point in the section scrap, else at the middle of the cut
      thpoint * xst = NULL;
      for (th2ddataobject * x2 = xs->fs2doptr; x2 != NULL; x2 = x2->nscrapoptr)
        if ((x2->get_class_id() == TT_POINT_CMD) && (dynamic_cast<thpoint*>(x2)->type == TT_POINT_TYPE_STATION) &&
            (dynamic_cast<thpoint*>(x2)->station_name.id != 0)) xst = dynamic_cast<thpoint*>(x2);
      if (xst != NULL) {
        d.station = station_of(xst);
        d.linked = "station point in the section";
      }
      if (cut != NULL) {
        double mx = 0.5 * (cut->first_point->point->x + cut->last_point->point->x),
          my = 0.5 * (cut->first_point->point->y + cut->last_point->point->y), best = 1e99;
        thpoint * bp = NULL;
        for (th2ddataobject * x2 = cs->fs2doptr; x2 != NULL; x2 = x2->nscrapoptr) {
          if (x2->get_class_id() != TT_POINT_CMD) continue;
          thpoint * pp = dynamic_cast<thpoint*>(x2);
          if ((pp->type != TT_POINT_TYPE_STATION) || (pp->station_name.id == 0)) continue;
          double dd = std::hypot(pp->point->x - mx, pp->point->y - my);
          if (dd < best) { best = dd; bp = pp; }
        }
        if ((bp != NULL) && (mp > 0.0) && (best * mp < 0.5)) {
          long k = station_of(bp);
          if (d.station < 0) { d.station = k; d.linked = "station at the middle of the cut"; }
          else if (d.station != k) d.linked += " (the cut is at another station)";
          else d.linked += " and at the middle of the cut";
        }
        lxVec a(cut->first_point->point->xt, cut->first_point->point->yt, 0.0), b(cut->last_point->point->xt, cut->last_point->point->yt, 0.0);
        d.cut = b - a;
        if (d.cut.Length() > 1e-9) d.cut.Normalize();
      }
      // the outline, from the station's place in the drawing: the wall lines
      // (of any subtype, as a debris floor) joined end to end, from the longest;
      // a gap left is closed straight across
      std::vector<std::vector<std::array<double, 2>>> wls;
      for (th2ddataobject * x2 = xs->fs2doptr; x2 != NULL; x2 = x2->nscrapoptr)
        if ((x2->get_class_id() == TT_LINE_CMD) && (dynamic_cast<thline*>(x2)->type == TT_LINE_TYPE_WALL)) {
          std::vector<std::array<double, 2>> lp = line_points(dynamic_cast<thline*>(x2));
          if (lp.size() > 1) wls.push_back(lp);
        }
      size_t nwalls = wls.size(), njoined = 0;
      std::vector<std::array<double, 2>> raw;
      if (!wls.empty()) {
        std::sort(wls.begin(), wls.end(), [](const std::vector<std::array<double, 2>> & a, const std::vector<std::array<double, 2>> & b) { return a.size() > b.size(); });
        raw = wls[0];
        std::vector<bool> done(wls.size(), false);
        done[0] = true;
        njoined = 1;
        double tol = (xm > 0.0) ? XSECTION_JOIN / xm : 0.0;
        for (bool more = true; more;) {
          more = false;
          for (size_t w = 1; w < wls.size(); w++) {
            if (done[w]) continue;
            const auto & e = raw.back();
            const auto & a = wls[w].front(), & b = wls[w].back();
            if (std::hypot(a[0] - e[0], a[1] - e[1]) <= tol) raw.insert(raw.end(), wls[w].begin() + 1, wls[w].end());
            else if (std::hypot(b[0] - e[0], b[1] - e[1]) <= tol) raw.insert(raw.end(), wls[w].rbegin() + 1, wls[w].rend());
            else {
              const auto & f = raw.front();
              if (std::hypot(b[0] - f[0], b[1] - f[1]) <= tol) raw.insert(raw.begin(), wls[w].begin(), wls[w].end() - 1);
              else if (std::hypot(a[0] - f[0], a[1] - f[1]) <= tol) raw.insert(raw.begin(), wls[w].rbegin(), wls[w].rend() - 1);
              else continue;
            }
            done[w] = true;
            njoined++;
            more = true;
          }
        }
      }
      double ox = 0.0, oy = 0.0;
      if (xst != NULL) {
        ox = xst->point->x; oy = xst->point->y; d.placed = true; d.at = "station point";
      } else if (!raw.empty()) {
        double x0 = 1e99, x1 = -1e99, y0 = 1e99, y1 = -1e99;
        for (auto & q : raw) { x0 = std::min(x0, q[0]); x1 = std::max(x1, q[0]); y0 = std::min(y0, q[1]); y1 = std::max(y1, q[1]); }
        // a station on the wall can be just outside the drawing
        double px = sp->point->x, py = sp->point->y, tol = (xm > 0.0) ? XSECTION_EDGE / xm : 0.0;
        if ((px >= x0 - tol) && (px <= x1 + tol) && (py >= y0 - tol) && (py <= y1 + tol)) { ox = px; oy = py; d.placed = true; d.at = "point section"; }
        else d.at = "not known (point section outside the drawing)";
      }
      if ((xm <= 0.0) || raw.empty()) {
        d.placed = false;
        d.at = raw.empty() ? "no wall line" : "section scale not known";
      }
      for (auto & q : raw) d.pts.push_back({(q[0] - ox) * xm, (q[1] - oy) * xm});
      if (nwalls > 1) d.at += fmt::format(" ({} of {} wall lines joined)", njoined, nwalls);
      if (cut == NULL) d.linked += " (no cut before it)";
      out.push_back(d);
      cut = NULL;
    }
  }
}

}  // namespace

void thsplaytube_build(thdatabase * db, thdb3ddata & walls, thdb3ddata & section_splays,
  thdb3ddata & polygons, thdb3ddata & junctions, thdb3ddata & virtual_splays,
  thdb3ddata & map_tubes, thdb3ddata & map_junctions, thdb3ddata & map_fills, thdb3ddata & drawn_sections)
{
  thdb1d & d1 = db->db1d;
  thdb1d_tree_node * nodes = d1.get_tree_nodes();
  size_t nst = d1.station_vec.size();
  if ((nodes == NULL) || (nst == 0)) return;

  const lxVec up(0.0, 0.0, 1.0);
  std::vector<section> sec(nst);
  std::vector<std::vector<size_t>> nbrs(nst);
  size_t nsections = 0, nnosection = 0, njunctions = 0, nused = 0, nfitted = 0, nunordered = 0,
    nvirtual = 0;
  std::string report, outside_report;
  drawn_walls drawn;
  drawn.collect(db);
  drawn_elevation elev;
  elev.collect(db);

  // centreline neighbours of each station (cave legs only)
  for (size_t i = 0; i < nst; i++) {
    thdb1ds & st = d1.station_vec[i];
    if ((st.uid != i + 1) || st.is_temporary() || !st.survey->is_selected()) continue;
    std::set<size_t> nb;
    for (thdb1d_tree_arrow * a = nodes[i].first_arrow; a != NULL; a = a->next_arrow) {
      thdataleg * l = a->leg->leg;
      if ((l->flags & (TT_LEGFLAG_SPLAY | TT_LEGFLAG_SURFACE)) != 0) continue;
      size_t j = a->end_node->uid - 1;
      if ((j != i) && !d1.station_vec[j].is_temporary()) nb.insert(j);
    }
    nbrs[i].assign(nb.begin(), nb.end());
  }

  // sections
  // Junctions with a cross section: for each pair of legs bending less than
  // JUNCTION_STRAIGHT from straight on, the section plane is found as at a
  // station with those two legs; the junction has a cross section across that
  // pair when the plane is fitted (not the bisector), at least JUNCTION_SPLAYS
  // section splays are found, and they go all round the passage (right wall,
  // roof, left wall and floor). Of several such pairs, the one whose bisector
  // the plane is turned least from is taken. The junction then gets that
  // section for the through passage like a station with those two legs; its
  // other legs join it as branches (see side()).
  std::string jsreport;
  size_t njsections = 0;
  std::vector<std::pair<size_t, size_t>> thru(nst, std::make_pair(nst, nst));
  // a leg from junction j (with a cross section) to i that is not part of its
  // through passage: it joins the junction as a branch
  auto side = [&](size_t j, size_t i) {
    return (thru[j].first < nst) && (i != thru[j].first) && (i != thru[j].second);
  };
  auto junction_section = [&](size_t i, size_t & ta, size_t & tb) {
    thdb1ds & st = d1.station_vec[i];
    lxVec c(st.x, st.y, st.z);
    std::vector<lxVec> splays;
    for (thdb1d_tree_arrow * a = nodes[i].first_arrow; a != NULL; a = a->next_arrow) {
      if ((a->leg->leg->flags & TT_LEGFLAG_SPLAY) == 0) continue;
      thdb1ds & e = d1.station_vec[a->end_node->uid - 1];
      lxVec w = lxVec(e.x, e.y, e.z) - c;
      if (w.Length() >= 0.05) splays.push_back(w);
    }
    if (splays.size() < JUNCTION_SPLAYS) return false;
    double bestturn = 1e9;
    size_t ba = 0, bb = 0, bused = 0;
    for (size_t ka = 0; ka < nbrs[i].size(); ka++)
      for (size_t kb = ka + 1; kb < nbrs[i].size(); kb++) {
        thdb1ds & sa = d1.station_vec[nbrs[i][ka]], & sb = d1.station_vec[nbrs[i][kb]];
        lxVec din = c - lxVec(sa.x, sa.y, sa.z), dout = lxVec(sb.x, sb.y, sb.z) - c;
        if ((din.Length() < 1e-6) || (dout.Length() < 1e-6)) continue;
        din.Normalize();
        dout.Normalize();
        if (deg(acos(std::max(-1.0, std::min(1.0, din * dout)))) > JUNCTION_STRAIGHT) continue;
        lxVec nref = din + dout;
        if (nref.Length() < 1e-3) continue;
        nref.Normalize();
        if (fabs(nref.z) < STEEP) {
          nref.z = 0.0;
          nref.Normalize();
        }
        lxVec n = nref;
        bool fitted = false;
        for (int iter = 0; iter < 2; iter++) {
          std::vector<lxVec> cand;
          for (const lxVec & w : splays)
            if (splay_angle(w, n) <= CANDIDATE_ANGLE) cand.push_back(w);
          if (cand.size() < 3) break;
          lxVec nf = fit_plane(cand);
          if (nf * nref < 0.0) nf = -1.0 * nf;
          if (deg(acos(std::min(1.0, nf * nref))) > MAX_PLANE_TURN) { n = nref; fitted = false; break; }
          n = nf;
          fitted = true;
        }
        if (!fitted) continue;
        lxVec h(n.x, n.y, 0.0);
        if (h.Length() < 0.2) h = lxVec(dout.x, dout.y, 0.0);
        if (h.Length() < 0.2) h = lxVec(0.0, 1.0, 0.0);
        lxVec u = h ^ up;
        u = u - (u * n) * n;
        u.Normalize();
        lxVec v = u ^ n;
        size_t nused = 0;
        int quarters = 0;
        for (const lxVec & w : splays) {
          if (splay_angle(w, n) > SECTION_ANGLE) continue;
          double a = deg(atan2(w * v, w * u));
          if (a < 0.0) a += 360.0;
          quarters |= 1 << (int(fmod(a + 45.0, 360.0) / 90.0) % 4);
          nused++;
        }
        double turn = deg(acos(std::min(1.0, fabs(n * nref))));
        if ((nused >= JUNCTION_SPLAYS) && (quarters == 15) && (turn < bestturn)) {
          bestturn = turn;
          ba = nbrs[i][ka];
          bb = nbrs[i][kb];
          bused = nused;
        }
      }
    if (bestturn > 1e8) return false;
    njsections++;
    ta = ba;
    tb = bb;
    jsreport += fmt::format("{:>20}  across {} - {}: {} section splays, plane turned {:.1f} deg\n",
      fmt::format("{}@{}", st.name, st.survey->full_name), d1.station_vec[ba].name, d1.station_vec[bb].name, bused, bestturn);
    return true;
  };
  for (size_t i = 0; i < nst; i++) {
    thdb1ds & st = d1.station_vec[i];
    if ((st.uid != i + 1) || st.is_temporary() || !st.survey->is_selected()) continue;
    std::vector<size_t> legs = nbrs[i];
    if (nbrs[i].size() > 2) {
      njunctions++;
      size_t ta, tb;
      if (!junction_section(i, ta, tb)) continue;
      thru[i] = std::make_pair(ta, tb);
      legs = {ta, tb};
    }
    if (legs.empty()) continue;
    section & s = sec[i];
    s.c = lxVec(st.x, st.y, st.z);
    thdb1ds & sa = d1.station_vec[legs[0]];
    lxVec din = s.c - lxVec(sa.x, sa.y, sa.z);
    din.Normalize();
    lxVec nref = din;
    if (legs.size() == 2) {
      thdb1ds & sb = d1.station_vec[legs[1]];
      lxVec dout = lxVec(sb.x, sb.y, sb.z) - s.c;
      dout.Normalize();
      lxVec bis = din + dout;
      nref = (bis.Length() > 1e-3) ? bis : dout;
    }
    nref.Normalize();
    // reference plane: leg bisector, kept vertical unless the passage is steep
    if (fabs(nref.z) < STEEP) {
      nref.z = 0.0;
      nref.Normalize();
    }

    // splay directions
    std::vector<lxVec> splays, cand;
    for (thdb1d_tree_arrow * a = nodes[i].first_arrow; a != NULL; a = a->next_arrow) {
      if ((a->leg->leg->flags & TT_LEGFLAG_SPLAY) == 0) continue;
      thdb1ds & e = d1.station_vec[a->end_node->uid - 1];
      lxVec w = lxVec(e.x, e.y, e.z) - s.c;
      if (w.Length() >= 0.05) splays.push_back(w);
    }

    // fit the plane of the section splays
    s.n = nref;
    const char * plane = "bisector";
    for (int iter = 0; iter < 2; iter++) {
      cand.clear();
      for (const lxVec & w : splays)
        if (splay_angle(w, s.n) <= CANDIDATE_ANGLE) cand.push_back(w);
      if (cand.size() < 3) break;
      lxVec nf = fit_plane(cand);
      if (nf * nref < 0.0) nf = -1.0 * nf;
      if (deg(acos(std::min(1.0, nf * nref))) > MAX_PLANE_TURN) {
        s.n = nref;
        plane = "bisector (fit turned too far)";
        break;
      }
      s.n = nf;
      plane = "fitted";
    }
    double turn = deg(acos(std::min(1.0, fabs(s.n * nref))));

    // right direction: horizontal, perpendicular to the passage; with up it
    // gives angles anticlockwise from the right looking along the passage
    lxVec h(s.n.x, s.n.y, 0.0);
    if (h.Length() < 0.2) h = lxVec(din.x, din.y, 0.0);
    if (h.Length() < 0.2) h = lxVec(0.0, 1.0, 0.0);
    s.u = h ^ up;
    s.u = s.u - (s.u * s.n) * s.n;
    s.u.Normalize();
    s.v = s.u ^ s.n;

    // section splays
    std::vector<corner> pts;
    double sum2 = 0.0, worst = 0.0;
    for (const lxVec & w : splays) {
      double pa = splay_angle(w, s.n);
      if (pa > SECTION_ANGLE) continue;
      double pu = w * s.u, pv = w * s.v, ang = atan2(pv, pu);
      if (ang < 0.0) ang += 2.0 * M_PI;
      pts.push_back({ang, hypot(pu, pv), w * s.n});
      sum2 += pa * pa;
      worst = std::max(worst, pa);
      thdb3dfc * fc = section_splays.insert_face(THDB3DFC_LINE_STRIP);
      fc->insert_vertex(section_splays.insert_vertex(s.c));
      fc->insert_vertex(section_splays.insert_vertex(s.c + w));
      nused++;
    }
    std::string sname = fmt::format("{}@{}", st.name, st.survey->full_name);
    if (pts.size() < 3) {
      report += fmt::format("{:>20} {:4} {:6} {:4}  {:>5} {:>6} {:>5} {:>6}  no section (too few splays)\n",
        sname, nbrs[i].size(), splays.size(), pts.size(), "-", "-", "-", "-");
      nnosection++;
      continue;
    }
    if (plane[0] == 'f') nfitted++;
    size_t nsplays = pts.size();
    std::sort(pts.begin(), pts.end());

    // the plan drawing: splay ends well outside the drawn walls are reported;
    // where a drawn wall is beyond the furthest point of the section on that
    // side, a virtual horizontal splay from the station to that wall is added
    std::string virt;
    std::vector<size_t> around = drawn.around(i, s.c.x, s.c.y, s.c.z);
    if (!around.empty()) {
      for (const corner & p : pts) {
        lxVec e = s.c + (p.r * cos(p.a)) * s.u + (p.r * sin(p.a)) * s.v + p.h * s.n;
        bool in = false;
        double wd = 1e99;
        for (size_t si : around) {
          if (drawn.inside(si, e.x, e.y)) in = true;
          wd = std::min(wd, drawn.wall_distance(si, e.x, e.y));
        }
        if (!in && (wd > OUTSIDE_REPORT))
          outside_report += fmt::format("{:>20}  splay end {:.2f} m outside the drawn walls\n", sname, wd);
      }
      lxVec hor(s.u.x, s.u.y, 0.0);
      if (hor.Length() > 0.2) {
        hor.Normalize();
        for (int side = 1; side >= -1; side -= 2) {
          double dist = -1.0;
          for (size_t si : around) {
            double t = drawn.ray(si, s.c.x, s.c.y, side * hor.x, side * hor.y);
            if ((t > 0.0) && ((dist < 0.0) || (t < dist))) dist = t;
          }
          if (dist <= 0.0) continue;
          lxVec w = (side * dist) * hor;
          double px = w * s.u, py = w * s.v, ang = atan2(py, px), rad = std::hypot(px, py);
          if (ang < 0.0) ang += 2.0 * M_PI;
          // how far the section reaches across on that side, at any height (as a plan drawing shows)
          double reach = 0.0;
          for (const corner & p : pts)
            reach = std::max(reach, side * p.r * cos(p.a));
          if (rad < reach + VIRTUAL_MIN) continue;
          pts.push_back({ang, rad, w * s.n});
          thdb3dfc * fc = virtual_splays.insert_face(THDB3DFC_LINE_STRIP);
          fc->insert_vertex(virtual_splays.insert_vertex(s.c));
          fc->insert_vertex(virtual_splays.insert_vertex(s.c + w));
          virt += fmt::format(", virtual splay {} {:.2f} m (section {:.2f} m)", (side > 0) ? "right" : "left", rad, reach);
          nvirtual++;
        }
        std::sort(pts.begin(), pts.end());
      }
    }
    // close large gaps at the station (station on a wall, floor or ceiling)
    double maxgap = 0.0, gapstart = 0.0;
    for (size_t k = 0; k < pts.size(); k++) {
      double a0 = pts[k].a, a1 = pts[(k + 1) % pts.size()].a;
      if (a1 <= a0) a1 += 2.0 * M_PI;
      if (a1 - a0 > maxgap) {
        maxgap = a1 - a0;
        gapstart = a0;
      }
    }
    if (maxgap > MAX_GAP / 180.0 * M_PI) {
      double am = fmod(gapstart + maxgap / 2.0, 2.0 * M_PI);
      pts.push_back({am, 0.05, 0.0});
      std::sort(pts.begin(), pts.end());
    }

    // polygon centre (area centroid in the plane)
    double area = 0.0, sx = 0.0, sy = 0.0, mx = 0.0, my = 0.0;
    size_t np = pts.size();
    for (size_t k = 0; k < np; k++) {
      const corner & p0 = pts[k], & p1 = pts[(k + 1) % np];
      double x0 = p0.r * cos(p0.a), y0 = p0.r * sin(p0.a), x1 = p1.r * cos(p1.a), y1 = p1.r * sin(p1.a),
        cr = cross2(x0, y0, x1, y1);
      area += cr;
      sx += (x0 + x1) * cr;
      sy += (y0 + y1) * cr;
      mx += x0 / np;
      my += y0 / np;
    }
    if (fabs(area) > 1e-6) {
      mx = sx / (3.0 * area);
      my = sy / (3.0 * area);
    }
    s.cx = mx;
    s.cy = my;
    for (const corner & p : pts)
      s.q.push_back(lxVec(p.r * cos(p.a), p.r * sin(p.a), p.h));
    std::string anchors = "";
    if (!s.set_positions()) {
      anchors = ", anchors out of order";
      nunordered++;
    }
    report += fmt::format("{:>20} {:4} {:6} {:4}  {:5.1f} {:6.1f} {:5.1f} {:6.2f}  {}\n",
      sname, nbrs[i].size(), splays.size(), nsplays, turn, sqrt(sum2 / nsplays), worst,
      hypot(s.cx, s.cy), std::string(plane) + anchors + s.moved + virt);

    thdb3dfc * pfc = polygons.insert_face(THDB3DFC_LINE_STRIP);
    np = s.q.size();  // anchors may have been listed twice
    for (size_t k = 0; k <= np; k++)
      pfc->insert_vertex(polygons.insert_vertex(s.point(k % np)));
    for (size_t k = 0; k < np; k++)
      s.vx.push_back(walls.insert_vertex(s.point(k)));
    s.ok = true;
    nsections++;
  }


  // Drawn cross-sections (see collect_xsections), shown beside the section
  // polygons from splays for comparison: upright through the station, across
  // the passage along the cut. Which way the drawing faces is not recorded
  // (the cuts are `-direction both`): of the drawing and its mirror image, the
  // one overlapping the station's splay polygon more is used; at stations
  // without one, the way most others face.
  std::vector<drawn_xsection> xsecs;
  std::string xsreport;
  size_t nxsplaced = 0, nxsused = 0;
  std::vector<bool> drawnsec(nst, false);  // the station's section is a drawn one
  collect_xsections(db, xsecs);
  {
    // overlap (intersection over union) of two polygons, on a grid
    auto overlap = [](const std::vector<std::array<double, 2>> & a, const std::vector<std::array<double, 2>> & b) {
      if ((a.size() < 3) || (b.size() < 3)) return -1.0;
      ring2 ra(a.begin(), a.end()), rb(b.begin(), b.end());
      double x0 = 1e99, x1 = -1e99, y0 = 1e99, y1 = -1e99;
      for (const ring2 * r : {&ra, &rb})
        for (auto & q : *r) { x0 = std::min(x0, q[0]); x1 = std::max(x1, q[0]); y0 = std::min(y0, q[1]); y1 = std::max(y1, q[1]); }
      int ni = 0, nu = 0;
      const int G = 80;
      for (int i = 0; i < G; i++)
        for (int j = 0; j < G; j++) {
          double x = x0 + (i + 0.5) * (x1 - x0) / G, y = y0 + (j + 0.5) * (y1 - y0) / G;
          bool ia = in_ring(ra, x, y), ib = in_ring(rb, x, y);
          if (ia && ib) ni++;
          if (ia || ib) nu++;
        }
      return (nu > 0) ? double(ni) / nu : 0.0;
    };
    // forward along the survey at a station: the leg arriving there (else the one leaving)
    auto forward = [&](size_t i) {
      lxVec f(0.0, 0.0, 0.0);
      thdb1ds & si = d1.station_vec[i];
      for (thdb1d_tree_arrow * ar = nodes[i].first_arrow; ar != NULL; ar = ar->next_arrow) {
        if ((ar->leg->leg->flags & (TT_LEGFLAG_SPLAY | TT_LEGFLAG_SURFACE)) != 0) continue;
        thdb1ds & sj = d1.station_vec[ar->end_node->uid - 1];
        lxVec d(si.x - sj.x, si.y - sj.y, 0.0);
        if (ar->is_reversed) return d;  // leg from j to i
        if (f.Length() < 1e-9) f = -1.0 * d;
      }
      return f;
    };
    std::vector<int> mirror(xsecs.size(), 0);
    std::vector<double> ov(xsecs.size(), -1.0), ovm(xsecs.size(), -1.0);
    int nfwd = 0, nback = 0;
    for (size_t k = 0; k < xsecs.size(); k++) {
      drawn_xsection & d = xsecs[k];
      if ((d.station < 0) || !d.placed || (d.cut.Length() < 1e-9) || !sec[d.station].ok) continue;
      section & s = sec[d.station];
      std::vector<std::array<double, 2>> sp, dm;
      for (size_t q = 0; q < s.q.size(); q++) {
        lxVec P = s.point(q) - s.c;
        sp.push_back({P * d.cut, P.z});
      }
      for (auto & q : d.pts) dm.push_back({-q[0], q[1]});
      ov[k] = overlap(d.pts, sp);
      ovm[k] = overlap(dm, sp);
      // only a clear difference decides it; otherwise the way most others face
      if (std::fabs(ov[k] - ovm[k]) < XSECTION_CLEAR) continue;
      mirror[k] = (ovm[k] > ov[k]) ? -1 : 1;
      lxVec f = forward(d.station), right(f.y, -f.x, 0.0);
      if ((mirror[k] * (d.cut * right)) > 0.0) nfwd++; else nback++;
    }
    for (size_t k = 0; k < xsecs.size(); k++) {
      drawn_xsection & d = xsecs[k];
      std::string sname = (d.station >= 0) ? fmt::format("{}@{}", d1.station_vec[d.station].name, d1.station_vec[d.station].survey->full_name) : "?";
      if ((d.station < 0) || !d.placed || (d.cut.Length() < 1e-9)) {
        xsreport += fmt::format("{:>20} {:>20}  station {}: not placed ({}; {})\n", d.scrap, d.plan, sname, d.linked, d.at);
        continue;
      }
      thdb1ds & st = d1.station_vec[d.station];
      lxVec c(st.x, st.y, st.z), f = forward(d.station), right(f.y, -f.x, 0.0);
      std::string facing = "by overlap with the splay polygon";
      if (mirror[k] == 0) {  // no splay polygon, or no clear difference: the way most others face
        facing = "as most others";
        int want = (nfwd >= nback) ? 1 : -1;
        mirror[k] = ((d.cut * right) * want > 0.0) ? 1 : -1;
      }
      lxVec across = double(mirror[k]) * d.cut;
      thdb3dfc * fc = drawn_sections.insert_face(THDB3DFC_LINE_STRIP);
      thdb3dvx * first = NULL;
      double a0 = 1e99, a1 = -1e99, z0 = 1e99, z1 = -1e99;
      for (auto & q : d.pts) {
        thdb3dvx * v = drawn_sections.insert_vertex(c + q[0] * across + lxVec(0.0, 0.0, q[1]));
        fc->insert_vertex(v);
        if (first == NULL) first = v;
        a0 = std::min(a0, q[0]); a1 = std::max(a1, q[0]); z0 = std::min(z0, q[1]); z1 = std::max(z1, q[1]);
      }
      fc->insert_vertex(first);
      nxsplaced++;
      std::string cmp = "no splay polygon";
      if (sec[d.station].ok) {
        section & s = sec[d.station];
        double b0 = 1e99, b1 = -1e99, w0 = 1e99, w1 = -1e99;
        for (size_t q = 0; q < s.q.size(); q++) {
          lxVec P = s.point(q) - s.c;
          double a = P * across;
          b0 = std::min(b0, a); b1 = std::max(b1, a); w0 = std::min(w0, P.z); w1 = std::max(w1, P.z);
        }
        cmp = fmt::format("drawn {:.2f} x {:.2f} m, splays {:.2f} x {:.2f} m; overlap {:.2f} (mirrored {:.2f})",
          a1 - a0, z1 - z0, b1 - b0, w1 - w0, std::max(ov[k], ovm[k]), std::min(ov[k], ovm[k]));
      }
      // The drawn section becomes the station's section (at passage stations and
      // junctions with a through passage): its plane, its outline as corners;
      // anchors, joining and tubes as for sections from splays. The section
      // polygon from splays is still shown.
      std::string used;
      if ((nbrs[d.station].size() > 2) && (thru[d.station].first >= nst)) used = "; at a junction without a through passage: shown only";
      else if (drawnsec[d.station]) used = "; the station already has a drawn section: shown only";
      else {
        std::vector<lxVec> q;
        for (auto & p : d.pts)
          if (q.empty() || (std::hypot(p[0] - q.back().x, p[1] - q.back().y) > 1e-3)) q.push_back(lxVec(p[0], p[1], 0.0));
        if ((q.size() > 1) && (std::hypot(q[0].x - q.back().x, q[0].y - q.back().y) <= 1e-3)) q.pop_back();
        double ar = 0.0, gx = 0.0, gy = 0.0;
        for (size_t m = 0; m < q.size(); m++) {
          const lxVec & q0 = q[m], & q1 = q[(m + 1) % q.size()];
          double cr = cross2(q0.x, q0.y, q1.x, q1.y);
          ar += cr; gx += (q0.x + q1.x) * cr; gy += (q0.y + q1.y) * cr;
        }
        if ((q.size() < 3) || (std::fabs(ar) < 1e-4)) used = "; no area: shown only";
        else {
          if (ar < 0.0) std::reverse(q.begin(), q.end());  // anticlockwise looking along the passage
          section & s = sec[d.station];
          bool was = s.ok;
          s.c = c;
          s.u = across;
          s.v = up;
          s.n = up ^ across;
          s.n.Normalize();
          s.cx = gx / (3.0 * ar);
          s.cy = gy / (3.0 * ar);
          s.q = q;
          s.pos.clear();
          s.moved.clear();
          if (!s.set_positions()) used = ", anchors out of order";
          s.vx.clear();
          for (size_t m = 0; m < s.q.size(); m++) s.vx.push_back(walls.insert_vertex(s.point(m)));
          s.ok = true;
          if (!was) nsections++;
          drawnsec[d.station] = true;
          nxsused++;
          used = "; used as the station's section" + used + s.moved;
        }
      }
      xsreport += fmt::format("{:>20} {:>20}  station {}: {}; station in the drawing at the {}; looking {} the survey ({}); {}{}\n",
        d.scrap, d.plan, sname, d.linked, d.at, ((across * right) > 0.0) ? "along" : "against", facing, cmp, used);
    }
    if (!xsecs.empty())
      xsreport = fmt::format("faces looking along the survey {}, against {} (by overlap with the splay polygon)\n", nfwd, nback) + xsreport;
  }

  // bodies at junctions and at stations without a section (not passage ends)
  size_t nbodies = 0;
  std::vector<bool> body(nst, false);
  for (size_t i = 0; i < nst; i++) {
    thdb1ds & st = d1.station_vec[i];
    if ((st.uid != i + 1) || st.is_temporary() || !st.survey->is_selected()) continue;
    if (sec[i].ok || (nbrs[i].size() < 2)) continue;
    lxVec c(st.x, st.y, st.z);
    std::vector<lxVec> splays;
    for (thdb1d_tree_arrow * a = nodes[i].first_arrow; a != NULL; a = a->next_arrow) {
      if ((a->leg->leg->flags & TT_LEGFLAG_SPLAY) == 0) continue;
      thdb1ds & e = d1.station_vec[a->end_node->uid - 1];
      lxVec w = lxVec(e.x, e.y, e.z) - c;
      if (w.Length() >= 0.05) splays.push_back(w);
    }
    body[i] = true;
    if (splay_body(c, splays, junctions)) nbodies++;
  }

  // Model from maps and sections: sections along each leg (see positions)
  // (each the tube between the two station sections at that point), stretched
  // sideways only, right and left separately, to the drawn plan walls on the
  // horizontal line across the passage there.
  bool usemap = !drawn.scraps.empty();
  std::string stretchreport;
  size_t nrings = 0, nstretched = 0, nwhole = 0;
  auto horizontal = [](const lxVec & u) {
    lxVec h(u.x, u.y, 0.0);
    if (h.Length() < 1e-6) h = lxVec(1.0, 0.0, 0.0);
    h.Normalize();
    return h;
  };
  auto stretch = [&](std::vector<lxVec> & P, const lxVec & c, const lxVec & h, const std::string & where) {
    double dr, dl;
    if (!drawn.width(c.x, c.y, c.z, h, dr, dl)) return;
    double R = 0.0, L = 0.0;
    lxVec pr = c, pl = c;
    for (const lxVec & p : P) {
      double x = (p - c) * h;
      if (x > R) { R = x; pr = p; }
      if (-x > L) { L = -x; pl = p; }
    }
    // the drawn wall on the line through each side's furthest corner, so that
    // the corner itself meets the wall (corners lie a little off the section
    // line where they are true splay ends)
    double r1, l1;
    lxVec q = pr - R * h;
    if ((R > 0.1) && drawn.width(q.x, q.y, c.z, h, r1, l1)) dr = r1;
    q = pl + L * h;
    if ((L > 0.1) && drawn.width(q.x, q.y, c.z, h, r1, l1)) dl = l1;
    // a side with (almost) nothing on it, as where the station is against a
    // wall, cannot be stretched by itself: the whole section is then scaled
    // and moved across to fill the drawn width
    if ((R + L > 0.1) && ((std::min(R, L) < 0.1) || (std::min(R, L) < WHOLE_STRETCH * std::max(R, L)))) {
      double a = (dr + dl) / (R + L), b = dr - a * R;
      for (lxVec & p : P) {
        double x = (p - c) * h;
        p += (x * (a - 1.0) + b) * h;
      }
      nstretched++;
      nwhole++;
      if ((a > STRETCH_REPORT) || (a < 1.0 / STRETCH_REPORT))
        stretchreport += fmt::format("{:>30}  whole width x{:.2f} ({:.2f} -> {:.2f} m), moved {:+.2f} m\n",
          where, a, R + L, dr + dl, b);
      return;
    }
    double fr = (R > 0.1) ? dr / R : 1.0, fl = (L > 0.1) ? dl / L : 1.0;
    for (lxVec & p : P) {
      double x = (p - c) * h;
      p += (x * (((x > 0.0) ? fr : fl) - 1.0)) * h;
    }
    nstretched++;
    if ((std::max(fr, fl) > STRETCH_REPORT) || (std::min(fr, fl) < 1.0 / STRETCH_REPORT))
      stretchreport += fmt::format("{:>30}  right x{:.2f} ({:.2f} -> {:.2f} m)  left x{:.2f} ({:.2f} -> {:.2f} m)\n",
        where, fr, R, dr, fl, L, dl);
  };
  // Roof and floor from the extended elevation at fraction f along the leg from
  // station i to j: the leg's position along the extended elevation, in the
  // scraps drawing both its stations (else either); steep or hidden legs are skipped.
  size_t nvscaled = 0;
  std::string vreport;
  auto steep_leg = [&](size_t a, size_t b) {
    thdb1ds & sa = d1.station_vec[a], & sb = d1.station_vec[b];
    return deg(atan2(std::fabs(sb.z - sa.z), std::hypot(sb.x - sa.x, sb.y - sa.y))) > ELEVATION_STEEP;
  };
  // the leg's position along the extended elevation at f, its height, the
  // scraps to use (both stations' scraps if any, else either's), and whether
  // it is steeper than ELEVATION_STEEP; hidden legs are skipped
  auto leg_ext = [&](size_t i, size_t j, double f, double & x, double & z,
                     std::vector<size_t> & both, std::vector<size_t> & either, bool & steep) {
    if (elev.rings.empty()) return false;
    for (thdb1d_tree_arrow * ar = nodes[i].first_arrow; ar != NULL; ar = ar->next_arrow) {
      thdataleg * l = ar->leg->leg;
      if (((l->flags & (TT_LEGFLAG_SPLAY | TT_LEGFLAG_SURFACE)) != 0) || (ar->end_node->uid - 1 != j)) continue;
      if ((l->extend & TT_EXTENDFLAG_HIDE) != 0) return false;
      thdb1ds & si = d1.station_vec[i];
      steep = steep_leg(i, j);
      double xi = ar->is_reversed ? l->txx : l->fxx, xj = ar->is_reversed ? l->fxx : l->txx;
      both.clear();
      either.clear();
      auto ii = elev.station_scraps.find(i), jj = elev.station_scraps.find(j);
      if (ii != elev.station_scraps.end())
        for (size_t s1 : ii->second) {
          either.push_back(s1);
          if ((jj != elev.station_scraps.end()) && (std::find(jj->second.begin(), jj->second.end(), s1) != jj->second.end()))
            both.push_back(s1);
        }
      if (jj != elev.station_scraps.end())
        for (size_t s2 : jj->second)
          if (std::find(either.begin(), either.end(), s2) == either.end()) either.push_back(s2);
      x = xi + f * (xj - xi);
      z = si.z + f * (d1.station_vec[j].z - si.z);
      return true;
    }
    return false;
  };
  // Which way the extended elevation is read across a leg: an angle from
  // vertical, leaning towards the leg's forward direction, taken at each end
  // station and interpolated along the leg. At a station with a section: the
  // section's own up direction (so the drawing is read in the section's plane;
  // both legs steeper than ELEVATION_STEEP give a level section, read across).
  // At a station without one, between two legs: vertical if neither leg is
  // steep, horizontal if both are, else square to the legs' bisector seen from
  // the side. At dead ends and junctions without a section, from the leg alone.
  // (To be reviewed with more vertical caves.)
  auto other_leg = [&](size_t s, size_t not_this) {
    return (nbrs[s].size() == 2) ? ((nbrs[s][0] == not_this) ? nbrs[s][1] : nbrs[s][0]) : nst;
  };
  // legs giving the leg i -> j's forward direction in plan: itself, else (about
  // vertical) the leg before it, else the one after it
  auto forward_legs = [&](size_t i, size_t j) {
    std::vector<std::pair<size_t, size_t>> r = {{i, j}};
    size_t k = other_leg(i, j), l = other_leg(j, i);
    if (k != nst) r.push_back({k, i});
    if (l != nst) r.push_back({j, l});
    return r;
  };
  auto plan_dir = [&](size_t a, size_t b, lxVec & h) {
    thdb1ds & sa = d1.station_vec[a], & sb = d1.station_vec[b];
    double hl = std::hypot(sb.x - sa.x, sb.y - sa.y);
    if (hl < ELEVATION_LEVEL) return false;
    h = lxVec((sb.x - sa.x) / hl, (sb.y - sa.y) / hl, 0.0);
    return true;
  };
  auto station_tilt = [&](size_t p, size_t s, size_t n, const lxVec & h) {  // legs p -> s -> n; p or n may be nst
    if (sec[s].ok) {
      lxVec v = sec[s].v;
      if (v.z < 0.0) v = -1.0 * v;
      return deg(atan2(v * h, v.z));
    }
    if (p == nst) return steep_leg(s, n) ? 90.0 : 0.0;
    if (n == nst) return steep_leg(p, s) ? 90.0 : 0.0;
    bool sp = steep_leg(p, s), sn = steep_leg(s, n);
    if (sp == sn) return sp ? 90.0 : 0.0;
    thdb1ds & a = d1.station_vec[p], & b = d1.station_vec[s], & c = d1.station_vec[n];
    double ia = std::hypot(b.x - a.x, b.y - a.y), iz = b.z - a.z, oa = std::hypot(c.x - b.x, c.y - b.y), oz = c.z - b.z,
      il = std::hypot(ia, iz), ol = std::hypot(oa, oz);
    if ((il < 1e-6) || (ol < 1e-6)) return 0.0;
    return deg(atan2(-(iz / il + oz / ol), ia / il + oa / ol));  // square to the bisector
  };
  auto leg_tilt = [&](size_t i, size_t j, double f) {
    lxVec h(0.0, 0.0, 0.0);
    for (auto & ab : forward_legs(i, j))
      if (plan_dir(ab.first, ab.second, h)) break;
    double ti = station_tilt(other_leg(i, j), i, j, h), tj = station_tilt(i, j, other_leg(j, i), h);
    if (tj - ti > 90.0) tj -= 180.0;
    if (ti - tj > 90.0) tj += 180.0;
    return ti + f * (tj - ti);
  };
  // The leg i -> j as drawn at f (elev.drawn_leg), the direction to read the
  // drawing across it ((dx, dy) in the scrap; vertical: floor and roof only),
  // and that direction in the cave (D), tilted towards the forward direction in
  // plan (forward_legs), which must run along the drawing
  auto drawn_across = [&](size_t i, size_t j, double f, bool vertical, size_t & rsi, double & rx, double & ry,
                          double & dx, double & dy, lxVec & D) {
    thdb1ds & si = d1.station_vec[i], & sj = d1.station_vec[j];
    if (!elev.drawn_leg(i, j, f, std::hypot(sj.x - si.x, sj.y - si.y), sj.z - si.z, rsi, rx, ry)) return false;
    double t = vertical ? 0.0 : rad(leg_tilt(i, j, f));
    dx = 0.0;
    dy = 1.0;
    D = lxVec(0.0, 0.0, 1.0);
    if (std::fabs(sin(t)) < 1e-6) return true;
    lxVec h;
    double sgn = 0.0;
    for (auto & ab : forward_legs(i, j)) {
      thdb1ds & sa = d1.station_vec[ab.first], & sb = d1.station_vec[ab.second];
      double al = 0.0, x0, y0;
      size_t s0;
      if (plan_dir(ab.first, ab.second, h)
          && elev.drawn_leg(ab.first, ab.second, 0.0, std::hypot(sb.x - sa.x, sb.y - sa.y), sb.z - sa.z, s0, x0, y0, NULL, &al)
          && (std::fabs(al) >= ELEVATION_LEVEL)) {
        sgn = (al > 0.0) ? 1.0 : -1.0;
        break;
      }
    }
    if (sgn == 0.0) return false;
    dx = sgn * sin(t);
    dy = cos(t);
    D = sin(t) * h + lxVec(0.0, 0.0, cos(t));
    return true;
  };
  // all drawn spans at f along the leg, as distances along D from the leg
  // (up, down); vertical: floor and roof only, steep legs skipped
  auto leg_spans = [&](size_t i, size_t j, double f, std::vector<std::pair<double, double>> & ud, lxVec & D,
                       bool vertical) {
    double x, z;
    bool steep;
    std::vector<size_t> both, either;
    ud.clear();
    D = lxVec(0.0, 0.0, 1.0);
    if (!leg_ext(i, j, f, x, z, both, either, steep) || (vertical && steep)) return false;
    size_t rsi;
    double rx, ry, dx, dy;
    if (drawn_across(i, j, f, vertical, rsi, rx, ry, dx, dy, D)) {  // as drawn
      for (auto & p : elev.raw_spans(rsi, rx, ry, dx, dy)) ud.push_back(std::make_pair(p.second, -p.first));
      return !ud.empty();
    }
    if (steep) return false;
    std::vector<std::pair<double, double>> sp = elev.spans(both, x);
    if (sp.empty()) sp = elev.spans(either, x);
    for (auto & p : sp) ud.push_back(std::make_pair(p.second - z, z - p.first));
    return !ud.empty();
  };
  // Up and down along a leg from the extended elevation, in the given order of
  // positions: the first takes the span containing the leg (or the nearest within
  // ELEVATION_NEAR); each next one the span continuing from the previous (most
  // overlap); where the height would change by more than BRANCH_OPENING times
  // from one position to the next, the previous up and down are kept.
  auto track_elevation = [&](size_t i, size_t j, const std::vector<double> & fs, const std::vector<size_t> & order,
                             std::vector<double> & up, std::vector<double> & down, std::vector<char> & ok,
                             std::vector<lxVec> & dirs, bool vertical) {
    up.assign(fs.size(), 0.0);
    down.assign(fs.size(), 0.0);
    ok.assign(fs.size(), 0);
    dirs.assign(fs.size(), lxVec(0.0, 0.0, 1.0));
    bool have = false;
    double pu = 0.0, pd = 0.0;
    for (size_t k : order) {
      std::vector<std::pair<double, double>> ud;
      if (!leg_spans(i, j, fs[k], ud, dirs[k], vertical)) continue;
      double bu = 0.0, bd = 0.0, best = -1e99;
      for (auto & c : ud) {
        double score;
        if (have) score = std::min(c.first, pu) + std::min(c.second, pd);  // overlap with the previous span
        else score = -((c.first < 0.0) ? -c.first : ((c.second < 0.0) ? -c.second : 0.0));  // contains the leg
        if (score > best) { best = score; bu = c.first; bd = c.second; }
      }
      if (!have && (best < -ELEVATION_NEAR)) continue;
      if (have && (best <= 0.0)) { bu = pu; bd = pd; }  // nothing continues: keep the previous
      if (have && ((bu + bd > BRANCH_OPENING * (pu + pd)) || (BRANCH_OPENING * (bu + bd) < pu + pd))) { bu = pu; bd = pd; }
      up[k] = bu;
      down[k] = bd;
      ok[k] = 1;
      pu = bu;
      pd = bd;
      have = true;
    }
  };
  std::set<std::pair<size_t, size_t>> elevlegs_drawn, elevlegs_placed;
  // the drawn span at f along the leg containing the leg (or the nearest
  // within ELEVATION_NEAR), as distances along D from the leg (lo, hi);
  // vertical: floor and roof only, steep legs skipped
  auto leg_elevation = [&](size_t i, size_t j, double f, double & lo, double & hi, lxVec & D, bool vertical = false) {
    double x, z;
    bool steep;
    std::vector<size_t> both, either;
    D = lxVec(0.0, 0.0, 1.0);
    if (!leg_ext(i, j, f, x, z, both, either, steep) || (vertical && steep)) return false;
    // as drawn: the leg's place in a scrap with both stations drawn in it
    size_t rsi;
    double rx, ry, dx, dy;
    if (drawn_across(i, j, f, vertical, rsi, rx, ry, dx, dy, D)) {
      double best = -1.0;
      for (auto & p : elev.raw_spans(rsi, rx, ry, dx, dy)) {
        double d = (p.first > 0.0) ? p.first : ((p.second < 0.0) ? -p.second : 0.0);
        if ((d <= ELEVATION_NEAR) && ((best < 0.0) || (d < best))) {
          best = d;
          lo = p.first;
          hi = p.second;
        }
      }
      elevlegs_drawn.insert(std::make_pair(std::min(i, j), std::max(i, j)));
      return best >= 0.0;
    }
    if (steep) return false;
    elevlegs_placed.insert(std::make_pair(std::min(i, j), std::max(i, j)));
    double floor, roof;
    if ((both.empty() || !elev.span(both, x, z, floor, roof)) && !elev.span(either, x, z, floor, roof)) return false;
    lo = floor - z;
    hi = roof - z;
    return true;
  };
  auto station_elevation = [&](size_t i, double & lo, double & hi, lxVec & D) {
    for (size_t k : nbrs[i])
      if (leg_elevation(i, k, 0.0, lo, hi, D)) return true;
    return false;
  };
  auto station_at = [&](size_t i) { return lxVec(d1.station_vec[i].x, d1.station_vec[i].y, d1.station_vec[i].z); };
  // scales a ring along D about c so that its extremes meet the drawn span
  // (lo, hi from c along D): roof and floor when D is vertical
  auto vstretch = [&](std::vector<lxVec> & P, const lxVec & c, const lxVec & D, double lo, double hi,
                      const std::string & where) {
    double top = -1e99, bot = 1e99;
    for (const lxVec & p : P) { top = std::max(top, (p - c) * D); bot = std::min(bot, (p - c) * D); }
    if ((top - bot < 0.1) || (hi - lo < 0.1)) return;
    double fac = (hi - lo) / (top - bot);
    for (lxVec & p : P) {
      double t = (p - c) * D;
      p = p + (lo + (t - bot) * fac - t) * D;
    }
    nvscaled++;
    if ((fac > STRETCH_REPORT) || (fac < 1.0 / STRETCH_REPORT))
      vreport += fmt::format("{:>30}  height x{:.2f} ({:.2f} -> {:.2f} m)\n", where, fac, top - bot, hi - lo);
  };
  // Fractions along a leg (up to fmax) for the sections: wherever the drawn
  // width or height has changed by more than spacing_change percent since the
  // last section, at least spacing_min and at most spacing_max apart. measure(f)
  // gives (width, height) at f, negative where not drawn.
  auto positions = [&](double len, double fmax, const std::function<lxVec(double)> & measure) {
    std::vector<double> fs = {0.0};
    double mn = std::max(0.05, thsplaytube_opt.spacing_min), mx = std::max(mn, thsplaytube_opt.spacing_max),
      chg = thsplaytube_opt.spacing_change / 100.0;
    if (len * fmax <= 1e-6) return fs;
    int ns = std::max(1, (int) ceil(fmax * len / mn));
    double step = fmax * len / ns, last = 0.0;
    lxVec ml = measure(0.0);
    auto changed = [chg](double a, double b) { return (a > 0.0) && (b > 0.0) && (std::fabs(b - a) > chg * a); };
    for (int k = 1; k <= ns; k++) {
      double f = fmax * k / ns, dist = (f - last) * len;
      lxVec mv = measure(f);
      if ((k == ns) || changed(ml.x, mv.x) || changed(ml.y, mv.y) || (dist + step > mx + 1e-9)) {
        fs.push_back(f);
        last = f;
        ml = mv;
      }
    }
    return fs;
  };
  size_t nlegspaced = 0;
  auto ring = [&](const std::vector<lxVec> & P) {
    std::vector<thdb3dvx *> r;
    for (const lxVec & p : P) r.push_back(map_tubes.insert_vertex(p));
    nrings++;
    return r;
  };
  auto loft = [&](const std::vector<thdb3dvx *> & ra, const std::vector<thdb3dvx *> & rb) {
    loft_rings(map_tubes.insert_face(THDB3DFC_TRIANGLES), ra, rb);
  };
  auto station_ring = [&](size_t i) {  // a station section stretched to the drawn walls
    std::vector<lxVec> P;
    for (size_t k = 0; k < sec[i].q.size(); k++) P.push_back(sec[i].point(k));
    if (drawnsec[i]) return P;  // a drawn section wins at its own station
    std::string where = fmt::format("{}@{}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name);
    stretch(P, sec[i].c, horizontal(sec[i].u), where);
    double ef, er;
    lxVec ed;
    if (station_elevation(i, ef, er, ed)) vstretch(P, station_at(i), ed, ef, er, where);
    return P;
  };
  // branches ending at junctions: end position, across direction, end ring
  // (points and vertices), whether it is joined to the junction, its section
  struct branch_end {
    lxVec c, h;
    std::vector<lxVec> P;
    std::vector<thdb3dvx *> vx;
    bool joined;
    size_t sec;  // the branch's station (a section, or the far end of an oval leg)
    bool oval;  // an oval leg: closed by a fan, not as a section
    bool steep;  // a steep leg (a pitch): not used for the junction footprint
  };
  std::map<size_t, std::vector<branch_end>> branch_ends;

  // tubes along legs between stations with sections
  size_t ntubes = 0, ntapers = 0, nextended = 0, ncaps = 0, nfans = 0, ndeadends = 0;
  std::string flipreport, tubereport;
  for (size_t i = 0; i < nst; i++) {
    if (!sec[i].ok) continue;
    for (size_t j : nbrs[i]) {
      if (side(i, j)) continue;  // built from the other end, as a branch into this junction
      // taper to a passage end without section
      if (!sec[j].ok && (nbrs[j].size() == 1)) {
        thdb1ds & e = d1.station_vec[j];
        thdb3dvx * ev = walls.insert_vertex(e.x, e.y, e.z);
        taper_ring(walls.insert_face(THDB3DFC_TRIANGLES), sec[i].vx, ev);
        if (usemap && drawn.around(j, e.x, e.y, e.z).empty()) {
          std::vector<thdb3dvx *> r = ring(station_ring(i));
          taper_ring(map_tubes.insert_face(THDB3DFC_TRIANGLES), r, map_tubes.insert_vertex(e.x, e.y, e.z));
        } else if (usemap) {
          // a drawn dead end: the section carried along the leg, stretched to the
          // drawn walls and the extended elevation, rounded off to the drawn end wall
          section & a = sec[i];
          lxVec h = horizontal(a.u), cj(e.x, e.y, e.z), out(cj.x - a.c.x, cj.y - a.c.y, 0.0);
          double len = out.Length();
          if (len > 1e-6) out = out / len;
          std::vector<double> fs = positions(len, 1.0, [&](double f) {
            lxVec c = a.c + f * (cj - a.c);
            double r1, l1, ef, er;
            lxVec ed;
            return lxVec(drawn.width(c.x, c.y, c.z, h, r1, l1) ? r1 + l1 : -1.0,
              leg_elevation(i, j, f, ef, er, ed) ? er - ef : -1.0, 0.0);
          });
          nlegspaced++;
          std::vector<size_t> dord;
          for (size_t m = 0; m < fs.size(); m++) dord.push_back(m);
          std::vector<double> dup, ddown;
          std::vector<char> dok;
          std::vector<lxVec> ddir;
          track_elevation(i, j, fs, dord, dup, ddown, dok, ddir, false);
          std::vector<lxVec> Plast = station_ring(i);
          std::vector<thdb3dvx *> prev = ring(Plast);
          for (size_t m = 1; m < fs.size(); m++) {
            double f = fs[m];
            lxVec c = a.c + f * (cj - a.c);
            std::vector<lxVec> P;
            for (size_t k = 0; k < a.q.size(); k++) P.push_back(a.point(k) + (c - a.c));
            std::string where = fmt::format("{}@{} -> {}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name, e.name);
            stretch(P, c, h, where);
            if (dok[m]) vstretch(P, station_at(i) + f * (lxVec(e.x, e.y, e.z) - station_at(i)), ddir[m], -ddown[m], dup[m], where);
            std::vector<thdb3dvx *> cur = ring(P);
            loft(prev, cur);
            prev = cur;
            Plast = P;
          }
          double wall = -1.0;
          for (size_t si2 = 0; si2 < drawn.scraps.size(); si2++) {
            if (!drawn.at_level(si2, cj.z) || !drawn.inside(si2, cj.x, cj.y)) continue;
            double t = drawn.ray(si2, cj.x, cj.y, out.x, out.y);
            if ((t > 0.0) && ((wall < 0.0) || (t < wall))) wall = t;
          }
          lxVec cc = centroid(prev);
          std::vector<thdb3dvx *> rr = prev;
          if (wall > 0.05) {  // half an ellipsoid out to the end wall
            for (int k = 1; k < 4; k++) {
              double t = k / 4.0, sc = std::sqrt(1.0 - t * t);
              std::vector<lxVec> Q;
              for (const lxVec & p : Plast) Q.push_back(cc + sc * (p - cc) + (t * wall) * out);
              std::vector<thdb3dvx *> nr = ring(Q);
              loft(rr, nr);
              rr = nr;
            }
            taper_ring(map_tubes.insert_face(THDB3DFC_TRIANGLES), rr, map_tubes.insert_vertex(cc + wall * out));
          } else {
            piece pc;
            thdb3dvx * cv = map_tubes.insert_vertex(cc);
            for (size_t k = 0; k < rr.size(); k++) pc.add(cv, rr[(k + 1) % rr.size()], rr[k], out);
            pc.emit(map_tubes.insert_face(THDB3DFC_TRIANGLES));
          }
          ndeadends++;
        }
        sec[i].ntubes++;
        ntapers++;
        continue;
      }
      // extend into a junction (or a station without section): the section
      // is carried unchanged to that station, where the junction body covers the join
      if (body[j] || side(j, i)) {
        section & a = sec[i];
        thdb1ds & e = d1.station_vec[j];
        lxVec shift = lxVec(e.x, e.y, e.z) - a.c;
        size_t nv = a.vx.size();
        std::vector<thdb3dvx *> rb;
        for (size_t k = 0; k < nv; k++)
          rb.push_back(walls.insert_vertex(a.point(k) + shift));
        loft_rings(walls.insert_face(THDB3DFC_TRIANGLES), a.vx, rb);
        // close the end inside the junction
        if (cap(a, rb, walls.insert_vertex(a.c + shift), walls, shift)) ncaps++; else nfans++;
        if (usemap) {
          // branch of the model from maps and sections: the section carried along
          // the leg, stretched to the drawn walls, up to where the passage opens
          // into the junction (or BRANCH_END before the junction station)
          lxVec h = horizontal(a.u), cj(e.x, e.y, e.z);
          double len = std::hypot(cj.x - a.c.x, cj.y - a.c.y);
          // on a steep leg (a pitch) the branch is measured along the leg, so that
          // it ends just above or below the junction, not at the far station
          bool steep = deg(atan2(std::fabs(cj.z - a.c.z), len)) > ELEVATION_STEEP;
          if (steep) len = (cj - a.c).Length();
          double fmax = (len > BRANCH_END) ? (len - BRANCH_END) / len : 0.0;
          double w0r, w0l, wr, wl;
          bool w0 = drawn.width(a.c.x, a.c.y, a.c.z, h, w0r, w0l);
          std::vector<lxVec> P0 = station_ring(i);
          std::vector<thdb3dvx *> prev = ring(P0);
          double fend = 0.0;
          std::vector<lxVec> Pend = P0;
          std::vector<double> fs = positions(len, fmax, [&](double f) {
            lxVec c = a.c + f * (cj - a.c);
            double r1, l1, ef, er;
            lxVec ed;
            return lxVec(drawn.width(c.x, c.y, c.z, h, r1, l1) ? r1 + l1 : -1.0,
              leg_elevation(i, j, f, ef, er, ed) ? er - ef : -1.0, 0.0);
          });
          nlegspaced++;
          // up and down followed from the station towards the junction
          std::vector<size_t> bord;
          for (size_t m = 0; m < fs.size(); m++) bord.push_back(m);
          std::vector<double> bup, bdown;
          std::vector<char> bok;
          std::vector<lxVec> bdir;
          track_elevation(i, j, fs, bord, bup, bdown, bok, bdir, false);
          for (size_t m = 1; (m < fs.size()) && (fmax > 0.0); m++) {
            double f = fs[m];
            lxVec c = a.c + f * (cj - a.c);
            if (w0 && drawn.width(c.x, c.y, c.z, h, wr, wl) && (wr + wl > BRANCH_OPENING * (w0r + w0l))) break;
            std::vector<lxVec> P;
            for (size_t k = 0; k < a.q.size(); k++) P.push_back(a.point(k) + (c - a.c));
            std::string where = fmt::format("{}@{} -> {}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name, e.name);
            stretch(P, c, h, where);
            if (bok[m]) vstretch(P, station_at(i) + f * (lxVec(e.x, e.y, e.z) - station_at(i)), bdir[m], -bdown[m], bup[m], where);
            std::vector<thdb3dvx *> cur = ring(P);
            loft(prev, cur);
            prev = cur;
            fend = f;
            Pend = P;
          }
          // the end ring is made vertical in the cut plane and left open: the
          // junction uses it as its opening to this branch
          if (fend == 0.0) {
            std::vector<thdb3dvx *> cur = ring(Pend);
            loft(prev, cur);
            prev = cur;
          }
          lxVec cend = a.c + fend * (cj - a.c), nh(-h.y, h.x, 0.0);
          for (size_t k = 0; k < prev.size(); k++) {
            lxVec v = prev[k]->get_vector();
            v = v - ((v - cend) * nh) * nh;
            prev[k]->x = v.x;
            prev[k]->y = v.y;
            prev[k]->z = v.z;
            Pend[k] = v;
          }
          branch_ends[j].push_back({cend, h, Pend, prev, false, i, false, steep});
        }
        a.ntubes++;
        nextended++;
        continue;
      }
      if ((j <= i) || !sec[j].ok) continue;
      section & a = sec[i], & b = sec[j];
      // match left and right of the two sections using the leg direction:
      // b's left is a's right if one normal points along the leg and the other against it
      lxVec d = b.centre() - a.centre();
      bool flip = ((a.n * d) < 0.0) != ((b.n * d) < 0.0);
      if (flip != ((a.n * b.n) < 0.0))
        flipreport += fmt::format("{:>20} - {:<20} normals {:.0f} deg apart, {:.0f} and {:.0f} deg to leg\n",
          fmt::format("{}@{}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name),
          fmt::format("{}@{}", d1.station_vec[j].name, d1.station_vec[j].survey->full_name),
          deg(acos(std::max(-1.0, std::min(1.0, a.n * b.n)))),
          deg(acos(std::max(-1.0, std::min(1.0, (a.n * d) / d.Length())))),
          deg(acos(std::max(-1.0, std::min(1.0, (b.n * d) / d.Length())))));
      // corner sequences on the 0-400 loop in a's frame; b mirrored if flipped
      size_t na = a.q.size(), nb = b.q.size();
      std::vector<size_t> bidx(nb);
      std::vector<double> bpos(nb);
      for (size_t k = 0; k < nb; k++) {
        bidx[k] = flip ? (nb - k) % nb : k;
        bpos[k] = flip ? fmod(600.0 - b.pos[bidx[k]], 400.0) : b.pos[k];
      }
      auto wrapdiff = [](double x, double y) {
        double dd = fmod(x - y + 600.0, 400.0) - 200.0;
        return dd;
      };
      // start at the closest pair of corners
      size_t i0 = 0, j0 = 0;
      for (size_t k = 0; k < na; k++)
        for (size_t m = 0; m < nb; m++)
          if (fabs(wrapdiff(a.pos[k], bpos[m])) < fabs(wrapdiff(a.pos[i0], bpos[j0]))) {
            i0 = k;
            j0 = m;
          }
      // unwrapped positions from the start pair, once round (+ closing corner)
      std::vector<double> ap(na + 1), bp(nb + 1);
      for (size_t k = 0; k <= na; k++) {
        ap[k] = a.pos[(i0 + k) % na] + (((i0 + k) >= na) ? 400.0 : 0.0);
        if (k == na) ap[k] = ap[0] + 400.0;
      }
      for (size_t k = 0; k <= nb; k++) {
        bp[k] = bpos[(j0 + k) % nb];
        if (k > 0)
          while (bp[k] < bp[k - 1] - 1e-9) bp[k] += 400.0;
        if (k == nb) bp[k] = bp[0] + 400.0;
      }
      double shift = ap[0] - bp[0] - wrapdiff(ap[0], bp[0]);
      for (double & x : bp) x += shift;
      // order-keeping pairing of the corners in between: most pairs, then smallest differences
      std::vector<std::vector<double>> dp(na, std::vector<double>(nb, -1e99));
      std::vector<std::vector<int>> how(na, std::vector<int>(nb, 0));
      dp[0][0] = 0.0;
      for (size_t k = 0; k < na; k++)
        for (size_t m = 0; m < nb; m++) {
          if ((k == 0) && (m == 0)) continue;
          double best = -1e99;
          int bh = 0;
          if ((k > 0) && (dp[k - 1][m] > best)) { best = dp[k - 1][m]; bh = 1; }
          if ((m > 0) && (dp[k][m - 1] > best)) { best = dp[k][m - 1]; bh = 2; }
          double dd = fabs(ap[k] - bp[m]);
          if ((k > 0) && (m > 0) && (dd <= PAIR_TOLERANCE) && (dp[k - 1][m - 1] + 1000.0 - dd > best)) {
            best = dp[k - 1][m - 1] + 1000.0 - dd;
            bh = 3;
          }
          dp[k][m] = best;
          how[k][m] = bh;
        }
      std::vector<std::pair<size_t, size_t>> pairs;
      pairs.push_back(std::make_pair(na, nb));
      for (size_t k = na - 1, m = nb - 1; (k > 0) || (m > 0);) {
        if (how[k][m] == 3) {
          pairs.push_back(std::make_pair(k, m));
          k--; m--;
        } else if (how[k][m] == 1) k--;
        else m--;
      }
      pairs.push_back(std::make_pair(0, 0));
      std::reverse(pairs.begin(), pairs.end());
      double maxdiff = fabs(ap[0] - bp[0]);
      for (auto & pr : pairs) maxdiff = std::max(maxdiff, fabs(ap[pr.first] - bp[pr.second]));
      tubereport += fmt::format("{:>20} - {:<20} {:3} of {:2}/{:2} corners joined directly, largest difference {:5.1f}\n",
        fmt::format("{}@{}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name),
        fmt::format("{}@{}", d1.station_vec[j].name, d1.station_vec[j].survey->full_name),
        pairs.size() - 1, na, nb, maxdiff);

      // points of a and b at a position within a span of corners
      auto at = [](const std::vector<double> & pp, size_t k0, size_t k1, double sfrac,
                   const std::function<lxVec(size_t)> & pt) {
        double span = pp[k1] - pp[k0];
        if (span <= 1e-9) {  // no position range: spread by corner count
          double fk = k0 + sfrac * (k1 - k0);
          size_t k = std::min(k1 - (k1 > k0 ? 1 : 0), (size_t) floor(fk));
          if (k1 == k0) return pt(k0);
          return pt(k) + (fk - k) * (pt(k + 1) - pt(k));
        }
        double target = pp[k0] + sfrac * span;
        size_t k = k0;
        while ((k + 1 < k1) && (pp[k + 1] <= target)) k++;
        double f = (pp[k + 1] > pp[k]) ? (target - pp[k]) / (pp[k + 1] - pp[k]) : 0.0;
        f = std::max(0.0, std::min(1.0, f));
        return pt(k) + f * (pt(k + 1) - pt(k));
      };
      std::function<lxVec(size_t)> apt = [&](size_t k) { return a.point((i0 + k) % na); };
      std::function<lxVec(size_t)> bpt = [&](size_t k) { return b.point(bidx[(j0 + k) % nb]); };
      std::vector<thdb3dvx *> ra, rb;
      for (size_t pi = 0; pi + 1 < pairs.size(); pi++) {
        size_t ka = pairs[pi].first, ka1 = pairs[pi + 1].first, kb = pairs[pi].second, kb1 = pairs[pi + 1].second;
        // fractions of the span at the corners of a and b in between
        std::vector<double> fr = {0.0};
        auto frac = [](const std::vector<double> & pp, size_t k0, size_t k1, size_t k) {
          double span = pp[k1] - pp[k0];
          return (span > 1e-9) ? (pp[k] - pp[k0]) / span : double(k - k0) / double(k1 - k0);
        };
        for (size_t k = ka + 1; k < ka1; k++) fr.push_back(frac(ap, ka, ka1, k));
        for (size_t k = kb + 1; k < kb1; k++) fr.push_back(frac(bp, kb, kb1, k));
        std::sort(fr.begin(), fr.end());
        for (double f : fr) {
          ra.push_back(walls.insert_vertex(at(ap, ka, ka1, f, apt)));
          rb.push_back(walls.insert_vertex(at(bp, kb, kb1, f, bpt)));
        }
      }
      loft_rings(walls.insert_face(THDB3DFC_TRIANGLES), ra, rb);
      if (usemap) {
        // sections every MAP_SPACING: the tube at that point, stretched to the drawn walls
        lxVec ha = horizontal(a.u), hb = horizontal(b.u);
        if (flip) hb = -1.0 * hb;
        double len = std::hypot(b.c.x - a.c.x, b.c.y - a.c.y);
        std::vector<double> fs = positions(len, 1.0, [&](double f) {
          lxVec c = (1.0 - f) * a.c + f * b.c, h = horizontal((1.0 - f) * ha + f * hb);
          double r1, l1, ef, er;
          lxVec ed;
          return lxVec(drawn.width(c.x, c.y, c.z, h, r1, l1) ? r1 + l1 : -1.0,
            leg_elevation(i, j, f, ef, er, ed) ? er - ef : -1.0, 0.0);
        });
        nlegspaced++;
        size_t nseg = fs.size() - 1;
        std::vector<thdb3dvx *> prev;
        for (size_t m = 0; m <= nseg; m++) {
          double f = fs[m];
          std::vector<lxVec> P;
          for (size_t k = 0; k < ra.size(); k++)
            P.push_back((1.0 - f) * ra[k]->get_vector() + f * rb[k]->get_vector());
          lxVec c = (1.0 - f) * a.c + f * b.c, h = (1.0 - f) * ha + f * hb;
          h = horizontal(h);
          std::string where = fmt::format("{}@{} - {} {:.0f}%", d1.station_vec[i].name,
            d1.station_vec[i].survey->full_name, d1.station_vec[j].name, 100.0 * f);
          // a drawn section wins at its own station
          if (!(((m == 0) && drawnsec[i]) || ((m == nseg) && drawnsec[j]))) {
            stretch(P, c, h, where);
            double ef, er;
            lxVec ed;
            bool ok = (m == 0) ? station_elevation(i, ef, er, ed) : ((m == nseg) ? station_elevation(j, ef, er, ed)
              : leg_elevation(i, j, f, ef, er, ed));
            if (ok) vstretch(P, (m == 0) ? station_at(i) : ((m == nseg) ? station_at(j)
              : station_at(i) + f * (station_at(j) - station_at(i))), ed, ef, er, where);
          }
          std::vector<thdb3dvx *> cur = ring(P);
          if (m > 0) loft(prev, cur);
          prev = cur;
        }
      }
      a.ntubes++;
      b.ntubes++;
      ntubes++;
    }
  }

  // Legs with no section at either end (dead ends without splays, chains of
  // stations without sections), where drawn: at each position the LRUD from the
  // drawings (left and right to the plan walls, up and down to the extended
  // elevation's roof and floor) gives an oval section, each quarter an ellipse;
  // with only one of the drawings, a circle. From a junction (or station without
  // a section) the oval leg starts OVAL_END away and joins its footprint; at a
  // dead end it is rounded off to the drawn end wall.
  size_t novallegs = 0, novalrings = 0;
  std::vector<double> ovup, ovdown;
  std::vector<char> ovok;
  for (size_t i = 0; (i < nst) && usemap; i++) {
    thdb1ds & si = d1.station_vec[i];
    if ((si.uid != i + 1) || si.is_temporary() || !si.survey->is_selected() || sec[i].ok) continue;
    for (size_t j : nbrs[i]) {
      if (sec[j].ok ? !side(j, i) : (j <= i)) continue;
      bool bi = body[i], bj = body[j] || side(j, i);  // ends joining a junction
      thdb1ds & sj = d1.station_vec[j];
      lxVec ci(si.x, si.y, si.z), cj(sj.x, sj.y, sj.z), dir(cj.x - ci.x, cj.y - ci.y, 0.0);
      double len = dir.Length();
      if (len < 0.1) continue;
      dir = dir / len;
      lxVec h(dir.y, -dir.x, 0.0);  // right, looking from i to j
      // where the passage has settled, from a junction end: the first position
      // (at least spacing_min from it, at most halfway) where the drawn width is
      // no more than BRANCH_OPENING times the typical (median) width of the other
      // half of the leg; without plan widths, OVAL_END from the junction
      double fmin = std::min(std::max(0.05, thsplaytube_opt.spacing_min) / len, 0.5);
      int nsmp = std::max(2, (int) ceil(len / std::max(0.05, thsplaytube_opt.spacing_min)));
      std::vector<double> wsmp(nsmp + 1, -1.0);
      for (int k = 0; k <= nsmp; k++) {
        lxVec c = ci + (double(k) / nsmp) * (cj - ci);
        double r1, l1;
        if (drawn.width(c.x, c.y, c.z, h, r1, l1)) wsmp[k] = r1 + l1;
      }
      auto settle = [&](bool from_i) {
        std::vector<double> far;
        for (int k = 0; k <= nsmp; k++)
          if ((wsmp[k] > 0.0) && (from_i ? (2 * k >= nsmp) : (2 * k <= nsmp))) far.push_back(wsmp[k]);
        if (far.empty()) return std::min(OVAL_END / len, 0.5);
        std::sort(far.begin(), far.end());
        double ref = far[far.size() / 2];
        for (int m = 0; 2 * m <= nsmp; m++) {
          int k = from_i ? m : nsmp - m;
          double fr = double(m) / nsmp;
          if (fr < fmin) continue;
          if ((wsmp[k] > 0.0) && (wsmp[k] <= BRANCH_OPENING * ref)) return fr;
        }
        return 0.5;
      };
      double f0 = bi ? settle(true) : 0.0, f1 = bj ? 1.0 - settle(false) : 1.0;
      if (f1 - f0 < 1e-3) continue;
      std::vector<double> fs = positions(len * (f1 - f0), 1.0, [&](double g) {
        double f = f0 + g * (f1 - f0), r1, l1, ef, er;
        lxVec c = ci + f * (cj - ci), ed;
        return lxVec(drawn.width(c.x, c.y, c.z, h, r1, l1) ? r1 + l1 : -1.0,
          leg_elevation(i, j, f, ef, er, ed, true) ? er - ef : -1.0, 0.0);
      });
      for (double & g : fs) g = f0 + g * (f1 - f0);
      // up and down followed along the leg from the end that is not a junction
      // (from the middle if both or neither are)
      size_t nf = fs.size();
      std::vector<size_t> fwd, back;
      if (bi && !bj) {
        for (size_t k = nf; k > 0; k--) back.push_back(k - 1);
      } else if (!bi && bj) {
        for (size_t k = 0; k < nf; k++) fwd.push_back(k);
      } else {
        for (size_t k = nf / 2; k < nf; k++) fwd.push_back(k);
        for (size_t k = nf / 2 + 1; k > 0; k--) back.push_back(k - 1);
      }
      ovup.assign(nf, 0.0);
      ovdown.assign(nf, 0.0);
      ovok.assign(nf, 0);
      for (const std::vector<size_t> * od : {&back, &fwd}) {
        if (od->empty()) continue;
        std::vector<double> tu, td;
        std::vector<char> tok;
        std::vector<lxVec> tdir;
        track_elevation(i, j, fs, *od, tu, td, tok, tdir, true);
        for (size_t k : *od) { ovup[k] = tu[k]; ovdown[k] = td[k]; ovok[k] = tok[k]; }
      }
      // oval ring at position k, or empty where nothing is drawn
      auto oval = [&](size_t k) {
        std::vector<lxVec> P;
        double f = fs[k], l = 0, r = 0, u = ovup[k], d = ovdown[k];
        lxVec c = ci + f * (cj - ci);
        bool horiz = drawn.width(c.x, c.y, c.z, h, r, l), vert = ovok[k];
        if (!horiz && !vert) return P;
        u = std::max(0.05, u);
        d = std::max(0.05, d);
        double ox = 0.0, oz = 0.0;
        if (!horiz) { double rr = 0.5 * (u + d); oz = 0.5 * (u - d); l = r = u = d = rr; }
        else if (!vert) { double rr = 0.5 * (l + r); ox = 0.5 * (r - l); l = r = u = d = rr; }
        for (int q = 0; q < OVAL_POINTS; q++) {
          double a = 2.0 * M_PI * q / OVAL_POINTS, ca = cos(a), sa = sin(a);
          P.push_back(c + (ox + ca * ((ca > 0.0) ? r : l)) * h + lxVec(0.0, 0.0, oz + sa * ((sa > 0.0) ? u : d)));
        }
        return P;
      };
      std::vector<thdb3dvx *> prev, first;
      std::vector<lxVec> Pfirst, Plast;
      double ffirst = -1.0, flast = -1.0;
      for (size_t k = 0; k < nf; k++) {
        double f = fs[k];
        std::vector<lxVec> P = oval(k);
        if (P.empty()) continue;  // not drawn here
        std::vector<thdb3dvx *> cur = ring(P);
        novalrings++;
        if (!prev.empty()) loft(prev, cur);
        else { first = cur; Pfirst = P; ffirst = f; }
        prev = cur;
        Plast = P;
        flast = f;
      }
      if (prev.empty()) continue;
      novallegs++;
      // ends: join a junction footprint, or round off a dead end to the drawn end wall
      auto close_end = [&](const std::vector<thdb3dvx *> & r, const std::vector<lxVec> & P, double f, bool at_j) {
        size_t st_end = at_j ? j : i, st_other = at_j ? i : j;
        lxVec c = ci + f * (cj - ci), out = at_j ? dir : -1.0 * dir;
        if (at_j ? bj : bi) {
          branch_ends[st_end].push_back({c, at_j ? h : -1.0 * h, P, r, false, st_other, true,
            deg(atan2(std::fabs(cj.z - ci.z), len)) > ELEVATION_STEEP});
          return;
        }
        double wall = -1.0;
        for (size_t si2 = 0; si2 < drawn.scraps.size(); si2++) {
          if (!drawn.at_level(si2, c.z) || !drawn.inside(si2, c.x, c.y)) continue;
          double t = drawn.ray(si2, c.x, c.y, out.x, out.y);
          if ((t > 0.0) && ((wall < 0.0) || (t < wall))) wall = t;
        }
        lxVec cc = centroid(r);
        std::vector<thdb3dvx *> rr = r;
        if (wall > 0.05) {  // half an ellipsoid out to the end wall
          for (int k = 1; k < 4; k++) {
            double t = k / 4.0, sc = std::sqrt(1.0 - t * t);
            std::vector<lxVec> Q;
            for (const lxVec & p : P) Q.push_back(cc + sc * (p - cc) + (t * wall) * out);
            std::vector<thdb3dvx *> nr = ring(Q);
            loft(rr, nr);
            rr = nr;
          }
          taper_ring(map_tubes.insert_face(THDB3DFC_TRIANGLES), rr, map_tubes.insert_vertex(cc + wall * out));
        } else {
          piece pc;
          thdb3dvx * cv = map_tubes.insert_vertex(cc);
          for (size_t k = 0; k < rr.size(); k++) pc.add(cv, rr[(k + 1) % rr.size()], rr[k], out);
          pc.emit(map_tubes.insert_face(THDB3DFC_TRIANGLES));
        }
      };
      close_end(first, Pfirst, ffirst, false);
      close_end(prev, Plast, flast, true);
    }
  }

  // close open tube ends
  for (size_t i = 0; i < nst; i++) {
    section & s = sec[i];
    if (!s.ok || (s.ntubes >= 2)) continue;
    // outwards: away from the neighbour the tube goes to
    lxVec out = s.n;
    for (size_t j : nbrs[i])
      if (sec[j].ok || body[j] || (nbrs[j].size() == 1)) {
        thdb1ds & nb = d1.station_vec[j];
        out = s.c - lxVec(nb.x, nb.y, nb.z);
        break;
      }
    if (cap(s, s.vx, walls.insert_vertex(s.c), walls, out)) ncaps++; else nfans++;
    if (usemap) {
      if (cap(s, ring(station_ring(i)), map_tubes.insert_vertex(s.c), map_tubes, out)) ncaps++; else nfans++;
    }
  }

  // Junctions in the model from maps and sections: the drawn outline around the
  // junction, cut across each branch where the branch ends; roof and floor from
  // the shapes of the branch ends and the junction splays.
  std::string junctionreport;
  // height of the top and bottom of a branch end ring at a position across it
  auto ring_extent = [](const branch_end & be, double a, double & top, double & bot, bool clamp) {
    size_t n = be.P.size();
    std::vector<double> ac(n);
    double amin = 1e99, amax = -1e99;
    for (size_t k = 0; k < n; k++) {
      ac[k] = (be.P[k].x - be.c.x) * be.h.x + (be.P[k].y - be.c.y) * be.h.y;
      amin = std::min(amin, ac[k]);
      amax = std::max(amax, ac[k]);
    }
    if (clamp) a = std::max(amin, std::min(amax, a));
    top = -1e99;
    bot = 1e99;
    for (size_t k = 0; k < n; k++) {
      size_t k1 = (k + 1) % n;
      double a0 = ac[k], a1 = ac[k1];
      if ((a0 - a) * (a1 - a) > 0.0) continue;
      if (std::fabs(a1 - a0) > 1e-9) {
        double z = be.P[k].z + (a - a0) / (a1 - a0) * (be.P[k1].z - be.P[k].z);
        top = std::max(top, z);
        bot = std::min(bot, z);
      } else {
        top = std::max(top, std::max(be.P[k].z, be.P[k1].z));
        bot = std::min(bot, std::min(be.P[k].z, be.P[k1].z));
      }
    }
    return top >= bot;
  };
  for (size_t j = 0; (j < nst) && usemap; j++) {
    if (!body[j]) continue;
    std::vector<branch_end> & bends = branch_ends[j];
    thdb1ds & js = d1.station_vec[j];
    std::string jname = fmt::format("{}@{}", js.name, js.survey->full_name);
    // the outline of the scrap the junction is drawn in (not one of a passage
    // above or below that only covers it in plan)
    const ring2 * outline = NULL;
    for (size_t si : drawn.around(j, js.x, js.y, js.z)) {
      for (const ring2 & r : drawn.rings[si])
        if (in_ring(r, js.x, js.y)) { outline = &r; break; }
      if (outline != NULL) break;
    }
    if (outline == NULL) {
      junctionreport += fmt::format("{:>20}  not inside a drawn outline\n", jname);
      continue;
    }
    ring2 fp = *outline;
    size_t ncut = 0;
    // boundary points on a cut: heights fixed to the branch end ring, and which cut
    std::map<std::pair<long, long>, std::array<double, 3>> fixed;
    auto key = [](double x, double y) { return std::make_pair(std::lround(x * 1e4), std::lround(y * 1e4)); };
    for (size_t bi = 0; bi < bends.size(); bi++) {
      branch_end & be = bends[bi];
      if (be.steep) continue;  // run into the junction and closed
      double tr, tl;
      size_t er, el;
      if (!ray_ring(fp, be.c.x, be.c.y, be.h.x, be.h.y, tr, er) ||
          !ray_ring(fp, be.c.x, be.c.y, -be.h.x, -be.h.y, tl, el) || (er == el)) continue;
      std::array<double, 2> pr = {be.c.x + tr * be.h.x, be.c.y + tr * be.h.y},
        pl = {be.c.x - tl * be.h.x, be.c.y - tl * be.h.y};
      // points along the cut: the positions across of the branch end ring
      std::vector<double> across;
      for (const lxVec & p : be.P) {
        double a = (p.x - be.c.x) * be.h.x + (p.y - be.c.y) * be.h.y;
        if ((a > -tl + 1e-3) && (a < tr - 1e-3)) across.push_back(a);
      }
      std::sort(across.begin(), across.end());
      across.erase(std::unique(across.begin(), across.end(), [](double x, double y) { return std::fabs(x - y) < 1e-3; }), across.end());
      // the two loops on either side of the cut; keep the one with the junction
      size_t n = fp.size();
      ring2 l1, l2;
      l1.push_back(pr);
      for (size_t k = (er + 1) % n; ; k = (k + 1) % n) { l1.push_back(fp[k]); if (k == el) break; }
      l1.push_back(pl);
      for (size_t k = 0; k < across.size(); k++)  // pl -> pr along the cut
        l1.push_back({be.c.x + across[k] * be.h.x, be.c.y + across[k] * be.h.y});
      l2.push_back(pl);
      for (size_t k = (el + 1) % n; ; k = (k + 1) % n) { l2.push_back(fp[k]); if (k == er) break; }
      l2.push_back(pr);
      for (size_t k = across.size(); k > 0; k--)  // pr -> pl along the cut
        l2.push_back({be.c.x + across[k - 1] * be.h.x, be.c.y + across[k - 1] * be.h.y});
      if (in_ring(l1, js.x, js.y)) fp = l1;
      else if (in_ring(l2, js.x, js.y)) fp = l2;
      else continue;
      std::vector<double> all = across;
      all.push_back(tr);
      all.push_back(-tl);
      for (double a : all) {
        double t, b;
        if (!ring_extent(be, a, t, b, true)) continue;
        fixed[key(be.c.x + a * be.h.x, be.c.y + a * be.h.y)] = {t, b, double(bi)};
      }
      be.joined = true;
      ncut++;
    }
    // legs to neighbours without a branch end (junctions, stations without a
    // section): cut across the leg halfway, so the footprint covers only its
    // own stretch of passage; passage ends without a section are kept whole
    size_t nhalf = 0;
    for (size_t k : nbrs[j]) {
      bool hasend = false;
      for (const branch_end & be : bends)
        if (be.sec == k) hasend = true;
      if (hasend || (!sec[k].ok && (nbrs[k].size() == 1))) continue;
      thdb1ds & ks = d1.station_vec[k];
      lxVec c(0.5 * (js.x + ks.x), 0.5 * (js.y + ks.y), 0.5 * (js.z + ks.z)), h(-(ks.y - js.y), ks.x - js.x, 0.0);
      if (h.Length() < 1e-6) continue;
      h.Normalize();
      double tr, tl;
      size_t er, el;
      if (!ray_ring(fp, c.x, c.y, h.x, h.y, tr, er) || !ray_ring(fp, c.x, c.y, -h.x, -h.y, tl, el) || (er == el)) continue;
      std::array<double, 2> pr = {c.x + tr * h.x, c.y + tr * h.y}, pl = {c.x - tl * h.x, c.y - tl * h.y};
      size_t n = fp.size();
      ring2 l1, l2;
      l1.push_back(pr);
      for (size_t q = (er + 1) % n; ; q = (q + 1) % n) { l1.push_back(fp[q]); if (q == el) break; }
      l1.push_back(pl);
      l2.push_back(pl);
      for (size_t q = (el + 1) % n; ; q = (q + 1) % n) { l2.push_back(fp[q]); if (q == er) break; }
      l2.push_back(pr);
      if (in_ring(l1, js.x, js.y)) fp = l1;
      else if (in_ring(l2, js.x, js.y)) fp = l2;
      else continue;
      nhalf++;
    }
    // its own stretch of passage only: where a passage of the same outline
    // continues without a leg from this junction (another passage drawn in the
    // same scrap, or one at the same level under or over it), the part nearer to
    // a point of another leg than to the junction's own legs (out to halfway) is
    // cut off, across the line halfway between them (nearest first); on legs at
    // the junction's neighbours only points well past the neighbour count, as
    // their branches may end at the neighbour
    size_t nlimit = 0;
    {
      double reach = 0.0;
      for (const auto & q : fp) reach = std::max(reach, std::hypot(q[0] - js.x, q[1] - js.y));
      std::vector<std::pair<double, std::array<double, 2>>> sites;
      auto isnbr = [&](size_t k) { return std::find(nbrs[j].begin(), nbrs[j].end(), k) != nbrs[j].end(); };
      for (size_t a = 0; a < nst; a++) {
        if (a == j) continue;
        thdb1ds & sa = d1.station_vec[a];
        for (size_t b : nbrs[a]) {
          if ((b <= a) || (b == j)) continue;
          thdb1ds & sb = d1.station_vec[b];
          int n = std::max(1, int(ceil(std::hypot(sb.x - sa.x, sb.y - sa.y) / JUNCTION_SITE)));
          for (int m = 0; m <= n; m++) {
            double f = double(m) / n, x = sa.x + f * (sb.x - sa.x), y = sa.y + f * (sb.y - sa.y), z = sa.z + f * (sb.z - sa.z);
            double d = std::hypot(x - js.x, y - js.y);
            if ((std::fabs(z - js.z) > FILL_LEVEL) || (d < 1e-3) || (d > reach)) continue;
            // beyond a neighbour, only well past it (its branch may end there)
            bool close = false;
            for (size_t k : {a, b})
              if (isnbr(k)) {
                thdb1ds & ks = d1.station_vec[k];
                if (std::hypot(x - ks.x, y - ks.y) < 0.5 * std::hypot(ks.x - js.x, ks.y - js.y)) close = true;
              }
            if (close) continue;
            sites.push_back(std::make_pair(d, std::array<double, 2>{x, y}));
          }
        }
      }
      std::sort(sites.begin(), sites.end(), [](const std::pair<double, std::array<double, 2>> & p, const std::pair<double, std::array<double, 2>> & q) { return p.first < q.first; });
      for (const auto & site : sites) {
        // the nearest point of the junction's own legs, out to halfway
        double ox = js.x, oy = js.y, best = 1e99;
        for (size_t k : nbrs[j]) {
          thdb1ds & ks = d1.station_vec[k];
          double dx = ks.x - js.x, dy = ks.y - js.y, l2 = dx * dx + dy * dy;
          double f = (l2 > 1e-9) ? std::max(0.0, std::min(0.5, ((site.second[0] - js.x) * dx + (site.second[1] - js.y) * dy) / l2)) : 0.0;
          double px = js.x + f * dx, py = js.y + f * dy, d = std::hypot(site.second[0] - px, site.second[1] - py);
          if (d < best) { best = d; ox = px; oy = py; }
        }
        if (best < 1e-3) continue;
        double mx = 0.5 * (ox + site.second[0]), my = 0.5 * (oy + site.second[1]);
        if (!in_ring(fp, mx, my)) continue;
        lxVec h(-(site.second[1] - oy), site.second[0] - ox, 0.0);
        h.Normalize();
        double tr, tl;
        size_t er, el;
        if (!ray_ring(fp, mx, my, h.x, h.y, tr, er) || !ray_ring(fp, mx, my, -h.x, -h.y, tl, el) || (er == el)) continue;
        std::array<double, 2> pr = {mx + tr * h.x, my + tr * h.y}, pl = {mx - tl * h.x, my - tl * h.y};
        size_t n = fp.size();
        ring2 l1, l2;
        l1.push_back(pr);
        for (size_t q = (er + 1) % n; ; q = (q + 1) % n) { l1.push_back(fp[q]); if (q == el) break; }
        l1.push_back(pl);
        l2.push_back(pl);
        for (size_t q = (el + 1) % n; ; q = (q + 1) % n) { l2.push_back(fp[q]); if (q == er) break; }
        l2.push_back(pr);
        if (in_ring(l1, js.x, js.y)) fp = l1;
        else if (in_ring(l2, js.x, js.y)) fp = l2;
        else continue;
        nlimit++;
      }
    }
    // Roof and floor: at a point of the footprint, the drawn extended elevation
    // of the nearest leg there. Where that is not drawn, each branch end gives
    // the heights of its section at the same position across the passage (its
    // outermost point beyond its walls); branches nearer the point count more,
    // and junction splay ends above the roof (below the floor) raise (lower) it nearby.
    auto branch_heights = [&](double x, double y, double & top, double & bot) {
      double sw = 0.0, st = 0.0, sb = 0.0;
      for (const branch_end & be : bends) {
        if (be.steep) continue;
        size_t n = be.P.size();
        std::vector<double> ac(n);
        double amin = 1e99, amax = -1e99;
        for (size_t k = 0; k < n; k++) {
          ac[k] = (be.P[k].x - be.c.x) * be.h.x + (be.P[k].y - be.c.y) * be.h.y;
          amin = std::min(amin, ac[k]);
          amax = std::max(amax, ac[k]);
        }
        double a = (x - be.c.x) * be.h.x + (y - be.c.y) * be.h.y;
        a = std::max(amin, std::min(amax, a));
        double t = -1e99, b = 1e99;
        for (size_t k = 0; k < n; k++) {
          size_t k1 = (k + 1) % n;
          double a0 = ac[k], a1 = ac[k1];
          if (((a0 - a) * (a1 - a) > 0.0)) continue;
          double z = (std::fabs(a1 - a0) > 1e-9) ? be.P[k].z + (a - a0) / (a1 - a0) * (be.P[k1].z - be.P[k].z)
            : std::max(be.P[k].z, be.P[k1].z);
          double zl = (std::fabs(a1 - a0) > 1e-9) ? z : std::min(be.P[k].z, be.P[k1].z);
          t = std::max(t, z);
          b = std::min(b, zl);
        }
        if (t < b) continue;
        double d2 = (x - be.c.x) * (x - be.c.x) + (y - be.c.y) * (y - be.c.y), w = 1.0 / (d2 + 0.01);
        sw += w;
        st += w * t;
        sb += w * b;
      }
      if (sw <= 0.0) {  // no branch sections: from the station's splays
        double rw = 0.0, rz = 0.0, fw = 0.0, fz = 0.0;
        for (thdb1d_tree_arrow * ar = nodes[j].first_arrow; ar != NULL; ar = ar->next_arrow) {
          if ((ar->leg->leg->flags & TT_LEGFLAG_SPLAY) == 0) continue;
          thdb1ds & e = d1.station_vec[ar->end_node->uid - 1];
          double w = 1.0 / ((e.x - x) * (e.x - x) + (e.y - y) * (e.y - y) + 0.01);
          if (e.z >= js.z) { rw += w; rz += w * e.z; } else { fw += w; fz += w * e.z; }
        }
        top = (rw > 0.0) ? rz / rw : js.z + 1.0;
        bot = (fw > 0.0) ? fz / fw : js.z - 0.5;
        return;
      }
      top = st / sw;
      bot = sb / sw;
    };
    std::vector<lxVec> jsplays;
    for (thdb1d_tree_arrow * ar = nodes[j].first_arrow; ar != NULL; ar = ar->next_arrow) {
      if ((ar->leg->leg->flags & TT_LEGFLAG_SPLAY) == 0) continue;
      thdb1ds & e = d1.station_vec[ar->end_node->uid - 1];
      if (in_ring(fp, e.x, e.y)) jsplays.push_back(lxVec(e.x, e.y, e.z));
    }
    // from the extended elevation: the drawn span of the nearest leg (its
    // stretch from the junction to halfway) at the point's position along it
    auto elevation_heights = [&](double x, double y, double & top, double & bot) {
      double best = 1e99, bf = 0.0;
      size_t bk = nst;
      for (size_t k : nbrs[j]) {
        thdb1ds & ks = d1.station_vec[k];
        double dx = ks.x - js.x, dy = ks.y - js.y, l2 = dx * dx + dy * dy;
        if (l2 < 1e-6) continue;
        double f = std::max(0.0, std::min(0.5, ((x - js.x) * dx + (y - js.y) * dy) / l2));
        double d = std::hypot(js.x + f * dx - x, js.y + f * dy - y);
        if (d < best) { best = d; bf = f; bk = k; }
      }
      double fl, rf;
      lxVec ed;
      if ((bk == nst) || !leg_elevation(j, bk, bf, fl, rf, ed, true) || (rf < fl + 0.1)) return false;
      double zb = js.z + bf * (d1.station_vec[bk].z - js.z);
      top = zb + rf;
      bot = zb + fl;
      return true;
    };
    size_t nelevheights = 0, nheightpts = 0;
    auto heights = [&](double x, double y, double & top, double & bot) {
      nheightpts++;
      if (elevation_heights(x, y, top, bot)) {
        nelevheights++;
        return;
      }
      branch_heights(x, y, top, bot);
      for (const lxVec & e : jsplays) {
        double d = std::hypot(e.x - x, e.y - y);
        if (d >= 1.0) continue;
        double w = (1.0 - d) * (1.0 - d);
        double et, eb;
        branch_heights(e.x, e.y, et, eb);
        if (e.z > et) top = std::max(top, top + w * (e.z - et));
        if (e.z < eb) bot = std::min(bot, bot + w * (e.z - eb));
      }
      if (top < bot + 0.1) {
        double zm = 0.5 * (top + bot);
        top = zm + 0.05;
        bot = zm - 0.05;
      }
    };
    // triangulate the footprint with points inside on a grid; roof, floor and walls round it
    std::vector<std::unique_ptr<p2t::Point>> holder;
    std::vector<p2t::Point *> line;
    std::map<p2t::Point *, size_t> idx;
    std::vector<p2t::Point *> allpts;
    for (size_t k = 0; k < fp.size(); k++) {
      if (!line.empty() && (std::hypot(fp[k][0] - line.back()->x, fp[k][1] - line.back()->y) < 1e-4)) continue;
      holder.push_back(std::make_unique<p2t::Point>(fp[k][0], fp[k][1]));
      line.push_back(holder.back().get());
    }
    if ((line.size() > 1) && (std::hypot(line[0]->x - line.back()->x, line[0]->y - line.back()->y) < 1e-4)) line.pop_back();
    if (line.size() < 3) continue;
    for (p2t::Point * pt : line) { idx[pt] = allpts.size(); allpts.push_back(pt); }
    size_t nboundary = allpts.size();
    double xmin = 1e99, xmax = -1e99, ymin = 1e99, ymax = -1e99;
    for (p2t::Point * pt : line) {
      xmin = std::min(xmin, pt->x); xmax = std::max(xmax, pt->x);
      ymin = std::min(ymin, pt->y); ymax = std::max(ymax, pt->y);
    }
    std::vector<p2t::Point *> steiner;
    const double grid = 0.5;
    for (double gx = xmin + grid / 2; gx < xmax; gx += grid)
      for (double gy = ymin + grid / 2; gy < ymax; gy += grid) {
        if (!in_ring(fp, gx, gy)) continue;
        bool near = false;
        for (size_t k = 0, n = fp.size(); (k < n) && !near; k++) {
          const auto & a = fp[k], & b = fp[(k + 1) % n];
          double dx = b[0] - a[0], dy = b[1] - a[1], l2 = dx * dx + dy * dy;
          double t = (l2 > 0.0) ? std::max(0.0, std::min(1.0, ((gx - a[0]) * dx + (gy - a[1]) * dy) / l2)) : 0.0;
          if (std::hypot(a[0] + t * dx - gx, a[1] + t * dy - gy) < 0.15) near = true;
        }
        if (near) continue;
        holder.push_back(std::make_unique<p2t::Point>(gx, gy));
        steiner.push_back(holder.back().get());
        idx[steiner.back()] = allpts.size();
        allpts.push_back(steiner.back());
      }
    std::vector<thdb3dvx *> top, bottom;
    std::vector<int> cutof(allpts.size(), -1);
    for (size_t k = 0; k < allpts.size(); k++) {
      p2t::Point * pt = allpts[k];
      double zr, zf;
      auto fx = fixed.find(key(pt->x, pt->y));
      if ((k < nboundary) && (fx != fixed.end())) {
        zr = fx->second[0];
        zf = fx->second[1];
        cutof[k] = int(fx->second[2]);
      } else heights(pt->x, pt->y, zr, zf);
      top.push_back(map_junctions.insert_vertex(pt->x, pt->y, zr));
      bottom.push_back(map_junctions.insert_vertex(pt->x, pt->y, zf));
    }
    thdb3dfc * fc = map_junctions.insert_face(THDB3DFC_TRIANGLES);
    try {
      p2t::CDT cdt(line);
      for (p2t::Point * sp : steiner) cdt.AddPoint(sp);
      cdt.Triangulate();
      for (p2t::Triangle * t : cdt.GetTriangles()) {
        size_t i0 = idx.at(t->GetPoint(0)), i1 = idx.at(t->GetPoint(1)), i2 = idx.at(t->GetPoint(2));
        tri(fc, top[i0], top[i1], top[i2], lxVec(0.0, 0.0, 1.0));
        tri(fc, bottom[i0], bottom[i1], bottom[i2], lxVec(0.0, 0.0, -1.0));
      }
    } catch (...) {
      junctionreport += fmt::format("{:>20}  footprint could not be triangulated\n", jname);
    }
    double sarea = 0.0;  // footprint orientation: > 0 anticlockwise
    for (size_t k = 0, n = nboundary; k < n; k++)
      sarea += allpts[k]->x * allpts[(k + 1) % n]->y - allpts[(k + 1) % n]->x * allpts[k]->y;
    for (size_t k = 0, n = nboundary; k < n; k++) {
      size_t k1 = (k + 1) % n;
      if ((cutof[k] >= 0) && (cutof[k] == cutof[k1])) continue;  // opening to a branch
      double dx = allpts[k1]->x - allpts[k]->x, dy = allpts[k1]->y - allpts[k]->y;
      lxVec out = (sarea > 0.0) ? lxVec(dy, -dx, 0.0) : lxVec(-dy, dx, 0.0);
      tri(fc, top[k], bottom[k], bottom[k1], out);
      tri(fc, top[k], bottom[k1], top[k1], out);
    }
    junctionreport += fmt::format("{:>20}  footprint from the drawn walls: {} corners, {} inside, {} of {} branches cut, "
      "{} legs cut halfway, {} cuts towards other legs, {} splays, heights from the extended elevation at {} of {} points\n",
      jname, line.size(), steiner.size(), ncut, bends.size(), nhalf, nlimit, jsplays.size(), nelevheights, nheightpts);
  }

  thlog(fmt::format("\n\n########################### splay tubes ###########################\n"
    "Section plane at each station: leg bisector (vertical unless steep), turned to the best\n"
    "fit of splays within {:.0f} deg of it; section splays are within {:.0f} deg of that plane.\n"
    "Splays inclined less than {:.0f} deg are measured by their direction seen from above\n"
    "(where the plane is within {:.0f} deg of vertical).\n"
    "turn = angle between fitted plane and bisector, spread = RMS angle of the section splays\n"
    "to the plane, worst = largest of these (all in degrees); centre = distance from the\n"
    "station to the centre of the section polygon (m), around which tubes join sections.\n\n"
    "{:>20} {:>4} {:>6} {:>4}  {:>5} {:>6} {:>5} {:>6}  {}\n{}",
    CANDIDATE_ANGLE, SECTION_ANGLE, LEVEL_SPLAY, TILTED_PLANE, "STATION", "LEGS", "SPLAYS", "USED", "TURN", "SPREAD", "WORST", "CENTRE", "PLANE",
    report));
  thlog(fmt::format("\nDrawn cross-sections (scraps with -projection none placed in the plan): {} found, {} placed, "
    "{} used as their station's section (not stretched to the drawn plan or elevation there).\n{}",
    xsecs.size(), nxsplaced, nxsused, xsreport));
  thlog(fmt::format("\nJunctions with a cross section (legs within {:.0f} deg of straight on, fitted plane,\n"
    "at least {} section splays, all round the passage): {} of {}.\n{}",
    JUNCTION_STRAIGHT, JUNCTION_SPLAYS, njsections, njunctions, jsreport));
  thlog(fmt::format("\nTube joins: corners of neighbouring sections joined directly when their\n"
    "positions (0-400 round the section) differ by at most {:.0f}:\n{}", PAIR_TOLERANCE, tubereport));
  if (!flipreport.empty())
    thlog("\nTubes where left and right are matched differently by leg direction and by normals:\n"
      + flipreport + "\n");
  // branches not joined to a junction run a little into it and are closed
  size_t noverlap = 0;
  for (auto & bj : branch_ends)
    for (branch_end & be : bj.second) {
      if (be.joined) continue;
      thdb1ds & js = d1.station_vec[bj.first];
      lxVec dir(js.x - be.c.x, js.y - be.c.y, 0.0);
      if (dir.Length() > 1e-6) dir.Normalize();
      std::vector<lxVec> P = be.P;
      for (lxVec & p : P) p += BRANCH_OVERLAP * dir;
      std::vector<thdb3dvx *> cur = ring(P);
      loft(be.vx, cur);
      if (be.oval) {
        piece pc;
        thdb3dvx * cv = map_tubes.insert_vertex(be.c + BRANCH_OVERLAP * dir);
        for (size_t k = 0; k < cur.size(); k++) pc.add(cv, cur[(k + 1) % cur.size()], cur[k], dir);
        pc.emit(map_tubes.insert_face(THDB3DFC_TRIANGLES));
      } else if (cap(sec[be.sec], cur, map_tubes.insert_vertex(be.c + BRANCH_OVERLAP * dir), map_tubes, dir)) ncaps++; else nfans++;
      noverlap++;
    }


  // Passages drawn in the plan without legs or splays. For each plan scrap, a
  // grid marks the cells inside the drawn walls and those covered by the model's
  // passages and junctions at about the scrap's level. A piece of drawn plan not
  // covered is a passage when it is deeper than the width of its opening onto
  // the model (and at least FILL_PASSAGE_DEPTH deep); it gets a tube along its
  // centre line, out to the drawn walls, with the parent passage's roof and
  // floor at the opening, run back into the parent and rounded off at its end.
  // Other pieces (slivers along the model, bays) are left.
  std::string fillreport;
  size_t nfills = 0, nfillcells = 0, nfillleft = 0;
  double fillleft = 0.0;
  if (usemap) {
    std::vector<std::array<thdb3dvx *, 3>> mtris;
    for (thdb3ddata * dd : {&map_tubes, &map_junctions})
      for (thdb3dfc * fc = dd->firstfc; fc != NULL; fc = fc->next) {
        if (fc->type != THDB3DFC_TRIANGLES) continue;
        for (thdb3dfx * fx = fc->firstfx; (fx != NULL) && (fx->next != NULL) && (fx->next->next != NULL); fx = fx->next->next->next)
          mtris.push_back({fx->vertex, fx->next->vertex, fx->next->next->vertex});
      }
    // stations drawn in each scrap: their height range
    std::vector<double> szmin(drawn.scraps.size(), 1e99), szmax(drawn.scraps.size(), -1e99);
    for (auto & ss : drawn.station_scraps)
      for (size_t si : ss.second) {
        szmin[si] = std::min(szmin[si], d1.station_vec[ss.first].z);
        szmax[si] = std::max(szmax[si], d1.station_vec[ss.first].z);
      }
    for (size_t si = 0; si < drawn.scraps.size(); si++) {
      if (szmin[si] > szmax[si]) continue;
      double x0 = 1e99, y0 = 1e99, x1 = -1e99, y1 = -1e99;
      for (const auto & e : drawn.scraps[si]) {
        x0 = std::min(x0, std::min(e[0], e[2])); x1 = std::max(x1, std::max(e[0], e[2]));
        y0 = std::min(y0, std::min(e[1], e[3])); y1 = std::max(y1, std::max(e[1], e[3]));
      }
      double cell = FILL_CELL;
      if ((x1 - x0) * (y1 - y0) / (cell * cell) > 4e6) cell = std::sqrt((x1 - x0) * (y1 - y0) / 4e6);
      int nx = int(ceil((x1 - x0) / cell)) + 1, ny = int(ceil((y1 - y0) / cell)) + 1;
      auto cx = [&](int i) { return x0 + (i + 0.5) * cell; };
      auto cy = [&](int j) { return y0 + (j + 0.5) * cell; };
      std::vector<char> drawnc(size_t(nx) * ny, 0), covered(size_t(nx) * ny, 0);
      std::vector<float> croof(size_t(nx) * ny, -1e30f), cfloor(size_t(nx) * ny, 1e30f);
      // drawn cells: scan lines
      for (int j = 0; j < ny; j++) {
        double y = cy(j);
        std::vector<double> xs;
        for (const auto & e : drawn.scraps[si])
          if ((e[1] > y) != (e[3] > y)) xs.push_back(e[0] + (y - e[1]) / (e[3] - e[1]) * (e[2] - e[0]));
        std::sort(xs.begin(), xs.end());
        for (size_t k = 0; k + 1 < xs.size(); k += 2)
          for (int i = std::max(0, int(floor((xs[k] - x0) / cell - 0.5)) + 1); (i < nx) && (cx(i) <= xs[k + 1]); i++)
            if (cx(i) >= xs[k]) drawnc[size_t(j) * nx + i] = 1;
      }
      // covered cells, with the model's roof and floor there
      for (const auto & t : mtris) {
        lxVec a = t[0]->get_vector(), b = t[1]->get_vector(), c = t[2]->get_vector();
        double zc = (a.z + b.z + c.z) / 3.0;
        if ((zc < szmin[si] - FILL_LEVEL) || (zc > szmax[si] + FILL_LEVEL)) continue;
        double den = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
        if (std::fabs(den) < 1e-12) continue;
        int i0 = std::max(0, int(floor((std::min(a.x, std::min(b.x, c.x)) - x0) / cell))),
          i1 = std::min(nx - 1, int(floor((std::max(a.x, std::max(b.x, c.x)) - x0) / cell))),
          j0 = std::max(0, int(floor((std::min(a.y, std::min(b.y, c.y)) - y0) / cell))),
          j1 = std::min(ny - 1, int(floor((std::max(a.y, std::max(b.y, c.y)) - y0) / cell)));
        for (int j = j0; j <= j1; j++)
          for (int i = i0; i <= i1; i++) {
            double px = cx(i), py = cy(j);
            double l1 = ((b.y - c.y) * (px - c.x) + (c.x - b.x) * (py - c.y)) / den,
              l2 = ((c.y - a.y) * (px - c.x) + (a.x - c.x) * (py - c.y)) / den, l3 = 1.0 - l1 - l2;
            if ((l1 < -1e-9) || (l2 < -1e-9) || (l3 < -1e-9)) continue;
            size_t q = size_t(j) * nx + i;
            float z = float(l1 * a.z + l2 * b.z + l3 * c.z);
            covered[q] = 1;
            croof[q] = std::max(croof[q], z);
            cfloor[q] = std::min(cfloor[q], z);
          }
      }
      // pieces of drawn plan not covered
      std::vector<int> comp(size_t(nx) * ny, -1);
      int ncomp = 0;
      for (size_t q0 = 0; q0 < comp.size(); q0++) {
        if (!drawnc[q0] || covered[q0] || (comp[q0] >= 0)) continue;
        std::vector<size_t> cells = {q0};
        comp[q0] = ncomp;
        for (size_t h = 0; h < cells.size(); h++) {
          int i = int(cells[h] % nx), j = int(cells[h] / nx);
          const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
          for (int k = 0; k < 4; k++) {
            int a2 = i + di[k], b2 = j + dj[k];
            if ((a2 < 0) || (b2 < 0) || (a2 >= nx) || (b2 >= ny)) continue;
            size_t q = size_t(b2) * nx + a2;
            if (drawnc[q] && !covered[q] && (comp[q] < 0)) { comp[q] = ncomp; cells.push_back(q); }
          }
        }
        ncomp++;
        if (int(cells.size()) < FILL_MIN_CELLS) continue;
        // heights: from the nearest model cell next to the piece (breadth first),
        // coming down to meet at the drawn wall
        std::map<size_t, size_t> src;  // cell -> covered cell giving its heights
        std::map<size_t, int> db, dw;
        std::vector<size_t> front, wfront;
        for (size_t q : cells) {
          int i = int(q % nx), j = int(q / nx);
          bool nearmodel = false, nearwall = false;
          size_t from = 0;
          for (int dj2 = -1; dj2 <= 1; dj2++)
            for (int di2 = -1; di2 <= 1; di2++) {
              int a2 = i + di2, b2 = j + dj2;
              if ((a2 < 0) || (b2 < 0) || (a2 >= nx) || (b2 >= ny)) { nearwall = true; continue; }
              size_t r = size_t(b2) * nx + a2;
              if (covered[r]) { nearmodel = true; from = r; }
              else if (!drawnc[r]) nearwall = true;
            }
          if (nearmodel) { src[q] = from; db[q] = 0; front.push_back(q); }
          if (nearwall) { dw[q] = 0; wfront.push_back(q); }
        }
        if (front.empty()) {
          fillreport += fmt::format("{:>20}  {:.1f} m2 of drawn plan away from the model, not filled\n",
            "", cells.size() * cell * cell);
          continue;
        }
        auto spread = [&](std::vector<size_t> fr, std::map<size_t, int> & dist, bool carry) {
          for (size_t h = 0; h < fr.size(); h++) {
            int i = int(fr[h] % nx), j = int(fr[h] / nx);
            const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
            for (int k = 0; k < 4; k++) {
              int a2 = i + di[k], b2 = j + dj[k];
              if ((a2 < 0) || (b2 < 0) || (a2 >= nx) || (b2 >= ny)) continue;
              size_t q = size_t(b2) * nx + a2;
              if ((comp[q] != comp[fr[h]]) || dist.count(q)) continue;
              dist[q] = dist[fr[h]] + 1;
              if (carry) src[q] = src[fr[h]];
              fr.push_back(q);
            }
          }
        };
        spread(front, db, true);
        spread(wfront, dw, false);
        // depth from the model and width of the opening onto it
        int dmax = 0;
        size_t qdeep = cells[0];
        for (auto & e : db) if (e.second > dmax) { dmax = e.second; qdeep = e.first; }
        double opening = 0.0;
        for (size_t k1 = 0; k1 < front.size(); k1++)
          for (size_t k2 = k1 + 1; k2 < front.size(); k2++)
            opening = std::max(opening, cell * std::hypot(double(int(front[k1] % nx) - int(front[k2] % nx)),
              double(int(front[k1] / nx) - int(front[k2] / nx))));
        opening += cell;
        double depth = (dmax + 1) * cell;
        // only passages drawn without legs or splays: deeper than the width of
        // their opening onto the model and at least FILL_PASSAGE_DEPTH deep;
        // slivers along the model and bays are left
        std::string near;
        double nearest = 1e99;
        {
          double px = cx(int(qdeep % nx)), py = cy(int(qdeep / nx));
          for (size_t k = 0; k < d1.station_vec.size(); k++) {
            thdb1ds & sk = d1.station_vec[k];
            if (sk.is_temporary() || (sk.z < szmin[si] - FILL_LEVEL) || (sk.z > szmax[si] + FILL_LEVEL)) continue;
            double dd = std::hypot(sk.x - px, sk.y - py);
            if (dd < nearest) { nearest = dd; near = fmt::format("{}@{}", sk.name, sk.survey->full_name); }
          }
        }
        if ((depth <= opening) || (depth < FILL_PASSAGE_DEPTH)) {
          nfillleft++;
          fillleft += cells.size() * cell * cell;
          continue;
        }
        // centre line: from the deepest cell back to the model, keeping as far
        // from the drawn walls as possible, then out to the middle of the opening
        std::vector<size_t> path = {qdeep};
        while (db[path.back()] > 0) {
          size_t qc = path.back(), next = qc;
          int i = int(qc % nx), j = int(qc / nx), bw = -1;
          const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
          for (int k = 0; k < 4; k++) {
            int a2 = i + di[k], b2 = j + dj[k];
            if ((a2 < 0) || (b2 < 0) || (a2 >= nx) || (b2 >= ny)) continue;
            size_t q = size_t(b2) * nx + a2;
            if ((comp[q] != comp[q0]) || !db.count(q) || (db[q] != db[qc] - 1)) continue;
            int w = dw.count(q) ? dw[q] : 0;
            if (w > bw) { bw = w; next = q; }
          }
          if (next == qc) break;
          path.push_back(next);
        }
        std::reverse(path.begin(), path.end());
        std::vector<lxVec> line;
        size_t r0 = src[path[0]];
        line.push_back(lxVec(0.5 * (cx(int(path[0] % nx)) + cx(int(r0 % nx))), 0.5 * (cy(int(path[0] / nx)) + cy(int(r0 / nx))), 0.0));
        for (size_t q : path) line.push_back(lxVec(cx(int(q % nx)), cy(int(q / nx)), 0.0));
        {
          std::vector<lxVec> sm = line;
          for (size_t k = 1; k + 1 < line.size(); k++) {
            lxVec a(0.0, 0.0, 0.0);
            int n = 0;
            for (size_t m = (k >= 2 ? k - 2 : 0); (m <= k + 2) && (m < line.size()); m++) { a += line[m]; n++; }
            sm[k] = (1.0 / n) * a;
          }
          line = sm;
        }
        // heights: the parent passage's roof and floor about the opening
        double top = -1e99, bot = 1e99, rad = 0.5 * opening + cell;
        for (int j = std::max(0, int(floor((line[0].y - rad - y0) / cell))); (j < ny) && (y0 + j * cell <= line[0].y + rad); j++)
          for (int i = std::max(0, int(floor((line[0].x - rad - x0) / cell))); (i < nx) && (x0 + i * cell <= line[0].x + rad); i++) {
            size_t q = size_t(j) * nx + i;
            if (!covered[q] || (std::hypot(cx(i) - line[0].x, cy(j) - line[0].y) > rad)) continue;
            top = std::max(top, double(croof[q]));
            bot = std::min(bot, double(cfloor[q]));
          }
        if (top < bot + 0.1) { top = croof[r0]; bot = cfloor[r0]; }
        if (top < bot + 0.1) { top = 0.5 * (top + bot) + 0.05; bot = top - 0.1; }
        double mid = 0.5 * (top + bot), up = top - mid, dn = mid - bot;
        // sections along the centre line: ovals out to the drawn walls (or to
        // where the piece meets the model), the parent's height
        std::vector<double> cum = {0.0};
        for (size_t k = 1; k < line.size(); k++) cum.push_back(cum.back() + std::hypot(line[k].x - line[k - 1].x, line[k].y - line[k - 1].y));
        double plen = cum.back();
        if (plen < cell) continue;
        auto at = [&](double t) {
          size_t k = 1;
          while ((k + 1 < cum.size()) && (cum[k] < t)) k++;
          double g = (cum[k] > cum[k - 1]) ? (t - cum[k - 1]) / (cum[k] - cum[k - 1]) : 0.0;
          g = std::max(0.0, std::min(1.0, g));
          return line[k - 1] + g * (line[k] - line[k - 1]);
        };
        auto dirat = [&](double t) {
          double e = std::min(0.3, 0.5 * plen);
          lxVec d = at(std::min(plen, t + e)) - at(std::max(0.0, t - e));
          d.z = 0.0;
          if (d.Length() < 1e-9) d = line.back() - line.front();
          d.Normalize();
          return d;
        };
        auto across = [&](const lxVec & p, const lxVec & d) {
          double t = 0.5 * cell;
          for (; t < 100.0; t += 0.5 * cell) {
            int i = int(floor((p.x + t * d.x - x0) / cell)), j = int(floor((p.y + t * d.y - y0) / cell));
            if ((i < 0) || (j < 0) || (i >= nx) || (j >= ny)) break;
            size_t q = size_t(j) * nx + i;
            if (comp[q] == comp[q0]) continue;
            if (!drawnc[q]) {
              double w = drawn.ray(si, p.x, p.y, d.x, d.y);
              if ((w > 0.0) && (w < t + cell)) return w;
            }
            break;
          }
          return t;
        };
        double step = std::max(thsplaytube_opt.spacing_min, std::min(thsplaytube_opt.spacing_max, 0.25 * plen));
        int nst = std::max(1, int(ceil((plen - cum[1]) / step)));
        // the sections; those near the end narrower than half the widest are
        // left for the rounded end
        struct psec { lxVec p, d, hx; double r, l; };
        std::vector<psec> ps;
        double wmax = 0.0;
        for (int k = 0; k <= nst; k++) {
          double t = cum[1] + (plen - cum[1]) * k / nst;
          lxVec p = at(t), d = dirat(t), hx(d.y, -d.x, 0.0);
          psec q = {p, d, hx, across(p, hx), across(p, -1.0 * hx)};
          ps.push_back(q);
          wmax = std::max(wmax, q.r + q.l);
        }
        while ((ps.size() > 1) && (ps.back().r + ps.back().l < 0.5 * wmax)) ps.pop_back();
        std::vector<thdb3dvx *> prev, first;
        std::vector<lxVec> Pfirst, Plast;
        lxVec dfirst, dlast, clast;
        for (const psec & q : ps) {
          lxVec c(q.p.x, q.p.y, mid);
          std::vector<lxVec> P;
          for (int m = 0; m < OVAL_POINTS; m++) {
            double a = 2.0 * M_PI * m / OVAL_POINTS, ca = cos(a), sa = sin(a);
            P.push_back(c + (ca * ((ca > 0.0) ? q.r : q.l)) * q.hx + lxVec(0.0, 0.0, sa * ((sa > 0.0) ? up : dn)));
          }
          std::vector<thdb3dvx *> cur;
          for (const lxVec & v : P) cur.push_back(map_fills.insert_vertex(v));
          if (prev.empty()) { first = cur; Pfirst = P; dfirst = q.d; }
          else loft_rings(map_fills.insert_face(THDB3DFC_TRIANGLES), prev, cur);
          prev = cur;
          Plast = P;
          dlast = q.d;
          clast = c;
        }
        // start: the first section run back into the parent passage and closed
        {
          double back = std::min(BRANCH_OVERLAP, 0.5 * opening) + cum[1];
          std::vector<thdb3dvx *> cur;
          for (const lxVec & v : Pfirst) cur.push_back(map_fills.insert_vertex(v - back * dfirst));
          loft_rings(map_fills.insert_face(THDB3DFC_TRIANGLES), cur, first);
          piece pc;
          thdb3dvx * cv = map_fills.insert_vertex(centroid(cur));
          for (size_t k = 0; k < cur.size(); k++) pc.add(cv, cur[k], cur[(k + 1) % cur.size()], -1.0 * dfirst);
          pc.emit(map_fills.insert_face(THDB3DFC_TRIANGLES));
        }
        // end: rounded off to the drawn end wall
        {
          double wall = drawn.ray(si, clast.x, clast.y, dlast.x, dlast.y);
          lxVec cc = centroid(prev);
          std::vector<thdb3dvx *> rr = prev;
          if (wall > 0.05) {
            for (int k = 1; k < 4; k++) {
              double t = k / 4.0, sc = std::sqrt(1.0 - t * t);
              std::vector<thdb3dvx *> nr;
              for (const lxVec & v : Plast) nr.push_back(map_fills.insert_vertex(cc + sc * (v - cc) + (t * wall) * dlast));
              loft_rings(map_fills.insert_face(THDB3DFC_TRIANGLES), rr, nr);
              rr = nr;
            }
            taper_ring(map_fills.insert_face(THDB3DFC_TRIANGLES), rr, map_fills.insert_vertex(cc + wall * dlast));
          } else {
            piece pc;
            thdb3dvx * cv = map_fills.insert_vertex(cc);
            for (size_t k = 0; k < rr.size(); k++) pc.add(cv, rr[(k + 1) % rr.size()], rr[k], dlast);
            pc.emit(map_fills.insert_face(THDB3DFC_TRIANGLES));
          }
        }
        fillreport += fmt::format("{:>30}  {:.1f} m away: {:.2f} m2, opening {:.2f} m, depth {:.2f} m, height {:.2f} m\n",
          near, nearest, cells.size() * cell * cell, opening, depth, top - bot);
        nfills++;
        nfillcells += cells.size();
      }
    }
  }

  if (usemap) {
    thlog(fmt::format("\nModel from maps and sections: {} sections stretched to the drawn walls ({} rings on {} legs;\n"
      "a section wherever the drawn width or height changes by more than {:.0f}%, {:.2f} to {:.2f} m apart).\n",
      nstretched, nrings, nlegspaced, thsplaytube_opt.spacing_change, thsplaytube_opt.spacing_min, thsplaytube_opt.spacing_max));
    thlog(fmt::format("Sections with (almost) nothing on one side, stretched as a whole across the drawn width: {}.\n", nwhole));
    thlog(fmt::format("Rings scaled to the extended elevation: {} ({} extended scraps).\n", nvscaled, elev.rings.size()));
    thlog(fmt::format("Legs read from the extended elevation as drawn (both stations in one scrap): {}; "
      "from the scraps as Therion places them: {}.\n", elevlegs_drawn.size(), elevlegs_placed.size()));
    thlog(fmt::format("Legs without sections drawn as ovals from the drawn LRUD: {} ({} rings).\n", novallegs, novalrings));
    thlog(fmt::format("Legs to drawn dead ends without a section, the section carried to the end wall: {}.\n", ndeadends));
    thlog(fmt::format("Passages drawn without legs or splays (deeper than their opening, at least {:.1f} m): {} ({:.1f} m2).\n{}"
      "Other drawn plan not covered by the model (slivers, bays), left: {} pieces, {:.1f} m2.\n",
      FILL_PASSAGE_DEPTH, nfills, nfillcells * FILL_CELL * FILL_CELL, fillreport, nfillleft, fillleft));
    if (!vreport.empty())
      thlog(fmt::format("Height scaled by more than {:.1f} times (or less than 1/{:.1f}):\n{}", STRETCH_REPORT, STRETCH_REPORT, vreport));
    thlog(fmt::format("Branches joined to junction footprints from the drawn walls: {}; run into the junction and closed: {}.\n",
      [&]() { size_t n = 0; for (auto & bj : branch_ends) for (auto & be : bj.second) if (be.joined) n++; return n; }(), noverlap));
    if (!stretchreport.empty())
      thlog(fmt::format("Stretched by more than {:.1f} times (or less than 1/{:.1f}):\n{}", STRETCH_REPORT, STRETCH_REPORT, stretchreport));
    if (!junctionreport.empty())
      thlog("Junctions:\n" + junctionreport);
  }
  if (!outside_report.empty())
    thlog(fmt::format("\nSplay ends more than {:.1f} m outside the drawn plan walls:\n{}", OUTSIDE_REPORT, outside_report));
  for (size_t i = 0; i < nst; i++)
    if ((nbrs[i].size() > 2) && (d1.station_vec[i].uid == i + 1) && d1.station_vec[i].survey->is_selected())
      thlog(fmt::format("{:>20} {:4}  junction\n",
        fmt::format("{}@{}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name), nbrs[i].size()));
  thlog(fmt::format("\n{} sections ({} fitted), {} tubes, {} tapered ends, {} tubes extended into "
    "junctions or stations without section, {} bodies from splays, {} section splays, "
    "{} stations with too few section splays, {} sections with anchors out of order, {} junctions, "
    "{} tube ends closed by triangulation ({} by a fan from the station), "
    "{} virtual splays to drawn walls ({} plan scraps)\n"
    "##################################################################\n",
    nsections, nfitted, ntubes, ntapers, nextended, nbodies, nused, nnosection, nunordered, njunctions, ncaps, nfans, nvirtual, drawn.scraps.size()));
}
