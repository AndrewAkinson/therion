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
#include "QuickHull.hpp"
#include <poly2tri/poly2tri.h>

#include <algorithm>
#include <functional>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

const double CANDIDATE_ANGLE = 45.0;  // max. angle of splays to the reference plane for plane fitting [deg]
const double SECTION_ANGLE = 20.0;  // max. angle of section splays to the fitted plane [deg]
const double MAX_PLANE_TURN = 35.0;  // max. angle between fitted and reference plane [deg]
const double STEEP = 0.7;  // |z| of the bisector above which the plane is not kept vertical
const double MAX_GAP = 150.0;  // larger angular gaps are closed at the station [deg]
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

// Angle between a direction and a plane given by its unit normal [deg].
double plane_angle(const lxVec & w, const lxVec & n) {
  double wl = std::hypot(std::hypot(w.x, w.y), w.z);
  return (wl > 0.0) ? deg(asin(std::min(1.0, fabs(w * n) / wl))) : 90.0;
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

// Closes a tube end: triangulates the section polygon (any shape) with the
// given vertices for its corners. Falls back to a fan from the station, which
// sees every corner of a polygon made of its splays.
bool cap(const section & s, const std::vector<thdb3dvx *> & ring, thdb3dvx * station, thdb3ddata & out) {
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
      for (p2t::Triangle * t : cdt.GetTriangles())
        for (int k = 0; k < 3; k++)
          fc->insert_vertex(vmap.at(t->GetPoint(k)));
      done = true;
    } catch (...) {
      done = false;
    }
  }
  if (!done) {
    size_t nv = ring.size();
    for (size_t k = 0; k < nv; k++) {
      fc->insert_vertex(station);
      fc->insert_vertex(ring[(k + 1) % nv]);
      fc->insert_vertex(ring[k]);
    }
  }
  return done;
}

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
  for (size_t k = 0; k < idx.size(); k++) {
    if (vx[idx[k]] == nullptr) vx[idx[k]] = out.insert_vertex(ends[idx[k]]);
    fc->insert_vertex(vx[idx[k]]);
  }
  return true;
}

}


void thsplaytube_build(thdatabase * db, thdb3ddata & walls, thdb3ddata & section_splays,
  thdb3ddata & polygons, thdb3ddata & junctions)
{
  thdb1d & d1 = db->db1d;
  thdb1d_tree_node * nodes = d1.get_tree_nodes();
  size_t nst = d1.station_vec.size();
  if ((nodes == NULL) || (nst == 0)) return;

  const lxVec up(0.0, 0.0, 1.0);
  std::vector<section> sec(nst);
  std::vector<std::vector<size_t>> nbrs(nst);
  size_t nsections = 0, nnosection = 0, njunctions = 0, nused = 0, nfitted = 0, nunordered = 0;
  std::string report;

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
  for (size_t i = 0; i < nst; i++) {
    thdb1ds & st = d1.station_vec[i];
    if ((st.uid != i + 1) || st.is_temporary() || !st.survey->is_selected()) continue;
    if ((nbrs[i].size() < 1) || (nbrs[i].size() > 2)) {
      if (nbrs[i].size() > 2) njunctions++;
      continue;
    }
    section & s = sec[i];
    s.c = lxVec(st.x, st.y, st.z);
    thdb1ds & sa = d1.station_vec[nbrs[i][0]];
    lxVec din = s.c - lxVec(sa.x, sa.y, sa.z);
    din.Normalize();
    lxVec nref = din;
    if (nbrs[i].size() == 2) {
      thdb1ds & sb = d1.station_vec[nbrs[i][1]];
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
        if (plane_angle(w, s.n) <= CANDIDATE_ANGLE) cand.push_back(w);
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
      double pa = plane_angle(w, s.n);
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
      hypot(s.cx, s.cy), std::string(plane) + anchors + s.moved);

    thdb3dfc * pfc = polygons.insert_face(THDB3DFC_LINE_STRIP);
    np = s.q.size();  // anchors may have been listed twice
    for (size_t k = 0; k <= np; k++)
      pfc->insert_vertex(polygons.insert_vertex(s.point(k % np)));
    for (size_t k = 0; k < np; k++)
      s.vx.push_back(walls.insert_vertex(s.point(k)));
    s.ok = true;
    nsections++;
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

  // tubes along legs between stations with sections
  size_t ntubes = 0, ntapers = 0, nextended = 0, ncaps = 0, nfans = 0;
  std::string flipreport, tubereport;
  for (size_t i = 0; i < nst; i++) {
    if (!sec[i].ok) continue;
    for (size_t j : nbrs[i]) {
      // taper to a passage end without section
      if (!sec[j].ok && (nbrs[j].size() == 1)) {
        thdb1ds & e = d1.station_vec[j];
        thdb3dvx * ev = walls.insert_vertex(e.x, e.y, e.z);
        thdb3dfc * fc = walls.insert_face(THDB3DFC_TRIANGLES);
        size_t nv = sec[i].vx.size();
        for (size_t k = 0; k < nv; k++) {
          fc->insert_vertex(sec[i].vx[k]);
          fc->insert_vertex(sec[i].vx[(k + 1) % nv]);
          fc->insert_vertex(ev);
        }
        sec[i].ntubes++;
        ntapers++;
        continue;
      }
      // extend into a junction (or a station without section): the section
      // is carried unchanged to that station, where the junction body covers the join
      if (body[j]) {
        section & a = sec[i];
        thdb1ds & e = d1.station_vec[j];
        lxVec shift = lxVec(e.x, e.y, e.z) - a.c;
        size_t nv = a.vx.size();
        std::vector<thdb3dvx *> rb;
        for (size_t k = 0; k < nv; k++)
          rb.push_back(walls.insert_vertex(a.point(k) + shift));
        thdb3dfc * fc = walls.insert_face(THDB3DFC_TRIANGLES);
        for (size_t k = 0; k < nv; k++) {
          size_t k1 = (k + 1) % nv;
          fc->insert_vertex(a.vx[k]);
          fc->insert_vertex(a.vx[k1]);
          fc->insert_vertex(rb[k1]);
          fc->insert_vertex(a.vx[k]);
          fc->insert_vertex(rb[k1]);
          fc->insert_vertex(rb[k]);
        }
        // close the end inside the junction
        if (cap(a, rb, walls.insert_vertex(a.c + shift), walls)) ncaps++; else nfans++;
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
      thdb3dfc * fc = walls.insert_face(THDB3DFC_TRIANGLES);
      size_t nv = ra.size();
      for (size_t k = 0; k < nv; k++) {
        size_t k1 = (k + 1) % nv;
        fc->insert_vertex(ra[k]);
        fc->insert_vertex(ra[k1]);
        fc->insert_vertex(rb[k1]);
        fc->insert_vertex(ra[k]);
        fc->insert_vertex(rb[k1]);
        fc->insert_vertex(rb[k]);
      }
      a.ntubes++;
      b.ntubes++;
      ntubes++;
    }
  }

  // close open tube ends
  for (size_t i = 0; i < nst; i++) {
    section & s = sec[i];
    if (!s.ok || (s.ntubes >= 2)) continue;
    if (cap(s, s.vx, walls.insert_vertex(s.c), walls)) ncaps++; else nfans++;
  }

  thlog(fmt::format("\n\n########################### splay tubes ###########################\n"
    "Section plane at each station: leg bisector (vertical unless steep), turned to the best\n"
    "fit of splays within {:.0f} deg of it; section splays are within {:.0f} deg of that plane.\n"
    "turn = angle between fitted plane and bisector, spread = RMS angle of the section splays\n"
    "to the plane, worst = largest of these (all in degrees); centre = distance from the\n"
    "station to the centre of the section polygon (m), around which tubes join sections.\n\n"
    "{:>20} {:>4} {:>6} {:>4}  {:>5} {:>6} {:>5} {:>6}  {}\n{}",
    CANDIDATE_ANGLE, SECTION_ANGLE, "STATION", "LEGS", "SPLAYS", "USED", "TURN", "SPREAD", "WORST", "CENTRE", "PLANE",
    report));
  thlog(fmt::format("\nTube joins: corners of neighbouring sections joined directly when their\n"
    "positions (0-400 round the section) differ by at most {:.0f}:\n{}", PAIR_TOLERANCE, tubereport));
  if (!flipreport.empty())
    thlog("\nTubes where left and right are matched differently by leg direction and by normals:\n"
      + flipreport + "\n");
  for (size_t i = 0; i < nst; i++)
    if ((nbrs[i].size() > 2) && (d1.station_vec[i].uid == i + 1) && d1.station_vec[i].survey->is_selected())
      thlog(fmt::format("{:>20} {:4}  junction\n",
        fmt::format("{}@{}", d1.station_vec[i].name, d1.station_vec[i].survey->full_name), nbrs[i].size()));
  thlog(fmt::format("\n{} sections ({} fitted), {} tubes, {} tapered ends, {} tubes extended into "
    "junctions or stations without section, {} bodies from splays, {} section splays, "
    "{} stations with too few section splays, {} sections with anchors out of order, {} junctions, "
    "{} tube ends closed by triangulation ({} by a fan from the station)\n"
    "##################################################################\n",
    nsections, nfitted, ntubes, ntapers, nextended, nbodies, nused, nnosection, nunordered, njunctions, ncaps, nfans));
}
