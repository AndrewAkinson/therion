# Plan for model wall sources and one wall list for every format

Branch `LochModel`. Written 2026-10-15 for work done while the user is away, so
every decision is recorded here with its reason. This file is a working plan and
is removed before any pull request upstream.

## Aim

The 3D model from maps and sections (built in `thsplaytube.cxx`) replaces the old
walls from maps (each plan scrap's outline made 3D). Every model format that has
walls should carry the same walls, so that Loch, VRML, 3DMF, DXF, glTF and HTML
differ only in how the file is written.

## Where things stand

`-wall-source` takes one value of `all`, `centreline`, `maps`, `scans`, `splays`
and `splay-tubes`, with `all` the default. Each writer has its own copy of the
code choosing walls, and they have drifted apart.

| Format | Walls written now |
|---|---|
| glTF, HTML | LRUD, splays, splay tubes and junction bodies, the model from maps and sections, working layers (old map walls removed in 9869b395) |
| VRML, 3DMF, DXF | LRUD, old map walls, splay tubes and junction bodies |
| Loch | LRUD, old map walls, scans, splays, splay tubes and junction bodies |
| therion (`.thm`) | LRUD, old map walls |
| Survex, Compass, KML, Esri | no walls |

Loch has four copies of the code turning a mesh into its triangles. The glTF and
HTML working layers were uncommitted test code until now.

## Decisions

1. **`maps` builds the model from maps and sections.** It is the walls from maps
   and sections, the junctions from maps and the passages drawn in the plan. The
   old outlines of plan scraps are no longer exported in any format. The user
   asked for the old walls from maps to be replaced, and keeping both under one
   name would put two overlapping walls in every model.
2. **`maps` falls back** when there is nothing drawn. With no selected plan scraps
   the model from maps is empty, so `maps` gives the splay tubes and junction
   bodies instead, and with no cross section splays either it gives the LRUD walls.
   A cave with no drawings still gets walls, as it did before.
3. **`splay-tubes`** gives the tubes from cross section splays and their junction
   bodies, ignoring the drawings. It is the quick check model before walls are
   drawn.
4. **`centreline`, `splays` and `scans`** keep their meaning.
5. **Several sources** may be given as one comma separated value, for example
   `-wall-source maps,centreline`. A single value works as before. Commas are used
   because the option takes one argument and `<a/b>` lists in the manual are
   already read as alternatives.
6. **`all`** still means every source (LRUD, model from maps, scans, splays and
   splay tubes), overlapping, as before.
7. **The default becomes `maps`**, not `all`. One clean wall is what most users
   want, and decision 2 means it never gives less than before. This changes the
   output of existing projects that do not set `-wall-source`, so it goes in
   CHANGES.
8. **`-model-debug on/off`** (default `off`) adds the working layers, which are the
   section splays, section polygons, virtual splays and drawn cross sections, and
   with `maps` the splay tubes and junction bodies underneath it. They are mostly
   lines, which only glTF and HTML can show beside walls, so the option affects
   those two formats only. A name of its own keeps it apart from `-wall-source`,
   which chooses what the walls are made from.
9. **One function lists the wall layers** for the chosen sources (mesh, name,
   group, colour, survey). Every writer loops over that list and keeps only its
   own file writing. Loch keeps one mesh per layer, attached to the top exported
   survey as the splay tubes are now (scraps, scans and station splays keep their
   own survey). The four Loch copies of the triangle code become one helper.
10. **Survex, Compass, KML and Esri** stay without walls. Their formats carry only
    the centreline here, and adding walls to them is a separate job.
11. **therion (`.thm`) format** uses the list too. It is read by xtherion's old 3D
    viewer, and keeping it in step costs nothing once the list exists.
12. **Documentation** follows CONTRIBUTING.md, a summary in CHANGES and the details
    in thbook (`ch03.tex`, model export options), including the undocumented
    `-section-spacing` option added earlier on this branch.
13. **Tests.** The parsing of `-wall-source` values becomes a small function with a
    unit test in `tests/`. Each stage is checked by building, running `utest`, and
    exporting Trident, FairyChimneys and LVS (in ~/Projects/TherionModels, local
    data, not cited in commits) in every model format, checking each file is
    written and the glTF meshes are the expected layers.
14. **Commits** are made after each stage and pushed to `LochModel` at the end.

## Stages

1. **Wall sources.** Decisions 1 to 7. The wall source value parser and its unit
   test, the shared wall layer list, glTF and HTML using it, the HTML layer names,
   CHANGES and thbook.
2. **Every format from the list.** Decisions 9 to 11. VRML, 3DMF, DXF, Loch and
   therion writers loop over the list; the Loch triangle helper.
3. **Debug layers.** Decision 8. `-model-debug`, the working layers in glTF and
   HTML (replacing the uncommitted test code), CHANGES and thbook.
4. **Wrap up.** Check every format on the test caves, update this plan with what
   was done and anything left, push.

## Progress

(Filled in as the stages are done.)
