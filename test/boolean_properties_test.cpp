// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "gtest/gtest.h"
#include "manifold/common.h"
#include "manifold/manifold.h"

using namespace manifold;

namespace {

// ---------------------------------------------------------------------------
// Local test material: a cube whose every edge is a property seam.
//
// Geometrically the 8 corners are shared, but each face owns its own four
// property vertices (24 in total), joined together with mergeFromVert /
// mergeToVert pairs. Two faces meeting at an edge therefore carry different
// values on either side of the seam, exactly like a texture or normal seam in
// an exported game asset.
//
// Seven extra property channels per vertex:
//   channel 0: per-face material tag (constant across the face)
//   channel 1..2: UV-like coordinates, an affine function of position with a
//              per-face offset. The offset makes them discontinuous across
//              every edge (a genuine texture seam), while being affine inside
//              each face so Boolean re-triangulation preserves them exactly.
//   channel 3: position-derived scalar, an affine function of position, whose
//              barycentric value is known independently at every point
//   channel 4..6: the geometric face normal, giving a normal seam on every
//              cube edge
// ---------------------------------------------------------------------------

constexpr int kNumProp = 7;

struct SeamedCube {
  MeshGL64 mesh;
  std::vector<vec3> corners;  // 8 geometric corner positions
};

SeamedCube MakeSeamedCube(vec3 origin, double tagBase, double uvBase) {
  SeamedCube out;
  const vec3 o = origin;
  out.corners = {
      o + vec3(-1, -1, -1),  // 0
      o + vec3(1, -1, -1),   // 1
      o + vec3(1, 1, -1),    // 2
      o + vec3(-1, 1, -1),   // 3
      o + vec3(-1, -1, 1),   // 4
      o + vec3(1, -1, 1),    // 5
      o + vec3(1, 1, 1),     // 6
      o + vec3(-1, 1, 1),    // 7
  };
  // CCW corner orders seen from outside; order also defines the normal.
  const int faces[6][4] = {
      {4, 5, 6, 7},  // +z
      {1, 0, 3, 2},  // -z
      {1, 2, 6, 5},  // +x
      {0, 4, 7, 3},  // -x
      {3, 7, 6, 2},  // +y
      {0, 1, 5, 4},  // -y
  };

  MeshGL64& m = out.mesh;
  m.numProp = 3 + kNumProp;
  m.runOriginalID.push_back(Manifold::ReserveIDs(1));

  // First property index created for each geometric corner becomes the
  // mergeToVert target; later duplicates merge into it.
  std::map<int, int> cornerToFirst;
  for (int face = 0; face < 6; ++face) {
    vec3 normal;
    {
      const vec3 a = out.corners[faces[face][0]];
      const vec3 b = out.corners[faces[face][1]];
      const vec3 c = out.corners[faces[face][2]];
      normal = la::normalize(la::cross(b - a, c - a));
    }
    const int triStart = static_cast<int>(m.NumTri());
    // Two triangles per face: (0,1,2) and (0,2,3).
    m.triVerts.insert(m.triVerts.end(),
                      {0, 0, 0, 0, 0, 0});  // indices patched below
    m.faceID.push_back(2 * face);
    m.faceID.push_back(2 * face + 1);
    for (int corner = 0; corner < 4; ++corner) {
      const int cornerId = faces[face][corner];
      const vec3 p = out.corners[cornerId];
      const uint64_t propIdx = m.NumVert();
      // xyz position
      for (int k : {0, 1, 2}) m.vertProperties.push_back(p[k]);
      // material tag
      m.vertProperties.push_back(tagBase + face);
      // UV-like seam coordinates: per-face offset plus a position-derived,
      // in-face affine term. Different faces meeting at an edge use different
      // offsets, so coincident geometric vertices stay distinct properties.
      m.vertProperties.push_back(uvBase + 10.0 * face + 0.3 * p.x +
                                 0.7 * p.y);
      m.vertProperties.push_back(uvBase + 10.0 * face + 0.5 * p.z -
                                 0.2 * p.x);
      // affine position scalar
      m.vertProperties.push_back(0.25 * p.x - 0.5 * p.y + 0.75 * p.z + 3.0);
      // face normal (seams on every edge)
      for (int k : {0, 1, 2}) m.vertProperties.push_back(normal[k]);

      auto it = cornerToFirst.find(cornerId);
      if (it == cornerToFirst.end()) {
        cornerToFirst[cornerId] = static_cast<int>(propIdx);
      } else {
        m.mergeFromVert.push_back(static_cast<uint64_t>(propIdx));
        m.mergeToVert.push_back(static_cast<uint64_t>(it->second));
      }
    }
    const int base = static_cast<int>(m.NumVert()) - 4;
    m.triVerts[3 * triStart + 0] = base;
    m.triVerts[3 * triStart + 1] = base + 1;
    m.triVerts[3 * triStart + 2] = base + 2;
    m.triVerts[3 * triStart + 3] = base;
    m.triVerts[3 * triStart + 4] = base + 2;
    m.triVerts[3 * triStart + 5] = base + 3;
  }
  m.runIndex = {0, static_cast<uint64_t>(m.triVerts.size())};
  return out;
}

double AffineScalar(vec3 p) {
  return 0.25 * p.x - 0.5 * p.y + 0.75 * p.z + 3.0;
}

// Every undirected geometric edge of a closed, manifold triangle soup has
// exactly two directed halfedges. Property seams intentionally duplicate
// vertices (same position, different properties), so the pairing must be
// checked on geometric positions rather than raw property indices. This does
// not depend on triangle ordering.
void ExpectClosedTopology(const MeshGL64& m) {
  ASSERT_EQ(m.triVerts.size() % 3, 0u);
  auto posKey = [&](uint64_t v) {
    std::array<int64_t, 3> key{};
    for (int k : {0, 1, 2}) {
      key[k] = std::llround(m.vertProperties[v * m.numProp + k] * 1e9);
    }
    return key;
  };
  std::map<std::array<int64_t, 3>, uint64_t> canonical;
  std::vector<uint64_t> geomVert(m.NumVert());
  for (uint64_t v = 0; v < m.NumVert(); ++v) {
    const auto key = posKey(v);
    auto it = canonical.find(key);
    if (it == canonical.end()) {
      geomVert[v] = canonical.size();
      canonical[key] = geomVert[v];
    } else {
      geomVert[v] = it->second;
    }
  }
  std::map<std::pair<uint64_t, uint64_t>, int> halfedges;
  for (size_t t = 0; t < m.NumTri(); ++t) {
    for (int i = 0; i < 3; ++i) {
      const uint64_t a = geomVert[m.triVerts[3 * t + i]];
      const uint64_t b = geomVert[m.triVerts[3 * t + ((i + 1) % 3)]];
      ASSERT_NE(a, b);
      ++halfedges[{a, b}];
    }
  }
  for (const auto& edge : halfedges) {
    auto reverse = halfedges.find({edge.first.second, edge.first.first});
    ASSERT_NE(reverse, halfedges.end())
        << "unpaired halfedge - mesh is not closed";
    EXPECT_EQ(edge.second, 1);
    EXPECT_EQ(reverse->second, 1);
  }
}

// The merge vectors of the output encode the property seams: every merge pair
// shares a position but carries a different property tuple.
void ExpectMergeVectorsEncodeSeams(const MeshGL64& m) {
  ASSERT_EQ(m.mergeFromVert.size(), m.mergeToVert.size());
  for (size_t i = 0; i < m.mergeFromVert.size(); ++i) {
    const uint64_t a = m.mergeFromVert[i];
    const uint64_t b = m.mergeToVert[i];
    ASSERT_LT(a, m.NumVert());
    ASSERT_LT(b, m.NumVert());
    for (int k : {0, 1, 2}) {
      EXPECT_NEAR(m.vertProperties[a * m.numProp + k],
                  m.vertProperties[b * m.numProp + k], 1e-12)
          << "merge pair must share a geometric position";
    }
    bool propsDiffer = false;
    for (size_t k = 3; k < m.numProp; ++k) {
      if (m.vertProperties[a * m.numProp + k] !=
          m.vertProperties[b * m.numProp + k]) {
        propsDiffer = true;
        break;
      }
    }
    EXPECT_TRUE(propsDiffer)
        << "a merge pair with identical properties means a seam was welded "
           "away by geometric welding";
  }
}

// Group property vertices by geometric position. A seam position shows up as a
// bucket with several distinct property tuples.
std::map<std::array<uint64_t, 3>, std::vector<uint64_t>> PositionBuckets(
    const MeshGL64& m, double tol = 1e-9) {
  std::map<std::array<uint64_t, 3>, std::vector<uint64_t>> buckets;
  for (uint64_t v = 0; v < m.NumVert(); ++v) {
    std::array<uint64_t, 3> key{};
    for (int k : {0, 1, 2}) {
      key[k] = static_cast<uint64_t>(
          std::llround(m.vertProperties[v * m.numProp + k] / tol));
    }
    buckets[key].push_back(v);
  }
  return buckets;
}

int CountSeamPositions(const MeshGL64& m) {
  int seams = 0;
  for (const auto& bucket : PositionBuckets(m)) {
    const auto& verts = bucket.second;
    bool distinct = false;
    for (size_t i = 1; i < verts.size(); ++i) {
      for (size_t k = 3; k < m.numProp; ++k) {
        if (m.vertProperties[verts[i] * m.numProp + k] !=
            m.vertProperties[verts[0] * m.numProp + k]) {
          distinct = true;
          break;
        }
      }
      if (distinct) break;
    }
    if (distinct) ++seams;
  }
  return seams;
}

vec3 TriPosition(const MeshGL64& m, uint64_t tri, int corner) {
  const uint64_t v = m.triVerts[3 * tri + corner];
  return m.GetVertPos(v);
}

std::array<double, kNumProp> TriProps(const MeshGL64& m, uint64_t tri,
                                      int corner) {
  const uint64_t v = m.triVerts[3 * tri + corner];
  std::array<double, kNumProp> p{};
  for (int k = 0; k < kNumProp; ++k)
    p[k] = m.vertProperties[v * m.numProp + 3 + k];
  return p;
}

// Barycentric coordinates of point p within triangle (a, b, c) by signed
// sub-triangle areas. Returns weights summing to one (including for points
// lying on an edge of the source triangle).
std::array<double, 3> Barycentric(vec3 p, vec3 a, vec3 b, vec3 c) {
  const vec3 n = la::cross(b - a, c - a);
  const double area2 = la::dot(n, n);
  std::array<double, 3> w{};
  // Weight at each corner equals the normalized area of the opposite
  // sub-triangle; all three are measured with the same orientation so their
  // signed sum is the full triangle area.
  w[0] = la::dot(la::cross(b - p, c - p), n) / area2;
  w[1] = la::dot(la::cross(c - p, a - p), n) / area2;
  w[2] = la::dot(la::cross(a - p, b - p), n) / area2;
  return w;
}

// Distance from point p to the plane through (a, b, c).
double PlaneDistance(vec3 p, vec3 a, vec3 b, vec3 c) {
  const vec3 normal = la::normalize(la::cross(b - a, c - a));
  return std::abs(la::dot(p - a, normal));
}

// For every triangle of `out`, locate its source run's input mesh
// geometrically: the source triangle is the unique one whose plane contains
// the output triangle's three vertices. That mapping is independent of any
// triangle reordering (faceID / triID bookkeeping). Then verify:
//  - every output vertex is an affine combination of the source corners,
//  - every extra property equals the barycentric interpolation of the source
//    triangle's corner properties.
// A wrong interpolation source, or a swapped corner correspondence, makes
// the property comparison fail even when the geometry alone lines up.
void ExpectBarycentricFromSources(
    const MeshGL64& out,
    const std::map<uint32_t, const MeshGL64*>& originals) {
  constexpr double kPlaneTol = 1e-7;
  for (size_t run = 0; run < out.NumRun(); ++run) {
    auto it = originals.find(out.runOriginalID[run]);
    ASSERT_NE(it, originals.end()) << "run with unknown OriginalID";
    const MeshGL64& in = *it->second;
    for (uint64_t tri = out.runIndex[run] / 3;
         tri < out.runIndex[run + 1] / 3; ++tri) {
      // The output triangle's own orientation selects the source face: an
      // intersection vertex can lie on an edge shared by two perpendicular
      // source faces, but the other two vertices of the output triangle pin
      // its plane to exactly one of them.
      const vec3 outNormal = la::normalize(
          la::cross(TriPosition(out, tri, 1) - TriPosition(out, tri, 0),
                    TriPosition(out, tri, 2) - TriPosition(out, tri, 0)));
      std::vector<uint64_t> candidates;
      for (uint64_t srcTri = 0; srcTri < in.NumTri(); ++srcTri) {
        const vec3 a = TriPosition(in, srcTri, 0);
        const vec3 b = TriPosition(in, srcTri, 1);
        const vec3 c = TriPosition(in, srcTri, 2);
        const vec3 srcNormal = la::normalize(la::cross(b - a, c - a));
        // Parallel (a reversed difference face is anti-parallel, hence the
        // absolute dot).
        if (std::abs(la::dot(srcNormal, outNormal)) < 0.9) continue;
        bool allOnPlane = true;
        for (int corner = 0; corner < 3; ++corner) {
          if (PlaneDistance(TriPosition(out, tri, corner), a, b, c) >
              kPlaneTol) {
            allOnPlane = false;
            break;
          }
        }
        if (allOnPlane) candidates.push_back(srcTri);
      }
      ASSERT_FALSE(candidates.empty())
          << "no source triangle plane contains an output face";
      // The two triangles of a coplanar quad share a plane; either is a valid
      // interpolation source because the properties are continuous inside one
      // input face. Every candidate must nevertheless reproduce the output
      // properties, so check each.
      for (const uint64_t srcTri : candidates) {
        const vec3 a = TriPosition(in, srcTri, 0);
        const vec3 b = TriPosition(in, srcTri, 1);
        const vec3 c = TriPosition(in, srcTri, 2);
        const auto pa = TriProps(in, srcTri, 0);
        const auto pb = TriProps(in, srcTri, 1);
        const auto pc = TriProps(in, srcTri, 2);
        for (int corner = 0; corner < 3; ++corner) {
          const vec3 p = TriPosition(out, tri, corner);
          const auto got = TriProps(out, tri, corner);
          const auto w = Barycentric(p, a, b, c);
          EXPECT_NEAR(w[0] + w[1] + w[2], 1.0, 1e-9);
          for (int k = 0; k < kNumProp; ++k) {
            const double expected =
                w[0] * pa[k] + w[1] * pb[k] + w[2] * pc[k];
            EXPECT_NEAR(got[k], expected, 1e-7)
                << "property channel " << k
                << " not barycentric in its source surface - wrong "
                   "interpolation source or flipped corner correspondence";
          }
        }
      }
    }
  }
}

// Material tags live in a disjoint value range per input, so they name the
// source surface unambiguously.
void ExpectTagsBelongToRunSource(
    const MeshGL64& out,
    const std::map<uint32_t, std::pair<double, double>>& tagRanges) {
  for (size_t run = 0; run < out.NumRun(); ++run) {
    auto range = tagRanges.find(out.runOriginalID[run]);
    ASSERT_NE(range, tagRanges.end());
    for (uint64_t tri = out.runIndex[run] / 3;
         tri < out.runIndex[run + 1] / 3; ++tri) {
      for (int corner = 0; corner < 3; ++corner) {
        const double tag = TriProps(out, tri, corner)[0];
        EXPECT_GE(tag, range->second.first - 1e-12);
        EXPECT_LE(tag, range->second.second + 1e-12);
      }
    }
  }
}

// Build the two overlapping seamed cubes used throughout:
//   A = [-1,1]^3, tags [0,5],   B = [0,2]^3, tags [10,15].
// The cut x = 0 / x = 1 lands on edge interiors, so every intersection
// vertex is newly created and must be barycentrically interpolated.
struct SeamedPair {
  SeamedCube a;
  SeamedCube b;
  Manifold ma;
  Manifold mb;
  uint32_t idA;
  uint32_t idB;
  std::map<uint32_t, const MeshGL64*> originals;
  std::map<uint32_t, std::pair<double, double>> tagRanges;
};

SeamedPair MakePair() {
  SeamedPair pair;
  pair.a = MakeSeamedCube(vec3(0.0), 0.0, 0.0);
  pair.b = MakeSeamedCube(vec3(1, 1, 1), 10.0, 100.0);
  pair.ma = Manifold(pair.a.mesh).AsOriginal();
  pair.mb = Manifold(pair.b.mesh).AsOriginal();
  pair.idA = static_cast<uint32_t>(pair.ma.OriginalID());
  pair.idB = static_cast<uint32_t>(pair.mb.OriginalID());
  pair.originals = {{pair.idA, &pair.a.mesh}, {pair.idB, &pair.b.mesh}};
  pair.tagRanges = {{pair.idA, {0.0, 5.0}}, {pair.idB, {10.0, 15.0}}};
  return pair;
}

}  // namespace

// ---------------------------------------------------------------------------
// Property seams survive geometric welding.
//
// The input cubes intentionally share geometric vertices that carry different
// properties on either side of every edge. A welding implementation that
// merges by position alone would erase those seams; the merge vectors must
// instead keep distinct property vertices for coincident positions.
// ---------------------------------------------------------------------------
TEST(BooleanProperties, SeamsSurviveInputWelding) {
  SeamedCube cube = MakeSeamedCube(vec3(0.0), 0.0, 0.0);
  const Manifold m(cube.mesh);
  ASSERT_EQ(m.Status(), Manifold::Error::NoError);
  // 24 property verts collapse to 8 geometric verts through the seam merges.
  EXPECT_EQ(m.NumVert(), 8u);
  EXPECT_EQ(m.NumPropVert(), 24u);

  const MeshGL64 out = m.GetMeshGL64();
  ExpectClosedTopology(out);
  ExpectMergeVectorsEncodeSeams(out);
  // All 8 geometric corners are seam positions (three faces meet there).
  EXPECT_EQ(CountSeamPositions(out), 8);
}

// ---------------------------------------------------------------------------
// All three Boolean operations keep the interpolation tied to the surface
// that actually contributed each face. Intersection vertices come from both
// input surfaces, so a run that interpolates from the wrong source triangle
// (or swaps two corners) cannot satisfy the barycentric check.
// ---------------------------------------------------------------------------
TEST(BooleanProperties, SourcesAndBarycentricInterpolation) {
  SeamedPair p = MakePair();

  for (const auto op :
       {std::pair<OpType, const char*>{OpType::Add, "union"},
        {OpType::Intersect, "intersect"},
        {OpType::Subtract, "subtract"}}) {
    SCOPED_TRACE(op.second);
    const Manifold result = p.ma.Boolean(p.mb, op.first);
    ASSERT_EQ(result.Status(), Manifold::Error::NoError);
    ASSERT_FALSE(result.IsEmpty());
    ASSERT_TRUE(result.MatchesTriNormals());

    const MeshGL64 out = result.GetMeshGL64();
    ExpectClosedTopology(out);
    ExpectMergeVectorsEncodeSeams(out);

    // Both input surfaces contribute retained faces for every operation.
    std::set<uint32_t> ids(out.runOriginalID.begin(),
                           out.runOriginalID.end());
    EXPECT_TRUE(ids.count(p.idA));
    EXPECT_TRUE(ids.count(p.idB));

    ExpectTagsBelongToRunSource(out, p.tagRanges);
    ExpectBarycentricFromSources(out, p.originals);

    // The affine position-derived scalar gives an independent, source-free
    // reference value at every vertex.
    for (uint64_t v = 0; v < out.NumVert(); ++v) {
      const vec3 pos = out.GetVertPos(v);
      EXPECT_NEAR(out.vertProperties[v * out.numProp + 3 + 3],
                  AffineScalar(pos), 1e-8);
    }
  }
}

// ---------------------------------------------------------------------------
// Difference semantics: the cut faces closing the hole are the tool's
// surface with reversed orientation, while the untouched region keeps the
// base surface's properties verbatim. The Backside run flag distinguishes the
// two, and either kind of face must still interpolate its properties from the
// surface it geometrically came from.
// ---------------------------------------------------------------------------
TEST(BooleanProperties, DifferenceCutFacesVersusUnchangedRegion) {
  SeamedPair p = MakePair();
  const Manifold diff = p.ma - p.mb;
  ASSERT_EQ(diff.Status(), Manifold::Error::NoError);
  const MeshGL64 out = diff.GetMeshGL64();
  ExpectClosedTopology(out);
  // A notched cube stays a genus-0 solid.
  EXPECT_EQ(diff.Genus(), 0);

  bool sawCut = false;
  bool sawUnchanged = false;
  for (size_t run = 0; run < out.NumRun(); ++run) {
    const uint32_t id = out.runOriginalID[run];
    const bool backside = out.Backside(run);
    const uint64_t triBegin = out.runIndex[run] / 3;
    const uint64_t triEnd = out.runIndex[run + 1] / 3;
    for (uint64_t tri = triBegin; tri < triEnd; ++tri) {
      // Every vertex of the difference result lies inside the base cube A.
      for (int corner = 0; corner < 3; ++corner) {
        const vec3 pos = TriPosition(out, tri, corner);
        for (int k : {0, 1, 2}) {
          EXPECT_GE(pos[k], -1.0 - 1e-9);
          EXPECT_LE(pos[k], 1.0 + 1e-9);
        }
      }
      if (id == p.idB) {
        // Tool surface only appears as reversed cut faces.
        EXPECT_TRUE(backside) << "tool faces of a difference must be backside";
        sawCut = true;
        // A cut face sits on a boundary plane of the tool cube.
        for (int corner = 0; corner < 3; ++corner) {
          const vec3 pos = TriPosition(out, tri, corner);
          const bool onToolBoundary =
              std::abs(pos.x - 0.0) < 1e-9 ||
              std::abs(pos.x - 2.0) < 1e-9 ||
              std::abs(pos.y - 0.0) < 1e-9 ||
              std::abs(pos.y - 2.0) < 1e-9 ||
              std::abs(pos.z - 0.0) < 1e-9 ||
              std::abs(pos.z - 2.0) < 1e-9;
          EXPECT_TRUE(onToolBoundary);
        }
      } else {
        EXPECT_EQ(id, p.idA);
        sawUnchanged = true;
        // Base-surface faces keep their original orientation.
        EXPECT_FALSE(backside);
        // The untouched region retains the base material tag exactly.
        for (int corner = 0; corner < 3; ++corner) {
          const double tag = TriProps(out, tri, corner)[0];
          EXPECT_GE(tag, 0.0 - 1e-12);
          EXPECT_LE(tag, 5.0 + 1e-12);
        }
      }
    }
  }
  EXPECT_TRUE(sawCut);
  EXPECT_TRUE(sawUnchanged);

  // Geometry check: a side-2 cube (volume 8) minus the 1x1x1 overlap with the
  // shifted cube leaves a volume of 7.
  EXPECT_NEAR(diff.Volume(), 7.0, 1e-9);
}

// ---------------------------------------------------------------------------
// Normal seams: per-face normals on coincident geometric positions stay
// distinct after a Boolean, including across runs from different inputs.
// ---------------------------------------------------------------------------
TEST(BooleanProperties, NormalSeamsStayDiscontinuous) {
  SeamedPair p = MakePair();
  const Manifold result = p.ma + p.mb;
  ASSERT_EQ(result.Status(), Manifold::Error::NoError);
  const MeshGL64 out = result.GetMeshGL64();

  int normalSeams = 0;
  for (const auto& bucket : PositionBuckets(out)) {
    const auto& verts = bucket.second;
    std::set<std::array<int64_t, 3>> normals;
    for (uint64_t v : verts) {
      std::array<int64_t, 3> key{};
      for (int k : {0, 1, 2}) {
        const double n = out.vertProperties[v * out.numProp + 3 + 4 + k];
        key[k] = std::llround(n * 1e9);
      }
      normals.insert(key);
    }
    if (normals.size() > 1) ++normalSeams;
    // Every stored normal remains unit length.
    for (uint64_t v : verts) {
      double len2 = 0;
      for (int k : {0, 1, 2}) {
        const double n = out.vertProperties[v * out.numProp + 3 + 4 + k];
        len2 += n * n;
      }
      EXPECT_NEAR(len2, 1.0, 1e-12);
    }
  }
  // The union keeps every original cube edge that survives; far more than one
  // position carries two different face normals.
  EXPECT_GT(normalSeams, 1);
}

// ---------------------------------------------------------------------------
// Chained Booleans: face provenance stays traceable back to the original
// inputs even after the property indices are reorganized several times, and a
// MeshGL round-trip preserves the surface semantics rather than the concrete
// index layout.
// ---------------------------------------------------------------------------
TEST(BooleanProperties, ChainedBooleansKeepProvenance) {
  SeamedPair p = MakePair();
  SeamedCube c = MakeSeamedCube(vec3(-4, 0, 0), 20.0, 200.0);
  const Manifold mc = Manifold(c.mesh).AsOriginal();
  const uint32_t idC = static_cast<uint32_t>(mc.OriginalID());
  std::map<uint32_t, const MeshGL64*> originals = p.originals;
  originals[idC] = &c.mesh;
  std::map<uint32_t, std::pair<double, double>> tagRanges = p.tagRanges;
  tagRanges[idC] = {20.0, 25.0};

  // A - B leaves a notch; union with a disjoint cube adds a second component.
  const Manifold chained = (p.ma - p.mb) + mc;
  ASSERT_EQ(chained.Status(), Manifold::Error::NoError);

  const MeshGL64 out = chained.GetMeshGL64();
  ExpectClosedTopology(out);
  ExpectMergeVectorsEncodeSeams(out);

  // All three originals remain traceable, and the component touching C keeps
  // the seam structure of its input (its vertices are never cut).
  std::set<uint32_t> ids(out.runOriginalID.begin(), out.runOriginalID.end());
  EXPECT_TRUE(ids.count(p.idA));
  EXPECT_TRUE(ids.count(p.idB));
  EXPECT_TRUE(ids.count(idC));
  ExpectTagsBelongToRunSource(out, tagRanges);
  ExpectBarycentricFromSources(out, originals);

  // Round-trip through MeshGL: indices get re-packed, but the per-face source
  // relation and the property discontinuities must survive.
  const Manifold rebuilt(out);
  ASSERT_EQ(rebuilt.Status(), Manifold::Error::NoError);
  const MeshGL64 out2 = rebuilt.GetMeshGL64();
  ExpectClosedTopology(out2);
  ExpectMergeVectorsEncodeSeams(out2);
  std::set<uint32_t> ids2(out2.runOriginalID.begin(),
                          out2.runOriginalID.end());
  EXPECT_EQ(ids, ids2);
  EXPECT_EQ(CountSeamPositions(out2), CountSeamPositions(out));
  ExpectTagsBelongToRunSource(out2, tagRanges);
  ExpectBarycentricFromSources(out2, originals);

  // A disjoint union preserves both volumes exactly.
  const double expectedVolume =
      (p.ma - p.mb).Volume() + mc.Volume();
  EXPECT_NEAR(chained.Volume(), expectedVolume, 1e-9);
}
