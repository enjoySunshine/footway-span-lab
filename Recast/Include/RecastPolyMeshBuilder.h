//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef RECASTPOLYMESHBUILDER_H
#define RECASTPOLYMESHBUILDER_H

#include "Recast.h"

#include <stddef.h>

/// @file RecastPolyMeshBuilder.h
///
/// The single contour-to-polygon-mesh pass used by both Recast and DetourTileCache.
///
/// Recast's rcBuildPolyMesh() and DetourTileCache's dtBuildTileCachePolyMesh() used to
/// carry separate copies of the triangulation, polygon merging, border-vertex removal and
/// adjacency code. The copies drifted. In particular the adjacency pass could annotate a
/// tile-border edge differently depending on which copy ran, so a tile rebuilt after an
/// obstacle moved could link its polygons differently than the original bake, and the same
/// start/end pair could straighten into a different path.
///
/// Both callers now build their mesh through rcBuildPolyMeshFromContours(). The algorithm
/// (and therefore the adjacency) is one implementation. Only the pieces that genuinely
/// differ between the two data formats are parameters:
///
///  - Which flag bit marks a removable build-border vertex.
///    Recast uses #RC_BORDER_VERTEX; DetourTileCache uses 0x80.
///  - How a tile-border edge is recorded as a portal to the neighboring tile.
///    Recast's contours carry no portal direction, so it tests the vertex coordinates
///    against the build border. DetourTileCache's layer format stores the portal direction
///    in each contour vertex, so it uses that.
///  - Whether triangulation may retry with relaxed diagonal tests when the strict pass
///    stalls. Recast does; DetourTileCache has never done so and keeps its behavior.
///  - The vertex hash size used when merging duplicate contour vertices, and the edge
///    budget for vertex removal. Both affect only degenerate inputs but are preserved so
///    each caller produces the exact mesh it produced before.
///  - The allocator used for temporary buffers. DetourTileCache passes its tile-scoped
///    linear allocator, so rebuilding the tiles under a group of obstacles uses memory
///    bounded by one tile rather than by the size of the whole navigation mesh.
///
/// Callers keep ownership of their own contour and output buffers.

/// Describes one simplified contour in voxel space.
struct rcPolyMeshBuilderContour
{
	/// Vertex data, 4 values per vertex: (x, y, z, flags).
	/// The flags carry border-vertex and portal information as interpreted by the caller.
	const int* verts;
	/// The number of vertices in the contour.
	int nverts;
	/// The region id assigned to polygons made from this contour.
	int reg;
	/// The area id assigned to polygons made from this contour.
	unsigned char area;

	rcPolyMeshBuilderContour() : verts(0), nverts(0), reg(0), area(0) {}
};

/// Selects how open contour edges are annotated as portals to a neighboring tile.
enum rcPolyMeshBuilderPortalMode
{
	/// Do not annotate any edge. (Used when rebuilding adjacency on an already-built mesh.)
	RC_POLYMESH_PORTAL_NONE = 0,
	/// Annotate edges whose endpoints both sit on the build border. (Recast mesh build.)
	RC_POLYMESH_PORTAL_BORDER = 1,
	/// Annotate edges using the direction stored in each contour vertex's flags.
	/// (DetourTileCache, where the direction is the portal bank of the layer.)
	RC_POLYMESH_PORTAL_BANK = 2
};

/// Optional allocator for the temporary buffers of a build.
///
/// When @p alloc is null the build uses rcAlloc()/rcFree(). DetourTileCache supplies an
/// adapter around its tile allocator so a tile rebuild does not allocate against the total
/// navigation mesh.
struct rcPolyMeshBuilderAllocator
{
	/// Allocates @p size bytes. Returns null on failure.
	void* (*alloc)(void* user, size_t size);
	/// Releases memory previously returned by @p alloc. May be null.
	void (*free)(void* user, void* ptr);
	/// Opaque pointer passed back to @p alloc and @p free.
	void* user;

	rcPolyMeshBuilderAllocator() : alloc(0), free(0), user(0) {}
};

/// Parameters for rcBuildPolyMeshFromContours().
struct rcPolyMeshBuilderParams
{
	/// The contours to polygonize.
	const rcPolyMeshBuilderContour* conts;
	/// The number of contours in @p conts.
	int nconts;
	/// The maximum number of vertices per polygon.
	int nvp;

	/// Flag bit that marks a contour vertex as a build-border vertex eligible for removal.
	/// Zero disables border-vertex removal.
	int borderVertexFlag;
	/// How to annotate portal edges. See rcPolyMeshBuilderPortalMode.
	int portalMode;
	/// Border width in cells, used only by RC_POLYMESH_PORTAL_BORDER.
	int width;
	/// Border height in cells, used only by RC_POLYMESH_PORTAL_BORDER.
	int height;

	/// Upper bound on the number of edges touched while removing a single vertex.
	/// Zero means unlimited. Applies to both the eligibility test and the temporary storage.
	int maxRemovalEdges;
	/// Whether triangulation may retry with relaxed diagonal tests when the strict pass stalls.
	bool allowLooseTriangulation;
	/// Number of hash buckets used when merging duplicate contour vertices.
	int vertexBucketCount;

	/// Optional allocator for temporary buffers.
	rcPolyMeshBuilderAllocator allocator;

	/// Output vertices. Must hold at least @p maxVerts * 3 values.
	unsigned short* outVerts;
	/// Capacity of @p outVerts, in vertices.
	int maxVerts;
	/// Output polygons. Must hold at least @p maxPolys * 2 * @p nvp values.
	unsigned short* outPolys;
	/// Optional output regions. Must hold at least @p maxPolys values. May be null.
	unsigned short* outRegs;
	/// Output areas. Must hold at least @p maxPolys values.
	unsigned char* outAreas;
	/// Capacity of @p outPolys, @p outRegs and @p outAreas, in polygons.
	int maxPolys;

	rcPolyMeshBuilderParams()
	: conts(0), nconts(0), nvp(0)
	, borderVertexFlag(0), portalMode(RC_POLYMESH_PORTAL_NONE), width(0), height(0)
	, maxRemovalEdges(0), allowLooseTriangulation(true), vertexBucketCount(4096)
	, allocator()
	, outVerts(0), maxVerts(0), outPolys(0), outRegs(0), outAreas(0), maxPolys(0)
	{
	}
};

/// Results of rcBuildPolyMeshFromContours().
struct rcPolyMeshBuilderResult
{
	/// The number of vertices written.
	int nverts;
	/// The number of polygons written.
	int npolys;

	rcPolyMeshBuilderResult() : nverts(0), npolys(0) {}
};

/// Builds a polygon mesh from simplified contours using the shared implementation.
///
/// Neighbor slots in @p outPolys are initialized to #RC_MESH_NULL_IDX and filled by the same
/// adjacency pass for every caller, so a mesh baked up front and a mesh rebuilt later link
/// their polygons identically.
///
/// @param[in]		ctx		The context to log to. May be null.
/// @param[in]		params	The contours, options and output buffers.
/// @param[out]		result	The number of vertices and polygons written.
/// @return True if the mesh was built.
/// @see rcBuildPolyMesh, dtBuildTileCachePolyMesh
bool rcBuildPolyMeshFromContours(
	rcContext* ctx,
	const rcPolyMeshBuilderParams& params,
	rcPolyMeshBuilderResult& result);

/// Rebuilds the neighbor links of an existing polygon mesh using the shared implementation.
///
/// This is the same adjacency pass used by rcBuildPolyMeshFromContours(). It is exposed so
/// that operations on already-built meshes (such as #rcMergePolyMeshes) use the same code.
///
/// Existing neighbor slots are not cleared, so portal markers a caller has already stored in
/// unconnected slots survive the rebuild.
///
/// @param[in,out]	polys	The polygon array. Neighbor slots (the second half of each polygon)
///							are rewritten for connected edges.
/// @param[in]		npolys	The number of polygons.
/// @param[in]		nverts	The number of vertices.
/// @param[in]		nvp		The maximum number of vertices per polygon.
/// @param[in]		verts	The mesh vertices, 3 values per vertex.
/// @param[in]		portalMode	How to annotate portal edges. See rcPolyMeshBuilderPortalMode.
/// @param[in]		width	Border width in cells, used only by RC_POLYMESH_PORTAL_BORDER.
/// @param[in]		height	Border height in cells, used only by RC_POLYMESH_PORTAL_BORDER.
/// @param[in]		conts	Contours used by RC_POLYMESH_PORTAL_BANK. May be null otherwise.
/// @param[in]		nconts	The number of contours in @p conts.
/// @param[in]		alloc	Optional temporary allocator. May be null.
/// @return True if adjacency was built.
bool rcBuildPolyMeshAdjacency(
	unsigned short* polys,
	const int npolys,
	const int nverts,
	const int nvp,
	const unsigned short* verts,
	const int portalMode,
	const int width,
	const int height,
	const rcPolyMeshBuilderContour* conts,
	const int nconts,
	const rcPolyMeshBuilderAllocator* alloc);

#endif  // RECASTPOLYMESHBUILDER_H
