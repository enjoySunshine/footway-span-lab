#include "catch2/catch_amalgamated.hpp"

#include "RecastPolyMeshBuilder.h"

#include <string.h>

namespace
{
// A unit square in cell coordinates, wound counter-clockwise the way Recast contours are
// (so that the interior lies to the left of each edge), with the given flag on every vertex.
void makeSquareContour(int verts[4][4], const int x0, const int z0, const int flags)
{
	const int cells[4][2] = {{x0, z0}, {x0, z0 + 1}, {x0 + 1, z0 + 1}, {x0 + 1, z0}};
	for (int i = 0; i < 4; ++i)
	{
		verts[i][0] = cells[i][0];
		verts[i][1] = 0;
		verts[i][2] = cells[i][1];
		verts[i][3] = flags;
	}
}

struct Buffers
{
	unsigned short verts[256 * 3];
	unsigned short polys[256 * 2 * 6];
	unsigned short regs[256];
	unsigned char areas[256];

	Buffers()
	{
		memset(verts, 0, sizeof(verts));
		memset(polys, 0xff, sizeof(polys));
		memset(regs, 0, sizeof(regs));
		memset(areas, 0, sizeof(areas));
	}
};

rcPolyMeshBuilderParams makeParams(const rcPolyMeshBuilderContour* conts, const int nconts, Buffers& buf)
{
	rcPolyMeshBuilderParams params;
	params.conts = conts;
	params.nconts = nconts;
	params.nvp = 6;
	params.borderVertexFlag = 0;
	params.portalMode = RC_POLYMESH_PORTAL_NONE;
	params.width = 0;
	params.height = 0;
	params.maxRemovalEdges = 0;
	params.allowLooseTriangulation = true;
	params.vertexBucketCount = 4096;
	params.outVerts = buf.verts;
	params.maxVerts = 256;
	params.outPolys = buf.polys;
	params.outRegs = buf.regs;
	params.outAreas = buf.areas;
	params.maxPolys = 256;
	return params;
}
}

TEST_CASE("rcBuildPolyMeshFromContours", "[recast, polymeshbuilder]")
{
	SECTION("Builds a single square polygon")
	{
		int square[4][4];
		makeSquareContour(square, 10, 20, 0);

		rcPolyMeshBuilderContour cont;
		cont.verts = &square[0][0];
		cont.nverts = 4;
		cont.reg = 7;
		cont.area = 63;

		Buffers buf;
		rcPolyMeshBuilderParams params = makeParams(&cont, 1, buf);
		rcPolyMeshBuilderResult result;

		REQUIRE(rcBuildPolyMeshFromContours(0, params, result));
		REQUIRE(result.nverts == 4);
		REQUIRE(result.npolys == 1);

		// The square has four edges and no neighbors.
		for (int i = 0; i < 4; ++i)
		{
			REQUIRE(buf.polys[6 + i] == RC_MESH_NULL_IDX);
		}
		REQUIRE(buf.areas[0] == 63);
		REQUIRE(buf.regs[0] == 7);
	}

	SECTION("Neighbor links are symmetric across a shared edge")
	{
		// Two squares side by side, sharing the edge x = 11.
		int left[4][4];
		int right[4][4];
		makeSquareContour(left, 10, 20, 0);
		makeSquareContour(right, 11, 20, 0);

		rcPolyMeshBuilderContour conts[2];
		conts[0].verts = &left[0][0];
		conts[0].nverts = 4;
		conts[0].reg = 1;
		conts[0].area = 63;
		conts[1].verts = &right[0][0];
		conts[1].nverts = 4;
		conts[1].reg = 2;
		conts[1].area = 63;

		Buffers buf;
		rcPolyMeshBuilderParams params = makeParams(conts, 2, buf);
		rcPolyMeshBuilderResult result;

		REQUIRE(rcBuildPolyMeshFromContours(0, params, result));
		REQUIRE(result.npolys == 2);

		// Find the one connected edge on each polygon and confirm it points at the other.
		int connections = 0;
		for (int poly = 0; poly < result.npolys; ++poly)
		{
			for (int edge = 0; edge < 6; ++edge)
			{
				const unsigned short nei = buf.polys[poly * 12 + 6 + edge];
				if (nei != RC_MESH_NULL_IDX)
				{
					REQUIRE(nei < result.npolys);
					REQUIRE(nei != poly);
					connections++;
				}
			}
		}
		REQUIRE(connections == 2);
	}

	SECTION("Adjacency rebuilt later matches the build-time adjacency")
	{
		// Build once, then rebuild only the adjacency pass on a copy. The two must agree.
		int left[4][4];
		int right[4][4];
		makeSquareContour(left, 10, 20, 0);
		makeSquareContour(right, 11, 20, 0);

		rcPolyMeshBuilderContour conts[2];
		conts[0].verts = &left[0][0];
		conts[0].nverts = 4;
		conts[0].reg = 1;
		conts[0].area = 63;
		conts[1].verts = &right[0][0];
		conts[1].nverts = 4;
		conts[1].reg = 2;
		conts[1].area = 63;

		Buffers built;
		rcPolyMeshBuilderParams params = makeParams(conts, 2, built);
		rcPolyMeshBuilderResult result;
		REQUIRE(rcBuildPolyMeshFromContours(0, params, result));

		unsigned short rebuilt[2 * 2 * 6];
		memcpy(rebuilt, built.polys, sizeof(rebuilt));
		for (int i = 0; i < result.npolys; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				rebuilt[i * 12 + 6 + j] = RC_MESH_NULL_IDX;
			}
		}

		REQUIRE(rcBuildPolyMeshAdjacency(rebuilt, result.npolys, result.nverts, 6, built.verts,
			RC_POLYMESH_PORTAL_NONE, 0, 0, 0, 0, 0));

		REQUIRE(memcmp(rebuilt, built.polys, sizeof(rebuilt)) == 0);
	}

	SECTION("Border portals are marked when requested")
	{
		// A square whose z = 0 edge lies on the border.
		int square[4][4];
		makeSquareContour(square, 10, 0, 0);

		rcPolyMeshBuilderContour cont;
		cont.verts = &square[0][0];
		cont.nverts = 4;
		cont.reg = 0;
		cont.area = 63;

		Buffers buf;
		rcPolyMeshBuilderParams params = makeParams(&cont, 1, buf);
		params.portalMode = RC_POLYMESH_PORTAL_BORDER;
		params.width = 48;
		params.height = 48;
		rcPolyMeshBuilderResult result;

		REQUIRE(rcBuildPolyMeshFromContours(0, params, result));

		int portals = 0;
		for (int edge = 0; edge < 6; ++edge)
		{
			const unsigned short nei = buf.polys[6 + edge];
			if (nei != RC_MESH_NULL_IDX)
			{
				REQUIRE((nei & 0x8000) != 0);
				portals++;
			}
		}
		REQUIRE(portals == 1);
	}

	SECTION("Fails cleanly when the polygon buffer is too small")
	{
		int left[4][4];
		int right[4][4];
		makeSquareContour(left, 10, 20, 0);
		makeSquareContour(right, 11, 20, 0);

		rcPolyMeshBuilderContour conts[2];
		conts[0].verts = &left[0][0];
		conts[0].nverts = 4;
		conts[0].reg = 0;
		conts[0].area = 63;
		conts[1].verts = &right[0][0];
		conts[1].nverts = 4;
		conts[1].reg = 0;
		conts[1].area = 63;

		Buffers buf;
		rcPolyMeshBuilderParams params = makeParams(conts, 2, buf);
		params.maxPolys = 1;
		rcPolyMeshBuilderResult result;

		REQUIRE_FALSE(rcBuildPolyMeshFromContours(0, params, result));
	}
}
