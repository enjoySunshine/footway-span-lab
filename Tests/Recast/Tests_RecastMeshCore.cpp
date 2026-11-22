#include "RecastMeshCore.h"

#include "catch2/catch_amalgamated.hpp"

#include <cstdlib>
#include <vector>

namespace
{
void* testAlloc(size_t size, void*)
{
	return std::malloc(size);
}

void testFree(void* ptr, void*)
{
	std::free(ptr);
}

rcMeshCore::Alloc makeTestAlloc()
{
	rcMeshCore::Alloc alloc;
	alloc.tempAlloc = testAlloc;
	alloc.tempFree = testFree;
	alloc.permAlloc = testAlloc;
	alloc.permFree = testFree;
	alloc.user = nullptr;
	return alloc;
}

// Builds a contour for an axis aligned square in cell coordinates.
std::vector<int> makeSquare(int minX, int minZ, int maxX, int maxZ)
{
	// Counter-clockwise winding in the x/z plane.
	return {
		minX, 0, minZ, 0,
		minX, 0, maxZ, 0,
		maxX, 0, maxZ, 0,
		maxX, 0, minZ, 0,
	};
}

rcMeshCore::Output build(const std::vector<int>& verts, unsigned short reg, unsigned char area)
{
	rcMeshCore::Contour cont;
	cont.verts = verts.data();
	cont.vertsByte = nullptr;
	cont.nverts = static_cast<int>(verts.size()) / 4;
	cont.reg = reg;
	cont.area = area;

	rcMeshCore::Input input;
	input.conts = &cont;
	input.nconts = 1;
	input.borderFlag = 0;
	input.width = 0;
	input.height = 0;
	input.portal = rcMeshCore::PORTAL_NONE;

	rcMeshCore::Alloc alloc = makeTestAlloc();
	rcMeshCore::Output out;
	REQUIRE(rcMeshCore::buildPolyMesh(alloc, input, 6, out, 0));
	return out;
}

bool hasNeighbour(const rcMeshCore::Output& out, int poly, int neighbour)
{
	for (int i = 0; i < out.nvp; ++i)
	{
		if (out.polys[poly * out.nvp * 2 + out.nvp + i] == neighbour)
		{
			return true;
		}
	}
	return false;
}
}

TEST_CASE("rcMeshCore::buildPolyMesh", "[recast][meshtopology]")
{
	SECTION("A single convex contour becomes one polygon")
	{
		const std::vector<int> square = makeSquare(0, 0, 10, 10);
		rcMeshCore::Output out = build(square, 7, 3);

		REQUIRE(out.nverts == 4);
		REQUIRE(out.npolys == 1);
		CHECK(out.nvp == 6);

		// Region and area are carried through from the contour.
		CHECK(out.regs[0] == 7);
		CHECK(out.areas[0] == 3);
	}

	SECTION("Two touching contours become two connected polygons")
	{
		// Two squares sharing the edge x == 10. Because contours are wound
		// consistently, the shared edge must be recognised as a neighbour.
		const std::vector<int> left = makeSquare(0, 0, 10, 10);
		const std::vector<int> right = makeSquare(10, 0, 20, 10);

		rcMeshCore::Contour conts[2];
		conts[0].verts = left.data();
		conts[0].vertsByte = nullptr;
		conts[0].nverts = 4;
		conts[0].reg = 0;
		conts[0].area = 1;
		conts[1].verts = right.data();
		conts[1].vertsByte = nullptr;
		conts[1].nverts = 4;
		conts[1].reg = 1;
		conts[1].area = 1;

		rcMeshCore::Input input;
		input.conts = conts;
		input.nconts = 2;
		input.borderFlag = 0;
		input.width = 0;
		input.height = 0;
		input.portal = rcMeshCore::PORTAL_NONE;

		rcMeshCore::Alloc alloc = makeTestAlloc();
		rcMeshCore::Output out;
		REQUIRE(rcMeshCore::buildPolyMesh(alloc, input, 6, out, 0));

		REQUIRE(out.npolys == 2);
		CHECK(hasNeighbour(out, 0, 1));
		CHECK(hasNeighbour(out, 1, 0));
	}

	SECTION("The portal pass marks only open edges")
	{
		const std::vector<int> square = makeSquare(0, 0, 10, 10);

		rcMeshCore::Contour cont;
		cont.verts = square.data();
		cont.vertsByte = nullptr;
		cont.nverts = 4;
		cont.reg = 0;
		cont.area = 1;

		rcMeshCore::Input input;
		input.conts = &cont;
		input.nconts = 1;
		input.borderFlag = 0;
		input.width = 10;
		input.height = 10;
		input.portal = rcMeshCore::PORTAL_COORD;

		rcMeshCore::Alloc alloc = makeTestAlloc();
		rcMeshCore::Output out;
		REQUIRE(rcMeshCore::buildPolyMesh(alloc, input, 6, out, 0));

		// Every edge of a lone square is on the cell boundary and should be
		// recorded as a portal, never as a polygon neighbour.
		REQUIRE(out.npolys == 1);
		for (int i = 0; i < out.nvp; ++i)
		{
			const unsigned short nei = out.polys[out.nvp + i];
			CHECK((nei & 0x8000) != 0);
		}
	}

	SECTION("Portal direction codes mark the matching open edge")
	{
		const std::vector<int> square = makeSquare(0, 0, 10, 10);

		// Direction codes are stored on the vertex that begins the segment they
		// describe (the tile-cache layout). Tag the segment along z == 0, which
		// runs from vertex 3 back to vertex 0, with direction code 3.
		std::vector<int> verts = square;
		verts[3 * 4 + 3] = 3;

		rcMeshCore::Contour cont;
		cont.verts = verts.data();
		cont.vertsByte = nullptr;
		cont.nverts = 4;
		cont.reg = 0;
		cont.area = 1;

		rcMeshCore::Input input;
		input.conts = &cont;
		input.nconts = 1;
		input.borderFlag = 0;
		input.width = 0;
		input.height = 0;
		input.portal = rcMeshCore::PORTAL_CONTOUR_DIR;

		rcMeshCore::Alloc alloc = makeTestAlloc();
		rcMeshCore::Output out;
		REQUIRE(rcMeshCore::buildPolyMesh(alloc, input, 6, out, 0));

		REQUIRE(out.npolys == 1);
		int portalled = 0;
		for (int i = 0; i < out.nvp; ++i)
		{
			const unsigned short nei = out.polys[out.nvp + i];
			if (nei == (0x8000 | 3))
			{
				portalled++;
			}
		}
		CHECK(portalled == 1);
	}
}
