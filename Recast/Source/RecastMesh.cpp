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

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"
#include "RecastPolyMeshBuilder.h"

// The contour-to-polygon-mesh pass (triangulation, polygon merging, border-vertex removal
// and adjacency) lives in RecastPolyMeshBuilder.cpp and is shared with DetourTileCache.
// This file keeps the vertex helper used by rcMergePolyMeshes.

static const int VERTEX_BUCKET_COUNT = (1<<12);

inline int computeVertexHash(int x, int y, int z)
{
	const unsigned int h1 = 0x8da6b343; // Large multiplicative constants;
	const unsigned int h2 = 0xd8163841; // here arbitrarily chosen primes
	const unsigned int h3 = 0xcb1ab31f;
	unsigned int n = h1 * x + h2 * y + h3 * z;
	return (int)(n & (VERTEX_BUCKET_COUNT-1));
}

/// Adds a vertex to the mesh, merging it with an existing vertex at the same cell.
/// Shared with rcBuildPolyMeshFromContours() through the same hashing rule.
static unsigned short addVertex(unsigned short x, unsigned short y, unsigned short z,
								unsigned short* verts, int* firstVert, int* nextVert, int& nv)
{
	int bucket = computeVertexHash(x, 0, z);
	int i = firstVert[bucket];

	while (i != -1)
	{
		const unsigned short* v = &verts[i*3];
		if (v[0] == x && (rcAbs(v[1] - y) <= 2) && v[2] == z)
			return (unsigned short)i;
		i = nextVert[i]; // next
	}

	// Could not find, create new.
	i = nv; nv++;
	unsigned short* v = &verts[i*3];
	v[0] = x;
	v[1] = y;
	v[2] = z;
	nextVert[i] = firstVert[bucket];
	firstVert[bucket] = i;

	return (unsigned short)i;
}


/// @par
///
/// @note If the mesh data is to be used to construct a Detour navigation mesh, then the upper 
/// limit must be restricted to <= #DT_VERTS_PER_POLYGON.
///
/// @see rcAllocPolyMesh, rcContourSet, rcPolyMesh, rcConfig
bool rcBuildPolyMesh(rcContext* ctx, const rcContourSet& cset, const int nvp, rcPolyMesh& mesh)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_BUILD_POLYMESH);

	rcVcopy(mesh.bmin, cset.bmin);
	rcVcopy(mesh.bmax, cset.bmax);
	mesh.cs = cset.cs;
	mesh.ch = cset.ch;
	mesh.borderSize = cset.borderSize;
	mesh.maxEdgeError = cset.maxError;

	int maxVertices = 0;
	int maxTris = 0;
	for (int i = 0; i < cset.nconts; ++i)
	{
		// Skip null contours.
		if (cset.conts[i].nverts < 3) continue;
		maxVertices += cset.conts[i].nverts;
		maxTris += cset.conts[i].nverts - 2;
	}

	if (maxVertices >= 0xfffe)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Too many vertices %d.", maxVertices);
		return false;
	}

	mesh.verts = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxVertices*3, RC_ALLOC_PERM);
	if (!mesh.verts)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.verts' (%d).", maxVertices);
		return false;
	}
	mesh.polys = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxTris*nvp*2, RC_ALLOC_PERM);
	if (!mesh.polys)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.polys' (%d).", maxTris*nvp*2);
		return false;
	}
	mesh.regs = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxTris, RC_ALLOC_PERM);
	if (!mesh.regs)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.regs' (%d).", maxTris);
		return false;
	}
	mesh.areas = (unsigned char*)rcAlloc(sizeof(unsigned char)*maxTris, RC_ALLOC_PERM);
	if (!mesh.areas)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.areas' (%d).", maxTris);
		return false;
	}

	mesh.nverts = 0;
	mesh.npolys = 0;
	mesh.nvp = nvp;
	mesh.maxpolys = maxTris;

	// Describe the contours for the shared polygonization pass.
	rcPolyMeshBuilderContour* builderContours = (rcPolyMeshBuilderContour*)rcAlloc(
		sizeof(rcPolyMeshBuilderContour)*cset.nconts, RC_ALLOC_TEMP);
	if (!builderContours)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'builderContours' (%d).", cset.nconts);
		return false;
	}
	for (int i = 0; i < cset.nconts; ++i)
	{
		builderContours[i].verts = cset.conts[i].verts;
		builderContours[i].nverts = cset.conts[i].nverts;
		builderContours[i].reg = cset.conts[i].reg;
		builderContours[i].area = cset.conts[i].area;
	}

	rcPolyMeshBuilderParams params;
	params.conts = builderContours;
	params.nconts = cset.nconts;
	params.nvp = nvp;
	params.borderVertexFlag = RC_BORDER_VERTEX;
	params.portalMode = (mesh.borderSize > 0) ? RC_POLYMESH_PORTAL_BORDER : RC_POLYMESH_PORTAL_NONE;
	params.width = cset.width;
	params.height = cset.height;
	params.maxRemovalEdges = 0;
	params.allowLooseTriangulation = true;
	params.vertexBucketCount = 4096;
	params.outVerts = mesh.verts;
	params.maxVerts = maxVertices;
	params.outPolys = mesh.polys;
	params.outRegs = mesh.regs;
	params.outAreas = mesh.areas;
	params.maxPolys = maxTris;

	rcPolyMeshBuilderResult result;
	const bool built = rcBuildPolyMeshFromContours(ctx, params, result);
	rcFree(builderContours);
	if (!built)
	{
		return false;
	}

	mesh.nverts = result.nverts;
	mesh.npolys = result.npolys;

	// Just allocate the mesh flags array. The user is resposible to fill it.
	mesh.flags = (unsigned short*)rcAlloc(sizeof(unsigned short)*mesh.npolys, RC_ALLOC_PERM);
	if (!mesh.flags)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: Out of memory 'mesh.flags' (%d).", mesh.npolys);
		return false;
	}
	memset(mesh.flags, 0, sizeof(unsigned short) * mesh.npolys);

	if (mesh.nverts > 0xffff)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.nverts, 0xffff);
	}
	if (mesh.npolys > 0xffff)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildPolyMesh: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.npolys, 0xffff);
	}

	return true;
}

/// @see rcAllocPolyMesh, rcPolyMesh
bool rcMergePolyMeshes(rcContext* ctx, rcPolyMesh** meshes, const int nmeshes, rcPolyMesh& mesh)
{
	rcAssert(ctx);
	
	if (!nmeshes || !meshes)
		return true;

	rcScopedTimer timer(ctx, RC_TIMER_MERGE_POLYMESH);

	mesh.nvp = meshes[0]->nvp;
	mesh.cs = meshes[0]->cs;
	mesh.ch = meshes[0]->ch;
	rcVcopy(mesh.bmin, meshes[0]->bmin);
	rcVcopy(mesh.bmax, meshes[0]->bmax);

	int maxVerts = 0;
	int maxPolys = 0;
	int maxVertsPerMesh = 0;
	for (int i = 0; i < nmeshes; ++i)
	{
		rcVmin(mesh.bmin, meshes[i]->bmin);
		rcVmax(mesh.bmax, meshes[i]->bmax);
		maxVertsPerMesh = rcMax(maxVertsPerMesh, meshes[i]->nverts);
		maxVerts += meshes[i]->nverts;
		maxPolys += meshes[i]->npolys;
	}
	
	mesh.nverts = 0;
	mesh.verts = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxVerts*3, RC_ALLOC_PERM);
	if (!mesh.verts)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.verts' (%d).", maxVerts*3);
		return false;
	}

	mesh.npolys = 0;
	mesh.polys = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxPolys*2*mesh.nvp, RC_ALLOC_PERM);
	if (!mesh.polys)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.polys' (%d).", maxPolys*2*mesh.nvp);
		return false;
	}
	memset(mesh.polys, 0xff, sizeof(unsigned short)*maxPolys*2*mesh.nvp);

	mesh.regs = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxPolys, RC_ALLOC_PERM);
	if (!mesh.regs)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.regs' (%d).", maxPolys);
		return false;
	}
	memset(mesh.regs, 0, sizeof(unsigned short)*maxPolys);

	mesh.areas = (unsigned char*)rcAlloc(sizeof(unsigned char)*maxPolys, RC_ALLOC_PERM);
	if (!mesh.areas)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.areas' (%d).", maxPolys);
		return false;
	}
	memset(mesh.areas, 0, sizeof(unsigned char)*maxPolys);

	mesh.flags = (unsigned short*)rcAlloc(sizeof(unsigned short)*maxPolys, RC_ALLOC_PERM);
	if (!mesh.flags)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'mesh.flags' (%d).", maxPolys);
		return false;
	}
	memset(mesh.flags, 0, sizeof(unsigned short)*maxPolys);
	
	rcScopedDelete<int> nextVert((int*)rcAlloc(sizeof(int)*maxVerts, RC_ALLOC_TEMP));
	if (!nextVert)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'nextVert' (%d).", maxVerts);
		return false;
	}
	memset(nextVert, 0, sizeof(int)*maxVerts);
	
	rcScopedDelete<int> firstVert((int*)rcAlloc(sizeof(int)*VERTEX_BUCKET_COUNT, RC_ALLOC_TEMP));
	if (!firstVert)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'firstVert' (%d).", VERTEX_BUCKET_COUNT);
		return false;
	}
	for (int i = 0; i < VERTEX_BUCKET_COUNT; ++i)
		firstVert[i] = -1;

	rcScopedDelete<unsigned short> vremap((unsigned short*)rcAlloc(sizeof(unsigned short)*maxVertsPerMesh, RC_ALLOC_PERM));
	if (!vremap)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Out of memory 'vremap' (%d).", maxVertsPerMesh);
		return false;
	}
	memset(vremap, 0, sizeof(unsigned short)*maxVertsPerMesh);
	
	for (int i = 0; i < nmeshes; ++i)
	{
		const rcPolyMesh* pmesh = meshes[i];
		
		const unsigned short ox = (unsigned short)floorf((pmesh->bmin[0]-mesh.bmin[0])/mesh.cs+0.5f);
		const unsigned short oz = (unsigned short)floorf((pmesh->bmin[2]-mesh.bmin[2])/mesh.cs+0.5f);
		
		bool isMinX = (ox == 0);
		bool isMinZ = (oz == 0);
		bool isMaxX = ((unsigned short)floorf((mesh.bmax[0] - pmesh->bmax[0]) / mesh.cs + 0.5f)) == 0;
		bool isMaxZ = ((unsigned short)floorf((mesh.bmax[2] - pmesh->bmax[2]) / mesh.cs + 0.5f)) == 0;
		bool isOnBorder = (isMinX || isMinZ || isMaxX || isMaxZ);

		for (int j = 0; j < pmesh->nverts; ++j)
		{
			unsigned short* v = &pmesh->verts[j*3];
			vremap[j] = addVertex(v[0]+ox, v[1], v[2]+oz,
								  mesh.verts, firstVert, nextVert, mesh.nverts);
		}
		
		for (int j = 0; j < pmesh->npolys; ++j)
		{
			unsigned short* tgt = &mesh.polys[mesh.npolys*2*mesh.nvp];
			unsigned short* src = &pmesh->polys[j*2*mesh.nvp];
			mesh.regs[mesh.npolys] = pmesh->regs[j];
			mesh.areas[mesh.npolys] = pmesh->areas[j];
			mesh.flags[mesh.npolys] = pmesh->flags[j];
			mesh.npolys++;
			for (int k = 0; k < mesh.nvp; ++k)
			{
				if (src[k] == RC_MESH_NULL_IDX) break;
				tgt[k] = vremap[src[k]];
			}

			if (isOnBorder)
			{
				for (int k = mesh.nvp; k < mesh.nvp * 2; ++k)
				{
					if (src[k] & 0x8000 && src[k] != 0xffff)
					{
						unsigned short dir = src[k] & 0xf;
						switch (dir)
						{
							case 0: // Portal x-
								if (isMinX)
									tgt[k] = src[k];
								break;
							case 1: // Portal z+
								if (isMaxZ)
									tgt[k] = src[k];
								break;
							case 2: // Portal x+
								if (isMaxX)
									tgt[k] = src[k];
								break;
							case 3: // Portal z-
								if (isMinZ)
									tgt[k] = src[k];
								break;
						}
					}
				}
			}
		}
	}

	// Calculate adjacency. Portal markers already copied into unconnected slots survive.
	if (!rcBuildPolyMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, mesh.nvp, mesh.verts,
		RC_POLYMESH_PORTAL_NONE, 0, 0, 0, 0, 0))
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: Adjacency failed.");
		return false;
	}

	if (mesh.nverts > 0xffff)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.nverts, 0xffff);
	}
	if (mesh.npolys > 0xffff)
	{
		ctx->log(RC_LOG_ERROR, "rcMergePolyMeshes: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.npolys, 0xffff);
	}
	
	return true;
}

bool rcCopyPolyMesh(rcContext* ctx, const rcPolyMesh& src, rcPolyMesh& dst)
{
	rcAssert(ctx);
	
	// Destination must be empty.
	rcAssert(dst.verts == 0);
	rcAssert(dst.polys == 0);
	rcAssert(dst.regs == 0);
	rcAssert(dst.areas == 0);
	rcAssert(dst.flags == 0);
	
	dst.nverts = src.nverts;
	dst.npolys = src.npolys;
	dst.maxpolys = src.npolys;
	dst.nvp = src.nvp;
	rcVcopy(dst.bmin, src.bmin);
	rcVcopy(dst.bmax, src.bmax);
	dst.cs = src.cs;
	dst.ch = src.ch;
	dst.borderSize = src.borderSize;
	dst.maxEdgeError = src.maxEdgeError;
	
	dst.verts = (unsigned short*)rcAlloc(sizeof(unsigned short)*src.nverts*3, RC_ALLOC_PERM);
	if (!dst.verts)
	{
		ctx->log(RC_LOG_ERROR, "rcCopyPolyMesh: Out of memory 'dst.verts' (%d).", src.nverts*3);
		return false;
	}
	memcpy(dst.verts, src.verts, sizeof(unsigned short)*src.nverts*3);
	
	dst.polys = (unsigned short*)rcAlloc(sizeof(unsigned short)*src.npolys*2*src.nvp, RC_ALLOC_PERM);
	if (!dst.polys)
	{
		ctx->log(RC_LOG_ERROR, "rcCopyPolyMesh: Out of memory 'dst.polys' (%d).", src.npolys*2*src.nvp);
		return false;
	}
	memcpy(dst.polys, src.polys, sizeof(unsigned short)*src.npolys*2*src.nvp);
	
	dst.regs = (unsigned short*)rcAlloc(sizeof(unsigned short)*src.npolys, RC_ALLOC_PERM);
	if (!dst.regs)
	{
		ctx->log(RC_LOG_ERROR, "rcCopyPolyMesh: Out of memory 'dst.regs' (%d).", src.npolys);
		return false;
	}
	memcpy(dst.regs, src.regs, sizeof(unsigned short)*src.npolys);
	
	dst.areas = (unsigned char*)rcAlloc(sizeof(unsigned char)*src.npolys, RC_ALLOC_PERM);
	if (!dst.areas)
	{
		ctx->log(RC_LOG_ERROR, "rcCopyPolyMesh: Out of memory 'dst.areas' (%d).", src.npolys);
		return false;
	}
	memcpy(dst.areas, src.areas, sizeof(unsigned char)*src.npolys);
	
	dst.flags = (unsigned short*)rcAlloc(sizeof(unsigned short)*src.npolys, RC_ALLOC_PERM);
	if (!dst.flags)
	{
		ctx->log(RC_LOG_ERROR, "rcCopyPolyMesh: Out of memory 'dst.flags' (%d).", src.npolys);
		return false;
	}
	memcpy(dst.flags, src.flags, sizeof(unsigned short)*src.npolys);
	
	return true;
}
