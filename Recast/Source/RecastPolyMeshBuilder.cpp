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

#include "RecastPolyMeshBuilder.h"

#include "RecastAlloc.h"

#include <math.h>
#include <string.h>

namespace
{

/// Temporary allocator that forwards to rcAlloc()/rcFree().
void* defaultAlloc(void* /*user*/, size_t size)
{
	return rcAlloc(size, RC_ALLOC_TEMP);
}

void defaultFree(void* /*user*/, void* ptr)
{
	rcFree(ptr);
}

/// A builder-local allocator handle.
struct Allocator
{
	void* (*allocFn)(void* user, size_t size);
	void (*freeFn)(void* user, void* ptr);
	void* user;

	Allocator() : allocFn(defaultAlloc), freeFn(defaultFree), user(0) {}

	void* alloc(size_t size) const { return allocFn(user, size); }
	void free(void* ptr) const
	{
		if (freeFn)
		{
			freeFn(user, ptr);
		}
	}
};

static const unsigned short NULL_IDX = RC_MESH_NULL_IDX;

static const int DIAGONAL_FLAG = (int)0x80000000u;
static const int DIAGONAL_INDEX_MASK = 0x0fffffff;

inline int prevIndex(int i, int n) { return i - 1 >= 0 ? i - 1 : n - 1; }
inline int nextIndex(int i, int n) { return i + 1 < n ? i + 1 : 0; }

inline int area2(const int* a, const int* b, const int* c)
{
	return (b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2]);
}

inline bool left(const int* a, const int* b, const int* c) { return area2(a, b, c) < 0; }
inline bool leftOn(const int* a, const int* b, const int* c) { return area2(a, b, c) <= 0; }
inline bool collinear(const int* a, const int* b, const int* c) { return area2(a, b, c) == 0; }

inline bool intersectProp(const int* a, const int* b, const int* c, const int* d)
{
	if (collinear(a, b, c) || collinear(a, b, d) || collinear(c, d, a) || collinear(c, d, b))
	{
		return false;
	}
	return (left(a, b, c) ^ left(a, b, d)) && (left(c, d, a) ^ left(c, d, b));
}

inline bool between(const int* a, const int* b, const int* c)
{
	if (!collinear(a, b, c))
	{
		return false;
	}
	if (a[0] != b[0])
	{
		return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
	}
	return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
}

inline bool intersect(const int* a, const int* b, const int* c, const int* d)
{
	if (intersectProp(a, b, c, d))
	{
		return true;
	}
	return between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b);
}

inline bool vequal(const int* a, const int* b) { return a[0] == b[0] && a[2] == b[2]; }

inline bool diagonalie(int i, int j, int n, const int* verts, int* indices)
{
	const int* d0 = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
	const int* d1 = &verts[(indices[j] & DIAGONAL_INDEX_MASK) * 4];
	for (int k = 0; k < n; k++)
	{
		const int k1 = nextIndex(k, n);
		if (!((k == i) || (k1 == i) || (k == j) || (k1 == j)))
		{
			const int* p0 = &verts[(indices[k] & DIAGONAL_INDEX_MASK) * 4];
			const int* p1 = &verts[(indices[k1] & DIAGONAL_INDEX_MASK) * 4];
			if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1))
			{
				continue;
			}
			if (intersect(d0, d1, p0, p1))
			{
				return false;
			}
		}
	}
	return true;
}

inline bool inCone(int i, int j, int n, const int* verts, int* indices)
{
	const int* pi = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
	const int* pj = &verts[(indices[j] & DIAGONAL_INDEX_MASK) * 4];
	const int* pi1 = &verts[(indices[nextIndex(i, n)] & DIAGONAL_INDEX_MASK) * 4];
	const int* pin1 = &verts[(indices[prevIndex(i, n)] & DIAGONAL_INDEX_MASK) * 4];
	if (leftOn(pin1, pi, pi1))
	{
		return left(pi, pj, pin1) && left(pj, pi, pi1);
	}
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}

inline bool diagonal(int i, int j, int n, const int* verts, int* indices)
{
	return inCone(i, j, n, verts, indices) && diagonalie(i, j, n, verts, indices);
}

inline bool diagonalieLoose(int i, int j, int n, const int* verts, int* indices)
{
	const int* d0 = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
	const int* d1 = &verts[(indices[j] & DIAGONAL_INDEX_MASK) * 4];
	for (int k = 0; k < n; k++)
	{
		const int k1 = nextIndex(k, n);
		if (!((k == i) || (k1 == i) || (k == j) || (k1 == j)))
		{
			const int* p0 = &verts[(indices[k] & DIAGONAL_INDEX_MASK) * 4];
			const int* p1 = &verts[(indices[k1] & DIAGONAL_INDEX_MASK) * 4];
			if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1))
			{
				continue;
			}
			if (intersectProp(d0, d1, p0, p1))
			{
				return false;
			}
		}
	}
	return true;
}

inline bool inConeLoose(int i, int j, int n, const int* verts, int* indices)
{
	const int* pi = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
	const int* pj = &verts[(indices[j] & DIAGONAL_INDEX_MASK) * 4];
	const int* pi1 = &verts[(indices[nextIndex(i, n)] & DIAGONAL_INDEX_MASK) * 4];
	const int* pin1 = &verts[(indices[prevIndex(i, n)] & DIAGONAL_INDEX_MASK) * 4];
	if (leftOn(pin1, pi, pi1))
	{
		return leftOn(pi, pj, pin1) && leftOn(pj, pi, pi1);
	}
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}

inline bool diagonalLoose(int i, int j, int n, const int* verts, int* indices)
{
	return inConeLoose(i, j, n, verts, indices) && diagonalieLoose(i, j, n, verts, indices);
}

int triangulate(int n, const int* verts, int* indices, int* tris, const bool allowLoose)
{
	int ntris = 0;
	int* dst = tris;

	for (int i = 0; i < n; i++)
	{
		const int i1 = nextIndex(i, n);
		const int i2 = nextIndex(i1, n);
		if (diagonal(i, i2, n, verts, indices))
		{
			indices[i1] |= DIAGONAL_FLAG;
		}
	}

	while (n > 3)
	{
		int minLen = -1;
		int mini = -1;
		for (int i = 0; i < n; i++)
		{
			const int i1 = nextIndex(i, n);
			if (indices[i1] & DIAGONAL_FLAG)
			{
				const int* p0 = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
				const int* p2 = &verts[(indices[nextIndex(i1, n)] & DIAGONAL_INDEX_MASK) * 4];
				const int dx = p2[0] - p0[0];
				const int dy = p2[2] - p0[2];
				const int len = dx * dx + dy * dy;
				if (minLen < 0 || len < minLen)
				{
					minLen = len;
					mini = i;
				}
			}
		}

		if (mini == -1)
		{
			if (allowLoose)
			{
				minLen = -1;
				mini = -1;
				for (int i = 0; i < n; i++)
				{
					const int i1 = nextIndex(i, n);
					const int i2 = nextIndex(i1, n);
					if (diagonalLoose(i, i2, n, verts, indices))
					{
						const int* p0 = &verts[(indices[i] & DIAGONAL_INDEX_MASK) * 4];
						const int* p2 = &verts[(indices[nextIndex(i2, n)] & DIAGONAL_INDEX_MASK) * 4];
						const int dx = p2[0] - p0[0];
						const int dy = p2[2] - p0[2];
						const int len = dx * dx + dy * dy;
						if (minLen < 0 || len < minLen)
						{
							minLen = len;
							mini = i;
						}
					}
				}
			}
			if (mini == -1)
			{
				return -ntris;
			}
		}

		const int i = mini;
		int i1 = nextIndex(i, n);
		const int i2 = nextIndex(i1, n);

		*dst++ = indices[i] & DIAGONAL_INDEX_MASK;
		*dst++ = indices[i1] & DIAGONAL_INDEX_MASK;
		*dst++ = indices[i2] & DIAGONAL_INDEX_MASK;
		ntris++;

		n--;
		for (int k = i1; k < n; k++)
		{
			indices[k] = indices[k + 1];
		}

		if (i1 >= n)
		{
			i1 = 0;
		}
		int prev = prevIndex(i1, n);
		if (diagonal(prevIndex(prev, n), i1, n, verts, indices))
		{
			indices[prev] |= DIAGONAL_FLAG;
		}
		else
		{
			indices[prev] &= DIAGONAL_INDEX_MASK;
		}

		if (diagonal(prev, nextIndex(i1, n), n, verts, indices))
		{
			indices[i1] |= DIAGONAL_FLAG;
		}
		else
		{
			indices[i1] &= DIAGONAL_INDEX_MASK;
		}
	}

	*dst++ = indices[0] & DIAGONAL_INDEX_MASK;
	*dst++ = indices[1] & DIAGONAL_INDEX_MASK;
	*dst++ = indices[2] & DIAGONAL_INDEX_MASK;
	ntris++;

	return ntris;
}

inline int computeVertexHash(int x, int y, int z, const int bucketCount)
{
	const unsigned int h1 = 0x8da6b343;
	const unsigned int h2 = 0xd8163841;
	const unsigned int h3 = 0xcb1ab31f;
	const unsigned int n = h1 * (unsigned int)x + h2 * (unsigned int)y + h3 * (unsigned int)z;
	return (int)(n & (unsigned int)(bucketCount - 1));
}

unsigned short addVertex(unsigned short x, unsigned short y, unsigned short z,
	unsigned short* verts, int* firstVert, int* nextVert, int& nv, const int bucketCount)
{
	const int bucket = computeVertexHash(x, 0, z, bucketCount);
	int i = firstVert[bucket];
	while (i != -1)
	{
		const unsigned short* v = &verts[i * 3];
		if (v[0] == x && abs((int)v[1] - (int)y) <= 2 && v[2] == z)
		{
			return (unsigned short)i;
		}
		i = nextVert[i];
	}

	i = nv;
	nv++;
	unsigned short* v = &verts[i * 3];
	v[0] = x;
	v[1] = y;
	v[2] = z;
	nextVert[i] = firstVert[bucket];
	firstVert[bucket] = i;
	return (unsigned short)i;
}

inline int countPolyVerts(const unsigned short* p, const int nvp)
{
	for (int i = 0; i < nvp; ++i)
	{
		if (p[i] == NULL_IDX)
		{
			return i;
		}
	}
	return nvp;
}

inline bool uleft(const unsigned short* a, const unsigned short* b, const unsigned short* c)
{
	return ((int)b[0] - (int)a[0]) * ((int)c[2] - (int)a[2]) -
	       ((int)c[0] - (int)a[0]) * ((int)b[2] - (int)a[2]) < 0;
}

int getPolyMergeValue(unsigned short* pa, unsigned short* pb, const unsigned short* verts,
	int& ea, int& eb, const int nvp)
{
	const int na = countPolyVerts(pa, nvp);
	const int nb = countPolyVerts(pb, nvp);

	if (na + nb - 2 > nvp)
	{
		return -1;
	}

	ea = -1;
	eb = -1;
	for (int i = 0; i < na; ++i)
	{
		unsigned short va0 = pa[i];
		unsigned short va1 = pa[(i + 1) % na];
		if (va0 > va1)
		{
			rcSwap(va0, va1);
		}
		for (int j = 0; j < nb; ++j)
		{
			unsigned short vb0 = pb[j];
			unsigned short vb1 = pb[(j + 1) % nb];
			if (vb0 > vb1)
			{
				rcSwap(vb0, vb1);
			}
			if (va0 == vb0 && va1 == vb1)
			{
				ea = i;
				eb = j;
				break;
			}
		}
	}

	if (ea == -1 || eb == -1)
	{
		return -1;
	}

	unsigned short va, vb, vc;
	va = pa[(ea + na - 1) % na];
	vb = pa[ea];
	vc = pb[(eb + 2) % nb];
	if (!uleft(&verts[va * 3], &verts[vb * 3], &verts[vc * 3]))
	{
		return -1;
	}

	va = pb[(eb + nb - 1) % nb];
	vb = pb[eb];
	vc = pa[(ea + 2) % na];
	if (!uleft(&verts[va * 3], &verts[vb * 3], &verts[vc * 3]))
	{
		return -1;
	}

	va = pa[ea];
	vb = pa[(ea + 1) % na];
	const int dx = (int)verts[va * 3 + 0] - (int)verts[vb * 3 + 0];
	const int dy = (int)verts[va * 3 + 2] - (int)verts[vb * 3 + 2];
	return dx * dx + dy * dy;
}

void mergePolyVerts(unsigned short* pa, unsigned short* pb, int ea, int eb,
	unsigned short* tmp, const int nvp)
{
	const int na = countPolyVerts(pa, nvp);
	const int nb = countPolyVerts(pb, nvp);

	memset(tmp, 0xff, sizeof(unsigned short) * nvp);
	int n = 0;
	for (int i = 0; i < na - 1; ++i)
	{
		tmp[n++] = pa[(ea + 1 + i) % na];
	}
	for (int i = 0; i < nb - 1; ++i)
	{
		tmp[n++] = pb[(eb + 1 + i) % nb];
	}
	memcpy(pa, tmp, sizeof(unsigned short) * nvp);
}

void pushFront(int v, int* arr, int& an)
{
	an++;
	for (int i = an - 1; i > 0; --i)
	{
		arr[i] = arr[i - 1];
	}
	arr[0] = v;
}

void pushBack(int v, int* arr, int& an)
{
	arr[an] = v;
	an++;
}

bool canRemoveVertex(const unsigned short* polys, const int npolys, const int nvp,
	const unsigned short rem, const int maxRemovalEdges, const Allocator& alloc)
{
	int numTouchedVerts = 0;
	int numRemainingEdges = 0;
	for (int i = 0; i < npolys; ++i)
	{
		const unsigned short* p = &polys[i * nvp * 2];
		const int nv = countPolyVerts(p, nvp);
		int numRemoved = 0;
		int numVerts = 0;
		for (int j = 0; j < nv; ++j)
		{
			if (p[j] == rem)
			{
				numTouchedVerts++;
				numRemoved++;
			}
			numVerts++;
		}
		if (numRemoved)
		{
			numRemainingEdges += numVerts - (numRemoved + 1);
		}
	}

	if (numRemainingEdges <= 2)
	{
		return false;
	}

	const int maxEdges = numTouchedVerts * 2;
	if (maxRemovalEdges > 0 && maxEdges > maxRemovalEdges)
	{
		return false;
	}
	if (maxEdges == 0)
	{
		return true;
	}

	int* edges = (int*)alloc.alloc(sizeof(int) * maxEdges * 3);
	if (!edges)
	{
		return false;
	}
	int nedges = 0;

	for (int i = 0; i < npolys; ++i)
	{
		const unsigned short* p = &polys[i * nvp * 2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0, k = nv - 1; j < nv; k = j++)
		{
			if (p[j] == rem || p[k] == rem)
			{
				int a = p[j];
				int b = p[k];
				if (b == rem)
				{
					rcSwap(a, b);
				}
				bool exists = false;
				for (int m = 0; m < nedges; ++m)
				{
					int* e = &edges[m * 3];
					if (e[1] == b)
					{
						e[2]++;
						exists = true;
					}
				}
				if (!exists)
				{
					int* e = &edges[nedges * 3];
					e[0] = a;
					e[1] = b;
					e[2] = 1;
					nedges++;
				}
			}
		}
	}

	int numOpenEdges = 0;
	for (int i = 0; i < nedges; ++i)
	{
		if (edges[i * 3 + 2] < 2)
		{
			numOpenEdges++;
		}
	}

	alloc.free(edges);
	return numOpenEdges <= 2;
}

/// Removes a vertex and retriangulates the hole.
/// @return 0 on success, 1 on out of memory, 2 when the removal exceeded its edge budget.
int removeVertex(unsigned short* verts, int& nverts,
	unsigned short* polys, int& npolys, unsigned short* regs, unsigned char* areas,
	const int nvp, const unsigned short rem, const int maxPolys,
	const int maxRemovalEdges, const bool allowLoose, const Allocator& alloc)
{
	int numRemovedVerts = 0;
	for (int i = 0; i < npolys; ++i)
	{
		const unsigned short* p = &polys[i * nvp * 2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0; j < nv; ++j)
		{
			if (p[j] == rem)
			{
				numRemovedVerts++;
			}
		}
	}

	const int maxEdges = numRemovedVerts * nvp;

	int* edges = (int*)alloc.alloc(sizeof(int) * (maxEdges * 4 + 1));
	int* hole = (int*)alloc.alloc(sizeof(int) * (maxEdges + 1));
	int* hreg = (int*)alloc.alloc(sizeof(int) * (maxEdges + 1));
	int* harea = (int*)alloc.alloc(sizeof(int) * (maxEdges + 1));
	if (!edges || !hole || !hreg || !harea)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		return 1;
	}

	int nedges = 0;
	int nhole = 0;
	int nhreg = 0;
	int nharea = 0;

	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i * nvp * 2];
		const int nv = countPolyVerts(p, nvp);
		bool hasRem = false;
		for (int j = 0; j < nv; ++j)
		{
			if (p[j] == rem)
			{
				hasRem = true;
			}
		}
		if (hasRem)
		{
			for (int j = 0, k = nv - 1; j < nv; k = j++)
			{
				if (p[j] != rem && p[k] != rem)
				{
					if (maxRemovalEdges > 0 && nedges >= maxRemovalEdges)
					{
						alloc.free(edges);
						alloc.free(hole);
						alloc.free(hreg);
						alloc.free(harea);
						return 2;
					}
					int* e = &edges[nedges * 4];
					e[0] = p[k];
					e[1] = p[j];
					e[2] = regs ? regs[i] : 0;
					e[3] = areas ? areas[i] : 0;
					nedges++;
				}
			}
			unsigned short* p2 = &polys[(npolys - 1) * nvp * 2];
			if (p != p2)
			{
				memcpy(p, p2, sizeof(unsigned short) * nvp);
			}
			memset(p + nvp, 0xff, sizeof(unsigned short) * nvp);
			if (regs)
			{
				regs[i] = regs[npolys - 1];
			}
			if (areas)
			{
				areas[i] = areas[npolys - 1];
			}
			npolys--;
			--i;
		}
	}

	for (int i = (int)rem; i < nverts - 1; ++i)
	{
		verts[i * 3 + 0] = verts[(i + 1) * 3 + 0];
		verts[i * 3 + 1] = verts[(i + 1) * 3 + 1];
		verts[i * 3 + 2] = verts[(i + 1) * 3 + 2];
	}
	nverts--;

	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i * nvp * 2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0; j < nv; ++j)
		{
			if (p[j] > rem)
			{
				p[j]--;
			}
		}
	}
	for (int i = 0; i < nedges; ++i)
	{
		if (edges[i * 4 + 0] > (int)rem)
		{
			edges[i * 4 + 0]--;
		}
		if (edges[i * 4 + 1] > (int)rem)
		{
			edges[i * 4 + 1]--;
		}
	}

	if (nedges == 0)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		return 0;
	}

	pushBack(edges[0], hole, nhole);
	pushBack(edges[2], hreg, nhreg);
	pushBack(edges[3], harea, nharea);

	while (nedges)
	{
		bool match = false;
		for (int i = 0; i < nedges; ++i)
		{
			const int ea = edges[i * 4 + 0];
			const int eb = edges[i * 4 + 1];
			const int r = edges[i * 4 + 2];
			const int a = edges[i * 4 + 3];
			bool add = false;
			if (hole[0] == eb)
			{
				if (maxRemovalEdges > 0 && nhole >= maxRemovalEdges)
				{
					alloc.free(edges);
					alloc.free(hole);
					alloc.free(hreg);
					alloc.free(harea);
					return 2;
				}
				pushFront(ea, hole, nhole);
				pushFront(r, hreg, nhreg);
				pushFront(a, harea, nharea);
				add = true;
			}
			else if (hole[nhole - 1] == ea)
			{
				if (maxRemovalEdges > 0 && nhole >= maxRemovalEdges)
				{
					alloc.free(edges);
					alloc.free(hole);
					alloc.free(hreg);
					alloc.free(harea);
					return 2;
				}
				pushBack(eb, hole, nhole);
				pushBack(r, hreg, nhreg);
				pushBack(a, harea, nharea);
				add = true;
			}
			if (add)
			{
				edges[i * 4 + 0] = edges[(nedges - 1) * 4 + 0];
				edges[i * 4 + 1] = edges[(nedges - 1) * 4 + 1];
				edges[i * 4 + 2] = edges[(nedges - 1) * 4 + 2];
				edges[i * 4 + 3] = edges[(nedges - 1) * 4 + 3];
				--nedges;
				match = true;
				--i;
			}
		}
		if (!match)
		{
			break;
		}
	}

	int* tris = (int*)alloc.alloc(sizeof(int) * (nhole * 3 + 1));
	int* tverts = (int*)alloc.alloc(sizeof(int) * (nhole * 4 + 1));
	int* thole = (int*)alloc.alloc(sizeof(int) * (nhole + 1));
	if (!tris || !tverts || !thole)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		alloc.free(tris);
		alloc.free(tverts);
		alloc.free(thole);
		return 1;
	}

	for (int i = 0; i < nhole; ++i)
	{
		const int pi = hole[i];
		tverts[i * 4 + 0] = verts[pi * 3 + 0];
		tverts[i * 4 + 1] = verts[pi * 3 + 1];
		tverts[i * 4 + 2] = verts[pi * 3 + 2];
		tverts[i * 4 + 3] = 0;
		thole[i] = i;
	}

	int ntris = triangulate(nhole, tverts, thole, tris, allowLoose);
	if (ntris < 0)
	{
		ntris = -ntris;
	}

	if (maxRemovalEdges > 0 && ntris > maxRemovalEdges)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		alloc.free(tris);
		alloc.free(tverts);
		alloc.free(thole);
		return 2;
	}

	unsigned short* newPolys = (unsigned short*)alloc.alloc(sizeof(unsigned short) * (ntris + 1) * nvp);
	unsigned short* pregs = (unsigned short*)alloc.alloc(sizeof(unsigned short) * (ntris + 1));
	unsigned char* pareas = (unsigned char*)alloc.alloc(sizeof(unsigned char) * (ntris + 1));
	if (!newPolys || !pregs || !pareas)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		alloc.free(tris);
		alloc.free(tverts);
		alloc.free(thole);
		alloc.free(newPolys);
		alloc.free(pregs);
		alloc.free(pareas);
		return 1;
	}

	unsigned short* tmpPoly = &newPolys[ntris * nvp];

	int npolysNew = 0;
	memset(newPolys, 0xff, ntris * nvp * sizeof(unsigned short));
	for (int j = 0; j < ntris; ++j)
	{
		const int* t = &tris[j * 3];
		if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2])
		{
			newPolys[npolysNew * nvp + 0] = (unsigned short)hole[t[0]];
			newPolys[npolysNew * nvp + 1] = (unsigned short)hole[t[1]];
			newPolys[npolysNew * nvp + 2] = (unsigned short)hole[t[2]];

			if (regs)
			{
				if (hreg[t[0]] != hreg[t[1]] || hreg[t[1]] != hreg[t[2]])
				{
					pregs[npolysNew] = RC_MULTIPLE_REGS;
				}
				else
				{
					pregs[npolysNew] = (unsigned short)hreg[t[0]];
				}
			}
			else
			{
				pregs[npolysNew] = 0;
			}
			pareas[npolysNew] = (unsigned char)harea[t[0]];
			npolysNew++;
		}
	}
	if (!npolysNew)
	{
		alloc.free(edges);
		alloc.free(hole);
		alloc.free(hreg);
		alloc.free(harea);
		alloc.free(tris);
		alloc.free(tverts);
		alloc.free(thole);
		alloc.free(newPolys);
		alloc.free(pregs);
		alloc.free(pareas);
		return 0;
	}

	if (nvp > 3)
	{
		for (;;)
		{
			int bestMergeVal = 0;
			int bestPa = 0;
			int bestPb = 0;
			int bestEa = 0;
			int bestEb = 0;
			for (int j = 0; j < npolysNew - 1; ++j)
			{
				unsigned short* pj = &newPolys[j * nvp];
				for (int k = j + 1; k < npolysNew; ++k)
				{
					unsigned short* pk = &newPolys[k * nvp];
					int ea, eb;
					const int v = getPolyMergeValue(pj, pk, verts, ea, eb, nvp);
					if (v > bestMergeVal)
					{
						bestMergeVal = v;
						bestPa = j;
						bestPb = k;
						bestEa = ea;
						bestEb = eb;
					}
				}
			}

			if (bestMergeVal > 0)
			{
				unsigned short* pa = &newPolys[bestPa * nvp];
				unsigned short* pb = &newPolys[bestPb * nvp];
				mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp);
				if (regs && pregs[bestPa] != pregs[bestPb])
				{
					pregs[bestPa] = RC_MULTIPLE_REGS;
				}

				unsigned short* last = &newPolys[(npolysNew - 1) * nvp];
				if (pb != last)
				{
					memcpy(pb, last, sizeof(unsigned short) * nvp);
				}
				pregs[bestPb] = pregs[npolysNew - 1];
				pareas[bestPb] = pareas[npolysNew - 1];
				npolysNew--;
			}
			else
			{
				break;
			}
		}
	}

	int status = 0;
	for (int i = 0; i < npolysNew; ++i)
	{
		if (npolys >= maxPolys)
		{
			status = 2;
			break;
		}
		unsigned short* p = &polys[npolys * nvp * 2];
		memset(p, 0xff, sizeof(unsigned short) * nvp * 2);
		for (int j = 0; j < nvp; ++j)
		{
			p[j] = newPolys[i * nvp + j];
		}
		if (regs)
		{
			regs[npolys] = pregs[i];
		}
		if (areas)
		{
			areas[npolys] = pareas[i];
		}
		npolys++;
	}

	alloc.free(edges);
	alloc.free(hole);
	alloc.free(hreg);
	alloc.free(harea);
	alloc.free(tris);
	alloc.free(tverts);
	alloc.free(thole);
	alloc.free(newPolys);
	alloc.free(pregs);
	alloc.free(pareas);
	return status;
}

struct Edge
{
	unsigned short vert[2];
	unsigned short polyEdge[2];
	unsigned short poly[2];
	unsigned short portal;
};

inline bool overlapRangeExl(const int amin, const int amax, const int bmin, const int bmax)
{
	return (amin >= bmax || amax <= bmin) ? false : true;
}

/// The single adjacency implementation used by every caller.
///
/// @param[in,out]	polys	Polygon array; neighbor slots for connected edges are written.
/// @param[in]		npolys	The number of polygons.
/// @param[in]		nverts	The number of vertices.
/// @param[in]		nvp		The maximum number of vertices per polygon.
/// @param[in]		verts	The mesh vertices, 3 values per vertex.
/// @param[in]		portalMode	PORTAL_NONE, PORTAL_BORDER or PORTAL_BANK.
/// @param[in]		conts	Contours (PORTAL_BANK only).
/// @param[in]		nconts	The number of contours.
/// @param[in]		alloc	Scratch allocator.
/// @return True on success.
bool buildAdjacency(unsigned short* polys, const int npolys, const int nverts, const int nvp,
	const unsigned short* verts, const int portalMode,
	const rcPolyMeshBuilderContour* conts, const int nconts, const Allocator& alloc)
{
	const int maxEdgeCount = npolys * nvp;
	if (maxEdgeCount == 0)
	{
		return true;
	}

	unsigned short* firstEdge = (unsigned short*)alloc.alloc(sizeof(unsigned short) * (nverts + maxEdgeCount + 1));
	Edge* edges = (Edge*)alloc.alloc(sizeof(Edge) * (maxEdgeCount + 1));
	if (!firstEdge || !edges)
	{
		alloc.free(firstEdge);
		alloc.free(edges);
		return false;
	}
	unsigned short* nextEdge = firstEdge + nverts;
	int edgeCount = 0;

	for (int i = 0; i < nverts; i++)
	{
		firstEdge[i] = NULL_IDX;
	}

	// Collect half-edges in canonical (v0 < v1) order.
	for (int i = 0; i < npolys; ++i)
	{
		const unsigned short* t = &polys[i * nvp * 2];
		for (int j = 0; j < nvp; ++j)
		{
			if (t[j] == NULL_IDX)
			{
				break;
			}
			const unsigned short v0 = t[j];
			const unsigned short v1 = (j + 1 >= nvp || t[j + 1] == NULL_IDX) ? t[0] : t[j + 1];
			if (v0 < v1)
			{
				Edge& edge = edges[edgeCount];
				edge.vert[0] = v0;
				edge.vert[1] = v1;
				edge.poly[0] = (unsigned short)i;
				edge.polyEdge[0] = (unsigned short)j;
				edge.poly[1] = (unsigned short)i;
				edge.polyEdge[1] = NULL_IDX;
				edge.portal = NULL_IDX;
				nextEdge[edgeCount] = firstEdge[v0];
				firstEdge[v0] = (unsigned short)edgeCount;
				edgeCount++;
			}
		}
	}

	// Match opposite half-edges. DetourTileCache also records unmatched directed edges as
	// open edges so their portal direction can be attached below.
	const bool insertOpenEdges = (portalMode == RC_POLYMESH_PORTAL_BANK);
	for (int i = 0; i < npolys; ++i)
	{
		const unsigned short* t = &polys[i * nvp * 2];
		for (int j = 0; j < nvp; ++j)
		{
			if (t[j] == NULL_IDX)
			{
				break;
			}
			const unsigned short v0 = t[j];
			const unsigned short v1 = (j + 1 >= nvp || t[j + 1] == NULL_IDX) ? t[0] : t[j + 1];
			if (v0 > v1)
			{
				bool found = false;
				for (unsigned short e = firstEdge[v1]; e != NULL_IDX; e = nextEdge[e])
				{
					Edge& edge = edges[e];
					if (edge.vert[1] == v0 && edge.poly[0] == edge.poly[1])
					{
						edge.poly[1] = (unsigned short)i;
						edge.polyEdge[1] = (unsigned short)j;
						found = true;
						break;
					}
				}
				if (!found && insertOpenEdges)
				{
					Edge& edge = edges[edgeCount];
					edge.vert[0] = v1;
					edge.vert[1] = v0;
					edge.poly[0] = (unsigned short)i;
					edge.polyEdge[0] = (unsigned short)j;
					edge.poly[1] = (unsigned short)i;
					edge.polyEdge[1] = NULL_IDX;
					edge.portal = NULL_IDX;
					nextEdge[edgeCount] = firstEdge[v1];
					firstEdge[v1] = (unsigned short)edgeCount;
					edgeCount++;
				}
			}
		}
	}

	// DetourTileCache portal annotation: match each open edge against the contour segment
	// that generated it and reuse polyEdge[1] to carry the portal direction.
	if (portalMode == RC_POLYMESH_PORTAL_BANK && conts)
	{
		for (int i = 0; i < nconts; ++i)
		{
			const rcPolyMeshBuilderContour& cont = conts[i];
			if (cont.nverts < 3)
			{
				continue;
			}
			for (int j = 0, k = cont.nverts - 1; j < cont.nverts; k = j++)
			{
				const int* va = &cont.verts[k * 4];
				const int* vb = &cont.verts[j * 4];
				const int dir = va[3] & 0xf;
				if (dir == 0xf)
				{
					continue;
				}

				if (dir == 0 || dir == 2)
				{
					const int x = va[0];
					int zmin = va[2];
					int zmax = vb[2];
					if (zmin > zmax)
					{
						rcSwap(zmin, zmax);
					}
					for (int m = 0; m < edgeCount; ++m)
					{
						Edge& e = edges[m];
						if (e.poly[0] != e.poly[1])
						{
							continue;
						}
						const unsigned short* eva = &verts[e.vert[0] * 3];
						const unsigned short* evb = &verts[e.vert[1] * 3];
						if ((int)eva[0] == x && (int)evb[0] == x)
						{
							int ezmin = eva[2];
							int ezmax = evb[2];
							if (ezmin > ezmax)
							{
								rcSwap(ezmin, ezmax);
							}
							if (overlapRangeExl(zmin, zmax, ezmin, ezmax))
							{
								e.portal = (unsigned short)dir;
							}
						}
					}
				}
				else
				{
					const int z = va[2];
					int xmin = va[0];
					int xmax = vb[0];
					if (xmin > xmax)
					{
						rcSwap(xmin, xmax);
					}
					for (int m = 0; m < edgeCount; ++m)
					{
						Edge& e = edges[m];
						if (e.poly[0] != e.poly[1])
						{
							continue;
						}
						const unsigned short* eva = &verts[e.vert[0] * 3];
						const unsigned short* evb = &verts[e.vert[1] * 3];
						if ((int)eva[2] == z && (int)evb[2] == z)
						{
							int exmin = eva[0];
							int exmax = evb[0];
							if (exmin > exmax)
							{
								rcSwap(exmin, exmax);
							}
							if (overlapRangeExl(xmin, xmax, exmin, exmax))
							{
								e.portal = (unsigned short)dir;
							}
						}
					}
				}
			}
		}
	}

	// Store adjacency.
	for (int i = 0; i < edgeCount; ++i)
	{
		const Edge& e = edges[i];
		if (e.poly[0] != e.poly[1])
		{
			unsigned short* p0 = &polys[e.poly[0] * nvp * 2];
			unsigned short* p1 = &polys[e.poly[1] * nvp * 2];
			p0[nvp + e.polyEdge[0]] = e.poly[1];
			p1[nvp + e.polyEdge[1]] = e.poly[0];
		}
		else if (e.portal != NULL_IDX)
		{
			unsigned short* p0 = &polys[e.poly[0] * nvp * 2];
			p0[nvp + e.polyEdge[0]] = (unsigned short)(0x8000 | e.portal);
		}
	}

	alloc.free(firstEdge);
	alloc.free(edges);
	return true;
}

/// Recast portal annotation: an unconnected polygon edge whose endpoints both lie on the
/// build border is marked with the border direction.
void markBorderPortals(unsigned short* polys, const int npolys, const int nvp,
	const unsigned short* verts, const int width, const int height)
{
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i * 2 * nvp];
		for (int j = 0; j < nvp; ++j)
		{
			if (p[j] == NULL_IDX)
			{
				break;
			}
			if (p[nvp + j] != NULL_IDX)
			{
				continue;
			}
			int nj = j + 1;
			if (nj >= nvp || p[nj] == NULL_IDX)
			{
				nj = 0;
			}
			const unsigned short* va = &verts[p[j] * 3];
			const unsigned short* vb = &verts[p[nj] * 3];

			if ((int)va[0] == 0 && (int)vb[0] == 0)
			{
				p[nvp + j] = 0x8000 | 0;
			}
			else if ((int)va[2] == height && (int)vb[2] == height)
			{
				p[nvp + j] = 0x8000 | 1;
			}
			else if ((int)va[0] == width && (int)vb[0] == width)
			{
				p[nvp + j] = 0x8000 | 2;
			}
			else if ((int)va[2] == 0 && (int)vb[2] == 0)
			{
				p[nvp + j] = 0x8000 | 3;
			}
		}
	}
}

}  // namespace

bool rcBuildPolyMeshFromContours(
	rcContext* ctx,
	const rcPolyMeshBuilderParams& params,
	rcPolyMeshBuilderResult& result)
{
	result.nverts = 0;
	result.npolys = 0;

	if (!params.outVerts || !params.outPolys || !params.outAreas || params.nvp <= 0)
	{
		if (ctx)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Missing output buffers.");
		}
		return false;
	}

	Allocator alloc;
	if (params.allocator.alloc)
	{
		alloc.allocFn = params.allocator.alloc;
		alloc.freeFn = params.allocator.free;
		alloc.user = params.allocator.user;
	}

	int maxVertices = 0;
	int maxTris = 0;
	int maxVertsPerCont = 0;
	for (int i = 0; i < params.nconts; ++i)
	{
		const rcPolyMeshBuilderContour& cont = params.conts[i];
		if (cont.nverts < 3)
		{
			continue;
		}
		maxVertices += cont.nverts;
		maxTris += cont.nverts - 2;
		if (cont.nverts > maxVertsPerCont)
		{
			maxVertsPerCont = cont.nverts;
		}
	}

	if (maxVertices == 0)
	{
		return true;
	}

	unsigned short* meshVerts = params.outVerts;
	unsigned short* meshPolys = params.outPolys;
	unsigned short* meshRegs = params.outRegs;
	unsigned char* meshAreas = params.outAreas;
	const int nvp = params.nvp;
	const int bucketCount = params.vertexBucketCount > 0 ? params.vertexBucketCount : 4096;

	memset(meshPolys, 0xff, sizeof(unsigned short) * params.maxPolys * nvp * 2);
	memset(meshVerts, 0, sizeof(unsigned short) * params.maxVerts * 3);
	if (meshRegs)
	{
		memset(meshRegs, 0, sizeof(unsigned short) * params.maxPolys);
	}
	memset(meshAreas, 0, sizeof(unsigned char) * params.maxPolys);

	int nverts = 0;
	int npolys = 0;

	unsigned char* vflags = (unsigned char*)alloc.alloc(sizeof(unsigned char) * (maxVertices + 1));
	int* nextVert = (int*)alloc.alloc(sizeof(int) * (maxVertices + 1));
	int* firstVert = (int*)alloc.alloc(sizeof(int) * bucketCount);
	int* indices = (int*)alloc.alloc(sizeof(int) * (maxVertsPerCont + 1));
	int* tris = (int*)alloc.alloc(sizeof(int) * (maxVertsPerCont * 3 + 1));
	unsigned short* polys = (unsigned short*)alloc.alloc(sizeof(unsigned short) * (maxVertsPerCont + 1) * nvp);
	if (!vflags || !nextVert || !firstVert || !indices || !tris || !polys)
	{
		alloc.free(vflags);
		alloc.free(nextVert);
		alloc.free(firstVert);
		alloc.free(indices);
		alloc.free(tris);
		alloc.free(polys);
		if (ctx)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Out of memory.");
		}
		return false;
	}
	memset(vflags, 0, maxVertices);
	memset(nextVert, 0, sizeof(int) * maxVertices);
	for (int i = 0; i < bucketCount; ++i)
	{
		firstVert[i] = -1;
	}
	unsigned short* tmpPoly = &polys[maxVertsPerCont * nvp];

	bool ok = true;

	for (int i = 0; i < params.nconts && ok; ++i)
	{
		const rcPolyMeshBuilderContour& cont = params.conts[i];
		if (cont.nverts < 3)
		{
			continue;
		}

		for (int j = 0; j < cont.nverts; ++j)
		{
			indices[j] = j;
		}

		int ntris = triangulate(cont.nverts, cont.verts, indices, tris, params.allowLooseTriangulation);
		if (ntris <= 0)
		{
			if (ctx)
			{
				ctx->log(RC_LOG_WARNING, "rcBuildPolyMeshFromContours: Bad triangulation Contour %d.", i);
			}
			ntris = -ntris;
		}

		for (int j = 0; j < cont.nverts; ++j)
		{
			const int* v = &cont.verts[j * 4];
			indices[j] = addVertex((unsigned short)v[0], (unsigned short)v[1], (unsigned short)v[2],
				meshVerts, firstVert, nextVert, nverts, bucketCount);
			if (params.borderVertexFlag && (v[3] & params.borderVertexFlag))
			{
				vflags[indices[j]] = 1;
			}
		}

		int contPolys = 0;
		memset(polys, 0xff, maxVertsPerCont * nvp * sizeof(unsigned short));
		for (int j = 0; j < ntris; ++j)
		{
			const int* t = &tris[j * 3];
			if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2])
			{
				polys[contPolys * nvp + 0] = (unsigned short)indices[t[0]];
				polys[contPolys * nvp + 1] = (unsigned short)indices[t[1]];
				polys[contPolys * nvp + 2] = (unsigned short)indices[t[2]];
				contPolys++;
			}
		}
		if (!contPolys)
		{
			continue;
		}

		if (nvp > 3)
		{
			for (;;)
			{
				int bestMergeVal = 0;
				int bestPa = 0;
				int bestPb = 0;
				int bestEa = 0;
				int bestEb = 0;
				for (int j = 0; j < contPolys - 1; ++j)
				{
					unsigned short* pj = &polys[j * nvp];
					for (int k = j + 1; k < contPolys; ++k)
					{
						unsigned short* pk = &polys[k * nvp];
						int ea, eb;
						const int v = getPolyMergeValue(pj, pk, meshVerts, ea, eb, nvp);
						if (v > bestMergeVal)
						{
							bestMergeVal = v;
							bestPa = j;
							bestPb = k;
							bestEa = ea;
							bestEb = eb;
						}
					}
				}

				if (bestMergeVal > 0)
				{
					unsigned short* pa = &polys[bestPa * nvp];
					unsigned short* pb = &polys[bestPb * nvp];
					mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp);
					unsigned short* lastPoly = &polys[(contPolys - 1) * nvp];
					if (pb != lastPoly)
					{
						memcpy(pb, lastPoly, sizeof(unsigned short) * nvp);
					}
					contPolys--;
				}
				else
				{
					break;
				}
			}
		}

		for (int j = 0; j < contPolys; ++j)
		{
			if (npolys >= params.maxPolys)
			{
				ok = false;
				break;
			}
			unsigned short* p = &meshPolys[npolys * nvp * 2];
			const unsigned short* q = &polys[j * nvp];
			for (int k = 0; k < nvp; ++k)
			{
				p[k] = q[k];
			}
			if (meshRegs)
			{
				meshRegs[npolys] = (unsigned short)cont.reg;
			}
			meshAreas[npolys] = cont.area;
			npolys++;
		}
	}

	if (!ok)
	{
		alloc.free(vflags);
		alloc.free(nextVert);
		alloc.free(firstVert);
		alloc.free(indices);
		alloc.free(tris);
		alloc.free(polys);
		if (ctx)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Too many polygons.");
		}
		return false;
	}

	// Remove border vertices that are safe to drop.
	for (int i = 0; i < nverts; ++i)
	{
		if (!vflags[i])
		{
			continue;
		}
		if (!canRemoveVertex(meshPolys, npolys, nvp, (unsigned short)i, params.maxRemovalEdges, alloc))
		{
			continue;
		}
		const int status = removeVertex(meshVerts, nverts, meshPolys, npolys, meshRegs, meshAreas,
			nvp, (unsigned short)i, params.maxPolys, params.maxRemovalEdges,
			params.allowLooseTriangulation, alloc);
		if (status == 1)
		{
			alloc.free(vflags);
			alloc.free(nextVert);
			alloc.free(firstVert);
			alloc.free(indices);
			alloc.free(tris);
			alloc.free(polys);
			if (ctx)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Out of memory removing vertex.");
			}
			return false;
		}
		if (status == 2)
		{
			alloc.free(vflags);
			alloc.free(nextVert);
			alloc.free(firstVert);
			alloc.free(indices);
			alloc.free(tris);
			alloc.free(polys);
			if (ctx)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Vertex removal exceeded its edge budget.");
			}
			return false;
		}
		for (int j = i; j < nverts; ++j)
		{
			vflags[j] = vflags[j + 1];
		}
		--i;
	}

	const bool adjacencyOk = buildAdjacency(meshPolys, npolys, nverts, nvp, meshVerts,
		params.portalMode, params.conts, params.nconts, alloc);
	if (adjacencyOk && params.portalMode == RC_POLYMESH_PORTAL_BORDER)
	{
		markBorderPortals(meshPolys, npolys, nvp, meshVerts, params.width, params.height);
	}

	alloc.free(vflags);
	alloc.free(nextVert);
	alloc.free(firstVert);
	alloc.free(indices);
	alloc.free(tris);
	alloc.free(polys);

	if (!adjacencyOk)
	{
		if (ctx)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildPolyMeshFromContours: Adjacency failed.");
		}
		return false;
	}

	result.nverts = nverts;
	result.npolys = npolys;
	return true;
}

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
	const rcPolyMeshBuilderAllocator* allocator)
{
	Allocator alloc;
	if (allocator && allocator->alloc)
	{
		alloc.allocFn = allocator->alloc;
		alloc.freeFn = allocator->free;
		alloc.user = allocator->user;
	}

	const bool ok = buildAdjacency(polys, npolys, nverts, nvp, verts, portalMode,
		conts, nconts, alloc);
	if (ok && portalMode == RC_POLYMESH_PORTAL_BORDER)
	{
		markBorderPortals(polys, npolys, nvp, verts, width, height);
	}
	return ok;
}
