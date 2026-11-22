// Shared, dependency-free contour-set -> polygon-mesh implementation.
// Used by both Recast (rcBuildPolyMesh) and the Detour tile cache
// (dtBuildTileCachePolyMesh) so the two never disagree about topology.
#ifndef RECAST_MESH_CORE_H
#define RECAST_MESH_CORE_H

#include <stdlib.h>
#include <string.h>

namespace rcMeshCore
{

static const unsigned short NULL_IDX = 0xffff;

// Vertex dedup hash buckets used while merging contour vertices. Kept small
// so a tile-cache rebuild stays within its fixed per-tile scratch budget.
static const int VERTEX_BUCKET_COUNT = (1<<8);

// Allocator hooks. temp* is scratch (freed before return); perm* owns output.
struct Alloc
{
	void* (*tempAlloc)(size_t size, void* user);
	void (*tempFree)(void* ptr, void* user);
	void* (*permAlloc)(size_t size, void* user);
	void (*permFree)(void* ptr, void* user);
	void* user;
	void* t(size_t s) const { return tempAlloc(s, user); }
	void tf(void* p) const { if (p) tempFree(p, user); }
	void* p(size_t s) const { return permAlloc(s, user); }
	void pf(void* p) const { if (p) permFree(p, user); }
};

// One input contour: verts are 4 ints (x,y,z,flags) per vertex.
struct Contour
{
	const int* verts;             // used when vertsByte is null
	const unsigned char* vertsByte;  // alternative byte-packed layout
	int nverts;
	unsigned short reg;
	unsigned char area;

	inline int vx(int i) const { return vertsByte ? (int)vertsByte[i*4+0] : verts[i*4+0]; }
	inline int vy(int i) const { return vertsByte ? (int)vertsByte[i*4+1] : verts[i*4+1]; }
	inline int vz(int i) const { return vertsByte ? (int)vertsByte[i*4+2] : verts[i*4+2]; }
	inline int vflags(int i) const { return vertsByte ? (int)vertsByte[i*4+3] : verts[i*4+3]; }
};

enum PortalMode
{
	PORTAL_NONE,        // no portal edges are marked (offline solo mesh)
	PORTAL_COORD,       // mark open edges on the outer cell boundary (tiled Recast)
	PORTAL_CONTOUR_DIR  // use per-vertex direction codes from the contour (tile cache)
};

struct Input
{
	const Contour* conts;
	int nconts;
	int borderFlag;    // bit in a contour vertex's 4th component marking a border vertex
	int width;         // cell width, used by PORTAL_COORD
	int height;        // cell height, used by PORTAL_COORD
	PortalMode portal;
};

struct Output
{
	unsigned short* verts;  // 3 per vertex
	unsigned short* polys;  // 2*nvp per polygon
	unsigned short* regs;   // optional (may be null)
	unsigned char* areas;   // optional (may be null)
	int nverts;
	int npolys;
	int nvp;
	int maxpolys;
};

inline int computeVertexHash(int x, int y, int z)
{
	const unsigned int h1 = 0x8da6b343;
	const unsigned int h2 = 0xd8163841;
	const unsigned int h3 = 0xcb1ab31f;
	unsigned int n = h1 * x + h2 * y + h3 * z;
	return (int)(n & (VERTEX_BUCKET_COUNT-1));
}

inline unsigned short addVertex(unsigned short x, unsigned short y, unsigned short z,
								unsigned short* verts, int* firstVert, int* nextVert, int& nv)
{
	int bucket = computeVertexHash(x, 0, z);
	int i = firstVert[bucket];
	while (i != -1)
	{
		const unsigned short* v = &verts[i*3];
		if (v[0] == x && (int(v[1]) - int(y) <= 2) && (int(y) - int(v[1]) <= 2) && v[2] == z)
			return (unsigned short)i;
		i = nextVert[i];
	}
	i = nv; nv++;
	unsigned short* v = &verts[i*3];
	v[0] = x; v[1] = y; v[2] = z;
	nextVert[i] = firstVert[bucket];
	firstVert[bucket] = i;
	return (unsigned short)i;
}

inline int prevIdx(int i, int n) { return i-1 >= 0 ? i-1 : n-1; }
inline int nextIdx(int i, int n) { return i+1 < n ? i+1 : 0; }

inline int area2(const int* a, const int* b, const int* c)
{
	return (b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2]);
}
inline bool leftOf(const int* a, const int* b, const int* c) { return area2(a, b, c) < 0; }
inline bool leftOn(const int* a, const int* b, const int* c) { return area2(a, b, c) <= 0; }
inline bool collinear(const int* a, const int* b, const int* c) { return area2(a, b, c) == 0; }

inline bool intersectProp(const int* a, const int* b, const int* c, const int* d)
{
	if (collinear(a,b,c) || collinear(a,b,d) || collinear(c,d,a) || collinear(c,d,b))
		return false;
	return (leftOf(a,b,c) != leftOf(a,b,d)) && (leftOf(c,d,a) != leftOf(c,d,b));
}
inline bool between(const int* a, const int* b, const int* c)
{
	if (!collinear(a, b, c)) return false;
	if (a[0] != b[0])
		return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
	return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
}
inline bool intersectSeg(const int* a, const int* b, const int* c, const int* d)
{
	if (intersectProp(a, b, c, d)) return true;
	if (between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b)) return true;
	return false;
}
inline bool vequal(const int* a, const int* b) { return a[0] == b[0] && a[2] == b[2]; }

inline bool diagonalie(int i, int j, int n, const int* verts, int* indices)
{
	const int* d0 = &verts[(indices[i] & 0x0fffffff) * 4];
	const int* d1 = &verts[(indices[j] & 0x0fffffff) * 4];
	for (int k = 0; k < n; k++)
	{
		int k1 = nextIdx(k, n);
		if (!((k == i) || (k1 == i) || (k == j) || (k1 == j)))
		{
			const int* p0 = &verts[(indices[k] & 0x0fffffff) * 4];
			const int* p1 = &verts[(indices[k1] & 0x0fffffff) * 4];
			if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1)) continue;
			if (intersectSeg(d0, d1, p0, p1)) return false;
		}
	}
	return true;
}
inline bool inCone(int i, int j, int n, const int* verts, int* indices)
{
	const int* pi = &verts[(indices[i] & 0x0fffffff) * 4];
	const int* pj = &verts[(indices[j] & 0x0fffffff) * 4];
	const int* pi1 = &verts[(indices[nextIdx(i, n)] & 0x0fffffff) * 4];
	const int* pin1 = &verts[(indices[prevIdx(i, n)] & 0x0fffffff) * 4];
	if (leftOn(pin1, pi, pi1))
		return leftOf(pi, pj, pin1) && leftOf(pj, pi, pi1);
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}
inline bool diagonal(int i, int j, int n, const int* verts, int* indices)
{
	return inCone(i, j, n, verts, indices) && diagonalie(i, j, n, verts, indices);
}
inline bool diagonalieLoose(int i, int j, int n, const int* verts, int* indices)
{
	const int* d0 = &verts[(indices[i] & 0x0fffffff) * 4];
	const int* d1 = &verts[(indices[j] & 0x0fffffff) * 4];
	for (int k = 0; k < n; k++)
	{
		int k1 = nextIdx(k, n);
		if (!((k == i) || (k1 == i) || (k == j) || (k1 == j)))
		{
			const int* p0 = &verts[(indices[k] & 0x0fffffff) * 4];
			const int* p1 = &verts[(indices[k1] & 0x0fffffff) * 4];
			if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1)) continue;
			if (intersectProp(d0, d1, p0, p1)) return false;
		}
	}
	return true;
}
inline bool inConeLoose(int i, int j, int n, const int* verts, int* indices)
{
	const int* pi = &verts[(indices[i] & 0x0fffffff) * 4];
	const int* pj = &verts[(indices[j] & 0x0fffffff) * 4];
	const int* pi1 = &verts[(indices[nextIdx(i, n)] & 0x0fffffff) * 4];
	const int* pin1 = &verts[(indices[prevIdx(i, n)] & 0x0fffffff) * 4];
	if (leftOn(pin1, pi, pi1))
		return leftOn(pi, pj, pin1) && leftOn(pj, pi, pi1);
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}
inline bool diagonalLoose(int i, int j, int n, const int* verts, int* indices)
{
	return inConeLoose(i, j, n, verts, indices) && diagonalieLoose(i, j, n, verts, indices);
}

inline int triangulate(int n, const int* verts, int* indices, int* tris)
{
	int ntris = 0;
	int* dst = tris;
	for (int i = 0; i < n; i++)
	{
		int i1 = nextIdx(i, n);
		int i2 = nextIdx(i1, n);
		if (diagonal(i, i2, n, verts, indices)) indices[i1] |= 0x80000000;
	}
	while (n > 3)
	{
		int minLen = -1;
		int mini = -1;
		for (int i = 0; i < n; i++)
		{
			int i1 = nextIdx(i, n);
			if (indices[i1] & 0x80000000)
			{
				const int* p0 = &verts[(indices[i] & 0x0fffffff) * 4];
				const int* p2 = &verts[(indices[nextIdx(i1, n)] & 0x0fffffff) * 4];
				int dx = p2[0] - p0[0];
				int dy = p2[2] - p0[2];
				int len = dx*dx + dy*dy;
				if (minLen < 0 || len < minLen) { minLen = len; mini = i; }
			}
		}
		if (mini == -1)
		{
			minLen = -1;
			mini = -1;
			for (int i = 0; i < n; i++)
			{
				int i1 = nextIdx(i, n);
				int i2 = nextIdx(i1, n);
				if (diagonalLoose(i, i2, n, verts, indices))
				{
					const int* p0 = &verts[(indices[i] & 0x0fffffff) * 4];
					const int* p2 = &verts[(indices[nextIdx(i2, n)] & 0x0fffffff) * 4];
					int dx = p2[0] - p0[0];
					int dy = p2[2] - p0[2];
					int len = dx*dx + dy*dy;
					if (minLen < 0 || len < minLen) { minLen = len; mini = i; }
				}
			}
			if (mini == -1)
				return -ntris;
		}
		int i = mini;
		int i1 = nextIdx(i, n);
		int i2 = nextIdx(i1, n);
		*dst++ = indices[i] & 0x0fffffff;
		*dst++ = indices[i1] & 0x0fffffff;
		*dst++ = indices[i2] & 0x0fffffff;
		ntris++;
		n--;
		for (int k = i1; k < n; k++) indices[k] = indices[k+1];
		if (i1 >= n) i1 = 0;
		i = prevIdx(i1, n);
		if (diagonal(prevIdx(i, n), i1, n, verts, indices)) indices[i] |= 0x80000000;
		else indices[i] &= 0x0fffffff;
		if (diagonal(i, nextIdx(i1, n), n, verts, indices)) indices[i1] |= 0x80000000;
		else indices[i1] &= 0x0fffffff;
	}
	*dst++ = indices[0] & 0x0fffffff;
	*dst++ = indices[1] & 0x0fffffff;
	*dst++ = indices[2] & 0x0fffffff;
	ntris++;
	return ntris;
}

inline int countPolyVerts(const unsigned short* p, int nvp)
{
	for (int i = 0; i < nvp; ++i)
		if (p[i] == NULL_IDX) return i;
	return nvp;
}
inline bool uleft(const unsigned short* a, const unsigned short* b, const unsigned short* c)
{
	return ((int)b[0] - (int)a[0]) * ((int)c[2] - (int)a[2]) -
		   ((int)c[0] - (int)a[0]) * ((int)b[2] - (int)a[2]) < 0;
}
inline int getPolyMergeValue(unsigned short* pa, unsigned short* pb,
							 const unsigned short* verts, int& ea, int& eb, int nvp)
{
	const int na = countPolyVerts(pa, nvp);
	const int nb = countPolyVerts(pb, nvp);
	if (na+nb-2 > nvp) return -1;
	ea = -1; eb = -1;
	for (int i = 0; i < na; ++i)
	{
		unsigned short va0 = pa[i];
		unsigned short va1 = pa[(i+1) % na];
		if (va0 > va1) { unsigned short tmp = va0; va0 = va1; va1 = tmp; }
		for (int j = 0; j < nb; ++j)
		{
			unsigned short vb0 = pb[j];
			unsigned short vb1 = pb[(j+1) % nb];
			if (vb0 > vb1) { unsigned short tmp = vb0; vb0 = vb1; vb1 = tmp; }
			if (va0 == vb0 && va1 == vb1) { ea = i; eb = j; break; }
		}
	}
	if (ea == -1 || eb == -1) return -1;
	unsigned short va, vb, vc;
	va = pa[(ea+na-1) % na];
	vb = pa[ea];
	vc = pb[(eb+2) % nb];
	if (!uleft(&verts[va*3], &verts[vb*3], &verts[vc*3])) return -1;
	va = pb[(eb+nb-1) % nb];
	vb = pb[eb];
	vc = pa[(ea+2) % na];
	if (!uleft(&verts[va*3], &verts[vb*3], &verts[vc*3])) return -1;
	va = pa[ea];
	vb = pa[(ea+1)%na];
	int dx = (int)verts[va*3+0] - (int)verts[vb*3+0];
	int dy = (int)verts[va*3+2] - (int)verts[vb*3+2];
	return dx*dx + dy*dy;
}
inline void mergePolyVerts(unsigned short* pa, unsigned short* pb, int ea, int eb,
						   unsigned short* tmp, int nvp)
{
	const int na = countPolyVerts(pa, nvp);
	const int nb = countPolyVerts(pb, nvp);
	memset(tmp, 0xff, sizeof(unsigned short)*nvp);
	int n = 0;
	for (int i = 0; i < na-1; ++i) tmp[n++] = pa[(ea+1+i) % na];
	for (int i = 0; i < nb-1; ++i) tmp[n++] = pb[(eb+1+i) % nb];
	memcpy(pa, tmp, sizeof(unsigned short)*nvp);
}
inline void pushFront(int v, int* arr, int& an)
{
	an++;
	for (int i = an-1; i > 0; --i) arr[i] = arr[i-1];
	arr[0] = v;
}
inline void pushBack(int v, int* arr, int& an) { arr[an] = v; an++; }

struct Edge { unsigned short vert[2]; unsigned short polyEdge[2]; unsigned short poly[2]; };

inline bool buildMeshAdjacency(const Alloc& al, const Input& in,
							   unsigned short* polys, int npolys,
							   const unsigned short* verts, int nverts, int nvp)
{
	const int maxEdgeCount = npolys*nvp;
	int* firstEdge = (int*)al.t((size_t)(maxEdgeCount + nverts)*sizeof(int));
	if (!firstEdge) return false;
	int* nextEdge = firstEdge + nverts;
	int edgeCount = 0;
	Edge* edges = (Edge*)al.t(sizeof(Edge)*(size_t)maxEdgeCount);
	if (!edges) { al.tf(firstEdge); return false; }
	unsigned short* openPoly = (unsigned short*)al.t(sizeof(unsigned short)*(size_t)maxEdgeCount);
	unsigned short* openEdge = (unsigned short*)al.t(sizeof(unsigned short)*(size_t)maxEdgeCount);
	unsigned short* openVert0 = (unsigned short*)al.t(sizeof(unsigned short)*(size_t)maxEdgeCount);
	unsigned short* openVert1 = (unsigned short*)al.t(sizeof(unsigned short)*(size_t)maxEdgeCount);
	int openEdgeCount = 0;
	if (!openPoly || !openEdge || !openVert0 || !openVert1)
	{
		al.tf(firstEdge); al.tf(edges); al.tf(openPoly); al.tf(openEdge); al.tf(openVert0); al.tf(openVert1);
		return false;
	}
	for (int i = 0; i < nverts; i++) firstEdge[i] = -1;

	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* t = &polys[i*nvp*2];
		for (int j = 0; j < nvp; ++j)
		{
			if (t[j] == NULL_IDX) break;
			unsigned short v0 = t[j];
			unsigned short v1 = (j+1 >= nvp || t[j+1] == NULL_IDX) ? t[0] : t[j+1];
			if (v0 < v1)
			{
				Edge& edge = edges[edgeCount];
				edge.vert[0] = v0; edge.vert[1] = v1;
				edge.poly[0] = (unsigned short)i; edge.polyEdge[0] = (unsigned short)j;
				edge.poly[1] = (unsigned short)i; edge.polyEdge[1] = 0;
				nextEdge[edgeCount] = firstEdge[v0];
				firstEdge[v0] = edgeCount;
				edgeCount++;
			}
		}
	}
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* t = &polys[i*nvp*2];
		for (int j = 0; j < nvp; ++j)
		{
			if (t[j] == NULL_IDX) break;
			unsigned short v0 = t[j];
			unsigned short v1 = (j+1 >= nvp || t[j+1] == NULL_IDX) ? t[0] : t[j+1];
			if (v0 > v1)
			{
				bool found = false;
				for (int e = firstEdge[v1]; e != -1; e = nextEdge[e])
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
				if (!found)
				{
					// Open edge wound low->high; add the reverse so the edge is
					// represented exactly once regardless of winding direction.
					Edge& edge = edges[edgeCount];
					edge.vert[0] = v1; edge.vert[1] = v0;
					edge.poly[0] = (unsigned short)i; edge.polyEdge[0] = (unsigned short)j;
					edge.poly[1] = (unsigned short)i; edge.polyEdge[1] = 0;
					nextEdge[edgeCount] = firstEdge[v1];
					firstEdge[v1] = edgeCount;
					edgeCount++;
				}
			}
		}
	}

	// Store adjacency for connected polygons. Open edges (which only have a
	// least-index "owner" pass) carry the owning polygon and edge index so the
	// portal pass below can find them regardless of winding direction.
	for (int i = 0; i < edgeCount; ++i)
	{
		const Edge& e = edges[i];
		if (e.poly[0] != e.poly[1])
		{
			unsigned short* p0 = &polys[e.poly[0]*nvp*2];
			unsigned short* p1 = &polys[e.poly[1]*nvp*2];
			p0[nvp + e.polyEdge[0]] = e.poly[1];
			p1[nvp + e.polyEdge[1]] = e.poly[0];
		}
		else
		{
			openPoly[openEdgeCount] = e.poly[0];
			openEdge[openEdgeCount] = e.polyEdge[0];
			openVert0[openEdgeCount] = e.vert[0];
			openVert1[openEdgeCount] = e.vert[1];
			openEdgeCount++;
		}
	}

	if (in.portal == PORTAL_CONTOUR_DIR)
	{
		// Mark open edges that coincide with a contour segment that carries a
		// portal direction. The direction is stored in the segment end vertex.
		// A portal may be traversed in either winding direction relative to the
		// owner polygon, so the segment endpoints are matched ignoring order.
		for (int i = 0; i < in.nconts; ++i)
		{
			const Contour& cont = in.conts[i];
			if (cont.nverts < 3) continue;
			for (int j = 0, k = cont.nverts-1; j < cont.nverts; k=j++)
			{
				const int dir = cont.vflags(k) & 0xf;
				if (dir == 0xf) continue;
				const bool vertical = (dir == 0 || dir == 2);
				const unsigned short x = (unsigned short)(vertical ? cont.vx(k) : cont.vz(k));
				unsigned short amin = (unsigned short)(vertical ? cont.vz(k) : cont.vx(k));
				unsigned short amax = (unsigned short)(vertical ? cont.vz(j) : cont.vx(j));
				if (amin > amax) { unsigned short t = amin; amin = amax; amax = t; }
				for (int e = 0; e < openEdgeCount; ++e)
				{
					const unsigned short* eva = &verts[openVert0[e]*3];
					const unsigned short* evb = &verts[openVert1[e]*3];
					const unsigned short ex = (unsigned short)(vertical ? eva[0] : eva[2]);
					if (ex != x) continue;
					const unsigned short eb = (unsigned short)(vertical ? evb[0] : evb[2]);
					if (eb != x) continue;
					unsigned short emin = (unsigned short)(vertical ? eva[2] : eva[0]);
					unsigned short emax = (unsigned short)(vertical ? evb[2] : evb[0]);
					if (emin > emax) { unsigned short t = emin; emin = emax; emax = t; }
					if (!(amin >= emax || amax <= emin))
					{
						unsigned short* p0 = &polys[openPoly[e]*nvp*2];
						p0[nvp + openEdge[e]] = (unsigned short)(0x8000 | dir);
					}
				}
			}
		}
	}
	else if (in.portal == PORTAL_COORD)
	{
		const int w = in.width, h = in.height;
		for (int i = 0; i < npolys; ++i)
		{
			unsigned short* p = &polys[i*2*nvp];
			for (int j = 0; j < nvp; ++j)
			{
				if (p[j] == NULL_IDX) break;
				if (p[nvp+j] != NULL_IDX) continue;
				int nj = j+1;
				if (nj >= nvp || p[nj] == NULL_IDX) nj = 0;
				const unsigned short* va = &verts[p[j]*3];
				const unsigned short* vb = &verts[p[nj]*3];
				if ((int)va[0] == 0 && (int)vb[0] == 0) p[nvp+j] = 0x8000 | 0;
				else if ((int)va[2] == h && (int)vb[2] == h) p[nvp+j] = 0x8000 | 1;
				else if ((int)va[0] == w && (int)vb[0] == w) p[nvp+j] = 0x8000 | 2;
				else if ((int)va[2] == 0 && (int)vb[2] == 0) p[nvp+j] = 0x8000 | 3;
			}
		}
	}

	al.tf(firstEdge);
	al.tf(edges);
	al.tf(openPoly);
	al.tf(openEdge);
	al.tf(openVert0);
	al.tf(openVert1);
	return true;
}

inline bool canRemoveVertex(const Alloc& al, unsigned short* polys, int npolys, int nverts, int nvp,
							unsigned short rem, bool& result)
{
	(void)al; (void)nverts;
	int numTouchedVerts = 0;
	int numRemainingEdges = 0;
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i*nvp*2];
		const int nv = countPolyVerts(p, nvp);
		int numRemoved = 0;
		int numVerts = 0;
		for (int j = 0; j < nv; ++j)
		{
			if (p[j] == rem) { numTouchedVerts++; numRemoved++; }
			numVerts++;
		}
		if (numRemoved) numRemainingEdges += numVerts-(numRemoved+1);
	}
	if (numRemainingEdges <= 2) { result = false; return true; }
	const int maxEdges = numTouchedVerts*2;
	int nedges = 0;
	int* edges = (int*)al.t(sizeof(int)*maxEdges*3 + 1);
	if (!edges) return false;
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i*nvp*2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0, k = nv-1; j < nv; k = j++)
		{
			if (p[j] == rem || p[k] == rem)
			{
				int a = p[j], b = p[k];
				if (b == rem) { int t = a; a = b; b = t; }
				bool exists = false;
				for (int m = 0; m < nedges; ++m)
				{
					int* e = &edges[m*3];
					if (e[1] == b) { e[2]++; exists = true; }
				}
				if (!exists)
				{
					int* e = &edges[nedges*3];
					e[0] = a; e[1] = b; e[2] = 1;
					nedges++;
				}
			}
		}
	}
	int numOpenEdges = 0;
	for (int i = 0; i < nedges; ++i)
		if (edges[i*3+2] < 2) numOpenEdges++;
	al.tf(edges);
	result = (numOpenEdges <= 2);
	return true;
}

inline bool removeVertex(const Alloc& al, unsigned short* polys, int& npolys, unsigned short*& verts, int& nverts,
						 unsigned short* regs, unsigned char* areas, int nvp, unsigned short rem, int maxTris, int multipleRegs)
{
	int numRemovedVerts = 0;
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i*nvp*2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0; j < nv; ++j)
			if (p[j] == rem) numRemovedVerts++;
	}
	int nedges = 0;
	int* edges = (int*)al.t(sizeof(int)*(size_t)numRemovedVerts*nvp*4 + 4);
	if (!edges) return false;
	int nhole = 0;
	int* hole = (int*)al.t(sizeof(int)*(size_t)numRemovedVerts*nvp + 4);
	if (!hole) { al.tf(edges); return false; }
	int nhreg = 0;
	int* hreg = (int*)al.t(sizeof(int)*(size_t)numRemovedVerts*nvp + 4);
	if (!hreg) { al.tf(edges); al.tf(hole); return false; }
	int nharea = 0;
	int* harea = (int*)al.t(sizeof(int)*(size_t)numRemovedVerts*nvp + 4);
	if (!harea) { al.tf(edges); al.tf(hole); al.tf(hreg); return false; }

	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i*nvp*2];
		const int nv = countPolyVerts(p, nvp);
		bool hasRem = false;
		for (int j = 0; j < nv; ++j) if (p[j] == rem) hasRem = true;
		if (hasRem)
		{
			for (int j = 0, k = nv-1; j < nv; k = j++)
			{
				if (p[j] != rem && p[k] != rem)
				{
					int* e = &edges[nedges*4];
					e[0] = p[k]; e[1] = p[j];
					e[2] = regs ? regs[i] : 0;
					e[3] = areas ? areas[i] : 0;
					nedges++;
				}
			}
			unsigned short* p2 = &polys[(npolys-1)*nvp*2];
			if (p != p2) memcpy(p, p2, sizeof(unsigned short)*nvp);
			memset(p+nvp, 0xff, sizeof(unsigned short)*nvp);
			if (regs) regs[i] = regs[npolys-1];
			if (areas) areas[i] = areas[npolys-1];
			npolys--;
			--i;
		}
	}
	for (int i = (int)rem; i < nverts - 1; ++i)
	{
		verts[i*3+0] = verts[(i+1)*3+0];
		verts[i*3+1] = verts[(i+1)*3+1];
		verts[i*3+2] = verts[(i+1)*3+2];
	}
	nverts--;
	for (int i = 0; i < npolys; ++i)
	{
		unsigned short* p = &polys[i*nvp*2];
		const int nv = countPolyVerts(p, nvp);
		for (int j = 0; j < nv; ++j) if (p[j] > rem) p[j]--;
	}
	for (int i = 0; i < nedges; ++i)
	{
		if (edges[i*4+0] > rem) edges[i*4+0]--;
		if (edges[i*4+1] > rem) edges[i*4+1]--;
	}
	if (nedges == 0) { al.tf(edges); al.tf(hole); al.tf(hreg); al.tf(harea); return true; }
	pushBack(edges[0], hole, nhole);
	pushBack(edges[2], hreg, nhreg);
	pushBack(edges[3], harea, nharea);
	while (nedges)
	{
		bool match = false;
		for (int i = 0; i < nedges; ++i)
		{
			const int ea = edges[i*4+0];
			const int eb = edges[i*4+1];
			const int r = edges[i*4+2];
			const int a = edges[i*4+3];
			bool add = false;
			if (hole[0] == eb) { pushFront(ea, hole, nhole); pushFront(r, hreg, nhreg); pushFront(a, harea, nharea); add = true; }
			else if (hole[nhole-1] == ea) { pushBack(eb, hole, nhole); pushBack(r, hreg, nhreg); pushBack(a, harea, nharea); add = true; }
			if (add)
			{
				edges[i*4+0] = edges[(nedges-1)*4+0];
				edges[i*4+1] = edges[(nedges-1)*4+1];
				edges[i*4+2] = edges[(nedges-1)*4+2];
				edges[i*4+3] = edges[(nedges-1)*4+3];
				--nedges;
				match = true;
				--i;
			}
		}
		if (!match) break;
	}
	int* tris = (int*)al.t(sizeof(int)*(size_t)nhole*3 + 4);
	int* tverts = (int*)al.t(sizeof(int)*(size_t)nhole*4 + 4);
	int* thole = (int*)al.t(sizeof(int)*(size_t)nhole + 4);
	bool ok = tris && tverts && thole;
	if (ok)
	{
		for (int i = 0; i < nhole; ++i)
		{
			const int pi = hole[i];
			tverts[i*4+0] = verts[pi*3+0];
			tverts[i*4+1] = verts[pi*3+1];
			tverts[i*4+2] = verts[pi*3+2];
			tverts[i*4+3] = 0;
			thole[i] = i;
		}
		int ntris = triangulate(nhole, tverts, thole, tris);
		if (ntris < 0) ntris = -ntris;
		unsigned short* polys2 = (unsigned short*)al.t(sizeof(unsigned short)*(ntris+1)*nvp);
		unsigned short* pregs = regs ? (unsigned short*)al.t(sizeof(unsigned short)*ntris) : 0;
		unsigned char* pareas = areas ? (unsigned char*)al.t(sizeof(unsigned char)*ntris) : 0;
		ok = polys2 && (!regs || pregs) && (!areas || pareas);
		if (ok)
		{
			unsigned short* tmpPoly = &polys2[ntris*nvp];
			int npolys2 = 0;
			memset(polys2, 0xff, ntris*nvp*sizeof(unsigned short));
			for (int j = 0; j < ntris; ++j)
			{
				int* t = &tris[j*3];
				if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2])
				{
					polys2[npolys2*nvp+0] = (unsigned short)hole[t[0]];
					polys2[npolys2*nvp+1] = (unsigned short)hole[t[1]];
					polys2[npolys2*nvp+2] = (unsigned short)hole[t[2]];
					if (regs)
					{
						if (hreg[t[0]] != hreg[t[1]] || hreg[t[1]] != hreg[t[2]]) pregs[npolys2] = (unsigned short)multipleRegs;
						else pregs[npolys2] = (unsigned short)hreg[t[0]];
					}
					if (areas) pareas[npolys2] = (unsigned char)harea[t[0]];
					npolys2++;
				}
			}
			if (npolys2 && nvp > 3)
			{
				for (;;)
				{
					int bestMergeVal = 0;
					int bestPa = 0, bestPb = 0, bestEa = 0, bestEb = 0;
					for (int j = 0; j < npolys2-1; ++j)
					{
						unsigned short* pj = &polys2[j*nvp];
						for (int k = j+1; k < npolys2; ++k)
						{
							unsigned short* pk = &polys2[k*nvp];
							int ea, eb;
							int v = getPolyMergeValue(pj, pk, verts, ea, eb, nvp);
							if (v > bestMergeVal) { bestMergeVal = v; bestPa = j; bestPb = k; bestEa = ea; bestEb = eb; }
						}
					}
					if (bestMergeVal > 0)
					{
						unsigned short* pa = &polys2[bestPa*nvp];
						unsigned short* pb = &polys2[bestPb*nvp];
						mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp);
						if (regs && pregs[bestPa] != pregs[bestPb]) pregs[bestPa] = (unsigned short)multipleRegs;
						unsigned short* last = &polys2[(npolys2-1)*nvp];
						if (pb != last) memcpy(pb, last, sizeof(unsigned short)*nvp);
						if (regs) pregs[bestPb] = pregs[npolys2-1];
						if (areas) pareas[bestPb] = pareas[npolys2-1];
						npolys2--;
					}
					else break;
				}
			}
			for (int i = 0; i < npolys2; ++i)
			{
				if (npolys >= maxTris) break;
				unsigned short* p = &polys[npolys*nvp*2];
				memset(p, 0xff, sizeof(unsigned short)*nvp*2);
				for (int j = 0; j < nvp; ++j) p[j] = polys2[i*nvp+j];
				if (regs) regs[npolys] = pregs[i];
				if (areas) areas[npolys] = pareas[i];
				npolys++;
			}
			if (polys2) al.tf(polys2);
			if (pregs) al.tf(pregs);
			if (pareas) al.tf(pareas);
		}
	}
	if (tris) al.tf(tris);
	if (tverts) al.tf(tverts);
	if (thole) al.tf(thole);
	al.tf(edges); al.tf(hole); al.tf(hreg); al.tf(harea);
	return ok;
}

// Builds a polygon mesh from a contour set. On success the output arrays are
// owned by the caller (via the perm allocator) and must be freed with it.
// multipleRegs is the region id written when a merged polygon spans regions.
inline bool buildPolyMesh(const Alloc& al, const Input& in, int nvp, Output& out, int multipleRegs)
{
	out.verts = 0; out.polys = 0; out.regs = 0; out.areas = 0;
	out.nverts = 0; out.npolys = 0; out.nvp = nvp; out.maxpolys = 0;

	int maxVertices = 0;
	int maxTris = 0;
	int maxVertsPerCont = 0;
	for (int i = 0; i < in.nconts; ++i)
	{
		if (in.conts[i].nverts < 3) continue;
		maxVertices += in.conts[i].nverts;
		maxTris += in.conts[i].nverts - 2;
		if (in.conts[i].nverts > maxVertsPerCont) maxVertsPerCont = in.conts[i].nverts;
	}
	if (maxVertices == 0 || maxTris == 0) return true;
	if (maxVertices >= 0xfffe) return false;

	unsigned char* vflags = (unsigned char*)al.t(maxVertices);
	if (!vflags) return false;
	memset(vflags, 0, maxVertices);
	unsigned short* verts = (unsigned short*)al.p(sizeof(unsigned short)*maxVertices*3);
	unsigned short* polys = (unsigned short*)al.p(sizeof(unsigned short)*maxTris*nvp*2);
	unsigned short* regs = (unsigned short*)al.p(sizeof(unsigned short)*maxTris);
	unsigned char* areas = (unsigned char*)al.p(sizeof(unsigned char)*maxTris);
	if (!verts || !polys || !regs || !areas) return false;
	int nverts = 0, npolys = 0;
	memset(verts, 0, sizeof(unsigned short)*maxVertices*3);
	memset(polys, 0xff, sizeof(unsigned short)*maxTris*nvp*2);
	memset(regs, 0, sizeof(unsigned short)*maxTris);
	memset(areas, 0, sizeof(unsigned char)*maxTris);

	int* nextVert = (int*)al.t(sizeof(int)*maxVertices);
	if (!nextVert) return false;
	memset(nextVert, 0, sizeof(int)*maxVertices);
	int* firstVert = (int*)al.t(sizeof(int)*VERTEX_BUCKET_COUNT);
	if (!firstVert) return false;
	for (int i = 0; i < VERTEX_BUCKET_COUNT; ++i) firstVert[i] = -1;
	int* indices = (int*)al.t(sizeof(int)*maxVertsPerCont);
	int* tris = (int*)al.t(sizeof(int)*maxVertsPerCont*3);
	int* contVerts = (int*)al.t(sizeof(int)*maxVertsPerCont*4);
	unsigned short* polysTmp = (unsigned short*)al.t(sizeof(unsigned short)*(maxVertsPerCont+1)*nvp);
	if (!indices || !tris || !polysTmp) return false;
	unsigned short* tmpPoly = &polysTmp[maxVertsPerCont*nvp];

	for (int i = 0; i < in.nconts; ++i)
	{
		const Contour& cont = in.conts[i];
		if (cont.nverts < 3) continue;
		for (int j = 0; j < cont.nverts; ++j) indices[j] = j;
		for (int j = 0; j < cont.nverts; ++j)
		{
			contVerts[j*4+0] = cont.vx(j);
			contVerts[j*4+1] = cont.vy(j);
			contVerts[j*4+2] = cont.vz(j);
			contVerts[j*4+3] = cont.vflags(j);
		}
		// Triangulation reads int vertices; give it the normalized copy.
		int ntris = triangulate(cont.nverts, contVerts, &indices[0], &tris[0]);
		if (ntris <= 0) ntris = -ntris;
		for (int j = 0; j < cont.nverts; ++j)
		{
			const int* v = &contVerts[j*4];
			indices[j] = addVertex((unsigned short)v[0], (unsigned short)v[1], (unsigned short)v[2],
								   verts, firstVert, nextVert, nverts);
			if (in.borderFlag && (v[3] & in.borderFlag)) vflags[indices[j]] = 1;
		}
		int tnp = 0;
		memset(polysTmp, 0xff, maxVertsPerCont*nvp*sizeof(unsigned short));
		for (int j = 0; j < ntris; ++j)
		{
			int* t = &tris[j*3];
			if (t[0] != t[1] && t[0] != t[2] && t[1] != t[2])
			{
				polysTmp[tnp*nvp+0] = (unsigned short)indices[t[0]];
				polysTmp[tnp*nvp+1] = (unsigned short)indices[t[1]];
				polysTmp[tnp*nvp+2] = (unsigned short)indices[t[2]];
				tnp++;
			}
		}
		if (!tnp) continue;
		if (nvp > 3)
		{
			for (;;)
			{
				int bestMergeVal = 0;
				int bestPa = 0, bestPb = 0, bestEa = 0, bestEb = 0;
				for (int j = 0; j < tnp-1; ++j)
				{
					unsigned short* pj = &polysTmp[j*nvp];
					for (int k = j+1; k < tnp; ++k)
					{
						unsigned short* pk = &polysTmp[k*nvp];
						int ea, eb;
						int v = getPolyMergeValue(pj, pk, verts, ea, eb, nvp);
						if (v > bestMergeVal) { bestMergeVal = v; bestPa = j; bestPb = k; bestEa = ea; bestEb = eb; }
					}
				}
				if (bestMergeVal > 0)
				{
					unsigned short* pa = &polysTmp[bestPa*nvp];
					unsigned short* pb = &polysTmp[bestPb*nvp];
					mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp);
					unsigned short* lastPoly = &polysTmp[(tnp-1)*nvp];
					if (pb != lastPoly) memcpy(pb, lastPoly, sizeof(unsigned short)*nvp);
					tnp--;
				}
				else break;
			}
		}
		for (int j = 0; j < tnp; ++j)
		{
			unsigned short* p = &polys[npolys*nvp*2];
			unsigned short* q = &polysTmp[j*nvp];
			for (int k = 0; k < nvp; ++k) p[k] = q[k];
			regs[npolys] = cont.reg;
			areas[npolys] = cont.area;
			npolys++;
			if (npolys > maxTris) return false;
		}
	}

	for (int i = 0; i < nverts; ++i)
	{
		if (vflags[i])
		{
			// If the test cannot allocate its scratch, leave the vertex in place
			// (matching the historical Recast behaviour) rather than failing.
			bool can = false;
			if (!canRemoveVertex(al, polys, npolys, nverts, nvp, (unsigned short)i, can) || !can) continue;
			if (!removeVertex(al, polys, npolys, verts, nverts, regs, areas, nvp, (unsigned short)i, maxTris, multipleRegs)) return false;
			for (int j = i; j < nverts; ++j) vflags[j] = vflags[j+1];
			--i;
		}
	}

	if (!buildMeshAdjacency(al, in, polys, npolys, verts, nverts, nvp)) return false;

	out.verts = verts; out.polys = polys; out.regs = regs; out.areas = areas;
	out.nverts = nverts; out.npolys = npolys; out.nvp = nvp; out.maxpolys = maxTris;
	return true;
}

} // namespace rcMeshCore

#endif // RECAST_MESH_CORE_H
