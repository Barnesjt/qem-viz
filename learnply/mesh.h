#ifndef __MESH_H__
#define __MESH_H__

#include "icVector.H"
#include <stdio.h>
#include <stddef.h>

#include <vector>

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/matrix_inverse.hpp"


using std::vector;

class Face;
class Mesh;
class Edge;
class Pair;

__declspec(align(16)) class Vertex {
public:
	double x, y, z;
	int index=-1;

	glm::mat4x4 Q = glm::mat4x4(0.);
	Mesh* mesh=NULL;
	vector <Face*> faces;
	vector <Edge*> edges;
	vector <Pair*> pairs;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
	double distance(Vertex * v1);
	bool makesEdgeWith(Vertex* v1);

	virtual ~Vertex()
	{
	}

	void* operator new(size_t i)
	{
		return _mm_malloc(i, 16);
	}

	void operator delete(void* p)
	{
		_mm_free(p);
	}
};

__declspec(align(16)) class Edge {
public:
	Vertex* verts[2];
	vector <Face*> faces;
	Mesh* mesh;

	int index=-1;

public:
	Edge(Vertex* v0, Vertex* v1, Mesh* meshIn) { verts[0] = v0; verts[1] = v1; mesh = meshIn; };

	virtual ~Edge()
	{
	}

	void* operator new(size_t i)
	{
		return _mm_malloc(i, 16);
	}

	void operator delete(void* p)
	{
		_mm_free(p);
	}
};

__declspec(align(16)) class Face {
public:
	Vertex* verts[3];
	Edge* edges[3];
	Mesh* mesh;

	glm::mat4x4 Quadric = glm::mat4x4(0.);
	icVector3 normal;
	int index=-1;

public:
	Face(Vertex* v0, Vertex* v1, Vertex* v2, icVector3 norm) { verts[0] = v0; verts[1] = v1; verts[2] = v2; normal = norm; }
	void calcQuadric();
	void calcNormal();

	virtual ~Face()
	{
	}

	void* operator new(size_t i)
	{
		return _mm_malloc(i, 16);
	}

	void operator delete(void* p)
	{
		_mm_free(p);
	}
};

__declspec(align(16)) class Pair {
public:
	Vertex* v0;
	Vertex* v1;
	Edge* edge = NULL;
	int index=0;
	glm::mat4x4 Quadric = glm::mat4x4(0.);
	double Error=-1;

public:
	Pair(Vertex* v0In, Vertex* v1In) { v0 = v0In; v1 = v1In; }
	Pair(Vertex* v0In, Vertex* v1In, Edge* edgeIn) { v0 = v0In; v1 = v1In; edge = edgeIn; }

	void updateQuadric();
	void updateError();
	icVector3 Vector();
	icVector3 QuadricVector();
	double QuadricError(icVector3 v);

	virtual ~Pair()
	{
	}

	void* operator new(size_t i)
	{
		return _mm_malloc(i, 16);
	}

	void operator delete(void* p)
	{
		_mm_free(p);
	}

	bool operator < (const Pair& pair) const
	{
		return (Error < pair.Error);
	}

};

struct PairComp
{
	bool const operator()(Pair* lhs, Pair* rhs) const
	{
		return (*lhs) < (*rhs);
	}
};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
	vector <Edge*> elist;
	unsigned char orientation=0; //0=ccw, 1=cw
	icVector3 center;
	double radius=-1;
	unsigned int maxFaces;
	vector <Pair*> validPairs;

public:

	/*constructors*/
	Mesh(char*);

	/*utilties*/
	void seedValidPairs(double maxDistForPairs);
	void seedInitialQuadrics();
	void edge2fContract(Pair* target);
	void edge1fContract(Pair* target);
	void nonEdgeContract(Pair* target);

	/*initialization functions*/
	void initialize(double maxDistForPairs);
	void finalize();

	void simplify(unsigned int target);

};



#endif