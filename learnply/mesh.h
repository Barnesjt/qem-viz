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

#include <sstream>
#include <string>


using std::vector;

class Face;
class Mesh;
class Pair;
class PairSimple;

__declspec(align(16)) class Vertex {
public:
	double x, y, z;

	glm::mat4x4 Q = glm::mat4x4(0.);
	vector <Face*> faces;
	vector <Vertex*> edges;
	vector <Vertex*> pairs;

public:
	Vertex(double x0, double y0, double z0) { x = x0; y = y0; z = z0; }
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

__declspec(align(16)) class Face {
public:
	Vertex* verts[3];

	glm::mat4x4 Quadric = glm::mat4x4(0.);
	icVector3 normal;

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

class PairSimple {
public:
	Vertex* v0;
	Vertex* v1;
	bool edge;

public:
	PairSimple(Vertex* v0In, Vertex* v1In, bool edgeIn) { v0 = v0In; v1 = v1In; edge = edgeIn; }

	//silly comparison function to prune duplicate pairs with an ordered set
	//bool operator < (const PairSimple& pair) const
	//{
	//	std::string rightStr;
	//	std::string leftStr;
	//	std::stringstream ssv0, ssv1, ssv0p, ssv1p;

	//	const void* addv0 = static_cast<const void*>(v0);
	//	ssv0 << addv0;
	//	const void* addv1 = static_cast<const void*>(v1);
	//	ssv1 << addv1;
	//	const void* addv0p = static_cast<const void*>(pair.v0);
	//	ssv0p << addv0p;
	//	const void* addv1p = static_cast<const void*>(pair.v1);
	//	ssv1p << addv1p;

	//	if (ssv0.str() < ssv1.str()) leftStr = ssv0.str() + ssv1.str();
	//	else leftStr = ssv1.str() + ssv0.str();
	//	if (ssv0p.str() < ssv1p.str()) rightStr = ssv0p.str() + ssv1p.str();
	//	else rightStr = ssv1p.str() + ssv0p.str();

	//	return (leftStr < rightStr);
	//}

};

__declspec(align(16)) class Pair {
public:
	Vertex* v0;
	Vertex* v1;
	bool edge = NULL;
	glm::mat4x4 Quadric = glm::mat4x4(0.);
	double Error=-1;

public:
	Pair(Vertex* v0In, Vertex* v1In, bool edgeIn) { v0 = v0In; v1 = v1In; edge = edgeIn; }

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

	bool operator > (const Pair& pair) const
	{
		return (Error > pair.Error);
	}

};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
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