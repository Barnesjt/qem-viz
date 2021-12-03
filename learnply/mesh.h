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

class Vertex {
public:
	double x, y, z;
	int index;

	glm::dmat4x4 Q;
	Mesh* mesh;
	vector <Face*> faces;
	vector <Edge*> edges;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
	double distance(Vertex * v1);
	bool makesEdgeWith(Vertex* v1);
};

class Edge {
public:
	Vertex* verts[2];
	vector <Face*> faces;
	Mesh* mesh;

	int index;

public:
	Edge(Vertex* v0, Vertex* v1, Mesh* meshIn) { verts[0] = v0; verts[1] = v1; mesh = meshIn; };
};

class Face {
public:
	Vertex* verts[3];
	Edge* edges[3];
	Mesh* mesh;

	glm::dmat4x4 Quadric;
	icVector3 normal;
	int index;

public:
	Face(Vertex* v0, Vertex* v1, Vertex* v2, icVector3 norm) { verts[0] = v0; verts[1] = v1; verts[2] = v2; normal = norm; }
	void calcQuadric();
};

class Pair {
public:
	Vertex* v0;
	Vertex* v1;
	Edge* edge = NULL;
	int index;
	glm::dmat4x4 Quadric;
	double Error;

public:
	Pair(Vertex* v0In, Vertex* v1In) { v0 = v0In; v1 = v1In; }
	Pair(Vertex* v0In, Vertex* v1In, Edge* edgeIn) { v0 = v0In; v1 = v1In; edge = edgeIn; }

	void updateQuadric();
	void updateError();
	icVector3 Vector();
	icVector3 QuadricVector();
	double QuadricError(icVector3 v);

};

class PairComparison
{
	bool reverse;
public:
	PairComparison(const bool& revparam = false)
	{
		reverse = revparam;
	}
	bool operator() ( Pair* lhs, Pair* rhs)
	{
		if (reverse) return (lhs->Error > rhs->Error);
		else return (lhs->Error < rhs->Error);
	}
};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
	vector <Edge*> elist;
	unsigned char orientation; //0=ccw, 1=cw
	icVector3 center;
	double radius;
	int maxFaces;

	vector <Pair*> validPairs;

public:

	/*constructors*/
	Mesh();
	Mesh(char*);

	/*utilties*/
	void seedValidPairs(double maxDistForPairs);
	void seedInitialQuadrics();

	/*initialization functions*/
	void initialize(double maxDistForPairs);
	void finalize();

	void simplify(int target);

};



#endif