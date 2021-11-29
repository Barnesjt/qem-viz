#ifndef __MESH_H__
#define __MESH_H__

#include "icVector.H"
#include <stdio.h>
#include <stddef.h>

#include <vector>

using std::vector;

class Face;
class Mesh;
class Edge;

class Vertex {
public:
	double x, y, z;
	int index;

	Mesh* mesh;
	vector <Face*> faces;
	vector <Edge*> edges;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
};

class Edge {
public:
	Vertex* verts[2];
	vector <Face*> faces;
	Mesh* mesh;

public:
	Edge(Vertex* v0, Vertex* v1, Mesh* meshIn) { verts[0] = v0; verts[0] = v1; mesh = meshIn; };
};

class Face {
public:
	Vertex* verts[3];
	Edge* edges[3];
	Mesh* mesh;

	icVector3 normal;
	int index;

public:
	Face(Vertex* v0, Vertex* v1, Vertex* v2, icVector3 norm) { verts[0] = v0; verts[1] = v1; verts[2] = v2; normal = norm; }
};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
	vector <Edge*> elist;
	unsigned char orientation; //0=ccw, 1=cw
	icVector3 center;
	double radius;

public:

	/*constructors*/
	Mesh();
	Mesh(char*);

	/*initialization functions*/
	void initialize();
	void finalize();

	/*utilties*/
};



#endif