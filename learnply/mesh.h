#ifndef __MESH_H__
#define __MESH_H__

#include "icVector.H"
#include <stdio.h>
#include <stddef.h>

#include <vector>

using std::vector;

class Face;

class Vertex {
public:
	double x, y, z;
	int index;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz; }
};

class Face {
public:
	Vertex* verts[3];
	icVector3 normal;
	int index;

public:
	Face(Vertex* v0, Vertex* v1, Vertex* v2, icVector3 norm) { verts[0] = v0; verts[1] = v1; verts[2] = v2; normal = norm; }
};

class Mesh {
public:

	vector <Face*> flist;
	vector <Vertex*> vlist;
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