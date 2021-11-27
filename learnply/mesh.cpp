#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
#include <fstream>
#include "icVector.H"
#include "icMatrix.H"
#include "mesh.h"

#include <vector>

// delimiters for parsing the obj file:

#define OBJDELIMS		" \t"

char* ReadRestOfLine(FILE*);
void ReadObjVTN(char*, int*, int*, int*);
void Cross(double[3], double[3], double[3]);
double Unit(double[3]);
double Unit(double[3], double[3]);

struct face
{
	int v, n, t;
};


Mesh::Mesh(char* file) {

	char* cmd;		// the command string
	char* str;		// argument string

	std::vector <Vertex*> Vertices;
	std::vector <Face*> Faces;

	Vertices.clear();
	Faces.clear();

	// open the input file:

	FILE* fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open .obj file '%s'\n", file);
	}
	else {
		for (; ; ) {
			char* line = ReadRestOfLine(fp);
			if (line == NULL)
				break;

			// skip this line if it is a comment:
			if (line[0] == '#')
				continue;

			// skip this line if it is something we don't feel like handling today:
			if (line[0] == 'g' || line[0] == 'm' || line[0] == 's' || line[0] == 'u')
				continue;

			// get the command string:
			cmd = strtok(line, OBJDELIMS);


			// skip this line if it is empty:
			if (cmd == NULL) continue;

			// skip vertex normals
			if (strcmp(cmd, "vn") == 0) continue;

			// skip texture coords
			if (strcmp(cmd, "vt") == 0) continue;

			if (strcmp(cmd, "v") == 0) {
				str = strtok(NULL, OBJDELIMS);
				double x = atof(str);

				str = strtok(NULL, OBJDELIMS);
				double y = atof(str);

				str = strtok(NULL, OBJDELIMS);
				double z = atof(str);

				Vertices.push_back(new Vertex(x, y, z));

				continue;
			}

			if (strcmp(cmd, "f") == 0)
			{
				struct face vertices[10];
				for (int i = 0; i < 10; i++)
				{
					vertices[i].v = 0;
					vertices[i].n = 0;
					vertices[i].t = 0;
				}

				int sizev = (int)Vertices.size();

				int numVertices = 0;
				bool valid = true;
				int vtx = 0;
				char* str;
				while ((str = strtok(NULL, OBJDELIMS)) != NULL)
				{
					int v, n, t;
					ReadObjVTN(str, &v, &t, &n);

					// if v, n, or t are negative, they are wrt the end of their respective list:

					if (v < 0)
						v += (sizev + 1);

					// be sure we are not out-of-bounds (<vector> will abort):

					if (v > sizev)
					{
						if (v != 0)
							fprintf(stderr, "Read vertex coord %d, but only have %d so far\n", v, sizev);
						v = 0;
						valid = false;
					}

					vertices[vtx].v = v;
					//vertices[vtx].n = n;
					//vertices[vtx].t = t;
					vtx++;

					if (vtx >= 10)
						break;

					numVertices++;
				}

				// if vertices are invalid, don't draw anything this time:
				if (!valid || (numVertices < 3)) continue;

				// list the vertices:
				int numTriangles = numVertices - 2;

				for (int it = 0; it < numTriangles; it++) {
					int vv[3];
					vv[0] = 0;
					vv[1] = it + 1;
					vv[2] = it + 2;

					// get the planar normal, in case vertex normals are not defined:

					Vertex* v0 = Vertices[vertices[vv[0]].v - 1];
					Vertex* v1 = Vertices[vertices[vv[1]].v - 1];
					Vertex* v2 = Vertices[vertices[vv[2]].v - 1];

					double v01[3], v02[3], norm[3];
					v01[0] = v1->x - v0->x;
					v01[1] = v1->y - v0->y;
					v01[2] = v1->z - v0->z;
					v02[0] = v2->x - v0->x;
					v02[1] = v2->y - v0->y;
					v02[2] = v2->z - v0->z;
					Cross(v01, v02, norm);
					Unit(norm, norm);
					Faces.push_back(new Face(v0, v1, v2, icVector3(norm)));
				}
				continue;
			}


			if (strcmp(cmd, "s") == 0)
			{
				continue;
			}

		}

		orientation = 0;
		flist = Faces;
		vlist = Vertices;

		fclose(fp);

	}

}


Mesh::Mesh() {
	//default constructor
}

void Mesh::initialize() {
}

void Mesh::finalize() {
	for (auto f : flist) {
		free(f->verts);
		free(f);
		flist.clear();
	}
	for (auto v : vlist) {
		free(v);
		vlist.clear();
	}

}

char* ReadRestOfLine(FILE* fp){
	static char* line;
	std::vector<char> tmp(1000);
	tmp.clear();

	for (; ; )
	{
		int c = getc(fp);

		if (c == EOF && tmp.size() == 0)
		{
			return NULL;
		}

		if (c == EOF || c == '\n')
		{
			delete[] line;
			line = new char[tmp.size() + 1];
			for (int i = 0; i < (int)tmp.size(); i++)
			{
				line[i] = tmp[i];
			}
			line[tmp.size()] = '\0';	// terminating null
			return line;
		}
		else
		{
			tmp.push_back(c);
		}
	}

	return "";
}

void ReadObjVTN(char* str, int* v, int* t, int* n){
	// can be one of v, v//n, v/t, v/t/n:

	if (strstr(str, "//"))				// v//n
	{
		*t = 0;
		sscanf(str, "%d//%d", v, n);
		return;
	}
	else if (sscanf(str, "%d/%d/%d", v, t, n) == 3)	// v/t/n
	{
		return;
	}
	else
	{
		*n = 0;
		if (sscanf(str, "%d/%d", v, t) == 2)		// v/t
		{
			return;
		}
		else						// v
		{
			*n = *t = 0;
			sscanf(str, "%d", v);
		}
	}
}

void Cross(double v1[3], double v2[3], double vout[3]) {

	double tmp[3];

	tmp[0] = v1[1] * v2[2] - v2[1] * v1[2];
	tmp[1] = v2[0] * v1[2] - v1[0] * v2[2];
	tmp[2] = v1[0] * v2[1] - v2[0] * v1[1];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}

double Unit(double v[3]) {
	double dist = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

	if (dist > 0.0) {
		dist = sqrt(dist);
		v[0] /= dist;
		v[1] /= dist;
		v[2] /= dist;
	}

	return dist;
}

double Unit(double vin[3], double vout[3]) {
	double dist = vin[0] * vin[0] + vin[1] * vin[1] + vin[2] * vin[2];

	if (dist > 0.0) {
		dist = sqrt(dist);
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	} else {
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}