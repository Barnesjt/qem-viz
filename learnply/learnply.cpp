#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "icVector.H"
#include "icMatrix.H"
#include "trackball.h"
#include "tmatrix.h"

#include "drawUtil.h"

// Data structure to hold Mesh data, includes code adapted from: http://web.engr.oregonstate.edu/~mjb/cs550/loadobjfile.cpp
#include "mesh.h"

//Added include, specifically for reading directory contents
#include <Windows.h>

//Added natural sort from: https://github.com/scopeInfinity/NaturalSort
//For sorting file names into best order (e.g. r2.ply before r10.ply)
#include "natural_sort.hpp"


//Vectors and Strings
using std::vector;
using std::string;

//Vector of Meshes instead (Because why not just load everything at once)
//Individual vectors for different data, this could be rolled into a single class, maybe later
vector<string> all_obj_files;
vector<Mesh*> all_meshes;

//The marker for which ply we're looking at, used to index into the vectors and cycle through them
int curr_mesh = 0;
int face_target = 0;
double t_ratio = 100.;
bool view_error = false;

string OBJ_PATH = "../data/geometry/";

Mesh* mesh;

/*scene related variables*/
const float zoomspeed = 0.9;
int win_width = 1024;
int win_height = 1024;
float aspectRatio = win_width / win_height;
const int view_mode = 0;		// 0 = othogonal, 1=perspective
const double radius_factor = 0.9;

/*
Use keys 1 to 0 to switch among different display modes.
Each display mode can be designed to show one type 
visualization result.

Predefined ones: 
display mode 1: solid rendering
display mode 2: show wireframes
display mode 3: render each quad with colors of vertices
*/
int display_mode = 1;

/*User Interaction related variabes*/
float s_old, t_old;
float rotmat[4][4];
double zoom = 1.0;
double translation[2] = { 0, 0 };
int mouse_mode = -2;	// -1 = no action, 1 = translate y, 2 = rotate

/******************************************************************************
Forward declaration of functions
******************************************************************************/

void init(void);

/*glut attaching functions*/
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);
void mousewheel(int wheel, int direction, int x, int y);
void reshape(int width, int height);

/*display vis results*/
void display_mesh(Mesh* mesh);

/*added functions*/
string dispModeString(int mode);
void DoRasterString(float x, float y, float z, char* s);
std::vector<std::string> get_all_obj_in_folder(std::string folder);

/******************************************************************************
Main program.
******************************************************************************/

int main(int argc, char* argv[]) {
	//Get all obj files in directory
	all_obj_files = get_all_obj_in_folder(OBJ_PATH);

	//Natural sorting of obj files
	std::sort(all_obj_files.begin(), all_obj_files.end(), SI::natural::compare<string>);

	for (std::string i : all_obj_files) {
		all_meshes.push_back(new Mesh(const_cast<char*>((OBJ_PATH + i).c_str())));
		all_meshes.back()->initialize(all_meshes.back()->radius/t_ratio);
	}

	//if no ply files, then exit with error. otherwise set the first ply
	if (all_meshes.empty()) {
		return 1;
	} else {
		mesh = all_meshes.at(0);
		curr_mesh = 0;
		face_target = mesh->flist.size();
	}

	/*init glut and create window*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(20, 20);
	glutInitWindowSize(win_width, win_height);
	glutCreateWindow("CS553 - Group 1 - QEM w/ Ellipsoids");

	/*initialize openGL*/
	init();
	
	/*the render function and callback registration*/
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutMouseWheelFunc(mousewheel);
	
	/*event processing loop*/
	glutMainLoop();
	
	/*clear memory before exit*/
	mesh->finalize();	// finalize everything
	return 0;
}

/******************************************************************************
Get all .obj files

Function adapted from: https://stackoverflow.com/a/20847429
******************************************************************************/
// Function adapted from: https://stackoverflow.com/a/20847429
std::vector<std::string> get_all_obj_in_folder(std::string folder) {
	std::vector<std::string> names;
	std::string search_path = folder + "/*.obj";
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(search_path.c_str(), &fd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			// read all (real) files in current folder
			// , delete '!' read other 2 default folder . and ..
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}

/******************************************************************************
Use glut to display a string of characters using a raster font

Borrowed from CS450 Sample Code
******************************************************************************/
void
DoRasterString(float x, float y, float z, char* s) {
	glRasterPos3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
	char c;			// one character to print
	for (; (c = *s) != '\0'; s++) {
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, c);
	}
}

/******************************************************************************
Convert a display mode number to a string
******************************************************************************/
string dispModeString(int mode) {
	switch (mode) {
	case 1:
		return "Solid";
	case 2:
		return "Solid w/ Error Ellipsoids";
	default:
		return "Unrecognized";
	}
}

/******************************************************************************
Set projection mode
******************************************************************************/

void set_view(GLenum mode) {
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };

	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (aspectRatio >= 1.0) {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom * aspectRatio, radius_factor * zoom * aspectRatio, -radius_factor * zoom, radius_factor * zoom, 0.1, 1000);
	}
	else {
		if (view_mode == 0)
			glOrtho(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, -1000, 1000);
		else
			glFrustum(-radius_factor * zoom, radius_factor * zoom, -radius_factor * zoom / aspectRatio, radius_factor * zoom / aspectRatio, 0.1, 1000);
	}

	GLfloat light_position[3];
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

/******************************************************************************
Update the scene
******************************************************************************/

void set_scene(GLenum mode, Mesh* mesh){
	glTranslatef(translation[0], translation[1], -3.0);

	/*multiply rotmat to current mat*/
	{
		int i, j, index = 0;

		GLfloat mat[16];

		for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++)
				mat[index++] = rotmat[i][j];

		glMultMatrixf(mat);
	}

	//glScalef(0.9, 0.9, 0.9);
	//glTranslatef(0.0, 0.0, 0.0);

	glScalef(0.9 / mesh->radius, 0.9 / mesh->radius, 0.9 / mesh->radius);
	glTranslatef(-mesh->center.entry[0], -mesh->center.entry[1], -mesh->center.entry[2]);
}

/******************************************************************************
Init scene
******************************************************************************/

void init(void) {

	mat_ident(rotmat);

	/* select clearing color */
	glClearColor(0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	//set pixel storage modes
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glEnable(GL_NORMALIZE);
	if (mesh->orientation == 0)
		glFrontFace(GL_CW);
	else
		glFrontFace(GL_CCW);
}

/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {

	switch (key) {
	case 27:	// set excape key to exit program
		for (auto m : all_meshes) m->finalize(); // finalize_everything
		exit(0);
		break;

	case 'q':
		int newFaceTarget;
		std::cout << "Input New Face Target: ";
		std::cin >> newFaceTarget;
		if (newFaceTarget < 4) std::cout << "Too few faces, press 'q' to try again";
		else if (newFaceTarget > mesh->maxFaces) std::cout << "Too many faces, press 'q' to try again";
		else face_target = newFaceTarget;

		if (face_target > mesh->flist.size()) {
			all_meshes.at(curr_mesh) = new Mesh(const_cast<char*>((OBJ_PATH + all_obj_files.at(curr_mesh)).c_str()));
			mesh = all_meshes.at(curr_mesh);
			mesh->initialize(mesh->radius / t_ratio);
		}
		if (face_target < mesh->flist.size()){
			mesh->simplify(face_target);
		}
		glutPostRedisplay();
		break;

	case 'e':
		view_error = !view_error;
		glutPostRedisplay();
		break;

	case 'r':	// reset rotation and transformation
		mat_ident(rotmat);
		translation[0] = 0;
		translation[1] = 0;
		zoom = 1.0;
		glutPostRedisplay();
		break;

	case 't':
		//toggle through all ply files, increments and modulo on the vector size
		curr_mesh = ++curr_mesh % all_meshes.size();

		//reset the poly we are displaying, w/ min and max
		mesh = all_meshes.at(curr_mesh);

		face_target = mesh->flist.size();

		//recall keyboard function to make sure the correct display mode is set
		keyboard('0' + display_mode, x, y);
		glutPostRedisplay();
		break;
	}
}

/******************************************************************************
Callback function for dragging mouse
******************************************************************************/

void motion(int x, int y) {
	float r[4];
	float s, t;

	s = (2.0 * x - win_width) / win_width;
	t = (2.0 * (win_height - y) - win_height) / win_height;

	if ((s == s_old) && (t == t_old))
		return;

	switch (mouse_mode) {
	case 2:

		float rvec[4];

		mat_to_quat(rotmat, rvec);
		trackball(r, s_old, t_old, s, t);
		add_quats(r, rvec, rvec);
		quat_to_mat(rvec, rotmat);

		s_old = s;
		t_old = t;

		display();
		break;

	case 1:

		translation[0] += (s - s_old);
		translation[1] += (t - t_old);

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

/******************************************************************************
Callback function for mouse clicks
******************************************************************************/

void mouse(int button, int state, int x, int y) {

	int key = glutGetModifiers();

	if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
		
		if (state == GLUT_DOWN) {
			float xsize = (float)win_width;
			float ysize = (float)win_height;

			float s = (2.0 * x - win_width) / win_width;
			float t = (2.0 * (win_height - y) - win_height) / win_height;

			s_old = s;
			t_old = t;

			/*translate*/
			if (button == GLUT_LEFT_BUTTON)
			{
				mouse_mode = 1;
			}

			/*rotate*/
			if (button == GLUT_RIGHT_BUTTON)
			{
				mouse_mode = 2;
			}
		}
		else if (state == GLUT_UP) {

			mouse_mode = -1;
		}
	}
}

/******************************************************************************
Callback function for mouse wheel scroll
******************************************************************************/

void mousewheel(int wheel, int direction, int x, int y) {
	if (direction == 1) {
		zoom *= zoomspeed;
		glutPostRedisplay();
	}
	else if (direction == -1) {
		zoom /= zoomspeed;
		glutPostRedisplay();
	}
}

/******************************************************************************
Callback function for window reshaping
******************************************************************************/

void reshape(int width, int height) {
	win_width = width;
	win_height = height;

	aspectRatio = (float)width / (float)height;

	glViewport(0, 0, width, height);

	set_view(GL_RENDER);
}


/******************************************************************************
Callback function for scene display
******************************************************************************/

void display(void) {
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	set_view(GL_RENDER);
	set_scene(GL_RENDER, mesh);

	/*display the mesh*/
	display_mesh(mesh);

	glFlush();
	glutSwapBuffers();
	glFinish();

	CHECK_GL_ERROR();
}

/******************************************************************************
Diaplay the polygon with visualization results
******************************************************************************/

void display_mesh(Mesh* mesh) {

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

	glBegin(GL_TRIANGLES);
	glColor3f(0.5, 0.5, 0.5);
	for (auto f : mesh->flist) {
		glNormal3d(f->normal.x, f->normal.y, f->normal.z);
		for (int i = 0; i < 3; i++) {
			glVertex3d(f->verts[i]->x, f->verts[i]->y, f->verts[i]->z);
		}

	}

	glEnd();


	/*
	switch (display_mode)
	{
	case 1:	// solid color display with lighting
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat mat_diffuse[4] = { 0.24, 0.4, 0.47, 0.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT, GL_SHININESS, 50.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];
			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}
	}
	break;

	case 2:	// wireframe display
	{
		glDisable(GL_LIGHTING);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.0);

		for (int i = 0; i < poly->nquads; i++) {
			Quad* temp_q = poly->qlist[i];

			glBegin(GL_POLYGON);
			for (int j = 0; j < 4; j++) {
				Vertex* temp_v = temp_q->verts[j];
				glNormal3d(temp_q->normal.entry[0], temp_q->normal.entry[1], temp_q->normal.entry[2]);
				glColor3f(0.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
		}

		glDisable(GL_BLEND);
	}
	break;

	default:
	{
		// don't draw anything
	}

	}

	*/
}
