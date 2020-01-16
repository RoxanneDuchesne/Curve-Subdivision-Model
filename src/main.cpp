#include <array>
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <cassert>	//assert

 // glad beforw glfw
#include "glad/glad.h"
//
#include "GLFW/glfw3.h"

#include "common_matrices.hpp"
//#include "obj_mesh_file_io.hpp"
//#include "obj_mesh.hpp"
#include "mat4f.hpp"
#include "mat3f.hpp"
#include "shader.hpp"
#include "program.hpp"
#include "triangle.hpp"
#include "vec3f.hpp"
#include "vec2f.hpp"
#include "shader_file_io.hpp"
#include "buffer_object.hpp"
#include "vertex_array_object.hpp"
#include "vbo_tools.hpp"

//#include "texture.hpp"
//#include "image.hpp"

using namespace math;
using namespace geometry;
using namespace opengl;

// GLOBAL Variables
Mat4f g_M = Mat4f::identity();
Mat4f g_V = Mat4f::identity();
Mat4f g_P = Mat4f::identity();

GLuint g_width = 1000, g_height = 1000;

// 1 = Left Pressed
// 2 = Right Pressed
int pressed = 0;

// 1 = click
int click = 0; 

int closest_point = -1;
double closest_distance = -1;

double clicked_x = -1;
double clicked_y = -1;

double x_1 = -1;
double y_1 = -1;

int depth = 1;

std::vector<Vec3f> controlPoints;
std::vector<Vec3f> outCurve;

// function declaration
using namespace std;

struct OBJMesh {
    IndicesTriangles triangles;
    Vertices vertices;
    Normals normals;
};

void setFrameBufferSize(GLFWwindow *window, int width, int height) {
	g_width = width;
	g_height = height;
    glViewport(0, 0, g_width/2, g_height);
	g_P = perspectiveProjection(30, float(g_width) / g_height, 0.01, 100.f);
}

static void cursorPositionCallback(GLFWwindow* window, double xPos, double yPos) {

    x_1 = (xPos - ((float)g_width * (3.f/4.f)))/((float)g_width/4.f);
    y_1 = (yPos - ((float)g_height)/2.f)/(-(float)g_height/2.f);

    // Get Closest Point
    if (click == 1){

        clicked_x = x_1;
        clicked_y = y_1;

        if ((x_1 <= 1.f) && (x_1 >= -1.f) && (y_1 <= 1.f) && (y_1 >= -1.f)){

                // Find closest point
                closest_point = -1;
                closest_distance = -1;

                for (int i = 0; i < controlPoints.size(); i++){

                    double cur_distance = pow(pow((x_1 - controlPoints[i].x), 2) + pow((y_1 - controlPoints[i].y), 2), 0.5f);

                    if (closest_distance == -1){
                        closest_distance = cur_distance;
                        closest_point = i;
                    }
                    else if (closest_distance > cur_distance){
                        closest_distance = cur_distance;
                        closest_point = i;
                    }
                }
        }
        click = 0;
    }
    // Move Point
    if (pressed == 1){

        double movement_x = x_1 - clicked_x;
        double movement_y = y_1 - clicked_y;

        controlPoints[closest_point].x = x_1;
        controlPoints[closest_point].y = y_1;
    }


}


static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mod) {

    if ((button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS))
    {
        click = 1;
        pressed = 1;
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        click = 0;
        pressed = 0;
    }
    if ((button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS))
    {
        // Remove Point
            clicked_x = x_1;
            clicked_y = y_1;

            if ((x_1 <= 1.f) && (x_1 >= -1.f) && (y_1 <= 1.f) && (y_1 >= -1.f)){

                    // Find closest point
                    closest_point = -1;
                    closest_distance = -1;

                    for (int i = 0; i < controlPoints.size(); i++){

                        double cur_distance = pow(pow((x_1 - controlPoints[i].x), 2) + pow((y_1 - controlPoints[i].y), 2), 0.5f);

                        if (closest_distance == -1){
                            closest_distance = cur_distance;
                            closest_point = i;
                        }
                        else if (closest_distance > cur_distance){
                            closest_distance = cur_distance;
                            closest_point = i;
                        }
                    }
            }
            if (controlPoints.size() > 0){
                controlPoints.erase(controlPoints.begin() + closest_point);
            }
    }
}

void setKeyboard(GLFWwindow *window, int key, int scancode, int action,
	int mods) {

	// Students should try and fix:
	// g_M holds the transformation of the model
	// we keep multiplying on the left or right with new matrices
	// What effect does this have on translation/rotation/scale

	if (GLFW_KEY_LEFT == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = rotateAboutYMatrix(5) * g_M;
		}
	}
	else if (GLFW_KEY_RIGHT == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = rotateAboutYMatrix(-5) * g_M;
		}
	}
	else if (GLFW_KEY_UP == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = uniformScaleMatrix(1.1) * g_M;
		}
	}
	else if (GLFW_KEY_DOWN == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = uniformScaleMatrix(1. / 1.1) * g_M;
		}
	}
	else if (GLFW_KEY_W == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = translateMatrix(0, 0.1, 0) * g_M;
		}
	}
	else if (GLFW_KEY_S == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = translateMatrix(0, -0.1, 0) * g_M;
		}
	}
	else if (GLFW_KEY_D == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = translateMatrix(0.1, 0, 0) * g_M;
		}
	}
	else if (GLFW_KEY_A == key) {
		if (GLFW_REPEAT == action || GLFW_PRESS == action) {
			g_M = translateMatrix(-0.1, 0, 0) * g_M;
		}
	}
	else if (GLFW_KEY_ESCAPE == key) {
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
    // Add Point
    else if (GLFW_KEY_P == key) {
        if (GLFW_REPEAT == action || GLFW_PRESS == action) {
            controlPoints[controlPoints.size()];
            controlPoints.push_back({controlPoints[controlPoints.size() - 1].x - 0.05f, controlPoints[controlPoints.size() - 1].y + 0.1f, 0});
        }
    }
    else if (GLFW_KEY_PERIOD == key) {
        if (GLFW_REPEAT == action || GLFW_PRESS == action) {
            if (depth < 10){
                depth++;
            }
        }
    }
    else if (GLFW_KEY_COMMA == key) {
        if (GLFW_REPEAT == action || GLFW_PRESS == action) {
            if (depth > 0){
                depth--;
            }
        }
    }
}

// user defined alias
opengl::Program createShaderProgram(std::string const &vertexShaderFile,
	std::string const &fragmentShaderFile) {
	using namespace opengl;
	auto vertexShaderSource = loadShaderStringFromFile(vertexShaderFile);
	auto fragmentShaderSource = loadShaderStringFromFile(fragmentShaderFile);

	std::cout << "[Log] compiling program " << vertexShaderFile << ' '
		<< fragmentShaderFile << '\n';
	return opengl::makeProgram(vertexShaderSource, fragmentShaderSource);
}

std::string glfwVersion() {
	std::ostringstream s("GLFW version: ", std::ios::in | std::ios::ate);
	// print version
	int glfwMajor, glfwMinor, glfwRevision;
	glfwGetVersion(&glfwMajor, &glfwMinor, &glfwRevision);
	s << glfwMajor << '.' << glfwMinor << '.' << glfwRevision;
	return s.str();
}

bool loadGeometryToGPU(std::vector<Vec3f> const &vertices, GLuint vboID) {
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
	glBufferData(
		GL_ARRAY_BUFFER,                  // destination
		sizeof(Vec3f) * vertices.size(), // size (bytes) of memcopy to GPU
		vertices.data(),                  // pointer to data (contiguous)
		GL_STATIC_DRAW                    // usage patter of the GPU buffer
	);

	return true;
}

// Lane-Risenfield Subdivison
std::vector<Vec3f> subdivideOpenCurve(std::vector<Vec3f> const & points, int depth) {

    std::vector<Vec3f> out;

    if (depth == 0){
        return points;
    }

    for (int i=0; i < depth; i++){
        std::vector<Vec3f> ref_points;
        if (i == 0){
            for (int k=0; k < points.size(); k++){
                ref_points.push_back(points[k]);
            }
        }
        else{
            for (int k=0; k < out.size(); k++){
                ref_points.push_back(out[k]);
            }
        }

        out.clear();

        for (int j=0; j < ref_points.size(); j++){

            if (j == 0){
                out.push_back(ref_points[j]);
                out.push_back(((1.f/2.f) * ref_points[j]) + ((1.f/2.f) * ref_points[j+1]));

            }
            else if (j == (ref_points.size() - 1)){
              out.push_back(ref_points[(ref_points.size() - 1)]);
            }
            else{
                out.push_back(((1.f/8.f) * ref_points[j - 1]) + ((3.f/4.f) * ref_points[j]) + ((1.f/8.f) * ref_points[j+1]));
                out.push_back(((1.f/2.f) * ref_points[j]) + ((1.f/2.f) * ref_points[j+1]));
            }

        }
    }


	return out;
}

// Created Rotated Points Based on Subdivided Curve
std::vector<Vec3f> generateRotatedPoints(std::vector<Vec3f> const & curve){

    std::vector<Vec3f> rotatePoints;

    int size = curve.size();

    // Rotated Points
    for (int i=0; i < 72; i++){
        for (int j=0; j < size; j++){
            rotatePoints.push_back(rotateAroundAxis(curve[j], Vec3f(0.f,1.f,0.f), i * 5.f));
        }
    }

    return rotatePoints;
}

// Create Trianlges Based on Rotated Points
std::vector<Vec3f> generateTrianlges(std::vector<Vec3f> const & curve, std::vector<Vec3f> rotatePoints){
    int size = curve.size();

    std::vector<Vec3f> trianglePoints;


    // Create Triangle
    for (int i=0; i < 72; i++){

       for (int j=0; j < (size - 1); j++){
           if (i == 71){
               trianglePoints.push_back(rotatePoints[j + (i * size)]);
               trianglePoints.push_back(rotatePoints[j + (i * size) + 1]);
               trianglePoints.push_back(rotatePoints[j]);

               trianglePoints.push_back(rotatePoints[j]);
               trianglePoints.push_back(rotatePoints[j + (i * size) + 1]);
               trianglePoints.push_back(rotatePoints[j + 1]);
           }
           else{
             trianglePoints.push_back(rotatePoints[j + (i * size)]);
             trianglePoints.push_back(rotatePoints[j + (i * size) + 1]);
             trianglePoints.push_back(rotatePoints[j + ((i + 1) * size)]);

             trianglePoints.push_back(rotatePoints[j + ((i + 1) * size)]);
             trianglePoints.push_back(rotatePoints[j + (i * size) + 1]);
             trianglePoints.push_back(rotatePoints[j + ((i + 1) * size) + 1]);
           }
       }
    }
    return trianglePoints;
}

// Creates vector of IndicesTriangles out of Vector of Triangle Vertices
std::vector<IndicesTriangle> generateIndices(std::vector<Vec3f> const & tri){

    std::vector<IndicesTriangle> indtri;

    for (int i = 0; i < tri.size()-2; i += 3){

        geometry::Indices tmp_1;
        tmp_1.vertexID() = i;
        tmp_1.textureCoordID() = 0;
        tmp_1.normalID() = i;

        geometry::Indices tmp_2;
        tmp_2.vertexID() = i + 1;
        tmp_2.textureCoordID() = 0;
        tmp_2.normalID() = i + 1;

        geometry::Indices tmp_3;
        tmp_3.vertexID() = i + 2;
        tmp_3.textureCoordID() = 0;
        tmp_3.normalID() = i + 2;

        geometry::Triangle_<geometry::Indices> next({tmp_1, tmp_2, tmp_3});
        indtri.push_back(next);
    }
    return indtri;
}

void setupVAO(GLuint vaoID, GLuint vboID) {

	glBindVertexArray(vaoID);
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
	
	// set up position input into vertex shader
	glEnableVertexAttribArray(0);          // match layout # in shader
	glVertexAttribPointer(                 //
		0,                                 // attribute layout # (in shader)
		3,                                 // number of coordinates per vertex
		GL_FLOAT,                          // type
		GL_FALSE,                          // normalized?
		sizeof(Vec3f),                    // stride
		0									 // array buffer offset
	);

	glBindVertexArray(0);
}

GLFWwindow *initWindow() {
	GLFWwindow *window = nullptr;

	if (!glfwInit()) {
		exit(EXIT_FAILURE);
	}

	std::cout << glfwVersion() << '\n';

	// set opengl version
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
	glfwWindowHint(GLFW_SAMPLES, 4);

	window = glfwCreateWindow(1000,           // width
		1000,           // height
		"Mesh Loading", // title
		NULL,           // unused
		NULL);
	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window);
	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

	glfwSwapInterval(1); // vsync
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_DEPTH_TEST);
	// glEnable(GL_CULL_FACE);
	// glCullFace(GL_BACK);

	//Polygon fill mode
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// setup callbacks
    glfwSetFramebufferSizeCallback(window, setFrameBufferSize);
	glfwSetKeyCallback(window, setKeyboard);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

	return window;
}

int main() {
    GLFWwindow *window = initWindow();

    auto vao_control = makeVertexArrayObject();
    auto vbo_control = makeBufferObject();
    auto vao_curve = makeVertexArrayObject();
    auto vbo_curve = makeBufferObject();
    auto vbo_vert = makeBufferObject();

    GLuint totalIndices = 0;

	Vec3f viewPosition(0, 0, 3);
	g_V = lookAtMatrix(viewPosition,    // eye position
		{ 0.f, 0.f, 0.f }, // look at
		{ 0.f, 1.f, 0.f }  // up vector
	);
	g_P = orthographicProjection(-1, 1, 1, -1, 0.001f, 10);

    auto phongShader = createShaderProgram("./shaders/phong_vs.glsl",
                                          "./shaders/phong_fs.glsl");

    assert(phongShader);

    phongShader.use();
    setUniformVec3f(phongShader.uniformLocation("lightPosition"), {0, 2, 5});
    setUniformVec3f(phongShader.uniformLocation("viewPosition"), viewPosition);

    setupVAO(vao_control.id(), vbo_control.id());
    setupVAO(vao_curve.id(), vbo_curve.id());

    //Set to one shader program
    opengl::Program *program = &phongShader;

	//Load control points
    controlPoints.push_back({ 0.3,-0.5,0 });
    controlPoints.push_back({ 0.6,-0.25,0 });
    controlPoints.push_back({ 0.8,0.0,0 });
    controlPoints.push_back({ 0.6,0.25,0 });
    controlPoints.push_back({ 0.3,0.5,0 });

    glPointSize(10);


    while (!glfwWindowShouldClose(window)) {


        outCurve = subdivideOpenCurve(controlPoints, depth);
        // Points rotated around axis 5 degrees at a time
        std::vector<Vec3f> rotated = generateRotatedPoints(outCurve);
        // Create Triangle Vertices
        std::vector<Vec3f> tris = generateTrianlges(outCurve, rotated);

        geometry::OBJMesh meshData;
        // Put Trianlge Vertices into Proper Structure for MeshData
        meshData.triangles = generateIndices(tris);
        meshData.vertices = tris;
        auto normals = geometry::calculateVertexNormals(meshData.triangles, meshData.vertices);

        auto vboData = opengl::makeConsistentVertexNormalIndices(meshData, normals);

        totalIndices = opengl::setup_vao_and_buffers(vao_curve, vbo_curve, vbo_vert, vboData);

        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

        program->use();
        setUniformMat4f(program->uniformLocation("model"), g_M, true);
        setUniformMat4f(program->uniformLocation("view"), g_V, true);
        setUniformMat4f(program->uniformLocation("projection"), g_P, true);

        //Draw Vase
        glViewport(0, 0, g_width/2, g_height);
        vao_curve.bind();

        glDrawElements(GL_TRIANGLES,
                       totalIndices,
                       GL_UNSIGNED_INT,
                       (void *)0);


        //Draw control points
        glViewport(g_width/2, 0, g_width/2, g_height);
		
        loadGeometryToGPU(controlPoints, vbo_control.id());
		vao_control.bind();
		glDrawArrays(GL_LINE_STRIP,   // type of drawing (rendered to back buffer)
			0,						  // offset into buffer
			controlPoints.size()	// number of vertices in buffer
		);

		glDrawArrays(GL_POINTS,   // type of drawing (rendered to back buffer)
			0,						  // offset into buffer
			controlPoints.size()	// number of vertices in buffer
        );


		glfwSwapBuffers(window); // swaps back buffer to front for drawing to screen
		glfwPollEvents();        // will process event queue and carry on
	}

	// cleaup window, and glfw before exit
	glfwDestroyWindow(window);
	glfwTerminate();

	return EXIT_SUCCESS;
}
