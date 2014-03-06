/*
CSCI 480
Assignment 3 Raytracer

Name: Jeremy Chen
*/

// libraries
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

// maximum number of triangles, spheres and lights in scene
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex {
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle {
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere {
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light {
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

// camera/screen coordinates (z is -1)
double lookAtY = tan(0.5235981); // (pi/180) * fov/2 
double lookAtX = (WIDTH/(double)HEIGHT)*lookAtY;

// integer representation of shapes
int sphere = -100;
int triangle = -200;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// functions used for vector math
vector<double> add_together(vector<double> first, vector<double> second){ // add first and second vector
  vector<double> result;
  for (int i = 0; i<3; i++) // add for each component of the vectors
    result.push_back(first[i] + second[i]);
  return result; // return a new vector
}

vector<double> subtract(vector<double> first, vector<double> second){ // subtract second from first
  vector<double> result;
  for (int i = 0; i<3; i++) // subtract for each component
    result.push_back(first[i] - second[i]);
  return result;  // return a new vactor
}

vector<double> multiply(vector<double> vect, double x, double y, double z){
  vector<double> product;
  product.push_back(vect[0] * x);
  product.push_back(vect[1] * y);
  product.push_back(vect[2] * z);
  return product; // return a new vector
}

// dot product
double dot(vector<double> first, vector<double> second){
  double value = 0;
  for(int i = 0; i<3; i++)
    value += second[i] * first[i];
  return value;
}

// cross product
vector<double> cross(vector<double> first, vector<double> second){
  vector<double> result;
  for(int i = 0; i<3; i++)
    result.push_back(first[(i+1)%3] * second[(i+2)%3] - first[(i+2)%3] * second[(i+1)%3]);
  return result;
}

// distance formula, used for the unit vector
double distance(vector<double> points){
  double distance = 0;
  for (int i = 0; i < 3; i++)
    distance += pow(points[i], 2); // each component squared
  return sqrt(distance); // then return square root
}

// make into unit vector, using distance function as denominator
vector<double> make_into_unit_vector(vector<double> points){
  double denominator = distance(points); // call distance formula
  vector<double> unit_vector;

  if (denominator > 0)
    for(int i=0; i<3; i++)
      unit_vector.push_back(points[i]/denominator);
  return unit_vector;
}

// class used for rays!
class Ray {
  public: 
    // constructors
    Ray(vector<double>, vector<double>);
    Ray(double, double, double, double, double, double);

    // vectors within Ray
    vector<double> initial; // initial vector point
    vector<double> direction; // direction it will be going
    vector<double> point_of_intersect; // point of intersection with object
    vector<double> normal; // normal
    vector<double> color; // rgb color value

    // get/set functions for private variables
    void set_index_of_intersect(int);
    void set_object_intersect(int);  
    void set_intersect(double);
    int get_index_of_intersect();
    int get_object_intersect();    
    double get_intersect();
    void set_ray(); // initialization function

  private:
    int index_of_intersect;
    int object_intersect;
    double intersect;
};

// constructor when given individual points: initial point (xyz) and end point (x1y1z1)
Ray::Ray(double x, double y, double z, double x1, double y1, double z1){
  // computes direction vector
  this->direction.push_back(x1 - x);
  this->direction.push_back(y1 - y);
  this->direction.push_back(z1 - z);
  this->direction = make_into_unit_vector(this->direction);
  // sets initial vector
  this->initial.push_back(x);
  this->initial.push_back(y);
  this->initial.push_back(z);
  this->set_ray();
}

// constructor when given vectors: initial and end
Ray::Ray(vector<double> initial, vector<double> end){
  this->initial = initial; 
  for (int i=0; i<3; i++)
    this->direction.push_back(end[i] - initial[i]);
  this->direction = make_into_unit_vector(this->direction); 
  this->set_ray();
}

// get/set functions
void Ray::set_index_of_intersect(int index){
  this->index_of_intersect = index;
}

void Ray::set_object_intersect(int object){
  this->object_intersect = object;
}

void Ray::set_intersect(double intersect){
  this->intersect = intersect;
}

int Ray::get_index_of_intersect(){
  return this->index_of_intersect;
}

int Ray::get_object_intersect(){
  return this->object_intersect;
}  

double Ray::get_intersect(){
  return this->intersect;
}

// sets up ray
void Ray::set_ray(){
  this->intersect = 99999; // default set as no intersection
  this->object_intersect = -1; // default set as no object intersecting it
  this->index_of_intersect = -1; // same as above
}

// interpolation function needed to compute illumination for triangles 
vector<double> interpolate(double x, double y, double vec_0[3], double vec_1[3], double vec_2[3]){
  if (x==0){ // if x = 0 just copy vec_0 
    vector<double> interpolated_vector;
    interpolated_vector.assign(vec_0, vec_0+3); // copy vec_0 into interpolated_vector
     
    return interpolated_vector;
  } else { // if x is not equal to 0
    vector<double> temp;
    vector<double> temp1;
    temp.assign(vec_1, vec_1+3);
    temp1.assign(vec_2, vec_2+3);        

    double x_ = x-1;
    double y_ = 1-y;      

    temp1 = multiply(temp1, y, y, y);
    temp = multiply(temp, y_, y_, y_);
    temp = add_together(temp, temp1);
        
    vector<double> interpolated_vector;
    interpolated_vector.assign(vec_0, vec_0+3);
    // make equal to vec_0 but modify according to x and temp vectors (vec_1 and vec_2)
    interpolated_vector = multiply(interpolated_vector, x_, x_, x_);
    interpolated_vector = add_together(interpolated_vector, temp);
    interpolated_vector = multiply(interpolated_vector, 1/x, 1/x, 1/x);    
    
    return interpolated_vector;
  }
}

// find closest object intersection with current ray
double get_intersection(Ray* ray){  
  int sphereIndex = -1; // initialize indexes to nothing
  int triangleIndex = -1;    
  int index = -1;

  double nearest_intersect = 99999; // default: no intersection

  // get closest sphere intersection  
  for (int i = 0; i < num_spheres; i++) {
    vector<double> center; // center of sphere
    center.assign(spheres[i].position, spheres[i].position+3);
    vector<double> sRay = ray->initial; // sphere ray
    sRay = subtract(sRay, center); // subtract center from sphere ray
    
    double result1 = -1; // initialize results to nothing
    double result2 = -1;
    double quadA = -1 * pow(distance(ray->direction), 2);
    double quadB = -2 * dot(ray->direction, sRay);
    double quadC = -1 * (pow(distance(sRay), 2) - pow(spheres[i].radius, 2));
    // quadratic equation, get two results
    if((pow(quadB,2) - 4*quadA*quadC) >= 0){
      result1 = (-1*quadB - sqrt((pow(quadB,2) - 4*quadA*quadC)))/(2*quadA);
      result2 = (-1*quadB + sqrt((pow(quadB,2) - 4*quadA*quadC)))/(2*quadA);
    }     
    
    double possible_intersect = -1;
    if (result1 >= 0 && result2 >= 0){
      if(result1>result2) // if - result is greater than than the + result
        possible_intersect = result2; // make it equal to the + result
      else
        possible_intersect = result1; // otherwise make it equal to the - result

      if(nearest_intersect > possible_intersect) { // if intersect is closer than current closest,
          nearest_intersect = possible_intersect;  // then set as the closest intersect
          sphereIndex = i;         
      }
    }
  } 
  double intersectionWithSphere = nearest_intersect; // set sphere intersection

  // get closest triangle intersection
  nearest_intersect = 99999; // default

  for (int i = 0; i < num_triangles; i++) {
    vector<double> vec_0; // represents vertices of triangle
    vec_0.assign(triangles[i].v[0].position, triangles[i].v[0].position+3);
    vector<double> vec_1;
    vec_1.assign(triangles[i].v[1].position, triangles[i].v[1].position+3);
    vector<double> vec_2;
    vec_2.assign(triangles[i].v[2].position, triangles[i].v[2].position+3);
    
    vector<double> vec_01 = vec_1; // represents edges of triangle
    vector<double> vec_02 = vec_2;

    vec_01=subtract(vec_01, vec_0); // modify to get correct
    vec_02=subtract(vec_02, vec_0);
    
    vector<double> normal = vec_01; // normal vector
    normal=cross(normal, vec_02); 
    normal=make_into_unit_vector(normal);
            
    if (dot(normal, ray->direction)>0)
        normal = multiply(normal, -1, -1, -1); // ray and normal opposites

    double m = -1 * dot(normal, vec_0); // negative dot product of normal and vec_0 
    
    if (dot(normal, ray->direction) != 0) {
      double possible_intersect = -1 * (dot(normal, ray->initial) + m) / dot(normal, ray->direction);            
      
      if (possible_intersect > 0){
        vector<double> point = ray->direction;
        point = multiply(point, possible_intersect, possible_intersect, possible_intersect);
        point = add_together(point, ray->initial);

        vector<double> vec_point = point;

        vec_point=subtract(vec_point, vec_0);

        double a = 0;
        double b = 0;

        double denom = vec_01[0]*vec_02[1] - vec_02[0]*vec_01[1];
        if (vec_01[0]*vec_02[1]==vec_02[0]*vec_01[1]){
          a = -1;
          b = -1;
        } else {
          a = (vec_02[0]*(vec_point[1]*-1) + vec_point[0]*vec_02[1])/denom;
          b = (vec_01[1]*(vec_point[0]*-1) + vec_01[0]*vec_point[1])/denom;
        }
                                
        if (a > 0 && b > 0 && (a + b) < 1) { // after calculating, if a and b are both positive but sum is less than 1..
            if (nearest_intersect > possible_intersect) { // and if determined intersect is closer than current closest
                nearest_intersect = possible_intersect; // make it the closest intersect
                triangleIndex = i;
            }
        }
      }
    }
  }

  double intersectionWithTriangle = nearest_intersect; // intersection closest with a triangle object

  nearest_intersect = 99999; // reset back to default

  if (intersectionWithSphere < intersectionWithTriangle){ // use whichever one is closer
      nearest_intersect = intersectionWithSphere;
      index = sphereIndex; // set index accordingly
  } else {
      nearest_intersect = intersectionWithTriangle;
      index = triangleIndex;
  }

  if (index == sphereIndex || index == triangleIndex) { // if there is indeed an intersection..
    ray->set_index_of_intersect(index); // set ray intersection variables
    ray->set_intersect(nearest_intersect);

    ray->point_of_intersect = ray->direction; // get this intersection point
    ray->point_of_intersect = multiply(ray->point_of_intersect, nearest_intersect, nearest_intersect, nearest_intersect);
    ray->point_of_intersect = add_together(ray->point_of_intersect, ray->initial); // must add to initial ray point
    ray->color.assign(ambient_light, ambient_light+3);             
            
    if (index == sphereIndex){ // if intersecting with sphere
      ray->set_object_intersect(sphere);            

      ray->normal.assign(spheres[index].position, spheres[index].position+3); // get normal at the intersection point
      ray->normal = multiply(ray->normal, -1, -1, -1);
      ray->normal = add_together(ray->normal, ray->point_of_intersect); // point_of_intersect is intersection point
      ray->normal = make_into_unit_vector(ray->normal);

      // scale the color by the diffused color number
      ray->color = multiply(ray->color, spheres[index].color_diffuse[0], spheres[index].color_diffuse[1],
              spheres[index].color_diffuse[2]);                
    } else { // if intersecting with triangle
      ray->set_object_intersect(triangle);
      
      ray->normal = interpolate(0,0, triangles[index].v[0].normal, triangles[index].v[1].normal,
              triangles[index].v[2].normal); // determine triangle normal by interpolation

      vector<double> color_diffuse; // get diffused color coefficient by interpolation
      color_diffuse = interpolate(0,0, triangles[index].v[0].color_diffuse, triangles[index].v[1].color_diffuse, triangles[index].v[2].color_diffuse);
      
      // scale color by diffused color number
      ray->color = multiply(ray->color, color_diffuse[0], color_diffuse[1], color_diffuse[2]);
    }
  }    
  return nearest_intersect;
}

void tracer(Ray* ray){ // tracer function
  double shininess;
  vector<double> diffuse_color (3);
  vector<double> specular_color (3);
  vector<double> diffuse (3);
  vector<double> specular (3);  
  vector<double> ray_direction = ray->direction;
  vector<double> normal = ray->normal;
  
  if (ray->get_object_intersect() == triangle){ // if intersecting with a triangle
    vector<double> vec_0, vec_1, vec_2;
    int index = ray->get_index_of_intersect(); // index is the index of intersection
    // copy the triangle point components into the new created vectors
    vec_0.assign(triangles[index].v[0].position, triangles[index].v[0].position+3);
    vec_1.assign(triangles[index].v[1].position, triangles[index].v[1].position+3);
    vec_2.assign(triangles[index].v[2].position, triangles[index].v[2].position+3);

    vector<double> vec_12 = vec_2;
    vector<double> vec_01 = vec_1;

    vec_12 = subtract(vec_12, vec_1);
    vec_12 = multiply(vec_12, -1, -1, -1);
    vec_01 = subtract(vec_01, vec_0);

    vector<double> point_vec = ray->point_of_intersect;
    point_vec=subtract(point_vec, vec_0); 

    double denom = point_vec[0]*vec_12[1] - vec_12[0]*point_vec[1];
    double interp_a = (vec_01[0]*vec_12[1] - (vec_12[0]*(vec_1[1])))/denom;
    double interp_b = (point_vec[0]*vec_01[1] - point_vec[1]*(vec_01[0]))/denom;

    // interpolate to get the diffuse color and also shininess
    diffuse = interpolate(interp_a, interp_b, triangles[ray->get_index_of_intersect()].v[0].color_diffuse,
          triangles[ray->get_index_of_intersect()].v[1].color_diffuse,
          triangles[ray->get_index_of_intersect()].v[2].color_diffuse);

    shininess = ((interp_a-1) * triangles[ray->get_index_of_intersect()].v[0].shininess) / interp_a;
  
  } else if (ray->get_object_intersect() == sphere){ // if it's sphere... it's a lot simpler
    shininess = spheres[ray->get_index_of_intersect()].shininess; // copy shininess of sphere
    // set the specular and diffuse color
    specular.assign(spheres[ray->get_index_of_intersect()].color_specular, spheres[ray->get_index_of_intersect()].color_specular+3);
    diffuse.assign(spheres[ray->get_index_of_intersect()].color_diffuse, spheres[ray->get_index_of_intersect()].color_diffuse+3);
  }

  ray_direction = multiply(ray_direction, -1, -1, -1);
  ray_direction = make_into_unit_vector(ray_direction);
  normal = make_into_unit_vector(normal);

  for (int i = 0; i < num_lights; i++) { // light logic  
    vector<double> light_position;
    light_position.assign(lights[i].position, lights[i].position+3); // set the light position based on current light
    Ray* checkShadow = new Ray(ray->point_of_intersect, light_position); // set ray to check if it's in shadow
    bool in_shadow = false;   // default is it's not in a shadow
          
    double nearest_intersect = get_intersection(checkShadow); // call get_intersection
    
    delete(checkShadow);

    vector<double> d_to_l;
    d_to_l.assign(lights[i].position, lights[i].position+3);
    d_to_l=subtract(d_to_l, ray->point_of_intersect);
    double distance_to_light = distance(d_to_l); // use and check if it's larger than the nearest_intersect

    if (nearest_intersect<distance_to_light && nearest_intersect != 99999)
      in_shadow = true; // if it's intersecting something then it's in shadow   

    vector<double> light_point;
    light_point.assign(lights[i].position, lights[i].position+3); // copy light position
    light_point = subtract(light_point, ray->point_of_intersect); // subtract point of intersection
    light_point = make_into_unit_vector(light_point);      // unit vector       

    double n_dot_l = dot(normal, light_point); // get dot product of normal and light point
    double r_dot_v = dot(ray->normal, ray_direction); // also get dot product of normal and ray direction

    if (!in_shadow) {       // if it's not in shadow, then give it color based on diffuse color and specular color          
      for (int j = 0; j < 3; j++) {
        diffuse_color[j] += lights[i].color[j] * diffuse[j] * n_dot_l; // alter with n_dot_l
        specular_color[j] += lights[i].color[j] * specular[j] * pow(r_dot_v, shininess); // give it shininess
      }
    }
  }
  ray->color = add_together(ray->color, diffuse_color); // add to the ray color
  ray->color = add_together(ray->color, specular_color); // add to the ray color
}

//MODIFY THIS FUNCTION
void draw_scene() {
  int newWidth = WIDTH/2;
  int newHeight = HEIGHT/2;
  int white = 255;

  for(int x=-newWidth; x<newWidth; x++) {
    glPointSize(2);  
    glBegin(GL_POINTS);
    for(int y=-newHeight; y<newHeight; y++) {
      Ray* ray = new Ray(0, 0, 0, x*lookAtX/(newWidth), y*lookAtY/(newHeight), -1); // make a new ray
      if(get_intersection(ray) != 99999) // if there is an intersection
        tracer(ray); // trace!
      if(ray->get_intersect() == 99999) // if there is no intersection
        plot_pixel(x, y, white, white, white);
      else { // if there is an intersection
        ray->color = multiply(ray->color, white, white, white);
        plot_pixel(x, y, ray->color[0], ray->color[1], ray->color[2]);
      }
      delete(ray); // delete
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b) {
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x+320, y+240);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b) {
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b) {
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg() {
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      
}

void parse_check(char *expected,char *found) {
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }
}

void parse_doubles(FILE*file, char *check, double p[3]) {
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r) {
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi) {
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++) {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);

    if(strcasecmp(type,"triangle")==0){
      printf("found triangle\n");
      int j;

      for(j=0;j < 3;j++){
          parse_doubles(file,"pos:",t.v[j].position);
          parse_doubles(file,"nor:",t.v[j].normal);
          parse_doubles(file,"dif:",t.v[j].color_diffuse);
          parse_doubles(file,"spe:",t.v[j].color_specular);
          parse_shi(file,&t.v[j].shininess);
        }

      if(num_triangles == MAX_TRIANGLES) {
          printf("too many triangles, you should increase MAX_TRIANGLES!\n");
          exit(0);
        }
      triangles[num_triangles++] = t;

    } else if (strcasecmp(type,"sphere")==0) {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES) {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;

    } else if (strcasecmp(type,"light")==0) {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS) {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;

    } else {
    printf("unknown type in scene description:\n%s\n",type);
    exit(0);
    }
  }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();

      if(mode == MODE_JPEG)
  save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}