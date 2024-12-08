#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "libraries.h"
#include "ant.h"


using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Polygon = CGAL::Polygon_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
using Segment_2 = K::Segment_2;
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;
using std_string = std::string;
typedef K::FT FT;

 

//Calculate squared distance between two points
double squared_distance_1(const Point_2& a, const Point_2& b);

//Check if a triangle is obtuse
bool is_obtuse(const Face_handle& face);
bool is_obtuse2(const Point_2& p1, const Point_2& p2, const Point_2& p3);

//Read JSON file
void read_json(const std::string& filename, value& jv);

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt, const Polygon& polygon);

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

//Return true if approves the flip
bool is_it_worth_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

void start_the_flips(Custom_CDT& cdt, const Polygon& polygon);

Point_2 find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3);

//Steiner methods
//void insert_circumcenter_centroid(Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon, Point_2& circum_or_centroid);
void insert_projection(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_projection, Segment_2& opposide_edge);
void insert_midpoint(Custom_CDT& custom_cdt, const Face_handle& face, Polygon& polygon, Point_2& in_midpoint, Segment_2& longest_edge);
bool insert_adjacent_steiner(Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon, Point_2& adjacent_steiner);
void insert_adjacent_steiner_local_search(Custom_CDT& custom_cdt, const Face_handle& face1, const Polygon& polygon, Point_2& adjacent_steiner);
bool insert_circumcenter(Custom_CDT& circumcenter_cdt, const Face_handle& face, const Polygon& polygon, Point_2& circumcenter_steiner);
void insert_centroid(Custom_CDT& centroid_cdt, const Face_handle& face, const Polygon& polygon, Point_2& centroid_steiner);

//Algorithms
void local_search(Custom_CDT& custom_cdt, Polygon& polygon, int& L);
void simulated_annealing(Custom_CDT& custom_cdt, Polygon& polygon, int max_iterations, const double& alpha, const double& beta, const int& batch_size);
void ant_colony(Custom_CDT& custom_cdt, Polygon& polygon, const double& alpha, const double& beta, const double& chi, const double& psi, const double& lamda, const int& L, const int& kappa);

bool should_accept_transition(const double deltaE,const double T);
double calculate_energy(const int obtuse_faces, const int steiner_points, const double alpha, const double beta);

//If the steiner inserted in the boundary of polygon, update the polygon.
//Take the steiner, and the edge (point1, point2) that we insert the steiner in cdt
void update_polygon(Polygon& polygon, const Point_2& steiner_point, const Point_2& point1, const Point_2& point2);

bool is_point_inside_region(const Point_2& point, const Polygon& polygon);

bool is_face_inside_region(const Face_handle& face, const Polygon& polygon);

//If and edge of the face lies into the border of the boundary, => face == face of boundary
//bool is_face_on_boundary(const Custom_CDT& cdt, Face_handle face);

bool is_face_on_boundary(const Custom_CDT& cdt, Face_handle face, const Polygon& polygon);

bool is_edge_inside_region(const Point_2& point1, const Point_2& point2, const Polygon& polygon);

bool is_edge_on_boundary(const Point_2& p1, const Point_2& p2, const Polygon& polygon);

bool worth_insert_centroid(Custom_CDT& custom_cdt, const Point_2& centroid, const Polygon& polygon);

Segment_2 find_longest_edge(const Face_handle& face);

//JSON OUTPUT METHODS 
bool is_steiner_point(Vertex_handle vertex, const vector<Point_2>& original_points);

void output_file(value jv, Custom_CDT custom_cdt, vector<Point_2> points, int obtuse_count, std_string output_path);

std::string convert_to_string(const FT& coord);
std::string format_double(double value);

pair<Point_2, double> find_obtuse_vertex_and_angle(const Point_2& v1, const Point_2& v2, const Point_2& v3);
Point_2 compute_centroid(const vector<Point_2>& points);
int count_vertices(const Custom_CDT& cdt);
void print_polygon_edges(const Polygon& polygon);
bool is_circumcenter_in_neighbor(const Custom_CDT& cdt, const Face_handle& face, const Point_2& circumcenter);
double calculate_radius_to_height(const Face_handle& face, const Custom_CDT& cdt);
double hta_vertex_projection(double rho);
double hta_circumcenter(double rho);
double hta_midpoint(double rho);
double hta_mean_adjacent(bool has_obtuse_neighbors);
bool has_obtuse_neighbors(const Custom_CDT& custom_cdt, const Face_handle& face, const Polygon& polygon);
Face_handle give_random_obtuse(Custom_CDT& custom_cdt, Polygon& polygon);
SteinerMethod selectSteinerMethod(const double& ro, const vector<double>& taf, vector<double>& hta, double chi, double psi, bool obtuse_neighbors);
void updatePheromones(vector<double>& taf, vector<double>& delta_taf, vector<Ant> selected_ants, double lamda);
bool is_polygon_convex(const std::set<Point_2>& unique_points);
bool are_faces_equal(const Face_handle& face1, const Face_handle& face2);
vector<Ant> save_the_best(vector<Ant>& ants);
//Check for conflict between 2 ants
bool have_conflict(Ant& ant1, Ant& ant2);
void printAntDetails(vector<Ant>& ants);
void affected_faces(Custom_CDT& best_cdt, Ant& ant);

#endif