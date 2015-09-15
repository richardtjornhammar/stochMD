double pbc_dx( t_frame, int , int , vector );
int    gen_velocities( t_frame );
int    check_frame( t_frame, int );
int    grand_canonical(t_frame,t_files);
double calc_energy(t_frame, int);
double calc_angle( vector, vector, double* );
double spec_angle(vector a, vector b);
double iprod( vector, vector );
double cprod( vector, vector, vector );
int    vector_inc(vector, vector);
int    vector_dec(vector, vector);
int    vector_sub(vector, vector, vector);
int    vector_add(vector, vector, vector);
int    svmul(double, vector, vector);
double dih_angle(int, int, int, int, t_frame, vector, vector, vector, vector, vector, double* );
int    do_dih_fup(int, int ,int ,int ,double , vector, vector, vector, vector, vector, t_frame);
double dopdihs(double, double, double, double, int mult, double, double, double*, double* );
double idihs(t_frame, int, int);
double pdihs(t_frame, int, int);
int    print_molecule( t_molecule, t_files );
double InvSqrt( double );
int    quaternion( vector, vector, vector, double, vector );
int    write_data(t_frame,t_files, double);

#define CUBE(X)   (X)*(X)*(X)
#define SQUARE(X) (X)*(X)
#define FLOAT_EPS    5.96046448E-08
#define DOUBLE_EPS   1.11022302E-16
