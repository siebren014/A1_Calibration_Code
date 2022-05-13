#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;


bool isvalid(const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d){

    if (points_3d.size() >= 6 && points_3d.size() == points_2d.size()){
        return true;
    }
    else{
        return false;
    }
}

bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{


//    check whether the input data is valid
    if (!isvalid(points_3d, points_2d)){
        throw std::invalid_argument( "invalid input" );
    }

    std::vector<double> all_point;

    int j = 0;
    for (const auto & point: points_3d ){
        all_point.push_back(point.x());
        all_point.push_back(point.y());
        all_point.push_back(point.z());
        all_point.push_back(points_2d[j].x());
        all_point.push_back(points_2d[j].y());
        j = j +1;
    }

    Matrix mat = Matrix(points_3d.size(), 5, all_point);

    int m = mat.rows()*2;
    int n = 12;

    Matrix A = Matrix(m, n, 0.0);
    int ii = 0;
    for (int i = 0; i < mat.rows(); i++) {
        A[ii][0] = mat[i][0];
        A[ii][1] = mat[i][1];
        A[ii][2] = mat[i][2];
        A[ii][3] = 1.0;
        A[ii][4] = 0.0;
        A[ii][5] = 0.0;
        A[ii][6] = 0.0;
        A[ii][7] = 0.0;
        A[ii][8] = -1.0*mat[i][0]*mat[i][3];
        A[ii][9] = -1.0*mat[i][1]*mat[i][3];
        A[ii][10] = -1.0*mat[i][2]*mat[i][3];
        A[ii][11] = -1.0*mat[i][3];
        ii += 1;
        A[ii][0] = 0.0;
        A[ii][1] = 0.0;
        A[ii][2] = 0.0;
        A[ii][3] = 0.0;
        A[ii][4] = mat[i][0];
        A[ii][5] = mat[i][1];
        A[ii][6] = mat[i][2];
        A[ii][7] = 1.0;
        A[ii][8] = -1.0*mat[i][0]*mat[i][4];
        A[ii][9] = -1.0*mat[i][1]*mat[i][4];
        A[ii][10] = -1.0*mat[i][2]*mat[i][4];
        A[ii][11] = -1.0*mat[i][4];
        ii += 1;
    }

    Vector M = Vector(n, 0.0);
    Matrix U = Matrix(m,m,0.0);
    Matrix S = Matrix(m,n, 0.0);
    Matrix V = Matrix(n,n,0.0);

    svd_decompose(A, U, S, V);
    for (int i = 0; i < n; i++) {
        M[i] = V[i][n-1];
    }

    Vector N = Vector(m,5.0);
    N = mult(A,M);

    Matrix34 MM = Matrix(3,4,0.0);
    MM[0][0]=M[0];
    MM[0][1]=M[1];
    MM[0][2]=M[2];
    MM[0][3]=M[3];
    MM[1][0]=M[4];
    MM[1][1]=M[5];
    MM[1][2]=M[6];
    MM[1][3]=M[7];
    MM[2][0]=M[8];
    MM[2][1]=M[9];
    MM[2][2]=M[10];
    MM[2][3]=M[11];


    Vector a1 = Vector3D(MM[0][0], MM[0][1], MM[0][2]);
    Vector a2 = Vector3D(MM[1][0], MM[1][1], MM[1][2]);
    Vector a3 = Vector3D(MM[2][0], MM[2][1], MM[2][2]);
    Vector b = Vector3D(MM[0][3], MM[1][3], MM[2][3]);

    double ro;
    double u0;
    double v0;
    double alpha;
    double beta;
    ro = 1.0/(a3.norm());

    u0 = pow(ro,2)*dot(a1,a3);
    v0 = pow(ro,2)*dot(a2,a3);

    double costheta = -((dot(cross(a1,a3),(cross(a2,a3))))/(((cross(a1,a3)).norm())*((cross(a2,a3)).norm())));
    double theta = acos(costheta);
    double sintheta = sin(theta);
    alpha = (pow(ro,2))*(norm(cross(a1,a3)))*sintheta;
    beta = (pow(ro,2))*(norm(cross(a2,a3)))*sintheta;


    Matrix33 K = (Matrix33(3,3,0.0));
    K[0][0]= alpha;
    K[0][1]= -alpha*(costheta/sintheta);
    K[0][2]= u0;
    K[1][0]= 0;
    K[1][1]= beta/sintheta;
    K[1][2]= v0;
    K[2][0]= 0;
    K[2][1]= 0;
    K[2][2]= 1;


    Vector3D r1 = (cross(a2,a3))/(norm(cross(a2,a3)));
    Vector3D r3 = ro*(a3);
    Vector3D r2 = cross(r3,r1);

    R.set_row(0,{r1[0], r1[1], r1[2]});
    R.set_row(1,{r2[0], r2[1], r2[2]});
    R.set_row(2,{r3[0], r3[1], r3[2]});

    t = ro* mult(inverse(K),b);

    //calibration validation parameters
    double cotheta = costheta / sintheta;
    skew = -alpha * cotheta;
    cx = u0;
    cy = v0;
    fx = alpha;
    fy = beta;

    Vector4D RR1= Vector4D(r1[0], r1[1], r1[2],t[0]);
    Vector4D RR2= Vector4D(r2[0], r2[1], r2[2],t[1]);
    Vector4D RR3= Vector4D(r3[0], r3[1], r3[2],t[2]);

    Vector4D P1 = (Vector4D(8.0,7.0,0.0,1.0));
    Vector4D P2 = (Vector4D(6.0,0.0,5.0,1.0));
    Vector4D P3 = (Vector4D(0.0,6.0,8.0,1.0));
    Vector4D P4 = (Vector4D(2.0,5.0,0.0,1.0));
    Vector4D P5 = (Vector4D(0.0,3.0,7.0,1.0));
    Vector4D P6 = (Vector4D(1.0,0.0,8.0,1.0));

    Vector3D Uit1 = mult(MM,P1);
    Vector3D Uit2 = mult(MM,P2);
    Vector3D Uit3 = mult(MM,P3);
    Vector3D Uit4 = mult(MM,P4);
    Vector3D Uit5 = mult(MM,P5);
    Vector3D Uit6 = mult(MM,P6);

    double coord1 [2] = {Uit1[0]/Uit1[2], Uit1[1]/Uit1[2]};
    double coord2 [2] = {Uit2[0]/Uit2[2], Uit2[1]/Uit2[2]};
    double coord3 [2] = {Uit3[0]/Uit3[2], Uit3[1]/Uit3[2]};
    double coord4 [2]= {Uit4[0]/Uit4[2], Uit4[1]/Uit4[2]};
    double coord5 [2] = {Uit5[0]/Uit5[2], Uit5[1]/Uit5[2]};
    double coord6 [2] = {Uit6[0]/Uit6[2], Uit6[1]/Uit6[2]};

    // check accuracy by backcalcuting the original points by multiplying the transformation matrix with the original points.
    std::cout<<"Coord 1:"<<"u = " <<coord1[0] <<" v = "<<coord1[1]<< std::endl;
    std::cout<<"Coord 2:"<<"u = " <<coord2[0] <<" v = "<<coord2[1]<< std::endl;
    std::cout<<"Coord 3:"<<"u = " <<coord3[0] <<" v = "<<coord3[1]<< std::endl;
    std::cout<<"Coord 4:"<<"u = " <<coord4[0] <<" v = "<<coord4[1]<< std::endl;
    std::cout<<"Coord 5:"<<"u = " <<coord5[0] <<" v = "<<coord5[1]<< std::endl;
    std::cout<<"Coord 6:"<<"u = " <<coord6[0] <<" v = "<<coord6[1]<< std::endl;


    return true;
}

