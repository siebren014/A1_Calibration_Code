#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;

// Check to see if the input data is valid by checking if it has at least 6 points and the 3d and 2d vectors have the same size.
bool isvalid(const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d){

    if (points_3d.size() >= 6 && points_3d.size() == points_2d.size()){
        return true;
    }
    else{
        return false;
    }
}

//calibration function
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{

    //check if the input data is valid, esle throw an error
    if (!isvalid(points_3d, points_2d)){
        throw std::invalid_argument( "invalid input" );
    }

    //vector for all points
    std::vector<double> all_point;

    // push all the points in the all_point vector
    int j = 0;
    for (const auto & point: points_3d ){
        all_point.push_back(point.x());
        all_point.push_back(point.y());
        all_point.push_back(point.z());
        all_point.push_back(points_2d[j].x());
        all_point.push_back(points_2d[j].y());
        j = j +1;
    }

    //matrix with the 3d and 2d points
    Matrix mat = Matrix(points_3d.size(), 5, all_point);
    int m = mat.rows()*2;
    int n = 12;

    //creating and filling the A matrix
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

    //Check if the data is not too noisy, if it is, throw an error.
    for (int i = 0; i < N.size(); i++){
        if (N[i]> 0.0001){
            throw std::invalid_argument( "input data too noisy" );
        }
    }

    //creating and filling MM matrix
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

    double ro = 0.0;
    double u0 = 0.0;
    double v0 = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    ro = 1.0/(a3.norm());

    u0 = pow(ro,2)*dot(a1,a3);
    v0 = pow(ro,2)*dot(a2,a3);


    double costheta = -((dot(cross(a1,a3),(cross(a2,a3))))/(((cross(a1,a3)).norm())*((cross(a2,a3)).norm())));
    double theta = acos(costheta);
    double sintheta = sin(theta);
    alpha = pow(ro,2)*norm(cross(a1,a3))*sintheta;
    beta = pow(ro,2)*norm(cross(a1,a3))*sintheta;

    Matrix33 K = (3,3,0.0);
    K[0][0]= alpha;
    K[0][1]= -alpha*(costheta*sintheta);
    K[0][2]= u0;
    K[1][0]= 0;
    K[1][1]= beta/sintheta;
    K[1][2]= v0;
    K[2][0]= 0;
    K[2][1]= 0;
    K[2][2]= 1;

    Vector3D r1 = (cross(a2,a3))/(norm(cross(a2,a1)));
    Vector3D r3 = ro*(a3);
    Vector3D r2 = cross(r3,r1);

    R.set_row(0,{r1[0], r1[1], r1[2]});
    R.set_row(1,{r2[0], r2[1], r2[2]});
    R.set_row(2,{r3[0], r3[1], r3[2]});

    Vector3D transl = ro* mult(inverse(K),b);



    // In case you want to check some things: here are some outputs.

    //    std::ofstream stream_out1;
//    std::string fileOut1 = "c:\\tmp\\mmatrix.dat";
//    stream_out1.open(fileOut1);
//    if (stream_out1.is_open()) {
//        stream_out1 << M << std::endl;
//        stream_out1.close();
//    }
//    std::cout << N << "\n" << std::endl;
//    std::cout<<std::endl<<"M matrix"<< std::endl;
//    std::cout<<M<<std::endl;
//    std::cout<<"transl "<<transl<<std::endl;
//    std::cout<<"R "<<R<<std::endl;
//    std::cout<<"MM matrix"<<MM<<std::endl;
//    std::cout << "a1 "<< a1 << "\n" << std::endl;
//    std::cout << "a2 "<< a2 << "\n" << std::endl;
//    std::cout << "a3 "<< a3 << "\n" << std::endl;
//    std::cout << "b "<< b << "\n" << std::endl;
//    std::cout << "ro "<< ro << "\n" << std::endl;
//    std::cout << "u0 "<< u0 << "\n" << std::endl;
//    std::cout << "v0 "<< v0 << "\n" << std::endl;
//    std::cout << "costheta "<< costheta << "\n" << std::endl;
//    std::cout << "theta "<< theta << "\n" << std::endl;
//    std::cout << "sintheta "<< sintheta << "\n" << std::endl;
//    std::cout<<"K matrix"<<K<<std::endl;
//    std::cout<<"ro"<<ro<<std::endl;
//    std::cout<<"u0"<<u0<<std::endl;
//    std::cout<<"v0"<<v0<<std::endl;
//    std::cout<<"costheta"<<costheta<<std::endl;
//    std::cout<<"alpha"<<alpha<<std::endl;
//    std::cout<<"beta"<<beta<<std::endl;
    //    std::cout<<"r1 "<<r1<<std::endl;
//    std::cout<<"r2 "<<r2<<std::endl;
//    std::cout<<"r3 "<<r3<<std::endl;


    return true;
}










