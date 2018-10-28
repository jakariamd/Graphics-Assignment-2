#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stack>
#include <vector>
#include <iomanip>
#include<math.h>

#define PI 3.14159265

using namespace std;
//typedef vector<vector<double> > mat;

ifstream infile;
ofstream outfile1,outfile2,outfile3;


struct point
{
    double x,y,z;
};

struct point eye,look,up;
struct point l,r,u;
double fovX,fovY,aspectRatio,near,far;


class transMat
{
public:
    double mat[4][4];
    bool mark;
    void setIMat();
    void set0Mat();
    void normalize();
    void printMat();
};

void transMat::normalize(){
    if(mat[3][3] != 1){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                mat[i][j] = mat[i][j] /mat[3][3] ;
            }
        }
    }
}

void transMat::setIMat()
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i==j)
                mat[i][j] = 1;
            else
                mat[i][j] = 0;
        }
    }
}

void transMat::set0Mat()
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
                mat[i][j] = 0;
        }
    }
}

void transMat::printMat()
{
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            cout << mat[i][j] << " " ;
        }
        cout << "\n" ;
    }
}

stack<transMat> mystack;
transMat V,P;


void getEyePos(string line)
{
    istringstream iss;
    iss.str(line);
    iss >> eye.x >> eye.y >> eye.z; // cout << eye.x << " " << eye.y << " " << eye.z << endl;
}

void getLookPos(string line)
{
    istringstream iss;
    iss.str(line);
    iss >> look.x >> look.y >> look.z; //cout << look.x << " " << look.y << " " << look.z << endl;
}

void getUpPos(string line)
{
    istringstream iss;
    iss.str(line);
    iss >> up.x >> up.y >> up.z; //cout << up.x << " " << up.y << " " << up.z << endl;
}

void getPerspective(string line)
{
    istringstream iss;
    iss.str(line);
    iss >> fovY >> aspectRatio >> near >> far; //cout << fovY << " " << fovR << " " << near << " " << far << endl;
}

void transPoint(double p[]){
   transMat tmat = mystack.top();
    double temppoint[4];
    for(int i=0; i< 4; i++){
        temppoint[i] = 0;
        for(int j = 0; j < 4; j++){
            temppoint[i] +=  tmat.mat[i][j] * p[j] ;
        }
    }
    for(int i=0; i < 3; i++){
        outfile1 << std::fixed << std::setprecision(7) << temppoint[i] <<  " ";
    }
    outfile1 << endl;

    double viewpoint[4];
    for(int i=0; i< 4; i++){
        viewpoint[i] = 0;
        for(int j = 0; j < 4; j++){
            viewpoint[i] +=  V.mat[i][j] * temppoint[j] ;
        }
    }
    for(int i=0; i < 3; i++){
        outfile2 << std::fixed << std::setprecision(7) << viewpoint[i] <<  " ";
    }
    outfile2 << endl;


    double projecpoint[4];
    for(int i=0; i< 4; i++){
        projecpoint[i] = 0;
        for(int j = 0; j < 4; j++){
            projecpoint[i] +=  P.mat[i][j] * viewpoint[j] ;
        }
    }
    for(int i=0; i < 3; i++){
        outfile3 << std::fixed << std::setprecision(7) << projecpoint[i]/projecpoint[3] <<  " ";
    }
    outfile3 << endl;
}

struct point calcRot(double x[], double a[], double theta){
    point temp;
    float costheta =  cos( theta * PI / 180.0 );
    float sintheta =  sin( theta * PI / 180.0 );
    double dotproduct = x[0]*a[0]+x[1]*a[1]+x[2]*a[2];

    temp.x = costheta*x[0] + (1-costheta)*dotproduct*a[0] + sintheta*(a[1]*x[2]-a[2]*x[1]);
    temp.y = costheta*x[1] + (1-costheta)*dotproduct*a[1] + sintheta*(a[2]*x[0]-a[0]*x[2]);
    temp.z = costheta*x[2] + (1-costheta)*dotproduct*a[2] + sintheta*(a[0]*x[1]-a[1]*x[0]);

    //cout  << std::fixed << std::setprecision(1) << temp.x << " " << temp.y << " " << temp.z << " " << endl;
    return temp;

}

void viewTransform(){
    l.x = look.x - eye.x;
    l.y = look.y - eye.y;
    l.z = look.z - eye.z;
    double val = pow(l.x*l.x+l.y*l.y+l.z*l.z,0.5);
    l.x = l.x / val;
    l.y = l.y / val;
    l.z = l.z / val;

    r.x = l.y*up.z - l.z*up.y;
    r.y = l.z*up.x - l.x*up.z;
    r.z = l.x*up.y - l.y*up.x;
    val = pow(r.x*r.x+r.y*r.y+r.z*r.z,0.5);
    r.x = r.x / val;
    r.y = r.y / val;
    r.z = r.z / val;

    u.x = r.y*l.z - r.z*l.y;
    u.y = r.z*l.x - r.x*l.z;
    u.z = r.x*l.y - r.y*l.x;

    transMat R,T;
    T.setIMat();
    T.mat[0][3] = -eye.x;
    T.mat[1][3] = -eye.y;
    T.mat[2][3] = -eye.z;

    //T.printMat();
    //cout << endl;

    R.setIMat();

    R.mat[0][0] = r.x ;
    R.mat[0][1] = r.y ;
    R.mat[0][2] = r.z ;

    R.mat[1][0] = u.x ;
    R.mat[1][1] = u.y ;
    R.mat[1][2] = u.z ;

    R.mat[2][0] = -l.x ;
    R.mat[2][1] = -l.y ;
    R.mat[2][2] = -l.z ;

    //R.printMat();

    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            for(int k = 0; k < 4; ++k){
                V.mat[i][j] += R.mat[i][k] * T.mat[k][j];
            }
        }
    }
    //cout << endl;
    //V.printMat();
}

void projectionTransform(){
    fovX = fovY * aspectRatio;
    double t = near * tan((fovY/2.0)*PI/180.0);
    double r = near * tan((fovX /2.0)*PI/180.0);

    P.set0Mat();
    P.mat[0][0] = near/r;
    P.mat[1][1] = near/t;
    P.mat[2][2] = -(far+near)/(far-near);
    P.mat[2][3] = -(2*far*near)/(far-near);
    P.mat[3][2] = -1;

    //P.printMat();
}

int main()
{
    infile.open("scene.txt");
    outfile1.open("stage1.txt");
    outfile2.open("stage2.txt");
    outfile3.open("stage3.txt");
    string line;
    //double x,y,z,w;

    getline(infile, line);
    getEyePos(line);
    getline(infile, line);
    getLookPos(line);
    getline(infile, line);
    getUpPos(line);
    getline(infile, line);
    getPerspective(line);

    viewTransform();
    projectionTransform();


    transMat idt;
    idt.setIMat();
    idt.mark= false;
    mystack.push(idt);

    bool ispush = false;

    while(true)
    {
        getline(infile, line);
        if(line == "triangle")
        {
            for(int i=0; i < 3; i++){
                getline(infile, line);
                double endp[4];
                istringstream iss;
                iss.str(line);
                iss >> endp[0] >> endp[1] >> endp[2];
                endp[3] = 1;
                transPoint(endp);
            }
            outfile1 << endl;
            outfile2 << endl;
            outfile3 << endl;
        }
        else if(line == "translate")
        {
            getline(infile, line);
            double trxp[4];
            istringstream iss;
            iss.str(line);
            iss >> trxp[0] >> trxp[1] >> trxp[2];
            trxp[3] = 1;

            transMat traslate;
            traslate.mark = ispush;
            ispush = false;

            traslate.setIMat();
            for(int i = 0; i < 4; i++)
            {
                traslate.mat[i][3] = trxp[i];
            }
            //traslate.printMat();
            transMat mult,top;
            top =  mystack.top();
            mult.set0Mat();
            mult.mark = traslate.mark;
            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    for(int k = 0; k < 4; ++k){
                        mult.mat[i][j] += top.mat[i][k] * traslate.mat[k][j];
                    }
                }
            }
            mult.normalize();
            mystack.push(mult);
            //cout << "translate" << endl; mult.printMat(); cout << mult.mark << endl;
        }

        else if(line=="scale"){
            getline(infile, line);
            double trxp[4];
            istringstream iss;
            iss.str(line);
            iss >> trxp[0] >> trxp[1] >> trxp[2];
            trxp[3] = 1;

            transMat scal;
            scal.mark = ispush;
            ispush = false;

            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(i==j)
                        scal.mat[i][j] = trxp[i] ;
                    else
                        scal.mat[i][j] = 0;
                }
            }

            transMat mult,top;
            top =  mystack.top();
            mult.set0Mat();
            mult.mark = scal.mark;

            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    for(int k = 0; k < 4; ++k){
                        mult.mat[i][j] += top.mat[i][k] * scal.mat[k][j];
                    }
                }
            }
            mult.normalize();
            mystack.push(mult);
            //cout << "scale" << endl;scal.printMat();cout << scal.mark << endl;
        }

        else if(line == "rotate"){
            getline(infile, line);
            double angle, a[3];
            istringstream iss;
            iss.str(line);
            iss >> angle >> a[0] >> a[1] >> a[2];

            double val = pow(a[0]*a[0]+a[1]*a[1]+a[2]*a[2],0.5);
            a[0] = a[0]/val;
            a[1] = a[1]/val;
            a[2] = a[2]/val;

            transMat rot;
            rot.mark = ispush;
            ispush = false;

            struct point c[3];
            double x[] = {1.0,0.0,0.0};
            c[0] = calcRot(x,a,angle);
            double y[] = {0.0,1.0,0.0};
            c[1]= calcRot(y,a,angle);
            double z[] = {0.0,0.0,1.0};
            c[2] = calcRot(z,a,angle);

            for(int i = 0; i < 3; i++){
                rot.mat[0][i] = c[i].x;
                rot.mat[1][i] = c[i].y;
                rot.mat[2][i] = c[i].z;
                rot.mat[3][i] = 0;
                rot.mat[i][3] = 0;
                rot.mat[3][3] = 1;
            }

            transMat mult,top;
            top =  mystack.top();
            mult.set0Mat();
            mult.mark = rot.mark;

            for(int i = 0; i < 4; ++i){
                for(int j = 0; j < 4; ++j){
                    for(int k = 0; k < 4; ++k){
                        mult.mat[i][j] += top.mat[i][k] * rot.mat[k][j];
                    }
                }
            }
            mult.normalize();
            mystack.push(mult);
        }

        else if(line == "push"){
            ispush = true;
        }

        else if(line == "pop"){
            while(true){
                if(mystack.empty()){
                    cout << "Error!";
                    return 1;
                }
                transMat top = mystack.top();
                if(top.mark){
                    mystack.pop();
                    break;
                }
                mystack.pop();
            }
        }

        else if(line == "end")
        {
            break;
        }
    }

    infile.close();
    outfile1.close();
    outfile2.close();
    outfile3.close();
}
