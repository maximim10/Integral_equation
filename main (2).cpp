#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
const double E=0.001;
using namespace std;
double K_at_point(double x, double t);
double f_at_point(double x);
double get_u_x(double x, double a, double b, vector <double> &u){
    int n=u.size();
    int flag=0;
    vector <double> t;//={a, a+(b-a)/3.0, a+2*(b-a)/3.0, b};
    for (int i=0;i<n;i++){
        t.push_back(a+i*(b-a)*1.0/(n-1));
    }
    double u_1x=0;
    vector <double> c={3/8.0*(b-a)/(n-1), 9/8.0*(b-a)/(n-1), 9/8.0*(b-a)/(n-1), 3/8.0*(b-a)/(n-1)};
    for (int i=0;i<(n-1)/3*4;i++){
        u_1x+=u[i-i/4]*K_at_point(x,t[i-i/4])*c[i%4]/(n-1);
    }
    u_1x+=f_at_point(x);
    return u_1x;
}
double f_at_point(double x){
    return exp(-x);
}
double K_at_point(double x, double t){
    return -x*exp(t)*0.5;
}
vector<double> gauss(vector<vector<double>> a, vector<double> y, int n) {
    vector<double> x(n);
    double max;
    int k, index;
    const double eps1 = 1e-7;
    k = 0;
    while (k < n) {
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > max) {
                max = abs(a[i][k]);
                index = i;
            }
        }
        if (max < eps1) {
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            exit(-1);
        }
        for (int j = 0; j < n; j++) {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        for (int i = k; i < n; i++) {
            double temp1 = a[i][k];
            if (abs(temp1) < eps1) continue;
            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] / temp1;
            }
            y[i] = y[i] / temp1;
            if (i == k) continue;
            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] - a[k][j];
            }
            y[i] = y[i] - y[k];
        }
        k++;
    }
    for (k = n - 1; k >= 0; k--) {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}

void int_from_a_to_b(double a, double b, vector <double> &u, int n){

    vector <double> c={3/8.0*(b-a)/(n-1), 9/8.0*(b-a)/(n-1), 9/8.0*(b-a)/(n-1), 3/8.0*(b-a)/(n-1)};
    vector <double> t;//={a, a+(b-a)/3.0, a+2*(b-a)/3.0, b};
    for (int i=0;i<n;i++){
        t.push_back(a+i*(b-a)*1.0/((n-1)*1.0));
    }
    vector < vector <double> > SLAE (n, vector <double> (n,0));
    vector < vector <double> > SLAE_temp=SLAE;
    vector <double> f(n, 0);
    int flag=0;
    for (int i=0;i<n;i++){//cout<<i<<" ";
        for (int j=0;j<(n-1)/3*4;j++){
            if ((i==j-j/4)&&(flag==0)){
                //SLAE[i][j-j/4]+=1+c[j%4];
                SLAE[i][j-j/4]+=1+K_at_point(t[i],t[j-j/4])*c[j%4];
                flag=1;
            } else {
                //[i][j-j/4]+=c[j%4];
                SLAE[i][j-j/4]+=K_at_point(t[i],t[j-j/4])*c[j%4];
            }
            f[i]=f_at_point(t[i]);

        }flag=0;
    }
    u=gauss(SLAE, f, f.size());
    for (int j=0;j<n;j++){
        //cout<<f[j]<<" ";
    }
}
double Int(double a, double b, vector<double> &u_n1, vector<double> &u_n2) {
    vector <double> c;//={3/8.0*(b-a), 9/8.0*(b-a), 9/8.0*(b-a), 3/8.0*(b-a)};
    vector <double> t,t1,t2;//={a, a+(b-a)/3.0, a+2*(b-a)/3.0, b};
    for (int i=0;i<u_n1.size();i++){
        t.push_back(a+i*(b-a)/3.0);
    }
    int n = 4, flag=0;
    double I = 0;
    double r_b = a;
    vector <double> u1_temp(0), u2_temp(0);
    while (r_b + (b - a)/ n < b) {
        while (true) {
            double h = (b - a) / n;
            double u_1x = 0;
            int n_1 = u_n1.size();
            int n_2 = u_n2.size();
            for (int i=0;i<4;i++){
                t.push_back(r_b+i*h/3.0);
            }
            c={3/8.0*h/3.0, 9/8.0*h/3.0, 9/8.0*h/3.0, 3/8.0*h/3.0};
            for (int j=0;j<4;j++){
                u_1x+=pow(get_u_x(t[j],a,b,u_n1)-get_u_x(t[j],a,b,u_n2),2)*c[j];
            }
            double h2 = (b - a) / 2 * n;
            double u_2x = 0;
            for (int i=0;i<7;i++){
                t1.push_back(r_b+i*h/3.0);
            }
            for (int j=0;j<4;j++){
                u_2x+=pow(get_u_x(t1[j],a,b,u_n1)-get_u_x(t1[j],a,b,u_n2),2)*c[j]/2;
            }
            for (int j=0;j<4;j++){
                u_2x+=pow(get_u_x(t1[j+3],a,b,u_n1)-get_u_x(t1[j+3],a,b,u_n2),2)*c[j]/2;
            }
            double delta = u_2x - u_1x;
            if (abs(delta) >  E/n) {
                n*=2;
            } else {
                I += u_1x;
                break;
            }
        }
        r_b += (b - a)/ n;
        n /= 2;
    }
    while (r_b < b) {
        double h = b - r_b;
        while (true) {
            double u_1x = 0;
            int n_1 = u_n1.size();
            int n_2 = u_n2.size();
            for (int i=0;i<4;i++){
                t.push_back(r_b+i*h/3.0);
            }
            c={3/8.0*h/3.0, 9/8.0*h/3.0, 9/8.0*h/3.0, 3/8.0*h/3.0};
            for (int j=0;j<4;j++){
                u_1x+=pow(get_u_x(t[j],a,b,u_n1)-get_u_x(t[j],a,b,u_n2),2)*c[j];
            }
            double h2 = (b - a) / 2 * n;
            double u_2x = 0;
            for (int i=0;i<7;i++){
                t1.push_back(r_b+i*h2/3.0);
            }
            for (int j=0;j<4;j++){
                u_2x+=pow(get_u_x(t1[j],a,b,u_n1)-get_u_x(t1[j],a,b,u_n2),2)*c[j]/2;
            }
            for (int j=0;j<4;j++){
                u_2x+=pow(get_u_x(t1[j+3],a,b,u_n1)-get_u_x(t1[j+3],a,b,u_n2),2)*c[j]/2;
            }
            double delta = u_2x - u_1x;
            if (abs(delta) > E){
                h/=2;
            } else {
                I += u_1x;
                break;
            }
        }
        r_b += h;
    }
    I=pow(I,0.5);
    cout << "I: " << I << endl;
    return I;
}
int main()
{
    int N=4;
    int a=0, b=1;
    vector <double> u, u_prev;
    int_from_a_to_b(a,b,u,3*N+1);
    while (true) {
        int_from_a_to_b(a, b, u_prev, 3*N+1);
        int_from_a_to_b(a, b, u, 3*(N*2)+1);
        if (Int(a, b, u_prev, u) < E) {
            u = u_prev;
            break;
        } else {
            N*=2;
        }
    }
    for (int i=0; i<u.size(); i++){
        cout<<u[i]<<" ";
    }
    return 0;
}
