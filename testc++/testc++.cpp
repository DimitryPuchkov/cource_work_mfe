#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

class mesh		//сетка и все параметры задачи
{
public:
   long count_r, count_z, count_t;  //количество узлов по r,z, t
   double *r, *z, *t;			//массивы координат и времени
   mesh(string file_name);	//конструктор mesh
   ~mesh();
   double lambda();	//выдают значение лямбда
   //double gamma(double r, double z);
   double f(double r, double z, double t);
   double ug(double r, double z, double t);
   double u(double r, double z, double t);
   double hi();
   double sigma();
};

mesh::mesh(string file_name)
{
   // считываем количество узлов и их координаты
   fstream in;
   long i;

   in.open(file_name);

   in >> count_r >> count_z>>count_t;
   r = new double[count_r];
   z = new double[count_z];
   t = new double[count_t];

   for (i = 0; i < count_r; i++)
      in >> r[i];
   for (i = 0; i < count_z; i++)
      in >> z[i];
   for (i = 0; i < count_t; i++)
      in >> t[i];

   in.close();
}
mesh::~mesh()
{
   delete[] r;
   delete[] z;
}

double mesh::lambda()
{
   return 1.0;
}
//double mesh::gamma(double r, double z)
//{
//   return 1.0;
//}
double mesh::f(double r, double z, double t)
{
   return 0;
}  
double mesh::ug(double r, double z, double t)
{
   return z;
}
double mesh::u(double r, double z, double t)
{
   return z;
}
double mesh::hi()
{
   return 1.0;
}
double mesh::sigma()
{
   return 1.0;
}

class sparse_sistem  //разряженная система
{
private:
   long dim;  //разреженно-строчная матрица
   vector<long> ig, jg;
   vector<double> di, gl;
   vector<double> b;   //вектор правой части
public:
   sparse_sistem(class mesh &MESH); //конструктор
   ~sparse_sistem();
   void make_slau(class mesh &MESH, vector<double> qj_2, vector<double> qj_1, double tj_2, double tj_1, double tj);  //сборка матрицы
   void add(double value, long i, long j); //вклад в матрицу
   void add(double value, long i);   //вклад в вектор
   void boundary1(long i, double value);  //учёт первого краевого с симметризацией матрицы
   void matrixvectormult(vector<double> *a, vector<double> *b);   //A*a=b
   void inversediagonalvectormult(vector<double> *a, vector<double> *b);   //D^(-1)*a=b
   double dotproduct(vector<double> *a, vector<double> *b);   //(a,b)
   void MSG(vector<double> *result);    //сопряжённые градиенты
   long fact(long a); // a!
   double integrate_L(long nu[], double det); // вычисляет интеграл по базисным функциям nu - массив с степенями
   double calc_m(long i, long j, double r[], double z[], double det);
   double calc_b(long i, double r[], double z[], double det, double tj, double qj_2[], double qj_1[], double tMhi_2, double tMhi_1, double tMsig_2, double tMsig_1, mesh &MESH);
};
sparse_sistem::sparse_sistem(mesh &MESH)
{
   // строим портрет матрицы
   long i, j, k, num, *tmp, **tmp_n;
   // (?) tmp - массив с количеством внедиагональных элементов построчно
   // tmp_n - двмерный массив где первое измерение - номер строки, второе - номер столбцов с ненулевыми элементами

   dim = MESH.count_r * MESH.count_z;

   tmp = new long[dim];

   tmp_n = new long *[dim];
   for (i = 0; i < dim; i++)
      tmp_n[i] = new long[3];

   tmp[0] = 0;
   for (i = 1; i < MESH.count_r - 1; i++)
      for (j = 1; j < MESH.count_z - 1; j++)
      {
         num = i * MESH.count_z + j; // (?) num = номер узла
         tmp[num] = 3;
         tmp_n[num][0] = num - MESH.count_z;
         tmp_n[num][1] = tmp_n[num][0] + 1;
         tmp_n[num][2] = num - 1;
      }

   for (i = 1; i < MESH.count_z; i++)
   {
      tmp[i] = 1;
      tmp_n[i][0] = i - 1;
   }

   for (i = 1; i < MESH.count_z - 1; i++)
   {
      num = MESH.count_z * (MESH.count_r - 1) + i;
      tmp[num] = 3;
      tmp_n[num][0] = num - MESH.count_z;
      tmp_n[num][1] = tmp_n[num][0] + 1;
      tmp_n[num][2] = num - 1;
   }

   for (i = 1; i < MESH.count_r; i++)
   {
      num = MESH.count_z * i;
      tmp[num] = 2;
      tmp_n[num][0] = num - MESH.count_z;
      tmp_n[num][1] = tmp_n[num][0] + 1;
   }

   for (i = 1; i < MESH.count_r; i++)
   {
      num = MESH.count_z * (i + 1) - 1;
      tmp[num] = 2;
      tmp_n[num][0] = num - MESH.count_z;
      tmp_n[num][1] = num - 1;
   }

   for (i = 0, j = 0; i < dim; i++)
   {
      j += tmp[i];
   }

   di.resize(dim);
   gl.resize(j);
   b.resize(dim);
   ig.resize(dim + 1);
   jg.resize(j);

   ig[0] = 0;
   for (i = 0, j = 0; i < dim; i++)
   {
      di[i] = 0.0;
      b[i] = 0.0;

      ig[i + 1] = ig[i] + tmp[i];
      for (k = 0; k < tmp[i]; k++, j++)
      {
         jg[j] = tmp_n[i][k];
         gl[j] = 0.0;
      }
   }

   delete[] tmp;

   for (i = 0; i < dim; i++)
      delete[] tmp_n[i];

   delete[] tmp_n;

   ofstream out;
   out.open("matr.txt");
   out << "ig" << endl;
   for (i = 0; i < dim + 1; i++)
      out << i << " " << ig[i] << endl;
   out << endl;

   out << "jg" << endl;
   for (i = 0; i < ig[dim]; i++)
      out << i << " " << jg[i] << endl;
   out << endl;
   //out.clear();//????????
   out.close();
}
sparse_sistem::~sparse_sistem()
{
   ig.clear();
   jg.clear();
   di.clear();
   gl.clear();
   b.clear();
}
void sparse_sistem::make_slau(mesh &MESH, vector<double> qj_2, vector<double> qj_1, double tj_2, double tj_1, double tj)
{
   long i, j;
   double a11, a12, a21, a22, a31, a32;
   double det;
   double h, sumr;
   double b;
   double r[3] = { 0, 0, 0 }, z[3] = { 0, 0, 0 }, n[3] = { 0, 0, 0 }, qj1_loc[3] = {0, 0, 0}, qj2_loc[3]={0, 0, 0};
   double dt, dt0, dt1, tMhi, tMhi_1, tMhi_2, tMsig, tMsig_1, tMsig_2, M_coef;

   dt = tj-tj_2;
   dt1 = tj_1 - tj_2;
   dt0 = tj-tj_1;
   tMhi = 2 / (dt * dt0);
   tMhi_1 = 2/(dt1*dt0);
   tMhi_2 = -2/(dt1*dt);
   tMsig = (dt + dt0) / (dt * dt0);
   tMsig_1 = dt/(dt1*dt0);
   tMsig_2 = -dt0/(dt*dt1);
   // For M in left part
   M_coef = (tMhi * MESH.hi() + tMsig * MESH.sigma());

   for (i = 0; i < MESH.count_r - 1; i++)//генерация вкладов с элементов
      for (j = 0; j < MESH.count_z - 1; j++)
      {
         //вклад с нижнего треугольника
         n[0] = i * MESH.count_z + j;
         n[1] = n[0] + 1;
         n[2] = n[0] + MESH.count_z;
         r[0] = MESH.r[i];
         r[1] = MESH.r[i];
         r[2] = MESH.r[i + 1];
         z[0] = MESH.z[j];
         z[1] = MESH.z[j + 1];
         z[2] = MESH.z[j];

         qj1_loc[0] = qj_1[n[0]];
         qj1_loc[1] = qj_1[n[1]];
         qj1_loc[2] = qj_1[n[2]];
         qj2_loc[0] = qj_2[n[0]];
         qj2_loc[1] = qj_2[n[1]];
         qj2_loc[2] = qj_2[n[2]];

         det = (z[1] - z[0]) * (r[2] - r[0]) - (z[2] - z[0]) * (r[1] - r[0]); // |detD|
         a11 = (r[1] - r[2]) / det;
         a12 = (z[2] - z[1]) / det;
         a21 = (r[2] - r[0]) / det;
         a22 = (z[0] - z[2]) / det;
         a31 = (r[0] - r[1]) / det;
         a32 = (z[1] - z[0]) / det;

         sumr = MESH.lambda() * (r[0] + r[1] + r[2]) * fabs(det) / 12;
         
         add(sumr * (a11 * a11 + a12 * a12) + M_coef * calc_m(1, 1, r, z, det), n[0], n[0]);

         add(sumr * (a11 * a21 + a22 * a12) + M_coef * calc_m(2, 1, r, z, det), n[1], n[0]);
         add(sumr * (a21 * a21 + a22 * a22) + M_coef * calc_m(2, 2, r, z, det), n[1], n[1]);

         add(sumr * (a31 * a11 + a32 * a12) + M_coef * calc_m(3, 1, r, z, det), n[2], n[0]);
         add(sumr * (a31 * a21 + a32 * a22) + M_coef * calc_m(3, 2, r, z, det), n[2], n[1]);
         add(sumr * (a31 * a31 + a32 * a32) + M_coef * calc_m(3, 3, r, z, det), n[2], n[2]);

         add(calc_b(1, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[0]);
         add(calc_b(2, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[1]);
         add(calc_b(3, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[2]);



         //Вклад с верхнего треугольника
         n[0] = i * MESH.count_z + j + 1;
         n[2] = n[0] + MESH.count_z;
         n[1] = n[2] - 1;
         r[0] = MESH.r[i];
         r[1] = MESH.r[i + 1];
         r[2] = MESH.r[i + 1];
         z[0] = MESH.z[j + 1];
         z[1] = MESH.z[j];
         z[2] = MESH.z[j + 1];

         det = (z[1] - z[0]) * (r[2] - r[0]) - (z[2] - z[0]) * (r[1] - r[0]); // |detD|
         a11 = (r[1] - r[2]) / det;
         a12 = (z[2] - z[1]) / det;
         a21 = (r[2] - r[0]) / det;
         a22 = (z[0] - z[2]) / det;
         a31 = (r[0] - r[1]) / det;
         a32 = (z[1] - z[0]) / det;
         sumr = MESH.lambda() * (r[0] + r[1] + r[2]) * fabs(det) / 12;

         add(sumr * (a11 * a11 + a12 * a12) + M_coef * calc_m(1, 1, r, z, det), n[0], n[0]);

         add(sumr * (a11 * a21 + a22 * a12) + M_coef * calc_m(2, 1, r, z, det), n[1], n[0]);
         add(sumr * (a21 * a21 + a22 * a22) + M_coef * calc_m(2, 2, r, z, det), n[1], n[1]);

         add(sumr * (a31 * a11 + a32 * a12) + M_coef * calc_m(3, 1, r, z, det), n[2], n[0]);
         add(sumr * (a31 * a21 + a32 * a22) + M_coef * calc_m(3, 2, r, z, det), n[2], n[1]);
         add(sumr * (a31 * a31 + a32 * a32) + M_coef * calc_m(3, 3, r, z, det), n[2], n[2]);


         //правая часть
         add(calc_b(1, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[0]);
         add(calc_b(2, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[1]);
         add(calc_b(3, r, z, det, tj, qj2_loc, qj1_loc, tMhi_2, tMhi_1, tMsig_2, tMsig_1, MESH), n[2]);

      }


   //Учёт первого краевого
   // Верхняя граница
   for (i = 0; i < MESH.count_z; i++)
      boundary1((MESH.count_r - 1) * (MESH.count_z) + i, MESH.ug(MESH.r[MESH.count_r - 1], MESH.z[i], tj));

   // Нижняя граница
   for (i = 0; i < MESH.count_z; i++)
      boundary1(i, MESH.ug(MESH.r[0], MESH.z[i], tj));

   //Левая граница
   for (i = 0; i < MESH.count_r; i++)
      boundary1(MESH.count_z * i, MESH.ug(MESH.r[i], MESH.z[0], tj));

   // Правая граница
   for (i = 1; i < MESH.count_r + 1; i++)
      boundary1(MESH.count_z * i - 1, MESH.ug(MESH.r[i - 1], MESH.z[MESH.count_z - 1], tj));



   ofstream out;
   out.open("slau.txt");


   out << "IG\n";
   for (i = 0; i < this->dim + 1; i++)
      out << this->ig[i] << " ";
   out << "\nJG\n";
   for (i = 0; i < this->ig[this->dim]; i++)
      out << this->jg[i] << " ";
   out << "\nDI\n";
   for (i = 0; i < this->dim; i++)
      out << this->di[i] << " ";
   out << "\nGL\n";
   for (i = 0; i < this->ig[this->dim]; i++)
      out << this->gl[i] << " ";
   out << "\nB\n";
   for (i = 0; i < this->dim; i++)
      out << this->b[i] << " ";

   out.close();
}
void sparse_sistem::add(double value, long i, long j)
{
   bool flag;
   long M, m, k;

   flag = true;
   if (i == j)
   {
      this->di[i] += value;
      flag = false;
   }
   else
   {
      M = (i > j ? i : j);
      m = (i < j ? i : j);
      for (k = this->ig[M]; flag && (k < this->ig[M + 1]); k++)
         if (this->jg[k] == m)
         {
            this->gl[k] += value;
            flag = false;
         }
   }

   if (flag)
   {
      cout << "Error add" << i << "\t" << j << "\n";
      exit(2);
   }

}
void sparse_sistem::add(double value, long i)
{
   this->b[i] += value;
}
void sparse_sistem::boundary1(long i, double value)
{
   long j, k;

   this->di[i] = 1.0;
   this->b[i] = value;

   for (j = this->ig[i]; j < this->ig[i + 1]; j++)
   {
      this->b[this->jg[j]] -= this->gl[j] * value;
      this->gl[j] = 0.0;
   }

   for (j = i + 1; j < this->dim; j++)
   {
      for (k = this->ig[j]; k < this->ig[j + 1]; k++)
         if (this->jg[k] == i)
         {
            this->b[j] -= this->gl[k] * value;
            this->gl[k] = 0.0;
         }
   }
}
void sparse_sistem::matrixvectormult(vector<double> *a, vector<double> *b)
{
   long i, j;

   for (i = 0; i < this->dim; i++)
      (*b)[i] = 0.0;

   for (i = 0; i < this->dim; i++)
   {
      (*b)[i] += (*a)[i] * this->di[i];
      for (j = this->ig[i]; j < this->ig[i + 1]; j++)
      {
         (*b)[i] += (*a)[this->jg[j]] * this->gl[j];
         (*b)[this->jg[j]] += (*a)[i] * this->gl[j];
      }
   }
}
void sparse_sistem::inversediagonalvectormult(vector<double> *a, vector<double> *b)
{
   long i;
   for (i = 0; i < this->dim; i++)
      (*b)[i] = (*a)[i] / this->di[i];
}
double sparse_sistem::dotproduct(vector<double> *a, vector<double> *b)
{
   long i;
   double res;

   res = 0.0;
   for (i = 0; i < this->dim; i++)
      res += (*a)[i] * (*b)[i];

   return res;
}
long sparse_sistem::fact(long a)
{
   long f = 1;
   for (long i = 1; i <= a; i++) f *= i;
   return f;
}
void sparse_sistem::MSG(vector<double> *result)
{
   long i, j, max_iter;
   double eps, nn, nb, rr, rr1, alpfa, betta;
   //double *x, *r, *p, *Ap, *Dr;
   vector<double> x, r, p, Ap, Dr;
   eps = 1.e-9;
   max_iter = 20000;

   x.resize(dim);
   r.resize(dim);
   p.resize(dim);
   Ap.resize(dim);
   Dr.resize(dim);

   for (i = 0; i < this->dim; i++)
      x[i] = 0.0;

   nb = sqrt(this->dotproduct(&b, &b));

   this->matrixvectormult(&x, &r);
   for (i = 0; i < this->dim; i++)
      r[i] = this->b[i] - r[i];

   this->inversediagonalvectormult(&r, &p);
   this->inversediagonalvectormult(&r, &Dr);

   j = 0;
   rr = this->dotproduct(&Dr, &r);
   cout << rr << endl;
   nn = 1.0;
   while (nn > eps && j < max_iter)
   {
      this->matrixvectormult(&p, &Ap);
      alpfa = rr / this->dotproduct(&Ap, &p);

      for (i = 0; i < this->dim; i++)
      {
         x[i] = x[i] + alpfa * p[i];
         r[i] = r[i] - alpfa * Ap[i];
      }

      this->inversediagonalvectormult(&r, &Dr);
      rr1 = this->dotproduct(&Dr, &r);
      betta = rr1 / rr;

      for (i = 0; i < this->dim; i++)
         p[i] = Dr[i] + betta * p[i];

      rr = rr1;
      nn = sqrt(this->dotproduct(&r, &r)) / nb;
      cout << j++ << "\t" << alpfa << "\t" << betta << "\t" << nn << endl;
   }

   for (i = 0; i < this->dim; i++)
      (*result)[i] = x[i];


   x.clear();
   r.clear();
   p.clear();
   Ap.clear();
   Dr.clear();
}
double sparse_sistem::integrate_L(long nu[], double det)
{
   double chisl = 1, znamen = 2;
   for (long i = 0; i < 3; i++)
   {
      chisl *= fact(nu[i]);
      znamen += nu[i];
   }
   znamen = fact(znamen);
   return (chisl * fabs(det)) / (znamen * 2);
}
double sparse_sistem::calc_m(long i, long j, double r[], double z[], double det)
{
   long nu[3];

   double sum = 0;
   for (int m = 0; m < 3; m++)
   {
      nu[0] = 0;
      nu[1] = 0;
      nu[2] = 0;
      nu[i - 1]++;
      nu[j - 1]++;
      nu[m]++;

      sum += r[m] * integrate_L(nu, det);
   }
   return sum;
}
double sparse_sistem::calc_b(long i, double r[], double z[], double det, double tj, double qj_2[], double qj_1[], double tMhi_2, double tMhi_1, double tMsig_2, double tMsig_1, mesh &MESH)
{
   double sum = 0;


   //bi
   for (long k = 0; k < 3; k++)
      sum += MESH.f(r[k], z[k], tj) * calc_m(i, k + 1, r, z, det);

   // Mhi qj_2
   for (long k = 0; k < 3; k++)
      sum += MESH.hi() * tMhi_2 * qj_2[k] * calc_m(i, k + 1, r, z, det);

   // Msig qj_2
   for (long k = 0; k < 3; k++)
      sum += MESH.sigma() * tMsig_2 * qj_2[k] * calc_m(i, k + 1, r, z, det);

   // Mhi qj_1
   for (long k = 0; k < 3; k++)
      sum += MESH.hi() * tMhi_1 * qj_1[k] * calc_m(i, k + 1, r, z, det);

   // Msig qj_1
   for (long k = 0; k < 3; k++)
      sum += MESH.sigma() * tMsig_1 * qj_1[k] * calc_m(i, k + 1, r, z, det);
   return sum;
}


int main()
{
   long i, j;
   double tj_2, tj_1, tj;
   vector<double> res, anal, raz, qj_2, qj_1;
   class mesh MESH("mesh.txt");
   class sparse_sistem SISTEM(MESH);

   res.resize(MESH.count_r * MESH.count_z);
   anal.resize(MESH.count_r * MESH.count_z);
   raz.resize(MESH.count_r * MESH.count_z);
   qj_2.resize(MESH.count_r * MESH.count_z);
   qj_1.resize(MESH.count_r * MESH.count_z);

   //for (int i=0; i < qj_2.size(); i++)
   //{
   //   qj_1[i] = MESH.u(MESH.r[i], MESH.z[j], MESH.t[0]);
   //   qj_2[i] = MESH.u(MESH.r[i], MESH.z[j], MESH.t[1]);
   //}
   for(int r = 0; r<MESH.count_r; r++)
      for (int z = 0; z < MESH.count_z; z++)
      {
         qj_1[r*MESH.count_z+z] = MESH.u(MESH.r[r], MESH.z[z], MESH.t[0]);
         qj_2[r*MESH.count_z+z] = MESH.u(MESH.r[r], MESH.z[z], MESH.t[1]);
      }

   ofstream out;
   out.open("res.xls");
   out <<"t\t" <<"r\t" << "z\t" << "u(r,z)\t" << "u*(r,z)\t" << "(u*-u)(r,z)\n";

   for (int j = 2; j < MESH.count_t; j++)
   {
      tj_2 = MESH.t[j-2];
      tj_1 = MESH.t[j-1];
      tj = MESH.t[j];
      SISTEM.make_slau(MESH, qj_2, qj_1, tj_2, tj_1, tj);
      SISTEM.MSG(&res);
      qj_2 = qj_1;
      qj_1 = res;

      for (int r = 0; r < MESH.count_r; r++)
         for (int z = 0; z < MESH.count_z; z++)
         {
            anal[r * MESH.count_z + z] = MESH.u(MESH.r[r], MESH.z[z], tj);
            raz[r * MESH.count_z + z] = fabs(res[r * MESH.count_z + z] - anal[r * MESH.count_z + z]);
            out << scientific << tj<<"\t" << MESH.r[r] << "\t" << MESH.z[z] << "\t" << res[r * MESH.count_z + z] << "\t" << anal[r * MESH.count_z + z] << "\t" << raz[r * MESH.count_z + z] << endl;
         }
      out << endl<<"ERROR: \t" << sqrt(SISTEM.dotproduct(&raz, &raz)) << endl;
   }
   
   

   

   
   

   out.close();
   res.clear();
   anal.clear();
   raz.clear();
   return 0;
}
